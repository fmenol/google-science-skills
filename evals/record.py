# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Coordinated, resumable recording pass over every ESM skill eval.

The Biohub account has a hard allowance of 100 credits/day, which makes a naive
"run all 11 live evals" impossible to repeat and easy to exhaust. This script
solves that:

* It runs the evals SERIALLY. Running them in parallel makes them race for the
  same credit pool, so a burst exhausts the quota and most evals fail for a
  reason that has nothing to do with the skill under test.
* It records every API response into a cassette directory. Subsequent runs
  replay from disk for free, so the suite becomes repeatable and CI-able.
* It is RESUMABLE. Responses are written as they arrive, so if the quota dies
  mid-pass the work already paid for is kept. Re-running tomorrow replays
  everything cached and only spends credits on what is still missing. The pass
  therefore converges even if the full suite costs more than one day's
  allowance.

Usage
  uv run --no-project evals/record.py              # record what is missing, then verify
  uv run --no-project evals/record.py --replay     # verify from cassettes only (free)
  uv run --no-project evals/record.py --wait       # block until the quota resets, then record
"""

# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import sys
import time
import urllib.error
import urllib.request

EVALS = pathlib.Path(__file__).resolve().parent
REPO = EVALS.parent
CASSETTE = EVALS / (__import__('os').environ.get('CASSETTE_DIR') or 'cassettes')
STATE = EVALS / '_record_state.json'

# Cheapest first, so a truncated pass still banks as many green skills as
# possible. `esm_protein_tracks` needs no credits at all.
ORDER = [
    'esm_protein_tracks',            # 0 billed calls (local structure parsing)
    'esm3_function_prediction',      # a few generate() calls
    'esm3_secondary_structure_sasa',                       # generate() on ss + sasa tracks
    'esmc_protein_embeddings',       # ~12 logits calls
    'esmc_embedding_layer_sweep',    # ~16 logits calls
    'esmc_sae_feature_interpretation',
    'esm3_guided_generation',
    'esm3_inverse_folding',
    'esmfold2_structure_prediction',
    'esmfold2_binder_screening',
    'esm3_protein_design',
    'esmc_mutation_effect_scoring',  # ~160 logits calls (leave-one-out scans)
]


def quota_available() -> bool:
  """Cheap probe: a 429 costs nothing, so this is free to call."""
  key = _api_key()
  payload = json.dumps({
      'model': 'esm3-open-2024-03',
      'inputs': {'sequence': '____'},
      'track': 'sequence', 'invalid_ids': [], 'schedule': 'cosine',
      'num_steps': 1, 'temperature': 1.0, 'top_p': 1.0,
      'condition_on_coordinates_only': True, 'strategy': 'entropy',
      'temperature_annealing': False, 'only_compute_backbone_rmsd': False,
  }).encode()
  req = urllib.request.Request(
      'https://biohub.ai/api/v1/generate', data=payload, method='POST'
  )
  req.add_header('Authorization', f'Bearer {key}')
  req.add_header('Content-Type', 'application/json')
  try:
    with urllib.request.urlopen(req, timeout=60):
      return True
  except urllib.error.HTTPError as exc:
    if exc.code == 429:
      return False
    return True  # any other error is not a quota problem
  except OSError:
    return False


def _api_key() -> str:
  if os.environ.get('BIOHUB_API_KEY'):
    return os.environ['BIOHUB_API_KEY']
  for env_path in (REPO / '.env', pathlib.Path.home() / '.env'):
    if env_path.is_file():
      for line in env_path.read_text(encoding='utf-8').splitlines():
        if line.startswith('BIOHUB_API_KEY='):
          return line.split('=', 1)[1].strip().strip('\'"')
  raise SystemExit('BIOHUB_API_KEY not found')


def run_eval(skill: str, mode: str, timeout: int) -> tuple[bool, str, float]:
  env = {
      **os.environ,
      'BIOHUB_CASSETTE': str(CASSETTE),
      'BIOHUB_CASSETTE_MODE': mode,
      # Serial pass: keep the request rate gentle; nothing is competing.
      'BIOHUB_QPS': os.environ.get('BIOHUB_QPS', '4'),
  }
  started = time.time()
  proc = subprocess.run(
      ['uv', 'run', '--no-project', '--quiet', str(EVALS / f'eval_{skill}.py')],
      capture_output=True, text=True, cwd=REPO, env=env, timeout=timeout,
  )
  out = proc.stdout + ('\n' + proc.stderr if proc.stderr.strip() else '')
  return proc.returncode == 0, out.strip(), time.time() - started


def main() -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--replay', action='store_true',
                      help='Verify from cassettes only; never call the API.')
  parser.add_argument('--wait', action='store_true',
                      help='Block until the daily quota resets, then record.')
  parser.add_argument('--skill', action='append', help='Limit to these skills.')
  parser.add_argument('--timeout', type=int, default=3600)
  args = parser.parse_args()

  mode = 'replay' if args.replay else 'auto'
  targets = args.skill or ORDER
  CASSETTE.mkdir(parents=True, exist_ok=True)

  if args.wait and not args.replay:
    waited = 0
    while not quota_available():
      if waited == 0:
        print('Daily credit allowance is exhausted. Waiting for reset '
              '(00:00 UTC); probing every 5 min. A 429 probe is free.',
              flush=True)
      time.sleep(300)
      waited += 300
      if waited % 1800 == 0:
        print(f'  ...still exhausted after {waited // 60} min', flush=True)
    print('Credits available. Starting the recording pass.\n', flush=True)

  print('=' * 78)
  print(f'{"REPLAY (free)" if args.replay else "RECORD (spends credits)"} — '
        f'{len(targets)} skills, serial')
  print(f'cassette: {CASSETTE}')
  print('=' * 78, flush=True)

  results: dict[str, tuple[bool, float]] = {}
  quota_hit = False
  for skill in targets:
    n_before = len(list(CASSETTE.glob('*.json')))
    ok, out, secs = run_eval(skill, mode, args.timeout)
    n_after = len(list(CASSETTE.glob('*.json')))
    results[skill] = (ok, secs)

    tail = '\n'.join(out.splitlines()[-14:])
    print(f'\n{"-" * 78}\n{"PASS" if ok else "FAIL"}  {skill}  '
          f'({secs:.0f}s, +{n_after - n_before} responses recorded)\n{"-" * 78}')
    print(tail, flush=True)

    if 'daily credit allowance exhausted' in out or 'credit limit' in out:
      quota_hit = True
      print('\n>>> Daily credit allowance exhausted. Stopping here.', flush=True)
      print('>>> Everything recorded so far is kept. Re-run this script after '
            '00:00 UTC to continue from where it stopped.', flush=True)
      break

  STATE.write_text(json.dumps(
      {s: {'passed': ok, 'seconds': round(t, 1)} for s, (ok, t) in results.items()},
      indent=2), encoding='utf-8')

  print('\n' + '=' * 78)
  print('SCOREBOARD')
  print('=' * 78)
  for skill in targets:
    if skill in results:
      ok, secs = results[skill]
      print(f'  {"PASS" if ok else "FAIL"}  {skill:36s} {secs:6.0f}s')
    else:
      print(f'  ----  {skill:36s} not reached')
  passed = sum(ok for ok, _ in results.values())
  print('-' * 78)
  print(f'  {passed}/{len(targets)} skills passing   '
        f'({len(list(CASSETTE.glob("*.json")))} responses cached)')
  if quota_hit:
    print('  QUOTA EXHAUSTED — pass is incomplete but resumable.')
  print('=' * 78)
  return 0 if passed == len(targets) else 1


if __name__ == '__main__':
  sys.exit(main())
