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

"""Runs every ESM skill eval and prints a scoreboard.

  uv run --no-project evals/run_evals.py                # all skills, serially
  uv run --no-project evals/run_evals.py --skill esmc_protein_embeddings

Each eval lives at evals/eval_<skill>.py, exits 0 on success, and must exercise
the skill's real CLI end to end. Structural checks (files written, schemas
honoured) and scientific checks (known-good biology, with negative controls) are
both required.

Runs SERIALLY by default. The Biohub account has a single 100-credit/day pool,
so running the evals concurrently against the live API makes them race for it —
a burst exhausts the quota and skills then fail for reasons that have nothing to
do with the skill under test.

To re-verify for free, record once and replay:
  uv run --no-project evals/record.py            # record (spends credits once)
  BIOHUB_CASSETTE=evals/cassettes BIOHUB_CASSETTE_MODE=replay \\
    uv run --no-project evals/run_evals.py --parallel
"""

# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

from __future__ import annotations

import argparse
import concurrent.futures
import os
import pathlib
import subprocess
import sys
import time

EVALS = pathlib.Path(__file__).resolve().parent
REPO = EVALS.parent

SKILLS = [
    'esmc_protein_embeddings',
    'esmc_embedding_layer_sweep',
    'esmc_mutation_effect_scoring',
    'esmc_sae_feature_interpretation',
    'esmfold2_structure_prediction',
    'esmfold2_binder_screening',
    'esm3_protein_design',
    'esm3_inverse_folding',
    'esm3_function_prediction',
    'esm3_guided_generation',
    'esm_protein_tracks',
    'esm3_secondary_structure_sasa',
]


def run_one(skill: str, timeout: int) -> tuple[str, bool, float, str]:
  path = EVALS / f'eval_{skill}.py'
  started = time.time()
  if not path.is_file():
    return skill, False, 0.0, f'MISSING eval script: evals/eval_{skill}.py'
  try:
    proc = subprocess.run(
        ['uv', 'run', '--no-project', '--quiet', str(path)],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=REPO,
    )
  except subprocess.TimeoutExpired:
    return skill, False, time.time() - started, f'TIMEOUT after {timeout}s'
  elapsed = time.time() - started
  output = proc.stdout + ('\n' + proc.stderr if proc.stderr.strip() else '')
  return skill, proc.returncode == 0, elapsed, output.strip()


def main() -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--skill', action='append', help='Run only these skills.')
  parser.add_argument(
      '--parallel',
      action='store_true',
      help='Run skills concurrently. NOT recommended against the live API: the '
      'daily credit allowance is a single shared pool, so a parallel burst '
      'exhausts it and skills fail for reasons unrelated to the skill under '
      'test. Safe only when replaying from cassettes.',
  )
  parser.add_argument('--timeout', type=int, default=2400, help='Per-skill seconds.')
  parser.add_argument('--quiet', action='store_true', help='Only show failures.')
  args = parser.parse_args()
  # Serial by default. Parallelism here is what exhausted the quota during
  # development: eleven evals racing for one 100-credit/day pool.
  args.serial = not args.parallel

  replaying = os.environ.get('BIOHUB_CASSETTE_MODE') == 'replay'
  if args.parallel and not replaying:
    print('WARNING: --parallel against the live API will race for the shared '
          'daily credit pool. Prefer `evals/record.py`, or set '
          'BIOHUB_CASSETTE_MODE=replay.\n', file=sys.stderr)
  if not replaying and not os.environ.get('BIOHUB_CASSETTE'):
    print('NOTE: no BIOHUB_CASSETTE set — this run will spend API credits '
          '(allowance: 100/day). Use `evals/record.py` to record once and '
          'replay for free.\n', file=sys.stderr)

  targets = args.skill or SKILLS
  unknown = [s for s in targets if s not in SKILLS]
  if unknown:
    print(f'Unknown skill(s): {unknown}\nKnown: {SKILLS}', file=sys.stderr)
    return 2

  # The shared client must be identical everywhere before anything else runs.
  print('=' * 78)
  print('Checking vendored esm_biohub.py is in sync...')
  sync = subprocess.run(
      ['uv', 'run', '--no-project', '--quiet',
       str(EVALS / 'sync_common.py'), '--check'],
      capture_output=True, text=True, cwd=REPO,
  )
  print(sync.stdout.strip() or sync.stderr.strip())
  if sync.returncode != 0:
    return 1

  print('=' * 78)
  print(f'Running {len(targets)} skill eval(s)'
        f'{" serially" if args.serial else " in parallel"}...\n')

  results: list[tuple[str, bool, float, str]] = []
  if args.serial:
    for skill in targets:
      results.append(run_one(skill, args.timeout))
      _report(results[-1], args.quiet)
  else:
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
      futures = {pool.submit(run_one, s, args.timeout): s for s in targets}
      for future in concurrent.futures.as_completed(futures):
        results.append(future.result())
        _report(results[-1], args.quiet)

  order = {s: i for i, s in enumerate(targets)}
  results.sort(key=lambda r: order[r[0]])

  print('\n' + '=' * 78)
  print('SCOREBOARD')
  print('=' * 78)
  passed = 0
  for skill, ok, elapsed, _ in results:
    print(f'  {"PASS" if ok else "FAIL"}  {skill:36s} {elapsed:6.0f}s')
    passed += ok
  print('-' * 78)
  print(f'  {passed}/{len(results)} skills passing')
  print('=' * 78)
  return 0 if passed == len(results) else 1


def _report(result: tuple[str, bool, float, str], quiet: bool) -> None:
  skill, ok, elapsed, output = result
  if ok and quiet:
    print(f'PASS {skill} ({elapsed:.0f}s)', flush=True)
    return
  print('\n' + '-' * 78)
  print(f'{"PASS" if ok else "FAIL"}  {skill}  ({elapsed:.0f}s)')
  print('-' * 78)
  print(output, flush=True)


if __name__ == '__main__':
  sys.exit(main())
