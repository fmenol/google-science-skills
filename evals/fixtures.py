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

"""Shared fixtures and assertions for ESM skill evals.

Ground truth here is either a literal biological fact (the ubiquitin sequence)
or a value measured against the live API and reproduced repeatedly. Thresholds
are set well inside the observed margin so they are not flaky.
"""

from __future__ import annotations

import json
import pathlib
import random
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
SKILLS = REPO / 'skills'

# --------------------------------------------------------------------------
# Sequences
# --------------------------------------------------------------------------

# Human ubiquitin (P0CG48, 76 aa). Extremely well characterised; the C-terminal
# LRGG (73-76) is the conjugation motif and is absolutely conserved.
UBIQUITIN = (
    'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
)

# Hen egg-white lysozyme, mature chain (P00698, 129 aa). PDB 1LYZ / 6LYZ.
LYSOZYME = (
    'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLC'
    'NIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL'
)

# Barnase (P00648 mature, 110 aa) and barstar (P11540, 89 aa) — the canonical
# high-affinity protein-protein complex (Kd ~ 10^-14 M). A true-positive binder.
BARNASE = (
    'AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADI'
    'NYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIR'
)
BARSTAR = (
    'KKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFRE'
    'AKAEGCDITIILS'
)

# Human carbonic anhydrase II (P00918, 260 aa), PDB 2CBA.
CA2 = (
    'MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDK'
    'AVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSA'
    'KPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKL'
    'NFNGEGEPEELMVDNWRPAQPLKNRQIKASFK'
)

# avGFP (P42212, 238 aa).
GFP = (
    'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPD'
    'HMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYI'
    'MADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAG'
    'ITHGMDELYK'
)


def scramble(sequence: str, seed: int = 0) -> str:
  """Shuffles a sequence, preserving composition but destroying the fold.

  This is the negative control used throughout: a real protein must score better
  than its own scramble on every meaningful metric.
  """
  chars = list(sequence)
  random.Random(seed).shuffle(chars)
  return ''.join(chars)


# --------------------------------------------------------------------------
# Assertions
# --------------------------------------------------------------------------


class EvalFailure(AssertionError):
  """Raised when a skill's behaviour does not meet its scientific contract."""


class Eval:
  """Collects PASS/FAIL results for one skill and reports at the end."""

  def __init__(self, skill: str):
    self.skill = skill
    self.results: list[tuple[bool, str, str]] = []
    self.started = time.time()

  def check(self, name: str, condition: bool, detail: str = '') -> bool:
    self.results.append((bool(condition), name, detail))
    status = 'PASS' if condition else 'FAIL'
    print(f'  {status}  {name}' + (f' — {detail}' if detail else ''), flush=True)
    return bool(condition)

  def check_close(self, name, value, low, high, detail='') -> bool:
    ok = low <= value <= high
    return self.check(
        name, ok, detail or f'{value:.4g} (want {low:g}..{high:g})'
    )

  def finish(self) -> int:
    failed = [n for ok, n, _ in self.results if not ok]
    elapsed = time.time() - self.started
    total = len(self.results)
    print(
        f'\n[{self.skill}] {total - len(failed)}/{total} passed '
        f'in {elapsed:.0f}s'
    )
    if failed:
      print(f'[{self.skill}] FAILED: {", ".join(failed)}')
      return 1
    print(f'[{self.skill}] ALL PASS')
    return 0


def run_script(
    skill_dir: str,
    script: str,
    *args: str,
    timeout: int = 900,
    cassette_ns: str | None = None,
) -> subprocess.CompletedProcess:
  """Runs `uv run skills/<skill_dir>/scripts/<script> <args>`.

  Args:
    cassette_ns: Optional cassette namespace for this invocation. Pass a value
      unique to the invocation (e.g. its objective + seed) when the script is a
      stochastic multi-step routine whose replay must not cross-talk with
      another invocation that issues identical payloads. Leave unset otherwise,
      so evals benefit from cross-invocation cache reuse.

  Raises:
    EvalFailure: If the script exits non-zero, with stdout+stderr attached.
  """
  path = SKILLS / skill_dir / 'scripts' / script
  if not path.is_file():
    raise EvalFailure(f'Script not found: {path}')
  # --no-project is required: without it uv walks up the tree, finds an
  # unrelated pyproject.toml, and tries to build that project instead.
  cmd = ['uv', 'run', '--no-project', '--quiet', str(path), *args]
  env = None
  if cassette_ns is not None:
    import os
    env = {**os.environ, 'BIOHUB_CASSETTE_NS': cassette_ns}
  proc = subprocess.run(
      cmd, capture_output=True, text=True, timeout=timeout, cwd=REPO, env=env
  )
  if proc.returncode != 0:
    raise EvalFailure(
        f'`{" ".join(cmd[2:])}` exited {proc.returncode}\n'
        f'--- stdout ---\n{proc.stdout[-2500:]}\n'
        f'--- stderr ---\n{proc.stderr[-2500:]}'
    )
  return proc


def load_json(path: str | pathlib.Path):
  with open(path, encoding='utf-8') as handle:
    return json.load(handle)


def write_fasta(path: str | pathlib.Path, records: list[tuple[str, str]]) -> str:
  path = pathlib.Path(path)
  path.parent.mkdir(parents=True, exist_ok=True)
  with open(path, 'w', encoding='utf-8') as handle:
    for name, seq in records:
      handle.write(f'>{name}\n')
      for i in range(0, len(seq), 60):
        handle.write(seq[i : i + 60] + '\n')
  return str(path)


def seq_identity(a: str, b: str) -> float:
  """Fraction of identical positions over the shorter sequence."""
  if not a or not b:
    return 0.0
  n = min(len(a), len(b))
  return sum(x == y for x, y in zip(a[:n], b[:n])) / n


if __name__ == '__main__':
  print('Fixtures:')
  for name, seq in [
      ('UBIQUITIN', UBIQUITIN), ('LYSOZYME', LYSOZYME), ('BARNASE', BARNASE),
      ('BARSTAR', BARSTAR), ('CA2', CA2), ('GFP', GFP),
  ]:
    print(f'  {name:10s} len={len(seq)}')
  sys.exit(0)
