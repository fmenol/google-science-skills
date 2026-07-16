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

"""Vendors skills/esm_common/esm_biohub.py into every ESM skill.

Skills must be self-contained (they can be installed individually), so the shared
client is copied rather than imported across directories. This script is the
single source of truth for those copies.

  uv run evals/sync_common.py          # copy canonical -> all skills
  uv run evals/sync_common.py --check  # verify copies are identical (CI/eval)
"""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import shutil
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
CANONICAL = REPO / 'skills' / 'esm_common' / 'esm_biohub.py'

ESM_SKILLS = [
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


def digest(path: pathlib.Path) -> str:
  return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def main() -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
      '--check',
      action='store_true',
      help='Verify the vendored copies match; do not write.',
  )
  args = parser.parse_args()

  if not CANONICAL.is_file():
    print(f'Missing canonical library: {CANONICAL}', file=sys.stderr)
    return 1
  want = digest(CANONICAL)

  drift = []
  for skill in ESM_SKILLS:
    target = REPO / 'skills' / skill / 'scripts' / 'esm_biohub.py'
    if args.check:
      if not target.is_file():
        drift.append(f'{skill}: MISSING')
      elif digest(target) != want:
        drift.append(f'{skill}: DRIFTED ({digest(target)} != {want})')
    else:
      target.parent.mkdir(parents=True, exist_ok=True)
      shutil.copyfile(CANONICAL, target)
      print(f'  vendored -> skills/{skill}/scripts/esm_biohub.py')

  if args.check:
    if drift:
      print('esm_biohub.py is out of sync:', file=sys.stderr)
      for line in drift:
        print(f'  {line}', file=sys.stderr)
      print('\nRun: uv run evals/sync_common.py', file=sys.stderr)
      return 1
    print(f'All {len(ESM_SKILLS)} copies of esm_biohub.py match ({want}).')
    return 0

  print(f'\nSynced {len(ESM_SKILLS)} skills from {CANONICAL.name} ({want}).')
  return 0


if __name__ == '__main__':
  sys.exit(main())
