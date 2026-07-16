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

"""Eval for the esm3_inverse_folding skill.

Setup folds ubiquitin with ESMFold2 and writes the result as a PDB; that
backbone is the design target. The skill must then recover ubiquitin's own
sequence from its own backbone, and its designs must fold back to that backbone.

Scientific contract:
  * Native-sequence recovery at T=0.1 > 30% for the best of 3 designs.
    (Measured: ~100%. The floor is deliberately far below the observed value.)
  * Self-consistency: the best design re-folds to the target with scTM > 0.7.
  * Negative control (recovery): recovery must beat 4x an explicitly computed
    random amino-acid baseline (~5%).
  * Negative control (scTM): because a perfect design saturates scTM at 1.000,
    that check alone cannot catch a scorer that always returns 1.0. So a decoy
    structure -- the same atoms with the residue correspondence shuffled -- is
    pushed through the skill's own scoring function and must be rejected.
    (Measured: decoy scTM 0.08, scRMSD ~15 A. No API calls.)
  * Diversity: 3 samples at T=1.0 are not all byte-identical.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
#   "numpy",
# ]
# ///

from __future__ import annotations

import csv
import pathlib
import random
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(
    0,
    str(
        pathlib.Path(__file__).resolve().parents[1]
        / 'skills' / 'esm3_inverse_folding' / 'scripts'
    ),
)

import esm_biohub as eb  # noqa: E402  (vendored in the skill under test)
from inverse_fold import self_consistency_scores  # noqa: E402  (skill under test)
from fixtures import (  # noqa: E402
    UBIQUITIN,
    Eval,
    load_json,
    run_script,
    seq_identity,
)

SKILL = 'esm3_inverse_folding'
SCRIPT = 'inverse_fold.py'


def random_recovery_baseline(native: str, trials: int = 500, seed: int = 0) -> float:
  """Mean identity to `native` of a uniformly random amino-acid sequence.

  This is the chance level the designs must beat. It is computed here rather
  than hardcoded so the control moves with the fixture.
  """
  rng = random.Random(seed)
  scores = [
      seq_identity(''.join(rng.choices(eb.AA20, k=len(native))), native)
      for _ in range(trials)
  ]
  return float(np.mean(scores))


def main() -> int:
  ev = Eval(SKILL)
  tmp = tempfile.TemporaryDirectory()
  work = pathlib.Path(tmp.name)

  # ---- Setup: fold ubiquitin, write the backbone we will design against ----
  print('Setup: folding ubiquitin with ESMFold2 to build the target backbone...')
  client = eb.BiohubClient()
  folded = client.fold(UBIQUITIN)
  coords = eb.to_array(folded['coordinates'])
  plddt = eb.to_array(folded['plddt'])
  target_pdb = work / 'ubiquitin_backbone.pdb'
  target_pdb.write_text(eb.atom37_to_pdb(coords, UBIQUITIN, plddt))

  ev.check(
      'setup: target backbone folded and written',
      coords.shape == (len(UBIQUITIN), 37, 3) and target_pdb.is_file(),
      f'coords {coords.shape}, pLDDT {np.nanmean(plddt):.3f}, '
      f'pTM {float(folded["ptm"]):.3f}',
  )

  baseline = random_recovery_baseline(UBIQUITIN)
  print(f'Random amino-acid recovery baseline: {baseline:.4f} '
        f'({baseline:.1%}); designs must beat 4x = {4 * baseline:.1%}\n')

  # ---- 1. design ----------------------------------------------------------
  print('--- design (T=0.1, 3 samples) ---')
  design_json = work / 'design.json'
  run_script(
      SKILL, SCRIPT, 'design',
      '--pdb', str(target_pdb),
      '--num-samples', '3',
      '--temperature', '0.1',
      '--output', str(design_json),
  )
  design_fasta = design_json.with_suffix('.fasta')

  ev.check(
      'design: JSON and FASTA written',
      design_json.is_file() and design_fasta.is_file(),
      f'{design_json.name}, {design_fasta.name}',
  )

  report = load_json(design_json)
  designs = report['designs']
  sequences = [d['sequence'] for d in designs]

  ev.check(
      'design: exactly 3 sequences returned',
      len(sequences) == 3,
      f'got {len(sequences)}',
  )
  ev.check(
      'design: every sequence is 76 residues',
      all(len(s) == len(UBIQUITIN) for s in sequences),
      f'lengths {[len(s) for s in sequences]}',
  )
  non_canonical = sorted(set(''.join(sequences)) - set(eb.AA20))
  ev.check(
      'design: canonical amino acids only',
      not non_canonical,
      f'unexpected characters: {non_canonical}' if non_canonical else '20 AAs',
  )
  fasta_records = eb.read_fasta(str(design_fasta))
  ev.check(
      'design: FASTA holds the same 3 sequences',
      [s for _, s in fasta_records] == sequences,
      f'{len(fasta_records)} records',
  )

  # SCIENTIFIC: native-sequence recovery.
  recoveries = [seq_identity(s, UBIQUITIN) for s in sequences]
  best_recovery = max(recoveries)
  ev.check(
      'SCIENCE design: best-of-3 native recovery > 30% at T=0.1',
      best_recovery > 0.30,
      f'best {best_recovery:.1%}, mean {np.mean(recoveries):.1%}, '
      f'worst {min(recoveries):.1%}',
  )
  # The skill must report recovery itself, not just expose sequences.
  ev.check(
      'design: skill reports native recovery, and it agrees with ours',
      report.get('recovery') is not None
      and abs(report['recovery']['best'] - best_recovery) < 1e-6,
      f'reported {report.get("recovery", {}).get("best")}',
  )

  # SCIENTIFIC: negative control against an explicit random baseline.
  ev.check(
      'SCIENCE control: recovery > 4x the random amino-acid baseline',
      best_recovery > 4 * baseline,
      f'best {best_recovery:.1%} vs 4x baseline {4 * baseline:.1%} '
      f'(baseline {baseline:.1%})',
  )

  # ---- 2. self-consistency ------------------------------------------------
  print('\n--- self-consistency (T=0.1, 3 samples) ---')
  sc_json = work / 'sc.json'
  run_script(
      SKILL, SCRIPT, 'self-consistency',
      '--pdb', str(target_pdb),
      '--num-samples', '3',
      '--temperature', '0.1',
      '--output', str(sc_json),
  )
  sc_csv = sc_json.with_suffix('.csv')

  ev.check(
      'self-consistency: ranked JSON and CSV written',
      sc_json.is_file() and sc_csv.is_file(),
      f'{sc_json.name}, {sc_csv.name}',
  )

  sc = load_json(sc_json)
  sc_designs = sc['designs']
  required = {'sctm', 'scrmsd', 'plddt'}
  ev.check(
      'self-consistency: scTM, scRMSD and pLDDT present for every design',
      len(sc_designs) == 3
      and all(required <= set(d) for d in sc_designs)
      and all(
          isinstance(d[k], (int, float)) and np.isfinite(d[k])
          for d in sc_designs for k in required
      ),
      f'{len(sc_designs)} designs, keys {sorted(required & set(sc_designs[0]))}',
  )
  ev.check(
      'self-consistency: designs are ranked best-first by scTM',
      [d['rank'] for d in sc_designs] == [1, 2, 3]
      and all(
          sc_designs[i]['sctm'] >= sc_designs[i + 1]['sctm']
          for i in range(len(sc_designs) - 1)
      ),
      f'scTM order {[round(d["sctm"], 3) for d in sc_designs]}',
  )

  with open(sc_csv, encoding='utf-8') as handle:
    csv_rows = list(csv.DictReader(handle))
  ev.check(
      'self-consistency: CSV carries the ranked scores',
      len(csv_rows) == 3
      and {'rank', 'sctm', 'scrmsd', 'plddt', 'sequence'} <= set(csv_rows[0]),
      f'{len(csv_rows)} rows, columns {list(csv_rows[0])}',
  )

  # SCIENTIFIC: the best design must actually fold back to the target.
  best_sctm = max(d['sctm'] for d in sc_designs)
  best = min(sc_designs, key=lambda d: d['rank'])
  ev.check(
      'SCIENCE self-consistency: best design re-folds to the target, scTM > 0.7',
      best_sctm > 0.70,
      f'best scTM {best_sctm:.3f}, scRMSD {best["scrmsd"]:.2f} A, '
      f'pLDDT {best["plddt"]:.3f}',
  )

  # SCIENTIFIC: prove the scorer discriminates rather than saturating at 1.0.
  # At T=0.1 the design *is* ubiquitin, so it re-folds onto the target exactly
  # and scTM legitimately hits 1.000 -- but that alone cannot tell a working
  # scorer apart from one that always returns 1.0 (e.g. a bug comparing the
  # target against itself). So push a deliberately wrong structure -- the very
  # same atoms with the residue correspondence shuffled -- through the exact
  # function the skill uses. Costs no API calls.
  # Measured: decoy scTM 0.08-0.09 (scRMSD ~15 A) across seeds.
  perm = np.random.default_rng(0).permutation(coords.shape[0])
  decoy_sctm, decoy_rmsd = self_consistency_scores(coords[perm], coords)
  ev.check(
      'SCIENCE control: scTM rejects a decoy structure (scorer discriminates)',
      decoy_sctm < 0.5 and (best_sctm - decoy_sctm) > 0.4,
      f'decoy scTM {decoy_sctm:.3f} (scRMSD {decoy_rmsd:.1f} A) vs best design '
      f'scTM {best_sctm:.3f}',
  )

  # ---- 3. recovery --------------------------------------------------------
  print('\n--- recovery (T=0.1, 3 samples) ---')
  rec_json = work / 'recovery.json'
  run_script(
      SKILL, SCRIPT, 'recovery',
      '--pdb', str(target_pdb),
      '--num-samples', '3',
      '--temperature', '0.1',
      '--output', str(rec_json),
  )
  rec_png = rec_json.with_suffix('.png')

  ev.check(
      'recovery: report and per-position plot written',
      rec_json.is_file() and rec_png.is_file() and rec_png.stat().st_size > 5000,
      f'{rec_png.name} ({rec_png.stat().st_size // 1024} KB)'
      if rec_png.is_file() else 'plot missing',
  )
  rec = load_json(rec_json)
  ev.check(
      'recovery: per-position agreement covers all 76 residues',
      len(rec['per_position']) == len(UBIQUITIN)
      and all(
          p['native'] == UBIQUITIN[p['position'] - 1] for p in rec['per_position']
      )
      and 0.0 <= rec['overall_recovery'] <= 1.0,
      f'overall {rec["overall_recovery"]:.1%}, best design '
      f'{rec["best_design_recovery"]:.1%}',
  )

  # ---- 4. diversity -------------------------------------------------------
  print('\n--- design (T=1.0, 3 samples) — diversity control ---')
  div_json = work / 'diverse.json'
  run_script(
      SKILL, SCRIPT, 'design',
      '--pdb', str(target_pdb),
      '--num-samples', '3',
      '--temperature', '1.0',
      '--output', str(div_json),
  )
  div = load_json(div_json)
  div_sequences = [d['sequence'] for d in div['designs']]
  unique = len(set(div_sequences))

  # SCIENTIFIC: sampling must actually sample. All-identical means the
  # temperature/resampling path is being ignored.
  ev.check(
      'SCIENCE diversity: 3 samples at T=1.0 are not all identical',
      unique > 1,
      f'{unique}/3 unique, mean pairwise identity '
      f'{div["diversity"]["mean_pairwise_identity"]:.1%}',
  )

  tmp.cleanup()
  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
