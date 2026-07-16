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

"""Eval for the esm3_secondary_structure_sasa skill.

Ground truth is the known DSSP topology of two well-characterised proteins:

* Ubiquitin (1UBQ): a β-grasp fold — an N-terminal β-hairpin, ONE prominent
  α-helix (residues ~23-34), and a mixed β-sheet. So its SS8 must contain both
  helix and strand, with the helix falling in that window.
* Myoglobin (a globin): an all-α protein — essentially no β-strand.

The buried-core check (some residues near-zero SASA, others highly exposed) and
the hydrophobic-vs-charged burial contrast test the SASA track against basic
biophysics.
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

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(
    0,
    str(pathlib.Path(__file__).resolve().parents[1] / 'skills' / 'esm_common'),
)

import numpy as np  # noqa: E402

from fixtures import UBIQUITIN, Eval, load_json, run_script, write_fasta  # noqa: E402

SKILL = 'esm3_secondary_structure_sasa'
SCRIPT = 'predict_ss_sasa.py'

# Sperm-whale myoglobin (P02185, mature 153 aa), PDB 1MBN — the classic all-α
# globin fold: eight helices, no β-strand.
MYOGLOBIN = (
    'VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILK'
    'KKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG'
)


def main() -> int:
  ev = Eval(SKILL)
  work = pathlib.Path('evals/_artifacts') / SKILL
  work.mkdir(parents=True, exist_ok=True)

  # ----------------------------------------------------------------------
  # 1. STRUCTURAL: predict writes a well-formed record for ubiquitin.
  # ----------------------------------------------------------------------
  print('--- predict (ubiquitin) ---')
  ubq_json = work / 'ubiquitin.json'
  run_script(
      SKILL, SCRIPT, 'predict',
      '--sequence', UBIQUITIN, '--temperature', '0.1',
      '--output', str(ubq_json),
  )
  ubq = load_json(ubq_json)

  ev.check('predict writes a JSON record', ubq_json.is_file())
  ss8 = ubq['secondary_structure_ss8']
  sasa = np.array(ubq['sasa'], dtype=np.float64)

  ev.check('SS8 length matches the sequence', len(ss8) == len(UBIQUITIN),
           f'{len(ss8)} vs {len(UBIQUITIN)}')
  ev.check('SASA length matches the sequence', len(sasa) == len(UBIQUITIN),
           f'{len(sasa)} vs {len(UBIQUITIN)}')
  ev.check('SS8 uses only the DSSP 8-state alphabet',
           set(ss8) <= set('HGIEBTSC'), f'alphabet {sorted(set(ss8))}')
  ev.check('burial labels align with the sequence',
           len(ubq['burial']) == len(UBIQUITIN))

  # ----------------------------------------------------------------------
  # 2. SCIENTIFIC: ubiquitin's known β-grasp topology.
  # ----------------------------------------------------------------------
  print('\n--- SCIENCE: ubiquitin topology ---')
  ss3 = ubq['secondary_structure_ss3']
  ev.check('ubiquitin has helix (it has one prominent α-helix)',
           'H' in ss3, f'helix content {ss3.count("H")}/{len(ss3)}')
  ev.check('ubiquitin has strand (it is a β-grasp fold)',
           'E' in ss3, f'strand content {ss3.count("E")}/{len(ss3)}')

  # The α-helix runs ~residue 23-34 (0-indexed 22-33). Require the helix to be
  # centred there, not scattered — the single strongest topological fact.
  helix_positions = [i for i, c in enumerate(ss3) if c == 'H']
  in_window = [i for i in helix_positions if 20 <= i <= 36]
  ev.check(
      'the α-helix sits at residues ~23-34',
      len(in_window) >= 6 and len(in_window) >= 0.5 * len(helix_positions),
      f'{len(in_window)} of {len(helix_positions)} helix residues in 21-37',
  )
  # β-strands cluster in the N- and C-terminal halves, not inside the helix.
  strand_positions = [i for i, c in enumerate(ss3) if c == 'E']
  ev.check(
      'β-strands avoid the helix core',
      all(not (23 <= i <= 33) for i in strand_positions)
      or sum(23 <= i <= 33 for i in strand_positions) <= 1,
      f'{sum(23 <= i <= 33 for i in strand_positions)} strand residues in the '
      'helix window (want ~0)',
  )

  # ----------------------------------------------------------------------
  # 3. SCIENTIFIC: SASA describes a real buried core.
  # ----------------------------------------------------------------------
  print('\n--- SCIENCE: solvent accessibility ---')
  ev.check('all SASA values are finite and non-negative',
           np.all(np.isfinite(sasa)) and np.all(sasa >= 0))
  ev.check('the protein has a buried core (some near-zero SASA)',
           float(sasa.min()) < 10.0, f'min SASA {sasa.min():.1f} Å²')
  ev.check('the protein has an exposed surface (some high SASA)',
           float(sasa.max()) > 80.0, f'max SASA {sasa.max():.1f} Å²')
  ev.check('SASA spans a wide range (core to surface)',
           float(sasa.max() - sasa.min()) > 60.0,
           f'range {sasa.max() - sasa.min():.1f} Å²')

  # Hydrophobic residues should be, on average, more buried than charged ones.
  hydrophobic = set('AILMFWVC')
  charged = set('DEKR')
  hyd_sasa = np.array([sasa[i] for i, a in enumerate(UBIQUITIN)
                       if a in hydrophobic])
  chg_sasa = np.array([sasa[i] for i, a in enumerate(UBIQUITIN)
                       if a in charged])
  ev.check(
      'hydrophobic residues are more buried than charged ones',
      float(hyd_sasa.mean()) < float(chg_sasa.mean()),
      f'hydrophobic mean {hyd_sasa.mean():.1f} vs charged {chg_sasa.mean():.1f} '
      'Å²',
  )

  # ----------------------------------------------------------------------
  # 4. SCIENTIFIC / NEGATIVE CONTROL: an all-α protein has ~no strand, and
  #    clearly more helix than the β-grasp ubiquitin.
  # ----------------------------------------------------------------------
  print('\n--- SCIENCE: all-α myoglobin vs β-grasp ubiquitin ---')
  myo_json = work / 'myoglobin.json'
  run_script(
      SKILL, SCRIPT, 'predict',
      '--sequence', MYOGLOBIN, '--temperature', '0.1',
      '--output', str(myo_json),
  )
  myo = load_json(myo_json)
  myo_helix = myo['summary']['percent_helix']
  myo_strand = myo['summary']['percent_strand']
  ubq_helix = ubq['summary']['percent_helix']

  ev.check('myoglobin is helix-rich (all-α globin)',
           myo_helix > 55.0, f'{myo_helix}% helix')
  ev.check('myoglobin has essentially no β-strand',
           myo_strand < 10.0, f'{myo_strand}% strand')
  ev.check(
      'the all-α protein is more helical than the β-grasp protein',
      myo_helix > ubq_helix,
      f'myoglobin {myo_helix}% vs ubiquitin {ubq_helix}% helix',
  )

  # ----------------------------------------------------------------------
  # 5. STRUCTURAL: batch + plot.
  # ----------------------------------------------------------------------
  print('\n--- batch + plot ---')
  fasta = write_fasta(work / 'two.fasta',
                      [('ubiquitin', UBIQUITIN), ('myoglobin', MYOGLOBIN)])
  batch_json = work / 'batch.json'
  batch_csv = work / 'batch.csv'
  run_script(SKILL, SCRIPT, 'batch', '--fasta', fasta,
             '--output', str(batch_json), '--summary-csv', str(batch_csv),
             '--temperature', '0.1')
  batch = load_json(batch_json)
  ev.check('batch returns one record per input', len(batch) == 2,
           f'{len(batch)}')
  ev.check('batch preserves order',
           [r['id'] for r in batch] == ['ubiquitin', 'myoglobin'])
  ev.check('batch writes a summary CSV', batch_csv.is_file())

  png = work / 'ubiquitin.png'
  run_script(SKILL, SCRIPT, 'plot', '--input', str(ubq_json),
             '--output', str(png))
  ev.check('plot writes a non-empty PNG',
           png.is_file() and png.stat().st_size > 1000,
           f'{png.stat().st_size if png.is_file() else 0} bytes')

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
