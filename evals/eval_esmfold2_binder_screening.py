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

"""Eval for esmfold2_binder_screening: a positive/negative control screen.

Target: BARNASE. Candidates:
  * BARSTAR              -- the true natural inhibitor (Kd ~ 1e-14 M). POSITIVE.
  * scramble(BARSTAR)    -- same composition, fold destroyed. DECOY.
  * LYSOZYME             -- an unrelated, non-binding protein. NEGATIVE.

The skill is only worth anything if it puts barstar first and both decoys well
behind it. That is the whole contract.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

from __future__ import annotations

import csv
import pathlib
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from fixtures import (  # pylint: disable=g-import-not-at-top
    BARNASE,
    BARSTAR,
    LYSOZYME,
    Eval,
    load_json,
    run_script,
    scramble,
    write_fasta,
)

SKILL = 'esmfold2_binder_screening'
SCRIPT = 'screen_binders.py'

# Keep the eval quick. The controls separate cleanly at these settings; they
# also separate at the script's defaults (3 / 100), just more slowly.
NUM_LOOPS = '3'
NUM_SAMPLING_STEPS = '50'


def metrics_by_name(payload) -> dict[str, dict]:
  return {row['candidate']: row['metrics'] for row in payload['results']}


def main() -> int:
  ev = Eval(SKILL)
  tmp = pathlib.Path(tempfile.mkdtemp(prefix='binder_screen_'))

  candidates = write_fasta(
      tmp / 'candidates.fasta',
      [
          ('barstar', BARSTAR),
          ('barstar_scrambled', scramble(BARSTAR, seed=0)),
          ('lysozyme', LYSOZYME),
      ],
  )
  out_dir = tmp / 'screen'

  # ---- screen ----------------------------------------------------------
  run_script(
      SKILL, SCRIPT, 'screen',
      '--target', BARNASE,
      '--candidates', candidates,
      '--output', str(out_dir),
      '--top-n', '1',
      '--num-loops', NUM_LOOPS,
      '--num-sampling-steps', NUM_SAMPLING_STEPS,
  )

  results_json = out_dir / 'results.json'
  results_csv = out_dir / 'results.csv'
  ev.check('screen writes results.json', results_json.is_file())
  ev.check('screen writes results.csv', results_csv.is_file())

  payload = load_json(results_json)
  rows = payload['results']
  ev.check('all 3 candidates screened', len(rows) == 3, f'{len(rows)} rows')

  m = metrics_by_name(payload)
  barstar = m['barstar']['iptm']
  scrambled = m['barstar_scrambled']['iptm']
  lysozyme = m['lysozyme']['iptm']
  print(
      f'\n  iPTM: barstar={barstar:.3f}  scrambled={scrambled:.3f}  '
      f'lysozyme={lysozyme:.3f}\n'
  )

  # ---- SCIENTIFIC ------------------------------------------------------
  ev.check(
      'BARSTAR ranks #1 by iPTM',
      rows[0]['candidate'] == 'barstar',
      f'#1 = {rows[0]["candidate"]}',
  )
  ev.check('BARSTAR iPTM > 0.8', barstar > 0.8, f'{barstar:.3f}')
  ev.check(
      'BARSTAR iPTM > scrambled-barstar iPTM',
      barstar > scrambled,
      f'{barstar:.3f} vs {scrambled:.3f}',
  )
  ev.check(
      'BARSTAR iPTM > lysozyme iPTM',
      barstar > lysozyme,
      f'{barstar:.3f} vs {lysozyme:.3f}',
  )
  ev.check(
      'BARSTAR has more interface contacts than LYSOZYME',
      (
          m['barstar']['confident_interface_contacts']
          > m['lysozyme']['confident_interface_contacts']
      ),
      f'{m["barstar"]["confident_interface_contacts"]} vs '
      f'{m["lysozyme"]["confident_interface_contacts"]} (PAE-gated)',
  )
  ev.check(
      'both decoys land in a non-binding iPTM band',
      max(scrambled, lysozyme) < 0.5,
      f'max decoy iPTM = {max(scrambled, lysozyme):.3f}',
  )

  # ---- STRUCTURAL ------------------------------------------------------
  iptms = [row['metrics']['iptm'] for row in rows]
  ev.check(
      'results sorted descending by the ranking metric',
      iptms == sorted(iptms, reverse=True),
      f'{[round(v, 3) for v in iptms]}',
  )
  ev.check(
      'ranking metric recorded',
      payload['ranking_metric'] == 'iptm',
      payload.get('ranking_metric', ''),
  )
  ev.check(
      'shortlist honours --top-n 1',
      payload['shortlist'] == ['barstar'],
      str(payload['shortlist']),
  )

  with open(results_csv, encoding='utf-8') as handle:
    csv_rows = list(csv.DictReader(handle))
  ev.check('CSV has one row per candidate', len(csv_rows) == 3)
  ev.check(
      'CSV is ranked barstar-first',
      csv_rows[0]['candidate'] == 'barstar' and csv_rows[0]['rank'] == '1',
  )

  pdbs = sorted((out_dir / 'structures').glob('*.pdb'))
  ev.check('one PDB per candidate', len(pdbs) == 3, f'{len(pdbs)} PDBs')
  two_chain = []
  for pdb in pdbs:
    chains = {
        line[21]
        for line in pdb.read_text(encoding='utf-8').splitlines()
        if line.startswith('ATOM')
    }
    two_chain.append(len(chains) == 2)
  ev.check(
      'every PDB has 2 chains',
      all(two_chain),
      f'{sum(two_chain)}/{len(pdbs)}',
  )

  # ---- analyze-interface ----------------------------------------------
  iface_dir = tmp / 'iface'
  run_script(
      SKILL, SCRIPT, 'analyze-interface',
      '--target', BARNASE,
      '--binder', BARSTAR,
      '--output', str(iface_dir),
      '--num-loops', NUM_LOOPS,
      '--num-sampling-steps', NUM_SAMPLING_STEPS,
  )
  iface_json = iface_dir / 'target_binder_interface.json'
  iface_png = iface_dir / 'target_binder_pae.png'
  ev.check('analyze-interface writes JSON', iface_json.is_file())
  ev.check('analyze-interface writes a PAE heatmap PNG', iface_png.is_file())

  if iface_json.is_file():
    iface = load_json(iface_json)
    contacts = iface['metrics']['confident_interface_contacts']
    ev.check(
        'analyze-interface reports > 0 interface contacts',
        contacts > 0,
        f'{contacts} contacts',
    )
    ev.check(
        'interface residues listed on both chains',
        iface['interface_residues_target'] and iface['interface_residues_binder'],
        f'{len(iface["interface_residues_target"])} on A, '
        f'{len(iface["interface_residues_binder"])} on B',
    )

  # ---- rank ------------------------------------------------------------
  ranked = tmp / 'ranked.json'
  run_script(
      SKILL, SCRIPT, 'rank',
      '--results', str(results_json),
      '--output', str(ranked),
      '--min-iptm', '0.8',
      '--top-n', '1',
  )
  kept = [row['candidate'] for row in load_json(ranked)['results']]
  ev.check(
      'rank --min-iptm 0.8 keeps barstar and drops both decoys',
      kept == ['barstar'],
      f'kept {kept}',
  )

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
