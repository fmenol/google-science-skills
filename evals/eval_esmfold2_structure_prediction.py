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

"""Eval for the `esmfold2_structure_prediction` skill.

Structural checks: the CLI writes a well-formed PDB (right residues, right
chains, B-factors rescaled to 0-100), a metrics JSON with the documented keys,
HETATM records for a CCD ligand, and both plots plus a report from `analyze`.

Scientific checks (the ones that matter):
  * Ubiquitin folds well:            pLDDT > 0.7 and pTM > 0.6.
  * Its own scramble folds WORSE on BOTH pLDDT and pTM. Same composition, no
    fold. If a shuffled sequence scored as well as the real one, the pipeline
    would be wired to something that is not looking at the sequence.
  * Barnase + barstar, a real 10^-14 M complex, gives iPTM > 0.8.

Measured on the live API at num_loops=10, num_sampling_steps=50:
  ubiquitin           pLDDT 0.823  pTM 0.782
  scrambled ubiquitin pLDDT 0.483  pTM 0.238
  barnase + barstar   pLDDT 0.951  pTM 0.967  iPTM 0.964
"""

# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

from __future__ import annotations

import collections
import pathlib
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from fixtures import BARNASE  # pylint: disable=g-import-not-at-top
from fixtures import BARSTAR
from fixtures import Eval
from fixtures import load_json
from fixtures import run_script
from fixtures import scramble
from fixtures import UBIQUITIN

SKILL = 'esmfold2_structure_prediction'

# Fast but still scientifically valid; the margins above are enormous.
FAST = ['--num-loops', '10', '--num-sampling-steps', '50']

PNG_MAGIC = b'\x89PNG\r\n\x1a\n'


def parse_pdb(path: str) -> dict:
  """Plain-text PDB parse. No biotite, no numpy — keep the eval dependency-free."""
  ca_by_chain = collections.Counter()
  residues_by_chain = collections.defaultdict(set)
  bfactors = []
  n_hetatm = 0
  with open(path, encoding='utf-8') as handle:
    for line in handle:
      if line.startswith('HETATM'):
        n_hetatm += 1
        continue
      if not line.startswith('ATOM'):
        continue
      chain = line[21]
      res_seq = line[22:26].strip()
      atom = line[12:16]
      bfactors.append(float(line[60:66]))
      residues_by_chain[chain].add(res_seq)
      if atom == ' CA ':
        ca_by_chain[chain] += 1
  return {
      'ca_by_chain': ca_by_chain,
      'residues_by_chain': {k: len(v) for k, v in residues_by_chain.items()},
      'bfactors': bfactors,
      'n_hetatm': n_hetatm,
  }


def main() -> int:
  ev = Eval(SKILL)
  tmp = pathlib.Path(tempfile.mkdtemp(prefix='esmfold2_struct_'))
  print(f'[{SKILL}] artifacts in {tmp}\n')

  # ---------------------------------------------------------------------
  # 1. fold ubiquitin (the reference real fold; also feeds `analyze`)
  # ---------------------------------------------------------------------
  ubi_pdb = tmp / 'ubiquitin.pdb'
  ubi_json = tmp / 'ubiquitin.json'
  run_script(
      SKILL, 'fold.py', 'fold',
      '--sequence', UBIQUITIN,
      '--include-pae',
      '--output-pdb', str(ubi_pdb),
      '--output-metrics', str(ubi_json),
      *FAST,
  )

  ev.check('fold writes a PDB', ubi_pdb.is_file(), str(ubi_pdb))
  pdb = parse_pdb(str(ubi_pdb))
  ev.check(
      'ubiquitin PDB has 76 residues',
      pdb['residues_by_chain'].get('A') == 76,
      f'{pdb["residues_by_chain"]}',
  )
  ev.check(
      'ubiquitin PDB has exactly one CA per residue',
      pdb['ca_by_chain'].get('A') == 76,
      f'{pdb["ca_by_chain"]["A"]} CA atoms',
  )
  ev.check(
      'ubiquitin PDB is a single chain A',
      list(pdb['residues_by_chain']) == ['A'],
      f'chains={list(pdb["residues_by_chain"])}',
  )
  bmin, bmax = min(pdb['bfactors']), max(pdb['bfactors'])
  # pLDDT is 0-1 from this API; the B-factor column must be the 0-100 rescale.
  ev.check(
      'B-factors are pLDDT on the 0-100 scale',
      0.0 <= bmin and bmax <= 100.0 and bmax > 1.5,
      f'{bmin:.2f}..{bmax:.2f}',
  )

  metrics = load_json(ubi_json)
  ev.check(
      'metrics JSON has plddt_mean, ptm, sequence_length',
      all(k in metrics for k in ('plddt_mean', 'ptm', 'sequence_length')),
      f'keys={sorted(metrics)[:8]}...',
  )
  ev.check(
      'metrics sequence_length == 76',
      metrics.get('sequence_length') == 76,
      str(metrics.get('sequence_length')),
  )
  ev.check(
      'metrics declares the 0-1 pLDDT scale',
      metrics.get('plddt_scale') == '0-1'
      and 0.0 <= metrics['plddt_mean'] <= 1.0,
      f'plddt_scale={metrics.get("plddt_scale")!r} '
      f'plddt_mean={metrics["plddt_mean"]:.4f}',
  )

  real_plddt = float(metrics['plddt_mean'])
  real_ptm = float(metrics['ptm'])

  # SCIENTIFIC: a real, well-behaved small protein must fold confidently.
  ev.check(
      'SCIENCE ubiquitin folds well (pLDDT > 0.7)',
      real_plddt > 0.7,
      f'plddt_mean={real_plddt:.4f}',
  )
  ev.check(
      'SCIENCE ubiquitin folds well (pTM > 0.6)',
      real_ptm > 0.6,
      f'ptm={real_ptm:.4f}',
  )

  # ---------------------------------------------------------------------
  # 2. analyze -> report + both plots
  # ---------------------------------------------------------------------
  report_json = tmp / 'ubiquitin_report.json'
  plddt_png = tmp / 'ubiquitin_plddt.png'
  pae_png = tmp / 'ubiquitin_pae.png'
  run_script(
      SKILL, 'fold.py', 'analyze',
      '--metrics', str(ubi_json),
      '--output-report', str(report_json),
      '--output-plddt-plot', str(plddt_png),
      '--output-pae-plot', str(pae_png),
  )

  for label, png in (('pLDDT plot', plddt_png), ('PAE heatmap', pae_png)):
    ok = png.is_file() and png.stat().st_size > 1000
    if ok:
      with open(png, 'rb') as handle:
        ok = handle.read(8) == PNG_MAGIC
    ev.check(
        f'analyze writes a non-empty {label} PNG',
        ok,
        f'{png.stat().st_size if png.is_file() else 0} bytes',
    )

  ev.check('analyze writes a JSON report', report_json.is_file())
  report = load_json(report_json)
  ev.check(
      'report interprets the numbers for the agent',
      'confidence_bands' in report
      and 'low_confidence_regions' in report
      and report.get('ptm', {}).get('verdict'),
      f'ptm verdict={report.get("ptm", {}).get("verdict")!r}',
  )
  bands = report['confidence_bands']
  total = sum(b['count'] for b in bands.values())
  ev.check(
      'confidence bands partition every residue',
      total == 76,
      f'{total} residues across {len(bands)} bands',
  )

  # ---------------------------------------------------------------------
  # 3. NEGATIVE CONTROL: scrambled ubiquitin must fold worse
  # ---------------------------------------------------------------------
  scr = scramble(UBIQUITIN, seed=0)
  scr_pdb = tmp / 'scrambled.pdb'
  scr_json = tmp / 'scrambled.json'
  run_script(
      SKILL, 'fold.py', 'fold',
      '--sequence', scr,
      '--output-pdb', str(scr_pdb),
      '--output-metrics', str(scr_json),
      *FAST,
  )
  scr_metrics = load_json(scr_json)
  scr_plddt = float(scr_metrics['plddt_mean'])
  scr_ptm = float(scr_metrics['ptm'])

  ev.check(
      'SCIENCE scramble folds worse than real ubiquitin (pLDDT)',
      scr_plddt < real_plddt,
      f'scrambled {scr_plddt:.4f} < real {real_plddt:.4f}',
  )
  ev.check(
      'SCIENCE scramble folds worse than real ubiquitin (pTM)',
      scr_ptm < real_ptm,
      f'scrambled {scr_ptm:.4f} < real {real_ptm:.4f}',
  )

  # ---------------------------------------------------------------------
  # 4. fold-complex: barnase + barstar, a true high-affinity complex
  # ---------------------------------------------------------------------
  cx_pdb = tmp / 'barnase_barstar.pdb'
  cx_json = tmp / 'barnase_barstar.json'
  run_script(
      SKILL, 'fold.py', 'fold-complex',
      '--protein', f'A:{BARNASE}',
      '--protein', f'B:{BARSTAR}',
      '--include-pae',
      '--output-pdb', str(cx_pdb),
      '--output-metrics', str(cx_json),
      *FAST,
  )
  cx = parse_pdb(str(cx_pdb))
  ev.check(
      'complex PDB has exactly 2 chains',
      len(cx['residues_by_chain']) == 2,
      f'chains={sorted(cx["residues_by_chain"])}',
  )
  ev.check(
      'complex chains are 110 (barnase) and 89 (barstar) residues',
      cx['residues_by_chain'].get('A') == 110
      and cx['residues_by_chain'].get('B') == 89,
      f'{dict(sorted(cx["residues_by_chain"].items()))}',
  )
  ev.check(
      'complex has one CA per residue in both chains',
      cx['ca_by_chain'].get('A') == 110 and cx['ca_by_chain'].get('B') == 89,
      f'{dict(sorted(cx["ca_by_chain"].items()))}',
  )

  cx_metrics = load_json(cx_json)
  iptm = cx_metrics.get('interface_ptm')
  ev.check(
      'complex metrics carry interface_ptm',
      isinstance(iptm, float),
      f'interface_ptm={iptm!r}',
  )
  # SCIENTIFIC: barnase/barstar is one of the tightest complexes known.
  ev.check(
      'SCIENCE barnase+barstar interface is confident (iPTM > 0.8)',
      iptm is not None and iptm > 0.8,
      f'interface_ptm={iptm:.4f}' if iptm is not None else 'None',
  )

  # analyze must also cope with a multi-chain metrics file.
  cx_report = tmp / 'barnase_barstar_report.json'
  run_script(
      SKILL, 'fold.py', 'analyze',
      '--metrics', str(cx_json),
      '--output-report', str(cx_report),
      '--output-plddt-plot', str(tmp / 'cx_plddt.png'),
      '--output-pae-plot', str(tmp / 'cx_pae.png'),
  )
  cx_rep = load_json(cx_report)
  ev.check(
      'analyze gives an iPTM verdict for the complex',
      cx_rep.get('interface_ptm', {}).get('verdict') == 'confident interface',
      f'verdict={cx_rep.get("interface_ptm", {}).get("verdict")!r}',
  )

  # ---------------------------------------------------------------------
  # 5. fold-complex with a CCD ligand -> HETATM records
  #    (Plumbing check for the all-atom/ligand path, not a binding claim:
  #     ubiquitin does not bind SAH, and the low iPTM correctly says so.)
  # ---------------------------------------------------------------------
  lig_pdb = tmp / 'ubiquitin_sah.pdb'
  lig_json = tmp / 'ubiquitin_sah.json'
  run_script(
      SKILL, 'fold.py', 'fold-complex',
      '--protein', f'A:{UBIQUITIN}',
      '--ligand-ccd', 'L:SAH',
      '--output-pdb', str(lig_pdb),
      '--output-metrics', str(lig_json),
      *FAST,
  )
  lig = parse_pdb(str(lig_pdb))
  ev.check(
      'CCD ligand is written as HETATM records',
      lig['n_hetatm'] > 0,
      f'{lig["n_hetatm"]} HETATM lines',
  )
  ev.check(
      'ligand complex keeps the protein chain intact',
      lig['residues_by_chain'].get('A') == 76,
      f'{dict(sorted(lig["residues_by_chain"].items()))}',
  )

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
