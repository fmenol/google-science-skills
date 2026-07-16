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

"""Eval for the esm3_protein_design skill.

Structural checks: every subcommand runs, honours its length contract, and
writes the FASTA + metrics JSON it promises.

Scientific checks, which are the point of the skill:
  1. ESM3 designs are more foldable than random sequence. Two de novo designs
     are compared against two uniformly-random amino-acid strings of the same
     length. The designs must win on mean pTM. If they do not, ESM3 is adding
     nothing over a random string generator.
  2. A grafted structural motif survives scaffolding: after folding the design,
     the motif's CA trace must superpose on the original crystal motif to
     better than 3.0 A.
  3. An SS8 prompt is actually obeyed: the design's secondary structure, read
     back by ESM3, must match what was asked for.
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
import random
import sys
import tempfile
import urllib.request

import numpy as np

EVALS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(EVALS))
sys.path.insert(0, str(EVALS.parent / 'skills' / 'esm3_protein_design' / 'scripts'))

import esm_biohub as eb  # pylint: disable=wrong-import-position
from fixtures import Eval, EvalFailure, load_json, run_script  # pylint: disable=wrong-import-position

SKILL = 'esm3_protein_design'
SCRIPT = 'design.py'

# Renal dipeptidase, the motif-scaffolding case from esm3_generate.ipynb.
MOTIF_PDB_ID = '1ITU'
MOTIF_CHAIN = 'A'
MOTIF_RANGE = '124-146'  # author numbering == tutorial's positional 123..145
SCAFFOLD_LENGTH = 200
MOTIF_START = 73  # 1-indexed == tutorial's 0-indexed 72

DESIGN_LENGTH = 60
# Helix - loop - helix, 44 residues.
SS8 = 'CCHHHHHHHHHHHHHHHHCCCCCCHHHHHHHHHHHHHHHHHHCC'

MOTIF_RMSD_MAX = 3.0
SS8_AGREEMENT_MIN = 0.60

CA = eb.ATOM37.index('CA')


def fetch_motif_pdb(work: pathlib.Path) -> pathlib.Path:
  """Downloads the reference crystal structure used for motif scaffolding."""
  path = work / f'{MOTIF_PDB_ID}.pdb'
  if path.is_file():
    return path
  url = f'https://files.rcsb.org/download/{MOTIF_PDB_ID}.pdb'
  try:
    with urllib.request.urlopen(url, timeout=120) as response:
      path.write_bytes(response.read())
  except OSError as exc:
    raise EvalFailure(
        f'Could not download {url} ({exc}). The motif-scaffolding check needs '
        'a real crystal structure to graft from.'
    ) from exc
  return path


def random_sequences(n: int, length: int, seed: int = 0) -> list[str]:
  """Uniformly-random amino-acid strings: the negative control."""
  rng = random.Random(seed)
  return [
      ''.join(rng.choice(eb.AA20) for _ in range(length)) for _ in range(n)
  ]


def main() -> int:
  ev = Eval(SKILL)
  work = pathlib.Path(tempfile.mkdtemp(prefix='esm3_design_eval_'))
  client = eb.BiohubClient()
  print(f'[{SKILL}] workdir: {work}\n')

  # ----------------------------------------------------------------------
  # 1. generate: de novo design
  # ----------------------------------------------------------------------
  print('--- generate (de novo) ---')
  prefix = work / 'denovo'
  # ESM3's server-side sampling is not seedable, so a two-sample mean is a noisy
  # estimator. Random 60-mers fold to a tight, low pTM (~0.15); designs average
  # ~0.32 but with real spread, so an unlucky pair of low designs can dip below
  # the random mean. Five samples make the mean comparison reliable, and the
  # skill folds each design itself, so the extra samples cost little.
  n_samples = 5
  run_script(
      SKILL, SCRIPT, 'generate',
      '--length', str(DESIGN_LENGTH),
      '--num-samples', str(n_samples),
      '--output-prefix', str(prefix),
  )

  fasta = prefix.with_suffix('.fasta')
  meta = prefix.with_suffix('.json')
  ev.check('generate writes FASTA', fasta.is_file(), str(fasta))
  ev.check('generate writes metrics JSON', meta.is_file(), str(meta))

  records = eb.read_fasta(str(fasta))
  ev.check(f'generate returns {n_samples} designs',
           len(records) == n_samples, f'{len(records)}')

  designs = [seq for _, seq in records]
  ev.check(
      'each design is exactly 60 residues',
      all(len(s) == DESIGN_LENGTH for s in designs),
      f'{[len(s) for s in designs]}',
  )
  noncanonical = sorted(set(''.join(designs)) - set(eb.AA20))
  ev.check(
      'designs use only the 20 canonical amino acids',
      not noncanonical,
      f'stray characters: {noncanonical}' if noncanonical else 'clean',
  )

  data = load_json(meta)
  ev.check(
      'metrics JSON carries pTM + pLDDT per design',
      all('ptm' in d and 'plddt_mean' in d for d in data['designs']),
  )
  ev.check(
      'designs are ranked by pTM',
      all(
          data['designs'][i]['ptm'] >= data['designs'][i + 1]['ptm']
          for i in range(len(data['designs']) - 1)
      ),
  )
  ev.check(
      'skill used the only reachable ESM3 model',
      data['esm3_model'] == 'esm3-open-2024-03',
      data['esm3_model'],
  )

  # ----------------------------------------------------------------------
  # 2. SCIENTIFIC: designs fold, random strings do not
  # ----------------------------------------------------------------------
  print('\n--- SCIENCE: designed vs random foldability ---')
  design_ptm = [float(d['ptm']) for d in data['designs']]
  design_plddt = [float(d['plddt_mean']) for d in data['designs']]

  controls = random_sequences(n_samples, DESIGN_LENGTH, seed=0)
  rand_ptm, rand_plddt = [], []
  for seq in controls:
    folded = client.fold(seq, eb.DEFAULT_ESMFOLD2)
    rand_ptm.append(float(folded['ptm']))
    rand_plddt.append(float(np.nanmean(eb.to_array(folded['plddt']))))

  mean_design_ptm = float(np.mean(design_ptm))
  mean_rand_ptm = float(np.mean(rand_ptm))
  mean_design_plddt = float(np.mean(design_plddt))
  mean_rand_plddt = float(np.mean(rand_plddt))

  print(f'  designed pTM   : {[round(x, 3) for x in design_ptm]} '
        f'-> mean {mean_design_ptm:.3f}')
  print(f'  random   pTM   : {[round(x, 3) for x in rand_ptm]} '
        f'-> mean {mean_rand_ptm:.3f}')
  print(f'  designed pLDDT : {[round(x, 3) for x in design_plddt]} '
        f'-> mean {mean_design_plddt:.3f}')
  print(f'  random   pLDDT : {[round(x, 3) for x in rand_plddt]} '
        f'-> mean {mean_rand_plddt:.3f}')

  ev.check(
      'designed pTM > random pTM',
      mean_design_ptm > mean_rand_ptm,
      f'{mean_design_ptm:.3f} vs {mean_rand_ptm:.3f}',
  )
  ev.check(
      'designed pLDDT > random pLDDT',
      mean_design_plddt > mean_rand_plddt,
      f'{mean_design_plddt:.3f} vs {mean_rand_plddt:.3f}',
  )

  # ----------------------------------------------------------------------
  # 3. design-with-ss
  # ----------------------------------------------------------------------
  print('\n--- design-with-ss ---')
  prefix = work / 'ss8'
  run_script(
      SKILL, SCRIPT, 'design-with-ss',
      '--ss8', SS8,
      '--num-samples', '1',
      '--output-prefix', str(prefix),
  )
  ss_records = eb.read_fasta(str(prefix.with_suffix('.fasta')))
  ss_design = ss_records[0][1]
  ev.check(
      'SS8-conditioned design has length == len(SS8)',
      len(ss_design) == len(SS8),
      f'{len(ss_design)} vs {len(SS8)}',
  )

  # Does ESM3 actually honour the SS8 prompt? Read the design's structure back.
  readback = client.generate(
      'secondary_structure', eb.DEFAULT_ESM3, sequence=ss_design,
      num_steps=1, temperature=0.0,
  )['outputs']['secondary_structure']
  agreement = sum(a == b for a, b in zip(SS8, readback)) / len(SS8)
  print(f'  requested: {SS8}')
  print(f'  readback : {readback}')
  ev.check(
      'the SS8 prompt is obeyed',
      agreement >= SS8_AGREEMENT_MIN,
      f'{agreement:.0%} agreement (need >= {SS8_AGREEMENT_MIN:.0%})',
  )

  # ----------------------------------------------------------------------
  # 4. scaffold-motif  (+ SCIENTIFIC: the motif survives)
  # ----------------------------------------------------------------------
  print('\n--- scaffold-motif ---')
  motif_pdb = fetch_motif_pdb(work)
  prefix = work / 'scaffold'
  run_script(
      SKILL, SCRIPT, 'scaffold-motif',
      '--motif-pdb', str(motif_pdb),
      '--chain', MOTIF_CHAIN,
      '--motif-range', MOTIF_RANGE,
      '--scaffold-length', str(SCAFFOLD_LENGTH),
      '--motif-start', str(MOTIF_START),
      '--num-samples', '2',
      '--output-prefix', str(prefix),
  )
  scaf = load_json(prefix.with_suffix('.json'))
  motif_seq = scaf['motif_sequence']
  offset = MOTIF_START - 1

  scaf_designs = [d['sequence'] for d in scaf['designs']]
  ev.check(
      'scaffold has the requested length',
      all(len(s) == SCAFFOLD_LENGTH for s in scaf_designs),
      f'{[len(s) for s in scaf_designs]}',
  )
  grafted = [s[offset : offset + len(motif_seq)] for s in scaf_designs]
  ev.check(
      'motif appears verbatim at the requested offset',
      all(g == motif_seq for g in grafted),
      f'want {motif_seq}, got {grafted}',
  )

  # The design a user would actually report: rank 1 by pTM.
  best = scaf['designs'][0]
  print(f'  best scaffold: pTM={best["ptm"]:.3f} '
        f'pLDDT={best["plddt_mean"]:.3f} motif CA RMSD='
        f'{best["motif_rmsd_ca"]:.2f} A')
  print(f'  all motif RMSDs: '
        f'{[round(d["motif_rmsd_ca"], 2) for d in scaf["designs"]]} A')

  # Independently recompute the RMSD from the PDB the skill wrote, so this is a
  # check on the structure, not on the skill's own arithmetic.
  ref_seq, ref_coords, res_ids = eb.parse_pdb_atom37(str(motif_pdb), MOTIF_CHAIN)
  lo, hi = (int(x) for x in MOTIF_RANGE.split('-'))
  picked = [i for i, r in enumerate(res_ids) if lo <= r <= hi]
  ref_ca = ref_coords[picked][:, CA]
  ev.check(
      'motif sequence matches the crystal structure',
      ''.join(ref_seq[i] for i in picked) == motif_seq,
      motif_seq,
  )

  _, design_coords, _ = eb.parse_pdb_atom37(best['pdb_file'])
  design_ca = design_coords[offset : offset + len(motif_seq), CA]
  rmsd, _ = eb.kabsch_rmsd(design_ca, ref_ca)
  print(f'  independently recomputed motif CA RMSD: {rmsd:.2f} A')
  ev.check(
      'motif is structurally preserved (CA RMSD < 3.0 A)',
      rmsd < MOTIF_RMSD_MAX,
      f'{rmsd:.2f} A (need < {MOTIF_RMSD_MAX} A)',
  )

  # ----------------------------------------------------------------------
  # 5. chain-of-thought
  # ----------------------------------------------------------------------
  print('\n--- chain-of-thought ---')
  prefix = work / 'cot'
  run_script(
      SKILL, SCRIPT, 'chain-of-thought',
      '--length', str(DESIGN_LENGTH),
      '--num-samples', '1',
      '--output-prefix', str(prefix),
  )
  cot = load_json(prefix.with_suffix('.json'))
  cot_design = cot['designs'][0]
  ev.check(
      'chain-of-thought produces a sequence of the requested length',
      len(cot_design['sequence']) == DESIGN_LENGTH,
      f'{len(cot_design["sequence"])}',
  )
  ev.check(
      'chain-of-thought decoded all three tracks',
      cot['tracks'] == ['secondary_structure', 'structure', 'sequence'],
      str(cot['tracks']),
  )
  ev.check(
      'chain-of-thought design is canonical',
      not set(cot_design['sequence']) - set(eb.AA20),
  )
  print(f'  cot design: pTM={cot_design["ptm"]:.3f} '
        f'pLDDT={cot_design["plddt_mean"]:.3f}')
  print(f'  {cot_design["sequence"]}')

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
