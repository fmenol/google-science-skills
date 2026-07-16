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

"""Eval for the esm3-guided-generation skill.

Budget: ~81 Biohub credits (1 credit == 1 API call). Free-tier keys are capped
at 100 credits PER DAY, so this eval can be run about once a day. Settings were
chosen with an offline Monte-Carlo of the decoding loop (a deliberately
pessimistic i.i.d. mock of ESM3) so that every scientific assertion has a
simulated pass probability >= 0.95; the real model has context effects that all
push in the favourable direction, so the true margin is larger.

  generate  no-cysteine  L=48, 4 steps x 3 samples ......... 22 credits
  compare   hydrophobicity L=48, 3 steps x 4 samples x 2 seeds  48 credits
  generate  ptm          L=48, 2 steps x 2 samples ......... 11 credits
"""

# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

from __future__ import annotations

import pathlib
import statistics
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from fixtures import Eval, load_json, run_script  # pylint: disable=g-import-not-at-top

SKILL = 'esm3_guided_generation'
SCRIPT = 'guided_design.py'
AA20 = set('ACDEFGHIKLMNPQRSTVWY')
LENGTH = 48


def main() -> int:
  ev = Eval(SKILL)
  out = pathlib.Path(tempfile.mkdtemp(prefix='eval_esm3_guided_'))
  print(f'  artifacts -> {out}\n')

  # ----------------------------------------------------------------------
  # STRUCTURAL: the objective registry is discoverable and annotates cost.
  # (0 credits)
  # ----------------------------------------------------------------------
  proc = run_script(SKILL, SCRIPT, 'list-objectives')
  listing = proc.stdout
  missing = [
      name for name in
      ('none', 'no-cysteine', 'hydrophobicity', 'isoelectric-point', 'ptm',
       'radius-of-gyration')
      if name not in listing
  ]
  ev.check('list-objectives prints the registry', not missing,
           f'missing: {missing}' if missing else f'{len(listing)} chars')
  ev.check(
      'list-objectives annotates API cost (cheap vs EXPENSIVE)',
      'EXPENSIVE' in listing and 'cheap' in listing,
  )

  # ----------------------------------------------------------------------
  # STRUCTURAL + SCIENTIFIC: no-cysteine.  (22 credits)
  # The killer deterministic test: guidance must drive the cysteine count to
  # exactly zero. Unguided ESM3 puts a mean of 1.67 cysteines in a 48-mer
  # (measured on this API; 83% of unguided 48-mers carry at least one), so this
  # is a real target, not a freebie.
  # ----------------------------------------------------------------------
  nocys_json = out / 'nocys.json'
  run_script(
      SKILL, SCRIPT, 'generate',
      '--length', str(LENGTH), '--objective', 'no-cysteine',
      '--num-decoding-steps', '4', '--num-samples-per-step', '3',
      '--seed', '0', '--output', str(nocys_json),
      cassette_ns='nocys',
  )
  report = load_json(nocys_json)
  sequence = report['sequence']

  ev.check('generate writes a JSON report', nocys_json.is_file())
  ev.check('generate writes a trajectory plot',
           (out / 'nocys.png').is_file())
  ev.check('sequence is exactly the requested length',
           len(sequence) == LENGTH, f'len={len(sequence)} want={LENGTH}')
  bad = sorted(set(sequence) - AA20)
  ev.check('sequence is all-canonical amino acids', not bad,
           f'non-canonical: {bad}' if bad else 'all 20-letter')
  ev.check(
      'report carries a per-step objective trajectory',
      len(report.get('objective_trajectory', [])) == 4,
      f'{report.get("objective_trajectory")}',
  )
  ptm = report['final_metrics'].get('ptm')
  ev.check(
      'generate folds the final design and reports pTM/pLDDT',
      ptm is not None and 0.0 <= ptm <= 1.0
      and 0.0 <= report['final_metrics']['mean_plddt'] <= 1.0,
      f'pTM={ptm} pLDDT={report["final_metrics"].get("mean_plddt")}',
  )

  n_cys = sequence.count('C')
  ev.check(
      'SCIENCE: no-cysteine guidance yields ZERO cysteines',
      n_cys == 0,
      f'{n_cys} cysteines in {sequence}',
  )

  # ----------------------------------------------------------------------
  # SCIENTIFIC: guidance actually moves the objective.  (48 credits)
  # Guided vs unguided through IDENTICAL machinery -- the unguided arm is the
  # same decoding loop with candidate selection switched off -- so this is an
  # exact ablation of the selection step, not a comparison of two different
  # algorithms. GRAVY is stochastic, so we average over 2 seeds.
  # ----------------------------------------------------------------------
  guided, unguided, plain_cys = [], [], []
  for seed in (0, 1):
    path = out / f'cmp_gravy_seed{seed}.json'
    run_script(
        SKILL, SCRIPT, 'compare',
        '--length', str(LENGTH), '--objective', 'hydrophobicity',
        '--num-decoding-steps', '3', '--num-samples-per-step', '4',
        '--seed', str(seed), '--output', str(path),
        cassette_ns=f'cmp_gravy_seed{seed}',
    )
    cmp_report = load_json(path)
    guided.append(cmp_report['guided']['objective_value'])
    unguided.append(cmp_report['unguided']['objective_value'])
    plain_cys.append(
        cmp_report['plain_generate_baseline']['metrics']['n_cysteine']
    )
    ev.check(
        f'compare(seed={seed}) reports both arms + a plain baseline',
        {'guided', 'unguided', 'plain_generate_baseline'} <= set(cmp_report),
    )

  mean_guided = statistics.mean(guided)
  mean_unguided = statistics.mean(unguided)
  ev.check(
      'SCIENCE: guidance improves GRAVY vs the unguided ablation (2-seed mean)',
      mean_guided > mean_unguided,
      f'guided {mean_guided:+.3f} vs unguided {mean_unguided:+.3f} '
      f'(per-seed guided={[round(x, 3) for x in guided]}, '
      f'unguided={[round(x, 3) for x in unguided]})',
  )

  # Context for the report: what an UNGUIDED design of the same length gives.
  print(
      f'\n  [context] unguided plain-generate 48-mers carried '
      f'{plain_cys} cysteines (guided no-cysteine run: {n_cys})\n'
  )

  # ----------------------------------------------------------------------
  # SCIENTIFIC: the ptm objective runs end to end.  (11 credits)
  # We deliberately do NOT assert guided-pTM > unguided-pTM: at 2 steps x 2
  # samples the selection sees four candidates total, and single-sequence pTM on
  # a 48-mer is noisy (unguided 48-mers measured 0.09-0.26 here). Asserting an
  # improvement on that budget would be a coin flip dressed up as science. We
  # assert the objective completes and returns a physically valid pTM.
  # The --constraint flag rides along for free: `gravy` is a cheap metric, so it
  # adds no fold calls, and it exercises parse -> violation -> feasibility-first
  # ranking -> report end to end.
  # ----------------------------------------------------------------------
  ptm_json = out / 'ptm.json'
  run_script(
      SKILL, SCRIPT, 'generate',
      '--length', str(LENGTH), '--objective', 'ptm',
      '--num-decoding-steps', '2', '--num-samples-per-step', '2',
      '--constraint', 'gravy>=-1.0',
      '--seed', '1', '--output', str(ptm_json),
      cassette_ns='ptm',
  )
  ptm_report = load_json(ptm_json)
  final_ptm = ptm_report['final_metrics']['ptm']
  ev.check(
      'SCIENCE: ptm objective completes and returns pTM in [0, 1]',
      isinstance(final_ptm, (int, float)) and 0.0 <= final_ptm <= 1.0,
      f'pTM={final_ptm}',
  )
  ev.check(
      'ptm run folds every candidate (trajectory carries pTM values)',
      all(t['best_value'] is not None for t in ptm_report['trajectory'])
      and len(ptm_report['trajectory']) == 2,
  )

  satisfied = ptm_report['constraint_satisfied']
  ev.check(
      'constraint is parsed, enforced and reported',
      ptm_report['constraint'] == 'gravy>=-1'
      and isinstance(satisfied, bool),
      f'constraint={ptm_report["constraint"]} satisfied={satisfied}',
  )
  # Self-consistency: a constraint reported as satisfied must actually hold.
  # This cannot flake -- it checks the code against itself, not against the model.
  ev.check(
      'constraint_satisfied agrees with the final metrics',
      (not satisfied) or ptm_report['final_metrics']['gravy'] >= -1.0,
      f'gravy={ptm_report["final_metrics"]["gravy"]} satisfied={satisfied}',
  )

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
