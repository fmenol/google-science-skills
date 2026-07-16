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

"""Eval for the esmc_mutation_effect_scoring skill.

Structural checks: the CLI writes the tensors and files it promises, and it
refuses to score a variant whose stated wild-type residue does not match the
sequence.

Scientific checks (the ones that actually matter):
  * The wild-type column of the LLR matrix is exactly 0 at every position. This
    is what proves the (L, 20) matrix rows are aligned to the sequence.
  * Ubiquitin's pseudo-perplexity is far below that of its own scramble. A real
    protein is vastly more predictable than a shuffle of its own composition;
    if the BOS offset were wrong this separation would collapse.
  * Substituting to proline is markedly worse than a conservative aliphatic
    swap. Proline breaks backbone hydrogen bonding in helices and strands, and
    ESM must know it.
  * Per-position entropy has real spread, which catches the failure mode where
    every position returns the same distribution.

Measured on esmc-600m-2024-12 (deterministic forward pass):
  pseudo-perplexity  ubiquitin 1.05  vs  scrambled 18.70   (17.8x)
  WT recovery        ubiquitin 100%  vs  scrambled 6.6%
  mean LLR -> Pro    -9.35  vs  aliphatic swaps -6.73      (gap 2.62)
  entropy            mean 0.23, std 0.34 bits
Thresholds below sit well inside those margins.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy",
# ]
# ///

from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from fixtures import (
    Eval, REPO, UBIQUITIN, load_json, run_script, scramble, write_fasta,
)

SKILL = 'esmc_mutation_effect_scoring'
SCRIPT = 'mutation_scoring.py'
MODEL = 'esmc-600m-2024-12'

# Leave-one-out scoring costs one API request per residue, and the polite-http
# rate limiter is shared across every process on the machine. When the sibling
# ESM evals run in parallel they contend for the same budget, so a scan that
# takes 30s alone can take many minutes under load. These timeouts are sized for
# the contended case; run_evals.py caps the whole skill at 2400s regardless.
SCAN_TIMEOUT = 1500  # 77 requests
PPPL_TIMEOUT = 1800  # 154 requests (ubiquitin + its scramble)
CHEAP_TIMEOUT = 600  # <= 7 requests, or none at all

AA20 = list('ACDEFGHIKLMNPQRSTVWY')
AA_COL = {aa: i for i, aa in enumerate(AA20)}
ALIPHATIC = ('I', 'L', 'V', 'M')


def run_expecting_failure(*args: str) -> subprocess.CompletedProcess:
  """Runs the CLI and returns the result WITHOUT raising on a non-zero exit."""
  path = REPO / 'skills' / SKILL / 'scripts' / SCRIPT
  return subprocess.run(
      ['uv', 'run', '--no-project', '--quiet', str(path), *args],
      capture_output=True, text=True, timeout=300, cwd=REPO, check=False,
  )


def main() -> int:
  ev = Eval(SKILL)
  work = pathlib.Path(tempfile.mkdtemp(prefix='eval_esmc_mut_'))
  scan_dir = work / 'scan'

  # ---------------------------------------------------------------- scan ---
  run_script(
      SKILL, SCRIPT, 'scan',
      '--sequence', UBIQUITIN,
      '--model', MODEL,
      '--output-dir', str(scan_dir),
      '--top', '10',
      timeout=SCAN_TIMEOUT,
  )

  written = {
      name: (scan_dir / name)
      for name in ('llr.npy', 'llr.csv', 'positions.csv', 'summary.json')
  }
  ev.check(
      'scan writes all four output files',
      all(p.is_file() and p.stat().st_size > 0 for p in written.values()),
      ', '.join(
          f'{n}={"ok" if p.is_file() and p.stat().st_size else "MISSING"}'
          for n, p in written.items()
      ),
  )

  llr = np.load(scan_dir / 'llr.npy')
  summary = load_json(scan_dir / 'summary.json')
  entropy = np.asarray(summary['entropy_bits'], dtype=np.float64)

  ev.check(
      'LLR matrix shape is (76, 20)',
      llr.shape == (len(UBIQUITIN), 20),
      f'{llr.shape}',
  )
  ev.check(
      'entropy vector shape is (76,)',
      entropy.shape == (len(UBIQUITIN),),
      f'{entropy.shape}',
  )
  ev.check(
      'summary reports the amino-acid column order',
      summary.get('amino_acid_columns') == AA20,
      f'{summary.get("amino_acid_columns")}',
  )

  # ----------------------------------------- SCIENTIFIC: LLR(WT) == 0 ------
  # Not merely an identity: it is what proves row i of the matrix really is
  # residue i. Any row misalignment makes these entries non-zero.
  wt_llr = np.array(
      [llr[i, AA_COL[aa]] for i, aa in enumerate(UBIQUITIN)]
  )
  ev.check(
      'LLR of the wild-type residue is 0 at every position',
      bool(np.all(np.abs(wt_llr) < 1e-5)),
      f'max |LLR(WT)| = {np.abs(wt_llr).max():.2e} over '
      f'{len(UBIQUITIN)} positions',
  )

  # ------------------------------- SCIENTIFIC: constraint is not uniform ---
  ev.check(
      'per-position entropy has real spread (std > 0.2 bits)',
      float(entropy.std()) > 0.2,
      f'std = {entropy.std():.3f} bits '
      f'(min {entropy.min():.3f}, max {entropy.max():.3f})',
  )

  # ------------------------------------ SCIENTIFIC: proline disruption -----
  # X->P averaged over every position that is not already proline, against
  # conservative swaps within the aliphatic set at aliphatic positions.
  pro_llr = np.array([
      llr[i, AA_COL['P']]
      for i, aa in enumerate(UBIQUITIN)
      if aa != 'P'
  ])
  conservative_llr = np.array([
      llr[i, AA_COL[sub]]
      for i, aa in enumerate(UBIQUITIN)
      if aa in ALIPHATIC
      for sub in ALIPHATIC
      if sub != aa
  ])
  gap = float(conservative_llr.mean() - pro_llr.mean())
  ev.check(
      'substitution to proline is markedly worse than a conservative '
      'aliphatic swap',
      gap > 1.0,
      f'mean LLR->Pro = {pro_llr.mean():.3f} (n={len(pro_llr)}) vs '
      f'aliphatic swap = {conservative_llr.mean():.3f} '
      f'(n={len(conservative_llr)}); gap = {gap:.3f}',
  )
  ev.check(
      'substitution to proline is strongly deleterious in absolute terms',
      float(pro_llr.mean()) < -3.0,
      f'mean LLR->Pro = {pro_llr.mean():.3f}',
  )

  # ------------------------------------------------------------- heatmap ---
  heatmap_png = work / 'heatmap.png'
  run_script(
      SKILL, SCRIPT, 'heatmap',
      '--scan-dir', str(scan_dir),
      '--output', str(heatmap_png),
      timeout=CHEAP_TIMEOUT,
  )
  ev.check(
      'heatmap writes a non-empty PNG',
      heatmap_png.is_file() and heatmap_png.stat().st_size > 5000,
      f'{heatmap_png.stat().st_size if heatmap_png.is_file() else 0} bytes',
  )

  # ------------------------------------------------------------- entropy ---
  entropy_png = work / 'entropy.png'
  run_script(
      SKILL, SCRIPT, 'entropy',
      '--scan-dir', str(scan_dir),
      '--output', str(entropy_png),
      timeout=CHEAP_TIMEOUT,
  )
  ev.check(
      'entropy writes a non-empty PNG and its CSV',
      (
          entropy_png.is_file() and entropy_png.stat().st_size > 5000
          and entropy_png.with_suffix('.csv').is_file()
      ),
      f'{entropy_png.stat().st_size if entropy_png.is_file() else 0} bytes',
  )

  # ------------------------------------------------------ score-variants ---
  # I44 and L8 line the hydrophobic patch; K48 is a polyubiquitin linkage site;
  # G76 is the C-terminal conjugation residue.
  variants_json = work / 'variants.json'
  run_script(
      SKILL, SCRIPT, 'score-variants',
      '--sequence', UBIQUITIN,
      '--model', MODEL,
      '--variants', 'I44A,K48R,L8A,G76A,M1P',
      '--output', str(variants_json),
      timeout=CHEAP_TIMEOUT,
  )
  scored = load_json(variants_json)
  records = scored['variants_ranked_most_to_least_deleterious']
  by_name = {r['variant']: r for r in records}
  ev.check(
      'score-variants scores every 1-indexed variant it was given',
      set(by_name) == {'I44A', 'K48R', 'L8A', 'G76A', 'M1P'},
      f'{sorted(by_name)}',
  )
  ev.check(
      'score-variants ranks most-deleterious first',
      all(
          records[i]['llr'] <= records[i + 1]['llr']
          for i in range(len(records) - 1)
      ),
      ', '.join(f'{r["variant"]}={r["llr"]:.2f}' for r in records),
  )
  # Cross-check the standalone path against the full scan: the same variant,
  # scored two different ways, must give the same number.
  ev.check(
      'score-variants agrees with the full scan on K48R',
      abs(
          by_name['K48R']['llr']
          - float(llr[47, AA_COL['R']])
      ) < 1e-3,
      f'score-variants={by_name["K48R"]["llr"]:.4f} vs '
      f'scan={float(llr[47, AA_COL["R"]]):.4f}',
  )

  # A stated wild-type residue that does not match the sequence must be fatal.
  # UBIQUITIN[47] is K, so 'A48R' is a lie and must not be silently scored.
  bad = run_expecting_failure(
      'score-variants',
      '--sequence', UBIQUITIN,
      '--model', MODEL,
      '--variants', 'A48R',
      '--output', str(work / 'never_written.json'),
  )
  message = (bad.stdout + bad.stderr).lower()
  ev.check(
      'a variant whose wild-type residue mismatches the sequence exits '
      'non-zero',
      bad.returncode != 0,
      f'exit={bad.returncode}',
  )
  ev.check(
      'the mismatch error names the position, the expected and the actual '
      'residue',
      all(token in message for token in ('a48r', '48', 'wild-type'))
      and 'k' in message,
      repr((bad.stderr.strip() or bad.stdout.strip())[:170]),
  )
  ev.check(
      'nothing is written when validation fails',
      not (work / 'never_written.json').exists(),
  )

  # ---------------- SCIENTIFIC: pseudo-perplexity, real << scrambled -------
  # Both sequences go in one FASTA so they are scored in a single invocation,
  # under identical conditions, and the subcommand's own ranking is exercised.
  scrambled_sequence = scramble(UBIQUITIN, seed=0)
  pppl_fasta = write_fasta(
      work / 'pppl.fasta',
      [('ubiquitin', UBIQUITIN), ('ubiquitin_scrambled', scrambled_sequence)],
  )
  pppl_json = work / 'pppl.json'
  run_script(
      SKILL, SCRIPT, 'pseudo-perplexity',
      '--fasta', pppl_fasta,
      '--model', MODEL,
      '--output', str(pppl_json),
      timeout=PPPL_TIMEOUT,
  )
  ranked = load_json(pppl_json)['results_ranked_best_to_worst']
  by_id = {record['id']: record for record in ranked}
  real, scrambled = by_id['ubiquitin'], by_id['ubiquitin_scrambled']
  real_pppl = float(real['pseudo_perplexity'])
  scrambled_pppl = float(scrambled['pseudo_perplexity'])

  ev.check(
      'pseudo-perplexity ranks the real protein above its scramble',
      ranked[0]['id'] == 'ubiquitin',
      f'best = {ranked[0]["id"]}',
  )

  ev.check(
      'pseudo-perplexity: ubiquitin is strictly lower than its own scramble',
      real_pppl < scrambled_pppl,
      f'ubiquitin {real_pppl:.3f} vs scrambled {scrambled_pppl:.3f}',
  )
  ev.check(
      'pseudo-perplexity: the gap is substantial (>3x), not marginal',
      scrambled_pppl > 3.0 * real_pppl,
      f'scrambled/real = {scrambled_pppl / real_pppl:.1f}x',
  )
  # A wrong BOS offset would read every position's distribution from its
  # neighbour, and wild-type recovery would collapse to chance (~5%).
  ev.check(
      'ESM C recovers ubiquitin under leave-one-out masking (BOS offset is '
      'right)',
      float(real['wt_recovery_fraction']) > 0.7,
      f'WT recovery: ubiquitin {real["wt_recovery_fraction"]:.1%} vs '
      f'scrambled {scrambled["wt_recovery_fraction"]:.1%}',
  )
  # The scan and the standalone subcommand must agree; they share one code path.
  ev.check(
      'scan and pseudo-perplexity report the same value for ubiquitin',
      abs(float(summary['pseudo_perplexity']) - real_pppl) < 1e-2,
      f'scan={summary["pseudo_perplexity"]:.4f} vs '
      f'pseudo-perplexity={real_pppl:.4f}',
  )

  print(
      f'\n  [numbers] pseudo-perplexity: ubiquitin {real_pppl:.4f}, '
      f'scrambled {scrambled_pppl:.4f} ({scrambled_pppl / real_pppl:.1f}x)'
  )
  print(
      f'  [numbers] mean LLR->Pro {pro_llr.mean():.3f}, '
      f'aliphatic swap {conservative_llr.mean():.3f}, gap {gap:.3f}'
  )
  print(f'  [artifacts] {work}')
  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
