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

"""Eval for the esmc_sae_feature_interpretation skill.

Structural checks
  * `extract` on ubiquitin gives (76, 64) indices/values, codebook 16384, and
    exactly 64 active features at every position (the SAE is top-k=64).
  * `top-features` emits both rankings with feature indices inside the codebook.
  * `describe` returns a non-empty `description` for each of the top-5 features
    (proves the unauthenticated feature-description API is wired up).
  * `plot` writes a non-empty PNG.

Scientific checks
  * Features are protein-specific. A single-point mutant of lysozyme must share
    far more feature machinery with lysozyme than ubiquitin does:
        jaccard(lysozyme, lysozyme_E35Q) > jaccard(lysozyme, ubiquitin)
    Guarded by an anti-degeneracy check: the mutant's activations must actually
    differ from wild-type (otherwise a Jaccard of 1.0 would only prove that the
    API ignored our input), and the largest perturbation must land on the
    mutated residue.
  * Descriptions are real annotations. Lysozyme's top features must mention
    terms diagnostic of *this* protein's biology (peptidoglycan / muramidase /
    lysozyme / glycoside / hydrolase / cell wall). To prove those terms are
    diagnostic rather than generic, ubiquitin is run as a negative control and
    must NOT hit them.

    Keywords were chosen only AFTER reading the live descriptions. 'glycan' and
    'disulfide' were considered and DELIBERATELY REJECTED: they leak into
    ubiquitin (feature 7865, 'ERAD-p97 UBL/UBX proteostasis', legitimately
    mentions N-glycanase), so they are not diagnostic.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy"]
# ///

from __future__ import annotations

import pathlib
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from fixtures import Eval, LYSOZYME, UBIQUITIN, load_json, run_script  # noqa: E402

SKILL = 'esmc_sae_feature_interpretation'
SCRIPT = 'sae_features.py'
CODEBOOK = 16384
TOP_K = 64

# Hen egg-white lysozyme (P00698). Residue 35 is the catalytic Glu; E35Q is the
# classic activity-killing point mutation. One residue out of 129 changes.
MUT_POS = 35  # 1-based
assert LYSOZYME[MUT_POS - 1] == 'E', 'fixture drift: residue 35 is not Glu'
LYSOZYME_E35Q = LYSOZYME[: MUT_POS - 1] + 'Q' + LYSOZYME[MUT_POS:]

# Terms genuinely diagnostic of lysozyme biology. Verified against the live
# feature-description API: 9/10 of lysozyme's top features hit these; 0/10 of
# ubiquitin's do.
DIAGNOSTIC = (
    'peptidoglycan',
    'lysozyme',
    'muramidase',
    'glycoside',
    'hydrolase',
    'cell wall',
    'cell-wall',
)


def describe_hits(payload: dict, limit: int) -> list[tuple[int, str, list[str]]]:
  """For the top-`limit` features by max activation, which diagnostic terms hit."""
  out = []
  descriptions = payload.get('descriptions', {})
  for item in payload['top_by_max_activation'][:limit]:
    idx = item['feature_index']
    record = descriptions.get(str(idx), {})
    blob = ' '.join(
        str(record.get(field, '') or '')
        for field in (
            'label',
            'summary',
            'description',
            'category',
            'activation_pattern',
            'exemplar_protein_families',
        )
    ).lower()
    hits = sorted({term for term in DIAGNOSTIC if term in blob})
    out.append((idx, str(record.get('label') or ''), hits))
  return out


def main() -> int:
  ev = Eval(SKILL)
  tmp = pathlib.Path(tempfile.mkdtemp(prefix='sae_eval_'))

  # ---------------------------------------------------------------- extract
  ubq_npz = tmp / 'ubiquitin.npz'
  run_script(
      SKILL, SCRIPT, 'extract',
      '--sequence', UBIQUITIN,
      '--output', str(ubq_npz),
      '--output-summary', str(tmp / 'ubiquitin_summary.json'),
  )
  ev.check('extract writes .npz', ubq_npz.is_file(), str(ubq_npz))

  with np.load(ubq_npz, allow_pickle=False) as data:
    indices = data['indices']
    values = data['values']
    codebook = int(data['codebook_size'])
    sequence = str(data['sequence'])

  ev.check(
      'indices shape == (76, 64)',
      indices.shape == (len(UBIQUITIN), TOP_K),
      f'{indices.shape}',
  )
  ev.check(
      'values shape == (76, 64)',
      values.shape == (len(UBIQUITIN), TOP_K),
      f'{values.shape}',
  )
  ev.check('codebook == 16384', codebook == CODEBOOK, f'{codebook}')
  ev.check(
      'BOS/EOS trimmed (rows == residues)',
      sequence == UBIQUITIN and indices.shape[0] == len(UBIQUITIN),
      f'{indices.shape[0]} rows for {len(UBIQUITIN)} residues',
  )

  # Sparsity: top-k means exactly 64 active features at every position.
  active = (values > 0).sum(axis=1)
  ev.check(
      'exactly 64 active features per position',
      bool((active == TOP_K).all()),
      f'min={int(active.min())} max={int(active.max())}',
  )
  distinct_per_row = [len(set(row.tolist())) for row in indices]
  ev.check(
      '64 DISTINCT features per position',
      all(n == TOP_K for n in distinct_per_row),
      f'min={min(distinct_per_row)}',
  )
  ev.check(
      'feature indices within codebook',
      bool((indices >= 0).all() and (indices < CODEBOOK).all()),
      f'range {int(indices.min())}..{int(indices.max())}',
  )

  summary = load_json(tmp / 'ubiquitin_summary.json')
  ev.check(
      'summary reports top_k=64, codebook=16384',
      summary['top_k'] == TOP_K and summary['codebook_size'] == CODEBOOK,
      f'top_k={summary["top_k"]} codebook={summary["codebook_size"]}',
  )

  # ----------------------------------------------------------- top-features
  top_json = tmp / 'ubiquitin_top.json'
  run_script(
      SKILL, SCRIPT, 'top-features',
      '--features', str(ubq_npz),
      '--limit', '10',
      '--output', str(top_json),
  )
  ranking = load_json(top_json)
  by_max = ranking.get('top_by_max_activation', [])
  by_prev = ranking.get('top_by_prevalence', [])
  ev.check(
      'top-features emits BOTH rankings',
      len(by_max) == 10 and len(by_prev) == 10,
      f'{len(by_max)} by max, {len(by_prev)} by prevalence',
  )
  all_ids = [i['feature_index'] for i in by_max + by_prev]
  ev.check(
      'ranked indices in [0, 16384)',
      all(0 <= i < CODEBOOK for i in all_ids),
      f'{min(all_ids)}..{max(all_ids)}',
  )
  ev.check(
      'max-activation ranking is monotonically decreasing',
      all(
          by_max[i]['max_activation'] >= by_max[i + 1]['max_activation']
          for i in range(len(by_max) - 1)
      ),
  )
  ev.check(
      'prevalence ranking is monotonically decreasing',
      all(
          by_prev[i]['prevalence'] >= by_prev[i + 1]['prevalence']
          for i in range(len(by_prev) - 1)
      ),
  )

  # -------------------------------------------------------------- describe
  ubq_desc = tmp / 'ubiquitin_described.json'
  run_script(
      SKILL, SCRIPT, 'describe',
      '--features', str(ubq_npz),
      '--limit', '10',
      '--output', str(ubq_desc),
      '--report', str(tmp / 'ubiquitin_report.md'),
  )
  described = load_json(ubq_desc)
  top5 = described['top_by_max_activation'][:5]
  non_empty = [
      bool(str(item.get('description') or '').strip()) for item in top5
  ]
  ev.check(
      'describe: non-empty description for each of the top-5',
      len(top5) == 5 and all(non_empty),
      f'{sum(non_empty)}/5 populated',
  )
  ev.check(
      'describe: labels populated too',
      all(str(i.get('label') or '').strip() for i in top5),
      '; '.join(str(i.get('label'))[:28] for i in top5[:3]),
  )
  ev.check(
      'describe: markdown report written',
      (tmp / 'ubiquitin_report.md').stat().st_size > 500,
  )

  # Ubiquitin's #1 feature should be about ubiquitin. Not a pass/fail gate on
  # a specific index, but a strong signal the pipeline is sane.
  ubq_top_label = str(top5[0].get('label') or '')
  ev.check(
      "ubiquitin's top feature is ubiquitin-related",
      'ubiquitin' in ubq_top_label.lower(),
      f'feature {top5[0]["feature_index"]}: {ubq_top_label}',
  )

  # ------------------------------------------------------------------ plot
  png = tmp / 'tracks.png'
  plot_ids = ','.join(str(i['feature_index']) for i in top5[:3])
  run_script(
      SKILL, SCRIPT, 'plot',
      '--features', str(ubq_npz),
      '--feature-ids', plot_ids,
      '--output', str(png),
  )
  ev.check(
      'plot writes a non-empty PNG',
      png.is_file() and png.stat().st_size > 5000,
      f'{png.stat().st_size if png.is_file() else 0} bytes',
  )

  # ============================ SCIENTIFIC ================================
  # 1. Features are protein-specific, not generic.
  lyz_npz = tmp / 'lysozyme.npz'
  mut_npz = tmp / 'lysozyme_e35q.npz'
  run_script(
      SKILL, SCRIPT, 'extract',
      '--sequence', LYSOZYME, '--output', str(lyz_npz),
  )
  run_script(
      SKILL, SCRIPT, 'extract',
      '--sequence', LYSOZYME_E35Q, '--output', str(mut_npz),
  )

  def dense(path):
    with np.load(path, allow_pickle=False) as d:
      out = np.zeros((d['indices'].shape[0], int(d['codebook_size'])))
      np.put_along_axis(
          out, d['indices'].astype(np.int64), d['values'].astype(np.float64), 1
      )
      return out

  d_lyz, d_mut = dense(lyz_npz), dense(mut_npz)

  # ANTI-DEGENERACY GUARD. Without this, a Jaccard of 1.0 between wild-type and
  # mutant would be indistinguishable from the API silently ignoring the
  # mutation. The activations must actually move, and they must move most at
  # the residue we mutated.
  ev.check(
      'mutant activations differ from wild-type (API saw the mutation)',
      not np.array_equal(d_lyz, d_mut),
      f'total |delta| = {np.abs(d_lyz - d_mut).sum():.3f}',
  )
  per_residue_delta = np.abs(d_lyz - d_mut).sum(axis=1)
  ranked = [int(p) + 1 for p in np.argsort(-per_residue_delta)]
  # In practice residue 35 ranks #1. Allow the top-3 so that server-side float
  # jitter cannot flake the eval; this is still ~1-in-43 by chance and remains
  # impossible to pass if the API had ignored the substitution.
  ev.check(
      'perturbation is localised to the mutated residue (E35Q)',
      MUT_POS in ranked[:3],
      f'most-perturbed residues {ranked[:3]}, mutated residue is {MUT_POS}',
  )

  cmp_mut = tmp / 'cmp_lyz_mut.json'
  cmp_ubq = tmp / 'cmp_lyz_ubq.json'
  run_script(
      SKILL, SCRIPT, 'compare',
      '--sequence-a', LYSOZYME, '--sequence-b', LYSOZYME_E35Q,
      '--top-n', '50', '--output', str(cmp_mut),
  )
  run_script(
      SKILL, SCRIPT, 'compare',
      '--sequence-a', LYSOZYME, '--sequence-b', UBIQUITIN,
      '--top-n', '50', '--output', str(cmp_ubq),
  )
  j_mut = load_json(cmp_mut)['jaccard']
  j_ubq = load_json(cmp_ubq)['jaccard']

  ev.check(
      'jaccard(lysozyme, lysozyme_E35Q) > jaccard(lysozyme, ubiquitin)',
      j_mut > j_ubq,
      f'point mutant {j_mut:.3f} vs unrelated protein {j_ubq:.3f}',
  )
  ev.check(
      'point mutant shares most of its top-50 machinery',
      j_mut >= 0.7,
      f'jaccard = {j_mut:.3f}',
  )
  ev.check(
      'unrelated protein shares little machinery',
      j_ubq <= 0.15,
      f'jaccard = {j_ubq:.3f}',
  )

  # 2. Descriptions are real annotations, with a negative control.
  lyz_desc = tmp / 'lysozyme_described.json'
  run_script(
      SKILL, SCRIPT, 'describe',
      '--features', str(lyz_npz),
      '--limit', '10',
      '--output', str(lyz_desc),
      '--report', str(tmp / 'lysozyme_report.md'),
  )
  lyz_hits = describe_hits(load_json(lyz_desc), 10)
  ubq_hits = describe_hits(described, 10)

  n_lyz = sum(1 for _, _, hits in lyz_hits if hits)
  n_ubq = sum(1 for _, _, hits in ubq_hits if hits)

  print('\n  lysozyme top-10 feature labels:')
  for idx, label, hits in lyz_hits:
    mark = '*' if hits else ' '
    print(f'   {mark} {idx:5d}  {label[:52]:54s} {hits if hits else ""}')

  ev.check(
      'lysozyme descriptions mention its real biology',
      n_lyz >= 5,
      f'{n_lyz}/10 top features hit {list(DIAGNOSTIC)}',
  )
  # This is what makes the keyword test meaningful rather than a freebie: the
  # same terms must NOT fire on an unrelated protein.
  ev.check(
      'NEGATIVE CONTROL: ubiquitin does not hit lysozyme terms',
      n_ubq <= 2,
      f'{n_ubq}/10 ubiquitin features hit the lysozyme terms (want ~0)',
  )
  ev.check(
      'lysozyme is far more enriched for these terms than ubiquitin',
      n_lyz > n_ubq,
      f'lysozyme {n_lyz}/10 vs ubiquitin {n_ubq}/10',
  )

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
