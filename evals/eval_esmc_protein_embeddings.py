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

"""Eval for the esmc_protein_embeddings skill.

Structural checks: the CLI writes the shapes and sidecars it promises.
Scientific checks (these define "working"):
  * Determinism   — ESMC is a deterministic encoder; the same sequence embedded
                    twice must give the identical vector.
  * Semantics     — a point mutant of ubiquitin must embed closer to ubiquitin
                    than lysozyme does. If this fails the embedding carries no
                    biology and every downstream use is void.
  * Clustering    — KMeans over 5 lysozyme variants + 5 barnase variants must
                    recover the two families (adjusted Rand index > 0.8).

Measured on the live API (esmc-600m-2024-12, mean-pooled output embedding):
  cos(ubiquitin, ubiquitin)          = 1.000000  (bit-identical)
  cos(ubiquitin, ubiquitin R42F)     = 0.984
  cos(ubiquitin, lysozyme)           = 0.507
  cos(ubiquitin, scrambled ubiquitin)= -0.103
  adjusted Rand index (lysozyme vs barnase variants, k=2) = 1.00
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
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

import fixtures as fx  # pylint: disable=g-import-not-at-top

SKILL = 'esmc_protein_embeddings'
SCRIPT = 'embed.py'
AA20 = 'ACDEFGHIKLMNPQRSTVWY'


def point_mutant(sequence: str, seed: int) -> tuple[str, str]:
  """Substitutes one residue at random. Returns (sequence, 'W62A'-style tag)."""
  rng = random.Random(seed)
  i = rng.randrange(len(sequence))
  alt = rng.choice([a for a in AA20 if a != sequence[i]])
  mutated = sequence[:i] + alt + sequence[i + 1 :]
  return mutated, f'{sequence[i]}{i + 1}{alt}'


def cosine(a: np.ndarray, b: np.ndarray) -> float:
  a = np.asarray(a, dtype=np.float64).ravel()
  b = np.asarray(b, dtype=np.float64).ravel()
  return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def read_matrix_csv(path) -> tuple[list[str], np.ndarray]:
  """Reads the `similarity` CSV into (ids, matrix)."""
  with open(path, encoding='utf-8') as handle:
    rows = list(csv.reader(handle))
  ids = rows[0][1:]
  matrix = np.array([[float(v) for v in row[1:]] for row in rows[1:]])
  return ids, matrix


def main() -> int:
  ev = fx.Eval(SKILL)
  work = pathlib.Path(tempfile.mkdtemp(prefix='eval_esmc_embeddings_'))
  print(f'[{SKILL}] workdir: {work}\n')

  # ---------------------------------------------------------------- structural
  # 1. Mean-pooled embedding of ubiquitin with the default 600M model.
  fx.run_script(
      SKILL, SCRIPT, 'embed',
      '--sequence', fx.UBIQUITIN,
      '--model', 'esmc-600m-2024-12',
      '--output', str(work / 'ubi.npy'),
  )
  ubi = np.load(work / 'ubi.npy')
  ev.check(
      'embed: mean embedding is (1152,)',
      ubi.shape == (1152,),
      f'got {ubi.shape}',
  )
  meta = fx.load_json(work / 'ubi.json')
  ev.check(
      'embed: sidecar metadata is written',
      meta.get('model') == 'esmc-600m-2024-12'
      and meta.get('shape') == [1152]
      and meta.get('sequence_length') == 76,
      f'model={meta.get("model")} shape={meta.get("shape")}',
  )

  # 2. Per-residue embedding: BOS/EOS trimmed, so exactly L rows.
  fx.run_script(
      SKILL, SCRIPT, 'embed',
      '--sequence', fx.UBIQUITIN,
      '--model', 'esmc-600m-2024-12',
      '--per-residue',
      '--output', str(work / 'ubi_pr.npy'),
  )
  per_residue = np.load(work / 'ubi_pr.npy')
  ev.check(
      'embed --per-residue: is (76, 1152) — BOS/EOS trimmed',
      per_residue.shape == (76, 1152),
      f'got {per_residue.shape} for a {len(fx.UBIQUITIN)}-residue protein',
  )

  # 3. A named hidden layer of the 300M model: D = 960.
  fx.run_script(
      SKILL, SCRIPT, 'embed',
      '--sequence', fx.UBIQUITIN,
      '--model', 'esmc-300m-2024-12',
      '--layer', '12',
      '--output', str(work / 'ubi_300m_l12.npy'),
  )
  layer12 = np.load(work / 'ubi_300m_l12.npy')
  ev.check(
      'embed --model esmc-300m-2024-12 --layer 12: is (960,)',
      layer12.shape == (960,),
      f'got {layer12.shape}',
  )

  # 4. Batch over the six shared fixtures.
  six = [
      ('ubiquitin', fx.UBIQUITIN), ('lysozyme', fx.LYSOZYME),
      ('barnase', fx.BARNASE), ('barstar', fx.BARSTAR),
      ('ca2', fx.CA2), ('gfp', fx.GFP),
  ]
  six_fasta = fx.write_fasta(work / 'six.fasta', six)
  fx.run_script(
      SKILL, SCRIPT, 'batch',
      '--fasta', six_fasta,
      '--model', 'esmc-600m-2024-12',
      '--output', str(work / 'six.npy'),
      '--max-workers', '6',
  )
  matrix = np.load(work / 'six.npy')
  ids = fx.load_json(work / 'six.json')['ids']
  ev.check(
      'batch: 6-record FASTA -> (6, 1152) matrix',
      matrix.shape == (6, 1152),
      f'got {matrix.shape}',
  )
  ev.check(
      'batch: ids preserved in FASTA order',
      ids == [name for name, _ in six],
      f'got {ids}',
  )
  # Row order must match ids, not just the id list: row 0 must BE ubiquitin.
  ev.check(
      'batch: row i really is the embedding of ids[i]',
      cosine(matrix[0], ubi) > 0.999,
      f'cos(batch row 0, standalone ubiquitin embedding) = '
      f'{cosine(matrix[0], ubi):.6f}',
  )

  # ---------------------------------------------------------------- scientific
  # 5. Determinism. ESMC is a deterministic encoder; a second call on the same
  #    sequence must return the same vector, or nothing downstream is stable.
  fx.run_script(
      SKILL, SCRIPT, 'embed',
      '--sequence', fx.UBIQUITIN,
      '--model', 'esmc-600m-2024-12',
      '--output', str(work / 'ubi_again.npy'),
  )
  ubi_again = np.load(work / 'ubi_again.npy')
  self_cos = cosine(ubi, ubi_again)
  ev.check(
      'DETERMINISM: ubiquitin embedded twice is identical',
      abs(self_cos - 1.0) < 1e-5,
      f'cosine = {self_cos:.8f} (want 1.0 +/- 1e-5); '
      f'max|delta| = {np.abs(ubi - ubi_again).max():.3g}',
  )

  # 6. Semantics, via the `similarity` subcommand: a single point mutation must
  #    perturb the embedding far less than swapping in an unrelated protein.
  ubi_mut, tag = point_mutant(fx.UBIQUITIN, seed=7)
  sim_fasta = fx.write_fasta(
      work / 'semantics.fasta',
      [
          ('ubiquitin', fx.UBIQUITIN),
          (f'ubiquitin_{tag}', ubi_mut),
          ('lysozyme', fx.LYSOZYME),
          ('ubiquitin_scrambled', fx.scramble(fx.UBIQUITIN, seed=1)),
      ],
  )
  fx.run_script(
      SKILL, SCRIPT, 'similarity',
      '--fasta', sim_fasta,
      '--model', 'esmc-600m-2024-12',
      '--top', '3',
      '--output', str(work / 'sim.csv'),
  )
  ev.check(
      'similarity: writes the CSV matrix and the JSON sidecar',
      (work / 'sim.csv').is_file() and (work / 'sim.json').is_file(),
      'sim.csv + sim.json',
  )
  sim_ids, sim = read_matrix_csv(work / 'sim.csv')
  index = {name: i for i, name in enumerate(sim_ids)}
  near = sim[index['ubiquitin'], index[f'ubiquitin_{tag}']]
  far = sim[index['ubiquitin'], index['lysozyme']]
  scrambled = sim[index['ubiquitin'], index['ubiquitin_scrambled']]
  ev.check(
      'SEMANTICS: cos(ubiquitin, point mutant) > cos(ubiquitin, lysozyme)',
      near > far,
      f'{tag} mutant = {near:.4f} vs lysozyme = {far:.4f} '
      f'(scrambled ubiquitin = {scrambled:.4f})',
  )
  top = fx.load_json(work / 'sim.json')['top_pairs'][0]
  ev.check(
      'similarity: top-ranked pair is the near-identical mutant',
      {top['a'], top['b']} == {'ubiquitin', f'ubiquitin_{tag}'},
      f'top pair = {top["a"]} <-> {top["b"]} at {top["cosine"]:.4f}',
  )

  # 7. Clustering: 5 point mutants of lysozyme + 5 of barnase. Each variant is a
  #    single substitution away from its parent, so the only signal KMeans can
  #    use is protein family. It must recover it.
  variants: list[tuple[str, str]] = []
  truth: list[tuple[str, str]] = []
  for k in range(5):
    seq, tag_l = point_mutant(fx.LYSOZYME, seed=100 + k)
    variants.append((f'lysozyme_{tag_l}', seq))
    truth.append((f'lysozyme_{tag_l}', 'lysozyme'))
  for k in range(5):
    seq, tag_b = point_mutant(fx.BARNASE, seed=200 + k)
    variants.append((f'barnase_{tag_b}', seq))
    truth.append((f'barnase_{tag_b}', 'barnase'))

  families_fasta = fx.write_fasta(work / 'families.fasta', variants)
  labels_csv = work / 'families_labels.csv'
  with open(labels_csv, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.writer(handle)
    writer.writerow(['id', 'label'])
    writer.writerows(truth)

  fx.run_script(
      SKILL, SCRIPT, 'cluster',
      '--fasta', families_fasta,
      '--model', 'esmc-600m-2024-12',
      '--clusters', '2',
      '--labels', str(labels_csv),
      '--output', str(work / 'clusters.csv'),
      '--max-workers', '5',
  )
  ev.check(
      'cluster: writes assignments CSV, summary JSON and PCA scatter PNG',
      all(
          (work / f'clusters{suffix}').is_file()
          for suffix in ('.csv', '.json', '.png')
      ),
      'clusters.csv + clusters.json + clusters.png',
  )
  summary = fx.load_json(work / 'clusters.json')
  ari = summary['adjusted_rand_index']
  ev.check(
      'CLUSTERING: KMeans(k=2) recovers lysozyme vs barnase, ARI > 0.8',
      ari > 0.8,
      f'adjusted Rand index = {ari:.3f} over 5 lysozyme + 5 barnase variants',
  )
  ev.check(
      'cluster: every sequence is assigned and both clusters are used',
      len(summary['assignments']) == 10 and len(summary['cluster_sizes']) == 2,
      f'sizes = {summary["cluster_sizes"]}',
  )

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
