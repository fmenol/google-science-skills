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

"""Embeds protein sequences with ESMC and analyses the embeddings."""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
#   "numpy",
#   "scikit-learn",
#   "matplotlib",
# ]
# ///

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')  # Headless: must be set before pyplot is imported.
import matplotlib.pyplot as plt  # pylint: disable=g-import-not-at-top

from esm_biohub import BiohubClient, BiohubError  # vendored, same directory
import esm_biohub as eb

# The Biohub API rejects a negative `ith_hidden_layer` with HTTP 500 (except -1,
# which means "all layers"). Negative indices are therefore resolved locally,
# against the (n_layers + 1)-row hidden-state stack, before any request is made.
MAX_WORKERS_CAP = 8


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def resolve_layer(model: str, layer: int | None) -> int | None:
  """Maps a possibly negative layer index onto the model's hidden-state stack.

  Hidden states have `n_layers + 1` rows: row 0 is the embedding layer and row
  `n_layers` is the final block. `-1` is the last row, `-2` the one before it.

  Args:
    model: An ESMC model name.
    layer: Layer index, possibly negative, or None for the output embedding.

  Returns:
    A non-negative layer index, or None.

  Raises:
    BiohubError: If the index is out of range for the model.
  """
  if layer is None:
    return None
  n_rows = eb.ESMC_MODELS[model]['n_layers'] + 1
  resolved = layer + n_rows if layer < 0 else layer
  if not 0 <= resolved < n_rows:
    raise BiohubError(
        f'--layer {layer} is out of range for {model}, which has {n_rows} '
        f'hidden-state rows (0..{n_rows - 1}; 0 is the embedding layer).'
    )
  return resolved


def load_records(args) -> list[tuple[str, str]]:
  """Reads (id, sequence) pairs from --fasta or --sequence."""
  if getattr(args, 'fasta', None):
    records = eb.read_fasta(args.fasta)
  else:
    records = [(getattr(args, 'id', None) or 'query', args.sequence)]
  return [(name, eb.validate_sequence(seq)) for name, seq in records]


def embed_many(
    client: BiohubClient,
    records: list[tuple[str, str]],
    model: str,
    layer: int | None,
    max_workers: int,
) -> np.ndarray:
  """Mean-pools one embedding per record, in parallel, preserving input order.

  Returns:
    An (N, D) float32 matrix whose row i is the embedding of records[i].
  """
  workers = max(1, min(max_workers, MAX_WORKERS_CAP, len(records)))
  vectors: list[np.ndarray | None] = [None] * len(records)

  with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    futures = {
        pool.submit(client.embed, seq, model, layer=layer): i
        for i, (_, seq) in enumerate(records)
    }
    done = 0
    for future in concurrent.futures.as_completed(futures):
      index = futures[future]
      vectors[index] = future.result()  # Re-raises BiohubError in the caller.
      done += 1
      print(f'  embedded {done}/{len(records)}', end='\r', file=sys.stderr)

  print(' ' * 40, end='\r', file=sys.stderr)
  return np.stack(vectors).astype(np.float32)  # type: ignore[arg-type]


def cosine_matrix(matrix: np.ndarray) -> np.ndarray:
  """Returns the (N, N) pairwise cosine-similarity matrix of the rows."""
  data = np.asarray(matrix, dtype=np.float64)
  norms = np.linalg.norm(data, axis=1, keepdims=True)
  norms[norms == 0.0] = 1.0
  return np.clip((data / norms) @ (data / norms).T, -1.0, 1.0)


def require_suffix(path: str, suffix: str) -> pathlib.Path:
  """Validates the --output extension and returns it as a Path.

  np.save() silently appends '.npy' to a path that lacks it, which would break
  the sidecar naming convention, so the extension is enforced up front.
  """
  out = pathlib.Path(path).expanduser()
  if out.suffix != suffix:
    raise BiohubError(f'--output must end in {suffix}; got {out.name!r}.')
  out.parent.mkdir(parents=True, exist_ok=True)
  return out


def save_matrix(array: np.ndarray, path: pathlib.Path) -> None:
  """Writes a .npy array and prints a one-line status."""
  np.save(path, array)
  print(f'Success! Array {tuple(array.shape)} written to: {path}')


def read_labels(path: str, ids: list[str]) -> list[str]:
  """Reads an id,label CSV (header optional) and aligns it to `ids`."""
  mapping: dict[str, str] = {}
  with open(path, encoding='utf-8') as handle:
    for row in csv.reader(handle):
      if len(row) < 2 or not row[0].strip():
        continue
      key, value = row[0].strip(), row[1].strip()
      if key.lower() in ('id', 'name', 'sequence_id') and not mapping:
        continue  # Header row.
      mapping[key] = value
  missing = [i for i in ids if i not in mapping]
  if missing:
    raise BiohubError(
        f'--labels {path} is missing {len(missing)} of {len(ids)} FASTA ids, '
        f'e.g. {missing[:3]}. Expected a CSV of "id,label" rows.'
    )
  return [mapping[i] for i in ids]


# --------------------------------------------------------------------------
# Subcommands
# --------------------------------------------------------------------------


def cmd_embed(args) -> None:
  """Embeds a single sequence, mean-pooled or per-residue."""
  records = load_records(args)
  if len(records) != 1:
    raise BiohubError(
        f'`embed` takes exactly one sequence but {args.fasta} holds '
        f'{len(records)} records. Use the `batch` subcommand instead.'
    )
  name, sequence = records[0]
  layer = resolve_layer(args.model, args.layer)

  client = BiohubClient()
  array = client.embed(
      sequence, args.model, layer=layer, per_residue=args.per_residue
  )

  out = require_suffix(args.output, '.npy')
  save_matrix(array, out)
  eb.write_json(
      {
          'sequence_id': name,
          'sequence_length': len(sequence),
          'model': args.model,
          'layer_requested': args.layer,
          'layer': layer,
          'layer_meaning': (
              'final output embedding (post-norm), not a hidden state'
              if layer is None
              else f'hidden-state row {layer} of '
              f'{eb.ESMC_MODELS[args.model]["n_layers"] + 1} '
              '(row 0 = embedding layer)'
          ),
          'pooling': 'none (per-residue)' if args.per_residue else 'mean',
          'per_residue': bool(args.per_residue),
          'shape': list(array.shape),
          'embedding_dim': int(array.shape[-1]),
          'dtype': str(array.dtype),
          'note': 'BOS/EOS are trimmed; row i is residue i+1 of the sequence.',
      },
      str(out.with_suffix('.json')),
  )


def cmd_batch(args) -> None:
  """Embeds every record in a FASTA into one (N, D) matrix."""
  records = eb.read_fasta(args.fasta)
  records = [(name, eb.validate_sequence(seq)) for name, seq in records]
  layer = resolve_layer(args.model, args.layer)

  client = BiohubClient()
  matrix = embed_many(client, records, args.model, layer, args.max_workers)

  out = require_suffix(args.output, '.npy')
  save_matrix(matrix, out)
  eb.write_json(
      {
          'ids': [name for name, _ in records],
          'lengths': [len(seq) for _, seq in records],
          'model': args.model,
          'layer_requested': args.layer,
          'layer': layer,
          'pooling': 'mean',
          'shape': list(matrix.shape),
          'embedding_dim': int(matrix.shape[1]),
          'note': 'Row i of the .npy is the embedding of ids[i].',
      },
      str(out.with_suffix('.json')),
  )


def cmd_similarity(args) -> None:
  """Writes the pairwise cosine-similarity matrix over a FASTA."""
  records = eb.read_fasta(args.fasta)
  records = [(name, eb.validate_sequence(seq)) for name, seq in records]
  if len(records) < 2:
    raise BiohubError('`similarity` needs at least 2 FASTA records.')
  layer = resolve_layer(args.model, args.layer)

  client = BiohubClient()
  matrix = embed_many(client, records, args.model, layer, args.max_workers)
  ids = [name for name, _ in records]
  similarity = cosine_matrix(matrix)

  out = require_suffix(args.output, '.csv')
  with open(out, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.writer(handle)
    writer.writerow(['id', *ids])
    for i, name in enumerate(ids):
      writer.writerow([name, *[f'{v:.6f}' for v in similarity[i]]])
  print(f'Success! Cosine matrix ({len(ids)}x{len(ids)}) written to: {out}')

  upper = [
      {'a': ids[i], 'b': ids[j], 'cosine': float(similarity[i, j])}
      for i in range(len(ids))
      for j in range(i + 1, len(ids))
  ]
  ranked = sorted(upper, key=lambda p: p['cosine'], reverse=True)
  off_diagonal = np.array([p['cosine'] for p in upper])
  eb.write_json(
      {
          'ids': ids,
          'model': args.model,
          'layer_requested': args.layer,
          'layer': layer,
          'n_sequences': len(ids),
          'n_pairs': len(upper),
          'mean_cosine': float(off_diagonal.mean()),
          'min_cosine': float(off_diagonal.min()),
          'max_cosine': float(off_diagonal.max()),
          'top_pairs': ranked[: args.top],
          'bottom_pairs': ranked[-args.top:][::-1],
      },
      str(out.with_suffix('.json')),
  )

  print(f'\nTop {min(args.top, len(ranked))} most similar pairs:')
  for pair in ranked[: args.top]:
    print(f'  {pair["cosine"]:+.4f}  {pair["a"]} <-> {pair["b"]}')


def cmd_cluster(args) -> None:
  """PCA + KMeans over a FASTA; writes assignments, a summary and a scatter."""
  from sklearn.cluster import KMeans  # pylint: disable=g-import-not-at-top
  from sklearn.decomposition import PCA  # pylint: disable=g-import-not-at-top
  from sklearn.metrics import adjusted_rand_score, silhouette_score  # pylint: disable=g-import-not-at-top

  records = eb.read_fasta(args.fasta)
  records = [(name, eb.validate_sequence(seq)) for name, seq in records]
  ids = [name for name, _ in records]
  if len(records) <= args.clusters:
    raise BiohubError(
        f'--clusters {args.clusters} needs more than {args.clusters} '
        f'sequences; {args.fasta} holds {len(records)}.'
    )
  layer = resolve_layer(args.model, args.layer)

  client = BiohubClient()
  matrix = embed_many(client, records, args.model, layer, args.max_workers)

  components = min(args.pca_components, len(records), matrix.shape[1])
  pca = PCA(n_components=components, random_state=0)
  projected = pca.fit_transform(matrix)
  kmeans = KMeans(n_clusters=args.clusters, n_init=10, random_state=0)
  assignments = kmeans.fit_predict(projected)

  truth = read_labels(args.labels, ids) if args.labels else None
  rand_index = (
      float(adjusted_rand_score(truth, assignments)) if truth else None
  )
  silhouette = (
      float(silhouette_score(projected, assignments))
      if len(set(assignments)) > 1
      else None
  )

  out = require_suffix(args.output, '.csv')
  with open(out, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.writer(handle)
    header = ['id', 'cluster', 'pc1', 'pc2'] + (['label'] if truth else [])
    writer.writerow(header)
    for i, name in enumerate(ids):
      pc2 = projected[i, 1] if components > 1 else 0.0
      row = [name, int(assignments[i]), f'{projected[i, 0]:.6f}', f'{pc2:.6f}']
      if truth:
        row.append(truth[i])
      writer.writerow(row)
  print(f'Success! Cluster assignments written to: {out}')

  plot_path = out.with_suffix('.png')
  _plot_clusters(projected, assignments, truth, pca, rand_index, plot_path)
  print(f'Success! PCA scatter written to: {plot_path}')

  sizes = {int(k): int(v) for k, v in zip(*np.unique(assignments, return_counts=True))}
  eb.write_json(
      {
          'ids': ids,
          'model': args.model,
          'layer_requested': args.layer,
          'layer': layer,
          'n_sequences': len(ids),
          'n_clusters': args.clusters,
          'pca_components': components,
          'explained_variance_ratio': [
              float(v) for v in pca.explained_variance_ratio_
          ],
          'cluster_sizes': sizes,
          'silhouette': silhouette,
          'adjusted_rand_index': rand_index,
          'assignments': [
              {
                  'id': name,
                  'cluster': int(assignments[i]),
                  **({'label': truth[i]} if truth else {}),
              }
              for i, name in enumerate(ids)
          ],
          'plot': str(plot_path),
      },
      str(out.with_suffix('.json')),
  )

  print(f'\nClusters (k={args.clusters}): {sizes}')
  if silhouette is not None:
    print(f'Silhouette:          {silhouette:+.3f}')
  if rand_index is not None:
    print(f'Adjusted Rand index: {rand_index:.3f} (vs --labels)')


def _plot_clusters(projected, assignments, truth, pca, rand_index, path) -> None:
  """Draws the PCA scatter: colour = KMeans cluster, marker = true label."""
  figure, axes = plt.subplots(figsize=(5.5, 5))
  x = projected[:, 0]
  y = projected[:, 1] if projected.shape[1] > 1 else np.zeros_like(x)
  markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']

  groups = sorted(set(truth)) if truth else [None]
  for g_index, group in enumerate(groups):
    mask = (
        np.array([t == group for t in truth])
        if truth
        else np.ones(len(x), dtype=bool)
    )
    axes.scatter(
        x[mask],
        y[mask],
        c=assignments[mask],
        cmap='tab10',
        vmin=0,
        vmax=9,
        marker=markers[g_index % len(markers)],
        s=90,
        edgecolors='black',
        linewidths=0.5,
        label=group,
    )

  variance = pca.explained_variance_ratio_
  axes.set_xlabel(f'PC1 ({variance[0] * 100:.1f}% var)')
  axes.set_ylabel(
      f'PC2 ({variance[1] * 100:.1f}% var)' if len(variance) > 1 else 'PC2'
  )
  title = 'PCA of ESMC mean embeddings\ncolour = KMeans cluster'
  if truth:
    title += ', marker = label'
  if rand_index is not None:
    title += f'\nadjusted Rand index = {rand_index:.2f}'
  axes.set_title(title, fontsize=10)
  if truth:
    axes.legend(title='label', fontsize=8, loc='best')
  figure.tight_layout()
  figure.savefig(path, dpi=150)
  plt.close(figure)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _add_common(parser, *, fasta_only: bool) -> None:
  """Adds the model/layer/output flags shared by every subcommand."""
  if fasta_only:
    parser.add_argument('--fasta', required=True, help='Input FASTA file.')
  else:
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--sequence', help='A single amino-acid sequence.')
    source.add_argument('--fasta', help='FASTA holding exactly one record.')
    parser.add_argument('--id', help='Name for --sequence (default: query).')
  parser.add_argument(
      '--model',
      default=eb.DEFAULT_ESMC,
      choices=sorted(eb.ESMC_MODELS),
      help=f'ESMC model (default: {eb.DEFAULT_ESMC}).',
  )
  parser.add_argument(
      '--layer',
      type=int,
      help=(
          'Hidden layer to pool. 0 = embedding layer, -2 = second-to-last '
          '(a good blind default; intermediate layers often beat the last). '
          'Omit to use the final output embedding.'
      ),
  )
  parser.add_argument('--output', required=True, help='Output file path.')


def _add_workers(parser) -> None:
  parser.add_argument(
      '--max-workers',
      type=int,
      default=4,
      help=f'Parallel requests (capped at {MAX_WORKERS_CAP}).',
  )


def main() -> None:
  parser = argparse.ArgumentParser(
      description='Embed proteins with ESMC and analyse the embeddings.'
  )
  sub = parser.add_subparsers(dest='command', required=True)

  p_embed = sub.add_parser('embed', help='Embed one sequence.')
  _add_common(p_embed, fasta_only=False)
  p_embed.add_argument(
      '--per-residue',
      action='store_true',
      help='Write an (L, D) per-residue matrix instead of a (D,) mean vector.',
  )
  p_embed.set_defaults(func=cmd_embed)

  p_batch = sub.add_parser('batch', help='Embed a FASTA into an (N, D) matrix.')
  _add_common(p_batch, fasta_only=True)
  _add_workers(p_batch)
  p_batch.set_defaults(func=cmd_batch)

  p_sim = sub.add_parser('similarity', help='Pairwise cosine over a FASTA.')
  _add_common(p_sim, fasta_only=True)
  _add_workers(p_sim)
  p_sim.add_argument(
      '--top',
      type=int,
      required=True,
      help='How many top/bottom pairs to report (the full matrix is written).',
  )
  p_sim.set_defaults(func=cmd_similarity)

  p_cluster = sub.add_parser('cluster', help='PCA + KMeans over a FASTA.')
  _add_common(p_cluster, fasta_only=True)
  _add_workers(p_cluster)
  p_cluster.add_argument(
      '--clusters', type=int, required=True, help='Number of KMeans clusters.'
  )
  p_cluster.add_argument(
      '--pca-components',
      type=int,
      default=2,
      help='PCA dimensions to cluster in (default: 2, as in the ESM tutorial).',
  )
  p_cluster.add_argument(
      '--labels',
      help='Optional "id,label" CSV of ground truth -> adjusted Rand index.',
  )
  p_cluster.set_defaults(func=cmd_cluster)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
