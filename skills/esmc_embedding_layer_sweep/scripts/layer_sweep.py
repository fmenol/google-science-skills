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

"""Chooses the optimal ESMC layer for a downstream supervised task.

Trains a linear probe on every layer's mean-pooled embedding with stratified
cross-validation and reports which layer carries the most task-relevant,
linearly-decodable signal. The last layer is specialised to the masked-LM
objective and is frequently NOT the best choice.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
#   "numpy",
#   "scikit-learn",
#   "matplotlib",
#   "joblib",
# ]
# ///

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import json
import pathlib
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import log_loss
from sklearn.metrics import matthews_corrcoef
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from esm_biohub import BiohubClient, BiohubError  # vendored, same directory
import esm_biohub as eb

# The probe is deliberately linear and strongly regularised. A linear probe
# answers "is this information present and linearly decodable?", which is the
# question a layer sweep exists to answer. A stronger head would confound layer
# quality with head capacity.
DEFAULT_C = 0.1
DEFAULT_MAX_ITER = 1000
DEFAULT_RANDOM_STATE = 44


# --------------------------------------------------------------------------
# Dataset cache (.npy + .meta.json sidecar)
# --------------------------------------------------------------------------


def _npy_path(output: str) -> pathlib.Path:
  """Normalises an output path to end in .npy."""
  path = pathlib.Path(output).expanduser()
  if path.suffix != '.npy':
    path = path.with_name(path.name + '.npy')
  return path


def _meta_path(npy: pathlib.Path) -> pathlib.Path:
  """Sidecar path holding ids, labels and the model used to embed them."""
  return npy.with_suffix('.meta.json')


def read_labels(
    path: str, id_column: str, label_column: str
) -> dict[str, str]:
  """Reads a labels CSV into {sequence_id: label}."""
  with open(path, newline='', encoding='utf-8') as handle:
    reader = csv.DictReader(handle)
    columns = reader.fieldnames or []
    for needed in (id_column, label_column):
      if needed not in columns:
        raise BiohubError(
            f'Column {needed!r} not found in {path}. Available columns: '
            f'{columns}. Use --id-column / --label-column to point at the '
            'right ones.'
        )
    labels: dict[str, str] = {}
    for row in reader:
      key = (row[id_column] or '').strip()
      value = (row[label_column] or '').strip()
      if not key or not value:
        continue
      labels[key] = value
  if not labels:
    raise BiohubError(f'No usable rows in {path}.')
  return labels


def load_dataset(npy: str) -> tuple[np.ndarray, dict]:
  """Loads a cached embedding array plus its metadata sidecar."""
  path = _npy_path(npy)
  if not path.is_file():
    raise BiohubError(
        f'Embedding cache not found: {path}. Run `embed-dataset` first.'
    )
  meta_file = _meta_path(path)
  if not meta_file.is_file():
    raise BiohubError(
        f'Metadata sidecar not found: {meta_file}. It is written next to the '
        '.npy by `embed-dataset` and carries the labels; do not move one '
        'without the other.'
    )
  embeddings = np.load(path)
  with open(meta_file, encoding='utf-8') as handle:
    meta = json.load(handle)
  if embeddings.ndim != 3:
    raise BiohubError(
        f'Expected a 3-D (n_sequences, n_layer_rows, dim) array, got shape '
        f'{embeddings.shape}.'
    )
  if embeddings.shape[0] != len(meta.get('labels', [])):
    raise BiohubError(
        f'{path} has {embeddings.shape[0]} rows but {meta_file} lists '
        f'{len(meta.get("labels", []))} labels. The cache is inconsistent; '
        're-run `embed-dataset`.'
    )
  return embeddings, meta


def encode_labels(meta: dict) -> tuple[np.ndarray, list[str]]:
  """Maps string labels to integer class ids with a stable class order."""
  labels = list(meta['labels'])
  classes = sorted(set(labels))
  if len(classes) < 2:
    raise BiohubError(
        f'A probe needs at least 2 classes; the dataset has {len(classes)} '
        f'({classes}).'
    )
  return np.array([classes.index(x) for x in labels]), classes


def check_n_splits(n_splits: int, labels: list[str]) -> None:
  """Validates the fold count against the rarest class.

  StratifiedKFold needs at least one member of every class in every fold, so
  `n_splits` can never exceed the smallest class count.
  """
  if n_splits < 2:
    raise BiohubError(f'--n-splits must be at least 2, got {n_splits}.')
  counts = collections.Counter(labels)
  rarest, smallest = min(counts.items(), key=lambda kv: kv[1])
  if n_splits > smallest:
    raise BiohubError(
        f'--n-splits {n_splits} exceeds the smallest class count '
        f'({smallest}, class {rarest!r}). Stratified cross-validation needs at '
        f'least one member of every class in each fold. Use --n-splits '
        f'{smallest} or fewer, or add more {rarest!r} examples. '
        f'Class counts: {dict(sorted(counts.items()))}'
    )
  if len(labels) < 20:
    print(
        f'[warning] Only {len(labels)} sequences. Cross-validated scores will '
        'be noisy and may saturate at MCC = 1.0 if the classes are easy; '
        'treat the winning layer as indicative, not definitive.',
        file=sys.stderr,
    )


def make_probe(c: float, max_iter: int, random_state: int) -> Pipeline:
  """StandardScaler + L2 logistic regression, exactly as in the tutorial."""
  return Pipeline([
      ('scaler', StandardScaler()),
      (
          'clf',
          LogisticRegression(
              C=c,
              max_iter=max_iter,
              solver='lbfgs',
              random_state=random_state,
          ),
      ),
  ])


def select_best_layer(
    mcc_mean: np.ndarray, loss_mean: np.ndarray
) -> tuple[int, list[int]]:
  """Picks the winning layer, breaking MCC ties with cross-validated log-loss.

  MCC is bounded above by 1.0. On an easy or small dataset many layers reach
  that ceiling at once, and `argmax` would then silently hand back the lowest
  tied index -- which is layer 0, the embedding layer. That is the worst
  possible advice a layer sweep can give: layer 0 has seen no attention at all
  and ties only because the task has no headroom left.

  Log-loss is a strictly proper scoring rule that keeps improving after the
  classification is already perfect (it rewards a confident, well-separated
  decision boundary), so it still discriminates when MCC cannot.

  Returns:
    (best layer index, the list of layers tied at the maximum MCC)
  """
  top = float(np.max(mcc_mean))
  tied = [int(i) for i in np.flatnonzero(mcc_mean >= top - 1e-9)]
  if len(tied) == 1:
    return tied[0], tied
  best = int(min(tied, key=lambda i: float(loss_mean[i])))
  return best, tied


def resolve_layer(layer: int, n_rows: int) -> int:
  """Resolves a possibly-negative layer index against the layer-row count.

  Row 0 is the embedding layer; row `n_rows - 1` is the final transformer
  block. Python-style negative indices are accepted so that the `[-2]` rule of
  thumb can be written literally.
  """
  index = layer + n_rows if layer < 0 else layer
  if not 0 <= index < n_rows:
    raise BiohubError(
        f'--layer {layer} is out of range. This model has {n_rows} layer rows, '
        f'so valid indices are 0..{n_rows - 1} (or -1..-{n_rows}). '
        'Row 0 is the embedding layer.'
    )
  return index


def layer_name(index: int, n_rows: int) -> str:
  """Human-readable name for a layer row."""
  if index == 0:
    return 'embedding (no attention)'
  if index == n_rows - 1:
    return f'block {index} (final)'
  return f'block {index}'


# --------------------------------------------------------------------------
# embed-dataset
# --------------------------------------------------------------------------


def cmd_embed_dataset(args) -> None:
  """Embeds a labelled FASTA once, caching every layer for offline sweeps."""
  records = eb.read_fasta(args.fasta)
  labels = read_labels(args.labels, args.id_column, args.label_column)

  pairs: list[tuple[str, str, str]] = []
  unlabelled: list[str] = []
  for name, sequence in records:
    if name in labels:
      pairs.append((name, eb.validate_sequence(sequence), labels[name]))
    else:
      unlabelled.append(name)
  if unlabelled:
    print(
        f'[warning] {len(unlabelled)} FASTA record(s) had no label and were '
        f'skipped: {unlabelled[:5]}{"..." if len(unlabelled) > 5 else ""}',
        file=sys.stderr,
    )
  if not pairs:
    raise BiohubError(
        'No FASTA record ids matched the labels CSV. Check that the FASTA '
        f'header ids match the {args.id_column!r} column.'
    )

  spec = eb.ESMC_MODELS[args.model]
  expected_rows = spec['n_layers'] + 1
  print(
      f'Embedding {len(pairs)} sequences with {args.model} '
      f'({expected_rows} layer rows x {spec["dim"]} dims). '
      'All layers arrive in one request per sequence.',
      file=sys.stderr,
  )

  client = BiohubClient()
  results: list[np.ndarray | None] = [None] * len(pairs)
  failures: list[tuple[str, str]] = []
  completed = 0
  with concurrent.futures.ThreadPoolExecutor(
      max_workers=args.max_workers
  ) as pool:
    futures = {
        pool.submit(client.mean_hidden_states, sequence, args.model): index
        for index, (_, sequence, _) in enumerate(pairs)
    }
    for future in concurrent.futures.as_completed(futures):
      index = futures[future]
      try:
        results[index] = future.result()
      except BiohubError as exc:
        failures.append((pairs[index][0], str(exc)))
      completed += 1
      if completed % 10 == 0 or completed == len(pairs):
        print(f'  {completed}/{len(pairs)} embedded', file=sys.stderr,
              flush=True)

  for name, message in failures:
    print(f'[warning] dropped {name}: {message}', file=sys.stderr)
  kept = [(pair, arr) for pair, arr in zip(pairs, results) if arr is not None]
  if not kept:
    raise BiohubError('Every sequence failed to embed. See the errors above.')

  embeddings = np.stack([arr for _, arr in kept]).astype(np.float32)
  if embeddings.shape[1:] != (expected_rows, spec['dim']):
    raise BiohubError(
        f'{args.model} returned layer stacks of shape {embeddings.shape[1:]}, '
        f'expected ({expected_rows}, {spec["dim"]}). The model registry and '
        'the API disagree; do not trust the sweep.'
    )

  output = _npy_path(args.output)
  output.parent.mkdir(parents=True, exist_ok=True)
  np.save(output, embeddings)

  kept_labels = [pair[2] for pair, _ in kept]
  counts = collections.Counter(kept_labels)
  meta = {
      'model': args.model,
      'n_transformer_layers': spec['n_layers'],
      'n_layer_rows': int(embeddings.shape[1]),
      'hidden_dim': int(embeddings.shape[2]),
      'n_sequences': int(embeddings.shape[0]),
      'layer_index_convention': (
          'Row 0 is the embedding layer (token embeddings, no attention). '
          f'Row k is the output of transformer block k. Row '
          f'{embeddings.shape[1] - 1} is the final layer.'
      ),
      'pooling': (
          'Server-side mean over all tokens, BOS and EOS included '
          '(return_mean_hidden_states=True).'
      ),
      'embeddings_file': output.name,
      'ids': [pair[0] for pair, _ in kept],
      'labels': kept_labels,
      'class_counts': dict(sorted(counts.items())),
      'dropped': [name for name, _ in failures],
  }
  eb.write_json(meta, str(_meta_path(output)))
  print(
      f'Success! Cached embeddings {tuple(embeddings.shape)} '
      f'(n_sequences, n_layer_rows, dim) to: {output}'
  )
  print(f'  classes: {dict(sorted(counts.items()))}')


# --------------------------------------------------------------------------
# sweep
# --------------------------------------------------------------------------


def cmd_sweep(args) -> None:
  """Trains one linear probe per layer and ranks the layers by CV MCC."""
  embeddings, meta = load_dataset(args.embeddings)
  y, classes = encode_labels(meta)
  check_n_splits(args.n_splits, list(meta['labels']))

  n_rows = embeddings.shape[1]
  splitter = StratifiedKFold(
      n_splits=args.n_splits, shuffle=True, random_state=args.random_state
  )
  folds = list(splitter.split(embeddings[:, 0, :], y))

  class_ids = list(range(len(classes)))
  mcc = np.zeros((n_rows, args.n_splits))
  acc = np.zeros((n_rows, args.n_splits))
  loss = np.zeros((n_rows, args.n_splits))
  print(
      f'Sweeping {n_rows} layers x {args.n_splits} folds '
      f'({n_rows * args.n_splits} probe fits)...',
      file=sys.stderr,
  )
  for row in range(n_rows):
    features = embeddings[:, row, :]
    for fold, (train, test) in enumerate(folds):
      probe = make_probe(args.C, args.max_iter, args.random_state)
      probe.fit(features[train], y[train])
      predicted = probe.predict(features[test])
      mcc[row, fold] = matthews_corrcoef(y[test], predicted)
      acc[row, fold] = accuracy_score(y[test], predicted)
      loss[row, fold] = log_loss(
          y[test], probe.predict_proba(features[test]), labels=class_ids
      )

  mcc_mean, mcc_std = mcc.mean(axis=1), mcc.std(axis=1)
  acc_mean, acc_std = acc.mean(axis=1), acc.std(axis=1)
  loss_mean, loss_std = loss.mean(axis=1), loss.std(axis=1)
  if not np.all(np.isfinite(mcc_mean)):
    raise BiohubError('Non-finite MCC produced by the probe; refusing to emit.')

  best, tied = select_best_layer(mcc_mean, loss_mean)
  saturated = len(tied) > 1
  last = n_rows - 1
  rows = [
      {
          'layer': index,
          'layer_name': layer_name(index, n_rows),
          'mcc_mean': round(float(mcc_mean[index]), 6),
          'mcc_std': round(float(mcc_std[index]), 6),
          'accuracy_mean': round(float(acc_mean[index]), 6),
          'accuracy_std': round(float(acc_std[index]), 6),
          'log_loss_mean': round(float(loss_mean[index]), 6),
          'log_loss_std': round(float(loss_std[index]), 6),
          'is_best': index == best,
          'is_last': index == last,
      }
      for index in range(n_rows)
  ]

  ranked = sorted(rows, key=lambda r: r['mcc_mean'], reverse=True)
  summary = {
      'model': meta['model'],
      'n_sequences': int(embeddings.shape[0]),
      'hidden_dim': int(embeddings.shape[2]),
      'n_layer_rows': n_rows,
      'layer_index_convention': meta['layer_index_convention'],
      'classes': classes,
      'class_counts': meta.get('class_counts', {}),
      'probe': {
          'estimator': 'StandardScaler + LogisticRegression',
          'C': args.C,
          'max_iter': args.max_iter,
          'solver': 'lbfgs',
      },
      'cv': {
          'splitter': 'StratifiedKFold',
          'n_splits': args.n_splits,
          'shuffle': True,
          'random_state': args.random_state,
      },
      'metric': 'matthews_corrcoef',
      'tie_break_metric': 'log_loss',
      'best_layer': best,
      'best_mcc_mean': round(float(mcc_mean[best]), 6),
      'best_mcc_std': round(float(mcc_std[best]), 6),
      'best_log_loss_mean': round(float(loss_mean[best]), 6),
      'last_layer': last,
      'last_mcc_mean': round(float(mcc_mean[last]), 6),
      'last_mcc_std': round(float(mcc_std[last]), 6),
      'last_log_loss_mean': round(float(loss_mean[last]), 6),
      'best_minus_last_mcc': round(float(mcc_mean[best] - mcc_mean[last]), 6),
      'embedding_layer_mcc_mean': round(float(mcc_mean[0]), 6),
      'embedding_layer_log_loss_mean': round(float(loss_mean[0]), 6),
      'best_layer_depth_fraction': round(best / last, 3) if last else 0.0,
      'mcc_saturated': saturated,
      'n_layers_tied_at_best_mcc': len(tied),
      'selection_metric': (
          'log_loss (MCC tied across layers)' if saturated else 'mcc'
      ),
      'top5_layers': [r['layer'] for r in ranked[:5]],
      'layers': rows,
  }
  eb.write_json(summary, args.output_json)

  csv_path = pathlib.Path(args.output_csv).expanduser()
  csv_path.parent.mkdir(parents=True, exist_ok=True)
  with open(csv_path, 'w', newline='', encoding='utf-8') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
  print(f'Success! Data written to: {csv_path}')

  print(
      f'\nBest layer:  {best:>3d} ({layer_name(best, n_rows)})  '
      f'MCC = {mcc_mean[best]:.3f} +/- {mcc_std[best]:.3f}, '
      f'log-loss = {loss_mean[best]:.4f}'
  )
  print(
      f'Last layer:  {last:>3d} ({layer_name(last, n_rows)})  '
      f'MCC = {mcc_mean[last]:.3f} +/- {mcc_std[last]:.3f}, '
      f'log-loss = {loss_mean[last]:.4f}'
  )
  print(
      f'Embed layer:   0 ({layer_name(0, n_rows)})  '
      f'MCC = {mcc_mean[0]:.3f} +/- {mcc_std[0]:.3f}, '
      f'log-loss = {loss_mean[0]:.4f}'
  )

  if saturated:
    print(
        f'\n[!] MCC SATURATED: {len(tied)} of {n_rows} layers tie at '
        f'MCC = {mcc_mean[best]:.3f}. The task has no headroom left, so MCC '
        'cannot rank these layers. The winner below was chosen by '
        'cross-validated log-loss, which keeps discriminating after the '
        'classification is already perfect.\n'
        '    Treat it as indicative, not definitive: the honest reading is '
        'that this dataset cannot separate those layers. Add harder or more '
        'numerous examples before committing to a layer.'
    )

  delta = float(mcc_mean[best] - mcc_mean[last])
  if best == last:
    print(
        f'\nThe final layer is already optimal for this task (--layer {last}).'
    )
  elif delta > 0:
    print(
        f'\nLayer {best} beats the final layer by {delta:.3f} MCC '
        f'({100 * best / last:.0f}% of the way through the network). '
        f'Use --layer {best} for this task.'
    )
  else:
    print(
        f'\nLayer {best} ties the final layer on MCC but has the lower '
        f'log-loss ({loss_mean[best]:.4f} vs {loss_mean[last]:.4f}) — a more '
        f'confident, better-separated boundary '
        f'({100 * best / last:.0f}% of the way through the network). '
        f'Use --layer {best} for this task.'
    )


# --------------------------------------------------------------------------
# plot
# --------------------------------------------------------------------------


def cmd_plot(args) -> None:
  """Plots CV MCC against layer depth with a +/-1 std band."""
  import matplotlib  # pylint: disable=g-import-not-at-top
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt  # pylint: disable=g-import-not-at-top

  with open(args.sweep, encoding='utf-8') as handle:
    summary = json.load(handle)

  rows = summary['layers']
  layers = np.array([row['layer'] for row in rows])
  mean = np.array([row['mcc_mean'] for row in rows])
  std = np.array([row['mcc_std'] for row in rows])
  loss = np.array([row['log_loss_mean'] for row in rows])
  loss_std = np.array([row['log_loss_std'] for row in rows])
  best = int(summary['best_layer'])
  last = int(summary['last_layer'])
  saturated = bool(summary.get('mcc_saturated'))

  fig, (ax, ax2) = plt.subplots(
      2, 1, figsize=(12, 8), sharex=True,
      gridspec_kw={'height_ratios': [3, 2]},
  )

  ax.plot(layers, mean, 'o-', linewidth=2, color='steelblue', label='MCC')
  ax.fill_between(
      layers, mean - std, mean + std, alpha=0.2, color='steelblue',
      label='+/-1 std',
  )
  ax.axvline(
      best, color='crimson', linestyle=':', linewidth=1.5,
      label=f'Best layer = {best} (MCC {mean[best]:.3f})',
  )
  ax.axvline(
      last, color='gray', linestyle=':', linewidth=1.5,
      label=f'Last layer = {last} (MCC {mean[last]:.3f})',
  )
  ax.set_ylabel('Cross-validated MCC', fontsize=12)
  classes = ', '.join(summary.get('classes', []))
  title = (
      f'Layer sweep: {summary["model"]}\n'
      f'mean-pooled embeddings + logistic regression, '
      f'{summary["cv"]["n_splits"]}-fold CV, n={summary["n_sequences"]} '
      f'({classes})'
  )
  if saturated:
    title += (
        f'\n[!] MCC saturated: {summary["n_layers_tied_at_best_mcc"]} layers '
        'tie at the ceiling — layer chosen by log-loss'
    )
  ax.set_title(title, fontsize=12)
  ax.grid(axis='y', alpha=0.3)
  ax.legend(fontsize=9, loc='lower right')

  # Log-loss keeps discriminating after MCC hits its ceiling, so it is the
  # panel to read whenever the MCC curve is flat.
  ax2.plot(layers, loss, 'o-', linewidth=2, color='darkorange',
           label='log-loss (lower is better)')
  ax2.fill_between(
      layers, loss - loss_std, loss + loss_std, alpha=0.2, color='darkorange',
      label='+/-1 std',
  )
  ax2.axvline(best, color='crimson', linestyle=':', linewidth=1.5,
              label=f'Best layer = {best} (log-loss {loss[best]:.4f})')
  ax2.axvline(last, color='gray', linestyle=':', linewidth=1.5)
  ax2.set_xlabel('ESMC layer (0 = embedding layer, no attention)', fontsize=12)
  ax2.set_ylabel('Cross-validated log-loss', fontsize=12)
  ax2.grid(axis='y', alpha=0.3)
  ax2.legend(fontsize=9, loc='upper right')
  ax2.set_xlim(-0.5, last + 0.5)
  fig.tight_layout()

  output = pathlib.Path(args.output).expanduser()
  output.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(output, dpi=150, bbox_inches='tight')
  plt.close(fig)
  print(f'Success! Figure written to: {output}')


# --------------------------------------------------------------------------
# probe
# --------------------------------------------------------------------------


def cmd_probe(args) -> None:
  """Fits, scores and saves a probe at one chosen layer."""
  import joblib  # pylint: disable=g-import-not-at-top

  embeddings, meta = load_dataset(args.embeddings)
  y, classes = encode_labels(meta)
  check_n_splits(args.n_splits, list(meta['labels']))

  n_rows = embeddings.shape[1]
  layer = resolve_layer(args.layer, n_rows)
  features = embeddings[:, layer, :]

  splitter = StratifiedKFold(
      n_splits=args.n_splits, shuffle=True, random_state=args.random_state
  )
  fold_mcc, fold_acc = [], []
  predicted_all = np.zeros_like(y)
  for train, test in splitter.split(features, y):
    probe = make_probe(args.C, args.max_iter, args.random_state)
    probe.fit(features[train], y[train])
    predicted = probe.predict(features[test])
    predicted_all[test] = predicted
    fold_mcc.append(float(matthews_corrcoef(y[test], predicted)))
    fold_acc.append(float(accuracy_score(y[test], predicted)))

  final = make_probe(args.C, args.max_iter, args.random_state)
  final.fit(features, y)
  bundle = {
      'pipeline': final,
      'classes': classes,
      'layer': layer,
      'model': meta['model'],
      'n_layer_rows': n_rows,
      'hidden_dim': int(embeddings.shape[2]),
      'pooling': meta['pooling'],
  }
  model_path = pathlib.Path(args.output_model).expanduser()
  model_path.parent.mkdir(parents=True, exist_ok=True)
  joblib.dump(bundle, model_path)
  print(f'Success! Probe written to: {model_path}')

  metrics = {
      'model': meta['model'],
      'layer': layer,
      'layer_name': layer_name(layer, n_rows),
      'n_layer_rows': n_rows,
      'classes': classes,
      'class_counts': meta.get('class_counts', {}),
      'n_sequences': int(embeddings.shape[0]),
      'probe': {'estimator': 'StandardScaler + LogisticRegression',
                'C': args.C, 'max_iter': args.max_iter},
      'cv': {'splitter': 'StratifiedKFold', 'n_splits': args.n_splits,
             'random_state': args.random_state},
      'mcc_mean': round(float(np.mean(fold_mcc)), 6),
      'mcc_std': round(float(np.std(fold_mcc)), 6),
      'mcc_per_fold': [round(v, 6) for v in fold_mcc],
      'accuracy_mean': round(float(np.mean(fold_acc)), 6),
      'accuracy_std': round(float(np.std(fold_acc)), 6),
      'confusion_matrix': confusion_matrix(
          y, predicted_all, labels=list(range(len(classes)))
      ).tolist(),
      'confusion_matrix_axes': 'rows = true class, cols = predicted class',
      'probe_file': str(model_path),
  }

  if args.predict_fasta:
    if not args.output_predictions:
      raise BiohubError(
          '--predict-fasta requires --output-predictions.'
      )
    records = eb.read_fasta(args.predict_fasta)
    client = BiohubClient()
    print(
        f'Embedding {len(records)} query sequences at layer {layer}...',
        file=sys.stderr,
    )
    # NOTE: this must use mean_hidden_states (server-side pooling that INCLUDES
    # BOS/EOS), the same call `embed-dataset` used. client.embed(layer=...)
    # pools a DIFFERENT way (residues only) and would silently shift the
    # features away from what the probe was trained on.
    queries = np.stack([
        client.mean_hidden_states(eb.validate_sequence(seq), meta['model'])[
            layer
        ]
        for _, seq in records
    ]).astype(np.float32)
    labels_out = final.predict(queries)
    proba = final.predict_proba(queries)
    predictions = [
        {
            'id': name,
            'predicted_label': classes[int(label)],
            'confidence': round(float(proba[i].max()), 6),
            **{
                f'p_{cls}': round(float(proba[i][j]), 6)
                for j, cls in enumerate(classes)
            },
        }
        for i, ((name, _), label) in enumerate(zip(records, labels_out))
    ]
    csv_path = pathlib.Path(args.output_predictions).expanduser()
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, 'w', newline='', encoding='utf-8') as handle:
      writer = csv.DictWriter(handle, fieldnames=list(predictions[0].keys()))
      writer.writeheader()
      writer.writerows(predictions)
    print(f'Success! Data written to: {csv_path}')
    metrics['predictions_file'] = str(csv_path)
    metrics['n_predicted'] = len(predictions)

  eb.write_json(metrics, args.output_json)
  print(
      f'\nLayer {layer} ({layer_name(layer, n_rows)}):  '
      f'MCC = {metrics["mcc_mean"]:.3f} +/- {metrics["mcc_std"]:.3f}, '
      f'accuracy = {metrics["accuracy_mean"]:.3f}'
  )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main() -> None:
  parser = argparse.ArgumentParser(
      description='Pick the best ESMC layer for a supervised protein task.'
  )
  sub = parser.add_subparsers(dest='command', required=True)

  p = sub.add_parser(
      'embed-dataset',
      help='Embed a labelled FASTA once, caching every layer to .npy.',
  )
  p.add_argument('--fasta', required=True, help='Input FASTA file.')
  p.add_argument(
      '--labels', required=True, help='CSV with one row per sequence id.'
  )
  p.add_argument(
      '--model',
      required=True,
      choices=sorted(eb.ESMC_MODELS),
      help='ESMC model. Determines the layer count and dimension.',
  )
  p.add_argument(
      '--output', required=True, help='Output .npy path for the layer stack.'
  )
  p.add_argument('--id-column', default='id', help='Default: id')
  p.add_argument('--label-column', default='label', help='Default: label')
  p.add_argument(
      '--max-workers', type=int, default=8, help='Concurrent requests (<=8).'
  )
  p.set_defaults(func=cmd_embed_dataset)

  p = sub.add_parser('sweep', help='Train one probe per layer; rank layers.')
  p.add_argument('--embeddings', required=True, help='Cached .npy.')
  p.add_argument(
      '--n-splits',
      type=int,
      required=True,
      help='Stratified CV folds. Must be <= the smallest class count.',
  )
  p.add_argument('--output-json', required=True, help='Sweep summary JSON.')
  p.add_argument('--output-csv', required=True, help='Per-layer scores CSV.')
  p.add_argument('--C', type=float, default=DEFAULT_C, help='Inverse L2.')
  p.add_argument('--max-iter', type=int, default=DEFAULT_MAX_ITER)
  p.add_argument('--random-state', type=int, default=DEFAULT_RANDOM_STATE)
  p.set_defaults(func=cmd_sweep)

  p = sub.add_parser('plot', help='Plot MCC vs layer with a std band.')
  p.add_argument('--sweep', required=True, help='Sweep summary JSON.')
  p.add_argument('--output', required=True, help='Output .png path.')
  p.set_defaults(func=cmd_plot)

  p = sub.add_parser('probe', help='Fit, score and save a probe at one layer.')
  p.add_argument('--embeddings', required=True, help='Cached .npy.')
  p.add_argument(
      '--layer',
      type=int,
      required=True,
      help='Layer row. 0 = embedding layer; negatives allowed (-1 = last).',
  )
  p.add_argument('--n-splits', type=int, required=True, help='Stratified folds.')
  p.add_argument('--output-model', required=True, help='Output .joblib path.')
  p.add_argument('--output-json', required=True, help='Output metrics JSON.')
  p.add_argument('--predict-fasta', help='Optional FASTA of new sequences.')
  p.add_argument('--output-predictions', help='CSV for --predict-fasta.')
  p.add_argument('--C', type=float, default=DEFAULT_C, help='Inverse L2.')
  p.add_argument('--max-iter', type=int, default=DEFAULT_MAX_ITER)
  p.add_argument('--random-state', type=int, default=DEFAULT_RANDOM_STATE)
  p.set_defaults(func=cmd_probe)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
