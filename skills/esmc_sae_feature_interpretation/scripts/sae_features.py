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

"""Sparse-autoencoder (SAE) feature interpretation for ESMC representations.

Decomposes an ESMC residue embedding into a sparse, human-interpretable feature
basis (16,384 features, top-k=64 active per residue), ranks the features that
fire on a protein, and joins them to natural-language descriptions.

The feature descriptions are AUTO-GENERATED HYPOTHESES, not curated annotations.
Never present them as established fact.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
#   "numpy",
#   "matplotlib",
# ]
# ///

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any

import numpy as np
from polite_http import http_client
from polite_http.http_client import HttpError

import esm_biohub as eb
from esm_biohub import BiohubError

# The feature-description service is a different host path from the inference
# API and needs no auth. Keep our own polite client for it.
FEATURE_API_QPS = 5.0
DEFAULT_CACHE_DIR = pathlib.Path.home() / '.cache' / 'esm_sae_features'

# Fields worth surfacing from the feature-description API.
DESCRIPTION_FIELDS = (
    'label',
    'summary',
    'description',
    'category',
    'activation_pattern',
    'exemplar_protein_families',
    'uniref90_frequency',
    'threshold',
)

# A feature is counted as "on" at a residue when its activation exceeds this.
# Matches the ESM cookbook tutorial. Override with --threshold.
DEFAULT_ACTIVATION_THRESHOLD = 0.01


# --------------------------------------------------------------------------
# Sparse feature helpers
#
# The SAE is top-k: each residue has exactly k active features out of `codebook`
# (k=64, codebook=16384). We never build the (L, 16384) dense matrix -- it is
# 65 MB for a 1000-residue protein and entirely unnecessary. All statistics are
# computed by scatter ops over the (L, k) sparse arrays instead.
# --------------------------------------------------------------------------


def max_activation(
    indices: np.ndarray, values: np.ndarray, codebook: int
) -> np.ndarray:
  """Peak activation of every feature across the sequence. Shape (codebook,)."""
  out = np.zeros(codebook, dtype=np.float64)
  np.maximum.at(out, indices.ravel(), values.ravel())
  return out


def prevalence(
    indices: np.ndarray, values: np.ndarray, codebook: int, threshold: float
) -> np.ndarray:
  """Number of residues where each feature exceeds `threshold`. (codebook,)."""
  out = np.zeros(codebook, dtype=np.int64)
  flat_idx = indices.ravel()
  flat_val = values.ravel()
  keep = flat_val > threshold
  np.add.at(out, flat_idx[keep], 1)
  return out


def mean_when_on(
    indices: np.ndarray, values: np.ndarray, codebook: int, threshold: float
) -> np.ndarray:
  """Mean activation over the residues where a feature is on. (codebook,)."""
  total = np.zeros(codebook, dtype=np.float64)
  count = np.zeros(codebook, dtype=np.int64)
  flat_idx = indices.ravel()
  flat_val = values.ravel()
  keep = flat_val > threshold
  np.add.at(total, flat_idx[keep], flat_val[keep])
  np.add.at(count, flat_idx[keep], 1)
  with np.errstate(divide='ignore', invalid='ignore'):
    return np.where(count > 0, total / np.maximum(count, 1), 0.0)


def feature_track(
    indices: np.ndarray, values: np.ndarray, feature_id: int
) -> np.ndarray:
  """Per-residue activation of one feature. Shape (L,); 0 where not in top-k."""
  track = np.zeros(indices.shape[0], dtype=np.float64)
  rows, cols = np.nonzero(indices == int(feature_id))
  track[rows] = values[rows, cols]
  return track


def rank_features(
    indices: np.ndarray,
    values: np.ndarray,
    codebook: int,
    threshold: float,
    limit: int,
) -> dict[str, Any]:
  """Ranks features two ways: by peak activation and by prevalence.

  Peak activation surfaces motif-like / local features (a sharp spike at a few
  residues). Prevalence surfaces domain- or family-level features (broadly on
  across the chain). They answer different questions; always read both.
  """
  peak = max_activation(indices, values, codebook)
  prev = prevalence(indices, values, codebook, threshold)
  mean_on = mean_when_on(indices, values, codebook, threshold)
  length = int(indices.shape[0])

  def row(feature_id: int) -> dict[str, Any]:
    feature_id = int(feature_id)
    track = feature_track(indices, values, feature_id)
    return {
        'feature_index': feature_id,
        'max_activation': round(float(peak[feature_id]), 6),
        'prevalence': int(prev[feature_id]),
        'prevalence_fraction': round(float(prev[feature_id]) / length, 4),
        'mean_activation_when_on': round(float(mean_on[feature_id]), 6),
        'peak_residue': int(np.argmax(track)) + 1,  # 1-based
    }

  # argsort is ascending; reverse for descending. Ties broken by feature index,
  # which keeps the ordering deterministic across runs.
  by_max = np.argsort(-peak, kind='stable')[:limit]
  by_prev = np.argsort(-prev, kind='stable')[:limit]

  return {
      'sequence_length': length,
      'codebook_size': int(codebook),
      'activation_threshold': threshold,
      'top_by_max_activation': [row(f) for f in by_max],
      'top_by_prevalence': [row(f) for f in by_prev],
  }


# --------------------------------------------------------------------------
# Persistence
# --------------------------------------------------------------------------


def save_features(
    path: str,
    sequence: str,
    indices: np.ndarray,
    values: np.ndarray,
    codebook: int,
    sae_model: str,
    normalized: bool,
) -> pathlib.Path:
  """Writes the sparse top-k representation to a .npz."""
  out = pathlib.Path(path).expanduser()
  if out.suffix != '.npz':
    out = out.with_suffix('.npz')
  out.parent.mkdir(parents=True, exist_ok=True)
  np.savez_compressed(
      out,
      indices=indices.astype(np.int32),
      values=values.astype(np.float32),
      codebook_size=np.int64(codebook),
      sequence=np.array(sequence),
      sae_model=np.array(sae_model),
      normalized=np.array(bool(normalized)),
  )
  return out


def load_features(path: str) -> dict[str, Any]:
  """Reads a .npz written by `extract`."""
  file = pathlib.Path(path).expanduser()
  if not file.is_file():
    raise BiohubError(
        f'Feature file not found: {file}. Run `extract` first to create it.'
    )
  with np.load(file, allow_pickle=False) as data:
    return {
        'indices': data['indices'].astype(np.int64),
        'values': data['values'].astype(np.float64),
        'codebook_size': int(data['codebook_size']),
        'sequence': str(data['sequence']),
        'sae_model': str(data['sae_model']),
        'normalized': bool(data['normalized']),
    }


# --------------------------------------------------------------------------
# Feature descriptions (separate, unauthenticated host; cached on disk)
# --------------------------------------------------------------------------


class FeatureDescriptions:
  """Fetches SAE feature descriptions, with an on-disk cache.

  The descriptions are static per feature index, so we cache them forever and
  only ever hit the network on a miss.
  """

  def __init__(self, cache_dir: str | pathlib.Path = DEFAULT_CACHE_DIR):
    self.cache_dir = pathlib.Path(cache_dir).expanduser()
    self.cache_dir.mkdir(parents=True, exist_ok=True)
    self._http = http_client.HttpClient(
        eb.BASE_URL, qps=FEATURE_API_QPS, timeout=60.0, max_retries=5
    )
    self.hits = 0
    self.misses = 0

  def get(self, feature_index: int) -> dict[str, Any]:
    """Returns the description record for one feature index."""
    feature_index = int(feature_index)
    cached = self.cache_dir / f'{feature_index}.json'
    if cached.is_file():
      try:
        with open(cached, encoding='utf-8') as handle:
          self.hits += 1
          return json.load(handle)
      except (ValueError, OSError):
        cached.unlink(missing_ok=True)  # Corrupt cache entry; refetch.

    url = f'{eb.SAE_FEATURE_API}/{feature_index}'
    try:
      data = self._http.fetch_json(url)
    except HttpError as exc:
      raise BiohubError(
          f'Feature description lookup failed for feature {feature_index} '
          f'({url}): {exc}'
      ) from exc
    if not isinstance(data, dict):
      raise BiohubError(
          f'Unexpected response for feature {feature_index}: {type(data)}'
      )
    self.misses += 1
    with open(cached, 'w', encoding='utf-8') as handle:
      json.dump(data, handle)
    return data

  def get_many(self, feature_indices) -> dict[int, dict[str, Any]]:
    return {int(f): self.get(int(f)) for f in feature_indices}


def slim(record: dict[str, Any]) -> dict[str, Any]:
  """Keeps the interpretable fields and drops the bulky activation lists."""
  out = {k: record.get(k) for k in DESCRIPTION_FIELDS}
  swissprot = record.get('top_swissprot_activations') or []
  out['top_swissprot'] = [
      {
          'uniprot_id': entry.get('uniprot_id'),
          'activation': round(float(entry.get('activation', 0.0)), 3),
      }
      for entry in swissprot[:5]
  ]
  return out


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------


def resolve_sequence(args) -> str:
  """Takes the sequence from --sequence or the first record of --fasta."""
  if getattr(args, 'fasta', None):
    records = eb.read_fasta(args.fasta)
    return eb.validate_sequence(records[0][1])
  return eb.validate_sequence(args.sequence)


def extract_features(
    client: eb.BiohubClient, sequence: str, sae_model: str, normalize: bool
) -> tuple[np.ndarray, np.ndarray, int]:
  """Calls the API. BOS/EOS are trimmed by eb.sae_features -- residue i is row i."""
  # `eb.sae_features` derives the ESMC model as `sae_model.split('-sae-')[0]`.
  # A name without '-sae-' (e.g. the underscore convention used in some ESM
  # cookbook snippets, 'esmc-600m-2024-12_k64_codebook16384_layer27') therefore
  # silently becomes a bogus ESMC model name and fails downstream with an
  # opaque "model not found". Catch it here with an actionable message.
  if '-sae-' not in sae_model:
    raise BiohubError(
        f'Malformed SAE model name {sae_model!r}. Expected the format '
        "'{esmc_model}-sae-layer{L}-k{k}-codebook{C}', e.g. "
        f'{eb.DEFAULT_SAE_MODEL}. (The underscore style used in some ESM '
        'cookbook snippets is not what this API serves.)'
    )
  if normalize and '300m' in sae_model.lower():
    raise BiohubError(
        'normalize_features=True is rejected by the API for ESMC 300M SAE '
        'models. Re-run with --no-normalize.'
    )
  return client.sae_features(
      sequence, sae_model=sae_model, normalize_features=normalize
  )


def cmd_extract(args) -> None:
  """Sequence -> sparse SAE features (.npz) + a JSON summary."""
  sequence = resolve_sequence(args)
  client = eb.BiohubClient()
  normalize = not args.no_normalize

  indices, values, codebook = extract_features(
      client, sequence, args.sae_model, normalize
  )

  if indices.shape[0] != len(sequence):
    raise BiohubError(
        f'Feature rows ({indices.shape[0]}) != sequence length '
        f'({len(sequence)}). BOS/EOS trimming is inconsistent; aborting rather '
        'than silently mis-aligning residues.'
    )

  npz = save_features(
      args.output, sequence, indices, values, codebook, args.sae_model,
      normalize,
  )

  active_per_residue = (values > 0).sum(axis=1)
  distinct = np.unique(indices)
  summary = {
      'sequence_length': len(sequence),
      'sae_model': args.sae_model,
      'esmc_model': args.sae_model.split('-sae-')[0],
      'normalize_features': normalize,
      'codebook_size': int(codebook),
      'top_k': int(indices.shape[1]),
      'active_features_per_residue': {
          'min': int(active_per_residue.min()),
          'max': int(active_per_residue.max()),
          'mean': round(float(active_per_residue.mean()), 3),
      },
      'distinct_features_fired': int(distinct.size),
      'fraction_of_codebook_used': round(float(distinct.size) / codebook, 5),
      'sparsity': round(
          1.0 - float(indices.shape[1]) / float(codebook), 6
      ),
      'features_npz': str(npz),
      'note': (
          'BOS/EOS already trimmed: row i of the feature arrays is residue i '
          '(0-based). Feature descriptions are auto-generated hypotheses.'
      ),
  }
  eb.write_json(summary, args.output_summary or f'{npz.with_suffix("")}.json')
  print(
      f'Extracted {indices.shape[1]} active features/residue over '
      f'{len(sequence)} residues ({distinct.size} distinct features fired) '
      f'-> {npz}'
  )


def cmd_top_features(args) -> None:
  """Ranks features by peak activation and by prevalence."""
  data = load_features(args.features)
  ranking = rank_features(
      data['indices'],
      data['values'],
      data['codebook_size'],
      args.threshold,
      args.limit,
  )
  ranking['sae_model'] = data['sae_model']
  ranking['normalize_features'] = data['normalized']
  eb.write_json(ranking, args.output)

  print(f'\nTop {args.limit} by MAX ACTIVATION (motif-like / local):')
  for rank, item in enumerate(ranking['top_by_max_activation'], 1):
    print(
        f'  {rank:2d}. feature {item["feature_index"]:5d}  '
        f'max={item["max_activation"]:.3f}  '
        f'prevalence={item["prevalence"]:3d}  '
        f'peak@residue {item["peak_residue"]}'
    )
  print(f'\nTop {args.limit} by PREVALENCE (domain / family-like):')
  for rank, item in enumerate(ranking['top_by_prevalence'], 1):
    print(
        f'  {rank:2d}. feature {item["feature_index"]:5d}  '
        f'prevalence={item["prevalence"]:3d}/'
        f'{ranking["sequence_length"]}  '
        f'max={item["max_activation"]:.3f}'
    )


def markdown_report(
    ranking: dict[str, Any],
    descriptions: dict[int, dict[str, Any]],
    sae_model: str,
) -> str:
  """Renders the joined ranking + descriptions as markdown."""
  lines = [
      '# SAE feature interpretation',
      '',
      f'* SAE model: `{sae_model}`',
      f'* Sequence length: {ranking["sequence_length"]} residues',
      f'* Codebook: {ranking["codebook_size"]} features, '
      f'top-k active per residue',
      f'* Activation threshold for prevalence: '
      f'{ranking["activation_threshold"]}',
      '',
      '> **These descriptions are auto-generated hypotheses, not curated**',
      '> **annotations.** They were produced by an automated agent that read '
      'each',
      '> feature\'s activation pattern across large protein databases. A high '
      'activation',
      '> means the model associates the concept with this region -- it is NOT '
      'evidence',
      '> that the protein truly has that function. Treat every row below as a '
      'lead to',
      '> verify, never as an established fact.',
      '',
  ]

  sections = [
      (
          'Top features by peak activation (motif-like / local)',
          'top_by_max_activation',
          'These spike at a few residues: candidate catalytic residues, '
          'binding motifs, or short conserved patterns.',
      ),
      (
          'Top features by prevalence (domain / family-like)',
          'top_by_prevalence',
          'These stay on across much of the chain: candidate fold class, '
          'domain identity, or protein-family signals.',
      ),
  ]

  for title, key, blurb in sections:
    lines += [f'## {title}', '', blurb, '']
    lines += [
        '| # | Feature | Label | Category | Max act. | Prevalence | '
        'Peak residue |',
        '|---|---|---|---|---|---|---|',
    ]
    for rank, item in enumerate(ranking[key], 1):
      idx = item['feature_index']
      info = descriptions.get(idx, {})
      label = str(info.get('label') or '(no label)').replace('|', '/')
      category = str(info.get('category') or '-').replace('|', '/')
      lines.append(
          f'| {rank} | `{idx}` | {label} | {category} | '
          f'{item["max_activation"]:.3f} | '
          f'{item["prevalence"]}/{ranking["sequence_length"]} | '
          f'{item["peak_residue"]} |'
      )
    lines.append('')

    lines += [f'### Detail: {title.split(" (")[0].lower()}', '']
    for rank, item in enumerate(ranking[key], 1):
      idx = item['feature_index']
      info = descriptions.get(idx, {})
      lines += [
          f'**{rank}. Feature `{idx}` — {info.get("label") or "(no label)"}**',
          '',
          f'* Category: {info.get("category") or "-"}',
          f'* Max activation: {item["max_activation"]:.3f} '
          f'(peak at residue {item["peak_residue"]})',
          f'* Prevalence: {item["prevalence"]}/{ranking["sequence_length"]} '
          f'residues above threshold',
          f'* Hypothesised concept: {info.get("summary") or "-"}',
          f'* Activation pattern: {info.get("activation_pattern") or "-"}',
          f'* Exemplar families: {info.get("exemplar_protein_families") or "-"}',
      ]
      swissprot = info.get('top_swissprot') or []
      if swissprot:
        ids = ', '.join(
            f'{e["uniprot_id"]} ({e["activation"]:.1f})' for e in swissprot[:3]
        )
        lines.append(f'* Strongest SwissProt activations: {ids}')
      lines.append('')

  return '\n'.join(lines) + '\n'


def cmd_describe(args) -> None:
  """Fetches natural-language descriptions and joins them onto the rankings."""
  descriptions_api = FeatureDescriptions(args.cache_dir)

  if args.feature_ids:
    wanted = [int(x) for x in args.feature_ids.split(',') if x.strip()]
    records = descriptions_api.get_many(wanted)
    payload = {
        'sae_model': args.sae_model,
        'features': {str(k): slim(v) for k, v in records.items()},
        'disclaimer': (
            'Auto-generated hypotheses from activation patterns, not curated '
            'annotations. Do not state them as fact.'
        ),
    }
    eb.write_json(payload, args.output)
    for idx in wanted:
      info = records[idx]
      print(f'  feature {idx:5d}: {info.get("label")} [{info.get("category")}]')
    print(
        f'\n{descriptions_api.misses} fetched, {descriptions_api.hits} from '
        f'cache.'
    )
    return

  if not args.features:
    raise BiohubError('Pass either --features (a .npz) or --feature-ids.')

  data = load_features(args.features)
  ranking = rank_features(
      data['indices'],
      data['values'],
      data['codebook_size'],
      args.threshold,
      args.limit,
  )
  wanted = {
      item['feature_index']
      for key in ('top_by_max_activation', 'top_by_prevalence')
      for item in ranking[key]
  }
  records = descriptions_api.get_many(sorted(wanted))
  slimmed = {idx: slim(rec) for idx, rec in records.items()}

  for key in ('top_by_max_activation', 'top_by_prevalence'):
    for item in ranking[key]:
      info = slimmed.get(item['feature_index'], {})
      item['label'] = info.get('label')
      item['category'] = info.get('category')
      item['summary'] = info.get('summary')
      item['description'] = info.get('description')

  payload = dict(ranking)
  payload['sae_model'] = data['sae_model']
  payload['descriptions'] = {str(k): v for k, v in slimmed.items()}
  payload['disclaimer'] = (
      'SAE feature descriptions are auto-generated hypotheses derived from '
      'activation patterns across protein databases. They are NOT curated '
      'annotations and must not be reported as established fact.'
  )
  eb.write_json(payload, args.output)

  if args.report:
    report = markdown_report(ranking, slimmed, data['sae_model'])
    report_path = pathlib.Path(args.report).expanduser()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding='utf-8')
    print(f'Success! Report written to: {report_path}')

  print(f'\nTop {args.limit} by max activation:')
  for rank, item in enumerate(ranking['top_by_max_activation'], 1):
    print(
        f'  {rank:2d}. feature {item["feature_index"]:5d}  '
        f'max={item["max_activation"]:.3f}  {item.get("label")}'
    )
  print(f'\nTop {args.limit} by prevalence:')
  for rank, item in enumerate(ranking['top_by_prevalence'], 1):
    print(
        f'  {rank:2d}. feature {item["feature_index"]:5d}  '
        f'prev={item["prevalence"]:3d}  {item.get("label")}'
    )
  print(
      f'\n{descriptions_api.misses} descriptions fetched, '
      f'{descriptions_api.hits} served from cache.'
  )
  print(
      'REMINDER: these labels are auto-generated hypotheses, not curated '
      'annotations.'
  )


def cmd_plot(args) -> None:
  """Plots per-residue activation tracks for chosen features."""
  import matplotlib

  matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  data = load_features(args.features)
  sequence = data['sequence']
  wanted = [int(x) for x in args.feature_ids.split(',') if x.strip()]
  if not wanted:
    raise BiohubError('--feature-ids must name at least one feature.')

  labels: dict[int, str] = {}
  if not args.no_labels:
    try:
      api = FeatureDescriptions(args.cache_dir)
      for feature_id in wanted:
        labels[feature_id] = str(api.get(feature_id).get('label') or '')
    except BiohubError as exc:
      print(f'[warn] could not fetch labels ({exc}); plotting without.',
            file=sys.stderr)

  positions = np.arange(1, len(sequence) + 1)
  fig, axes = plt.subplots(
      len(wanted), 1, figsize=(13, 2.2 * len(wanted)), sharex=True, squeeze=False
  )
  axes = axes[:, 0]

  for ax, feature_id in zip(axes, wanted):
    track = feature_track(data['indices'], data['values'], feature_id)
    ax.fill_between(positions, track, color='steelblue', alpha=0.35)
    ax.plot(positions, track, color='steelblue', linewidth=1.2)
    title = f'Feature {feature_id}'
    if labels.get(feature_id):
      title += f' — {labels[feature_id]}'
    ax.set_title(title, fontsize=10, loc='left')
    ax.set_ylabel('activation', fontsize=8)
    ax.margins(x=0)
    # Headroom so the peak annotation cannot collide with the title above.
    ax.set_ylim(0, max(float(track.max()), 1e-6) * 1.25)

    if track.max() > 0:
      peak = int(np.argmax(track))
      ax.annotate(
          f'{sequence[peak]}{peak + 1}',
          (peak + 1, track[peak]),
          textcoords='offset points',
          xytext=(0, 5),
          ha='center',
          fontsize=8,
      )

  axes[-1].set_xlabel('Residue position (1-based)')
  fig.suptitle(
      'SAE feature activation tracks — labels are auto-generated hypotheses',
      fontsize=9,
      y=1.0,
  )
  fig.tight_layout()

  out = pathlib.Path(args.output).expanduser()
  out.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(out, dpi=150, bbox_inches='tight')
  plt.close(fig)
  print(f'Success! Plot written to: {out}')


def cmd_compare(args) -> None:
  """Compares the feature sets of two sequences (Jaccard overlap)."""
  seq_a = eb.validate_sequence(args.sequence_a)
  seq_b = eb.validate_sequence(args.sequence_b)
  client = eb.BiohubClient()
  normalize = not args.no_normalize

  idx_a, val_a, codebook = extract_features(
      client, seq_a, args.sae_model, normalize
  )
  idx_b, val_b, _ = extract_features(client, seq_b, args.sae_model, normalize)

  peak_a = max_activation(idx_a, val_a, codebook)
  peak_b = max_activation(idx_b, val_b, codebook)
  top_a = set(int(f) for f in np.argsort(-peak_a, kind='stable')[: args.top_n])
  top_b = set(int(f) for f in np.argsort(-peak_b, kind='stable')[: args.top_n])

  shared = top_a & top_b
  union = top_a | top_b
  jaccard = len(shared) / len(union) if union else 0.0

  api = FeatureDescriptions(args.cache_dir)
  shared_rows = []
  for feature_id in sorted(
      shared, key=lambda f: -min(peak_a[f], peak_b[f])
  )[: args.limit]:
    info = api.get(feature_id)
    shared_rows.append({
        'feature_index': feature_id,
        'label': info.get('label'),
        'category': info.get('category'),
        'max_activation_a': round(float(peak_a[feature_id]), 4),
        'max_activation_b': round(float(peak_b[feature_id]), 4),
    })

  def distinct_rows(own, other, peak):
    rows = []
    for feature_id in sorted(own - other, key=lambda f: -peak[f])[: args.limit]:
      info = api.get(feature_id)
      rows.append({
          'feature_index': feature_id,
          'label': info.get('label'),
          'category': info.get('category'),
          'max_activation': round(float(peak[feature_id]), 4),
      })
    return rows

  result = {
      'sae_model': args.sae_model,
      'top_n': args.top_n,
      'length_a': len(seq_a),
      'length_b': len(seq_b),
      'jaccard': round(jaccard, 4),
      'n_shared': len(shared),
      'n_union': len(union),
      'shared_features': shared_rows,
      'only_in_a': distinct_rows(top_a, top_b, peak_a),
      'only_in_b': distinct_rows(top_b, top_a, peak_b),
      'interpretation': _interpret_jaccard(jaccard),
      'disclaimer': (
          'Feature labels are auto-generated hypotheses, not curated '
          'annotations.'
      ),
  }
  eb.write_json(result, args.output)
  print(
      f'\nJaccard(top-{args.top_n}) = {jaccard:.3f} '
      f'({len(shared)}/{len(union)} features shared)'
  )
  print(f'  {result["interpretation"]}')
  if shared_rows:
    print('\nShared features (strongest first):')
    for row in shared_rows[:10]:
      print(f'  {row["feature_index"]:5d}  {row["label"]}')


def _interpret_jaccard(jaccard: float) -> str:
  """Plain-language reading of the overlap. Calibrated on this SAE."""
  if jaccard >= 0.7:
    return (
        'Very high overlap: near-identical feature machinery (point mutants, '
        'close homologs).'
    )
  if jaccard >= 0.3:
    return (
        'Substantial overlap: the model sees shared domains, fold, or family '
        'signals.'
    )
  if jaccard >= 0.1:
    return (
        'Modest overlap: some shared machinery, but the proteins are largely '
        'distinct.'
    )
  return (
      'Low overlap: the model represents these as unrelated proteins '
      '(unrelated proteins typically score < 0.05).'
  )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          'Interpret ESMC representations through a sparse autoencoder: which '
          'of the 16,384 learned features fire on your protein, and what do '
          'they mean.'
      )
  )
  sub = parser.add_subparsers(dest='command', required=True)

  def add_sae_model(p):
    p.add_argument(
        '--sae-model',
        default=eb.DEFAULT_SAE_MODEL,
        help=f'SAE codebook name (default: {eb.DEFAULT_SAE_MODEL}).',
    )
    p.add_argument(
        '--no-normalize',
        action='store_true',
        help=(
            'Disable TF-IDF feature normalization. REQUIRED for ESMC 300M SAE '
            'models, which reject normalize_features=True.'
        ),
    )

  def add_cache(p):
    p.add_argument(
        '--cache-dir',
        default=str(DEFAULT_CACHE_DIR),
        help=(
            'On-disk cache for feature descriptions (they are static). '
            f'Default: {DEFAULT_CACHE_DIR}'
        ),
    )

  # extract
  p = sub.add_parser(
      'extract', help='Sequence -> sparse SAE features (.npz + JSON summary).'
  )
  group = p.add_mutually_exclusive_group(required=True)
  group.add_argument('--sequence', help='Protein sequence (one-letter codes).')
  group.add_argument('--fasta', help='FASTA file; the first record is used.')
  p.add_argument('--output', required=True, help='Output .npz path.')
  p.add_argument(
      '--output-summary', help='Output JSON summary path (default: <output>.json).'
  )
  add_sae_model(p)
  p.set_defaults(func=cmd_extract)

  # top-features
  p = sub.add_parser(
      'top-features',
      help='Rank features by peak activation and by prevalence.',
  )
  p.add_argument('--features', required=True, help='.npz from `extract`.')
  p.add_argument(
      '--limit',
      type=int,
      required=True,
      help='How many features to return in each ranking.',
  )
  p.add_argument(
      '--threshold',
      type=float,
      default=DEFAULT_ACTIVATION_THRESHOLD,
      help=(
          'Activation above which a feature counts as "on" at a residue, for '
          f'the prevalence ranking (default: {DEFAULT_ACTIVATION_THRESHOLD}).'
      ),
  )
  p.add_argument('--output', required=True, help='Output JSON path.')
  p.set_defaults(func=cmd_top_features)

  # describe
  p = sub.add_parser(
      'describe',
      help='Fetch natural-language descriptions and join them to the rankings.',
  )
  p.add_argument('--features', help='.npz from `extract`.')
  p.add_argument(
      '--feature-ids',
      help='Comma-separated feature indices to describe directly (no .npz).',
  )
  p.add_argument(
      '--limit',
      type=int,
      default=10,
      help='Features per ranking to describe (ignored with --feature-ids).',
  )
  p.add_argument(
      '--threshold',
      type=float,
      default=DEFAULT_ACTIVATION_THRESHOLD,
      help='Activation threshold for the prevalence ranking.',
  )
  p.add_argument('--output', required=True, help='Output JSON path.')
  p.add_argument('--report', help='Optional markdown report path.')
  p.add_argument(
      '--sae-model',
      default=eb.DEFAULT_SAE_MODEL,
      help='Recorded in the output (descriptions are codebook-wide).',
  )
  add_cache(p)
  p.set_defaults(func=cmd_describe)

  # plot
  p = sub.add_parser(
      'plot', help='Per-residue activation track for chosen features (PNG).'
  )
  p.add_argument('--features', required=True, help='.npz from `extract`.')
  p.add_argument(
      '--feature-ids',
      required=True,
      help='Comma-separated feature indices to plot, e.g. 2773,8581.',
  )
  p.add_argument('--output', required=True, help='Output PNG path.')
  p.add_argument(
      '--no-labels',
      action='store_true',
      help='Skip the description lookup used for plot titles.',
  )
  add_cache(p)
  p.set_defaults(func=cmd_plot)

  # compare
  p = sub.add_parser(
      'compare',
      help='Jaccard overlap of the top-N feature sets of two sequences.',
  )
  p.add_argument('--sequence-a', required=True, help='First protein sequence.')
  p.add_argument('--sequence-b', required=True, help='Second protein sequence.')
  p.add_argument(
      '--top-n',
      type=int,
      required=True,
      help='Size of each top-feature set (by peak activation), e.g. 50.',
  )
  p.add_argument(
      '--limit',
      type=int,
      default=15,
      help='How many shared/distinct features to describe in the output.',
  )
  p.add_argument('--output', required=True, help='Output JSON path.')
  add_sae_model(p)
  add_cache(p)
  p.set_defaults(func=cmd_compare)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
