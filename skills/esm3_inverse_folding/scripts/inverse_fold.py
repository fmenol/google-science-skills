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

"""Inverse folding with ESM3: design sequences that fold into a given backbone.

The dedicated `/inverse_fold` endpoint is NOT supported by any model this API key
can reach (every model answers HTTP 422 "does not support the 'inverse_fold'
endpoint"). Inverse folding therefore goes through the ESM3 generate endpoint
with coordinates supplied and the sequence track left masked:

    client.generate('sequence', coordinates=coords_atom37, ...)

Designs are validated by self-consistency: re-fold each design with ESMFold2 and
compare the refolded structure back to the input backbone (scTM / scRMSD).
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
import concurrent.futures
import csv
import itertools
import pathlib
import sys

import numpy as np

from esm_biohub import BiohubClient, BiohubError  # vendored, same directory
import esm_biohub as eb

# ESM3 reads the residue frame from N, CA and C only. We additionally pass O
# because it costs nothing and keeps the written backbone chemically complete.
# VERIFIED against the live API on ubiquitin: sending all 37 atoms, N/CA/C/O, or
# N/CA/C all return byte-identical designs (100% native recovery at T=0.1), so
# stripping sidechains is a numeric no-op that *guarantees* no sidechain identity
# can leak into the design. Sending a CA-only trace does NOT error -- it silently
# returns a ~5%-identity random sequence. See `_check_backbone`.
BACKBONE_ATOMS = ('N', 'CA', 'C', 'O')
FRAME_ATOMS = ('N', 'CA', 'C')

CA_INDEX = eb.ATOM37.index('CA')

# Max decoding steps the API accepts (it is also capped at the sequence length).
MAX_NUM_STEPS = 100

MAX_WORKERS = 8


# --------------------------------------------------------------------------
# Backbone loading and validation
# --------------------------------------------------------------------------


def _expand_to_atom37(coords: np.ndarray) -> np.ndarray:
  """Accepts (L, 37, 3), (L, 4, 3) = N,CA,C,O or (L, 3, 3) = N,CA,C.

  The 3- and 4-atom layouts are what backbone generators (RFdiffusion and
  friends) emit, so they are accepted directly and padded out to atom37 with NaN.

  Raises:
    BiohubError: If the array is not one of the recognised layouts.
  """
  coords = np.asarray(coords, dtype=np.float64)
  if coords.ndim != 3 or coords.shape[2] != 3:
    raise BiohubError(
        f'Coordinates must have shape (L, n_atoms, 3); got {coords.shape}.'
    )
  n_atoms = coords.shape[1]
  if n_atoms == 37:
    return coords
  if n_atoms in (3, 4):
    names = FRAME_ATOMS if n_atoms == 3 else BACKBONE_ATOMS
    out = np.full((coords.shape[0], 37, 3), np.nan, dtype=np.float64)
    for slot, name in enumerate(names):
      out[:, eb.ATOM37.index(name)] = coords[:, slot]
    return out
  raise BiohubError(
      f'Coordinates have {n_atoms} atoms per residue. Expected 37 (atom37), '
      '4 (N, CA, C, O) or 3 (N, CA, C).'
  )


def _check_backbone(coords: np.ndarray) -> dict[str, float]:
  """Verifies the N/CA/C frame exists. A broken frame silently ruins designs.

  ESM3 builds its geometric representation from the N-CA-C frame of each
  residue. If those atoms are absent the API does *not* complain -- it returns a
  fluent-looking sequence with chance-level (~5%) identity to the native. Measured
  on ubiquitin: a CA-only trace scores 0.046 recovery, i.e. pure noise. So we
  refuse the input instead of designing garbage.

  Returns:
    Per-atom occupancy fractions for N, CA and C.

  Raises:
    BiohubError: If fewer than half the residues carry a complete N/CA/C frame.
  """
  occupancy = {
      name: float(
          np.isfinite(coords[:, eb.ATOM37.index(name)]).all(axis=-1).mean()
      )
      for name in FRAME_ATOMS
  }
  frame_ok = np.ones(coords.shape[0], dtype=bool)
  for name in FRAME_ATOMS:
    frame_ok &= np.isfinite(coords[:, eb.ATOM37.index(name)]).all(axis=-1)
  complete = float(frame_ok.mean())

  if complete < 0.5:
    missing = [n for n, f in occupancy.items() if f < 0.5]
    raise BiohubError(
        f'Only {complete:.0%} of residues have a complete N/CA/C backbone '
        f'frame (missing/sparse: {missing}). ESM3 derives its structure '
        'representation from that frame; with it absent the API still returns '
        'a sequence, but it is chance-level noise (measured: 4.6% identity for '
        'a CA-only trace). Supply a full backbone (N, CA, C, O), not a CA '
        'trace.'
    )
  if complete < 1.0:
    print(
        f'[warning] {1.0 - complete:.1%} of residues lack a complete N/CA/C '
        'frame; those positions are designed from context alone.',
        file=sys.stderr,
    )
  occupancy['complete_frame'] = complete
  return occupancy


def _strip_to_backbone(coords: np.ndarray) -> np.ndarray:
  """NaNs out every atom37 slot except N, CA, C, O."""
  out = np.full_like(coords, np.nan)
  for name in BACKBONE_ATOMS:
    idx = eb.ATOM37.index(name)
    out[:, idx] = coords[:, idx]
  return out


def load_backbone(args) -> tuple[np.ndarray, str | None]:
  """Loads the target backbone from --pdb or --coords.

  Returns:
    (coordinates (L, 37, 3), native sequence or None)
  """
  if getattr(args, 'pdb', None):
    native, coords, _ = eb.parse_pdb_atom37(args.pdb, getattr(args, 'chain', None))
  else:
    coords = _expand_to_atom37(np.load(args.coords))
    native = None

  _check_backbone(coords)
  if not getattr(args, 'all_atom', False):
    coords = _strip_to_backbone(coords)

  native_override = getattr(args, 'native', None)
  if native_override:
    native = eb.validate_sequence(native_override)
    if len(native) != coords.shape[0]:
      raise BiohubError(
          f'--native has {len(native)} residues but the backbone has '
          f'{coords.shape[0]}.'
      )
  return coords, native


def _resolve_num_steps(requested: int | None, length: int) -> int:
  """num_steps must be <= sequence length; the API caps it at 100."""
  if requested is None:
    return max(1, min(length, MAX_NUM_STEPS))
  if requested < 1:
    raise BiohubError('--num-steps must be >= 1.')
  if requested > length:
    raise BiohubError(
        f'--num-steps ({requested}) must be <= the backbone length ({length}).'
    )
  return min(requested, MAX_NUM_STEPS)


# --------------------------------------------------------------------------
# Design
# --------------------------------------------------------------------------


def design_sequences(
    client: BiohubClient,
    coords: np.ndarray,
    *,
    num_samples: int,
    temperature: float,
    num_steps: int,
    model: str,
) -> list[str]:
  """Draws `num_samples` independent designs for one backbone.

  Diversity comes from RESAMPLING, not from raising the temperature: each call
  is an independent stochastic decode. Cranking the temperature past ~1.5
  destroys the design (measured on ubiquitin: recovery 1.00 at T=0.1, 0.98 at
  T=1.0, 0.90 at T=1.5, 0.40 at T=2.0).
  """

  def one(_: int) -> str:
    out = client.generate(
        'sequence',
        model=model,
        coordinates=coords,
        num_steps=num_steps,
        temperature=temperature,
    )
    sequence = (out.get('outputs') or {}).get('sequence')
    if not sequence:
      raise BiohubError(
          'ESM3 generate returned no sequence. Response keys: '
          f'{list((out.get("outputs") or {}).keys())}'
      )
    return sequence

  with concurrent.futures.ThreadPoolExecutor(
      max_workers=min(MAX_WORKERS, num_samples)
  ) as pool:
    return list(pool.map(one, range(num_samples)))


def refold(
    client: BiohubClient, sequences: list[str], model: str
) -> list[dict]:
  """Re-folds each design with ESMFold2 (independent of ESM3). Threaded."""

  def one(sequence: str) -> dict:
    out = client.fold(sequence, model=model)
    return {
        'coordinates': eb.to_array(out['coordinates']),
        'plddt': float(np.nanmean(eb.to_array(out['plddt']))),
        'ptm': float(out['ptm']) if out.get('ptm') is not None else float('nan'),
    }

  with concurrent.futures.ThreadPoolExecutor(
      max_workers=min(MAX_WORKERS, len(sequences))
  ) as pool:
    return list(pool.map(one, sequences))


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------


def identity(a: str, b: str) -> float:
  """Fraction of identical positions over the shorter sequence."""
  n = min(len(a), len(b))
  if not n:
    return 0.0
  return sum(x == y for x, y in zip(a[:n], b[:n])) / n


def pairwise_identity(sequences: list[str]) -> tuple[list[list[float]], float]:
  """Returns (identity matrix, mean off-diagonal identity)."""
  n = len(sequences)
  matrix = np.eye(n)
  for i, j in itertools.combinations(range(n), 2):
    value = identity(sequences[i], sequences[j])
    matrix[i, j] = matrix[j, i] = value
  if n < 2:
    return matrix.tolist(), float('nan')
  mean = float(matrix[np.triu_indices(n, k=1)].mean())
  return matrix.tolist(), mean


def self_consistency_scores(
    design_coords: np.ndarray, target_coords: np.ndarray
) -> tuple[float, float]:
  """scTM and scRMSD of a refolded design against the target backbone.

  Both structures are the same length and already in residue correspondence
  (position i of the design occupies position i of the backbone), so the CA
  traces can be superposed directly.
  """
  design_ca = design_coords[:, CA_INDEX]
  target_ca = target_coords[:, CA_INDEX]
  n = min(len(design_ca), len(target_ca))
  design_ca, target_ca = design_ca[:n], target_ca[:n]
  sctm = eb.tm_score(design_ca, target_ca)
  scrmsd, _ = eb.kabsch_rmsd(design_ca, target_ca)
  return float(sctm), float(scrmsd)


def _verdict(sctm: float) -> str:
  if sctm > 0.8:
    return 'PASS: very likely adopts the target fold'
  if sctm >= 0.5:
    return 'UNCERTAIN: partially recovers the fold'
  return 'FAIL: does not fold back to the target'


def _write_fasta(path: pathlib.Path, records: list[tuple[str, str]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with open(path, 'w', encoding='utf-8') as handle:
    for name, sequence in records:
      handle.write(f'>{name}\n')
      for i in range(0, len(sequence), 60):
        handle.write(sequence[i : i + 60] + '\n')
  print(f'Success! Data written to: {path}')


def _sibling(output: str, suffix: str, override: str | None) -> pathlib.Path:
  if override:
    return pathlib.Path(override).expanduser()
  return pathlib.Path(output).expanduser().with_suffix(suffix)


# --------------------------------------------------------------------------
# Subcommands
# --------------------------------------------------------------------------


def cmd_design(args) -> None:
  """Designs N sequences for a backbone and reports recovery + diversity."""
  client = BiohubClient()
  coords, native = load_backbone(args)
  length = coords.shape[0]
  num_steps = _resolve_num_steps(args.num_steps, length)

  print(
      f'Designing {args.num_samples} sequence(s) for a {length}-residue '
      f'backbone (T={args.temperature}, num_steps={num_steps})...'
  )
  sequences = design_sequences(
      client,
      coords,
      num_samples=args.num_samples,
      temperature=args.temperature,
      num_steps=num_steps,
      model=args.model,
  )

  matrix, mean_pairwise = pairwise_identity(sequences)
  designs = []
  for i, sequence in enumerate(sequences, start=1):
    entry = {
        'id': f'design_{i}',
        'sequence': sequence,
        'length': len(sequence),
        'native_recovery': identity(sequence, native) if native else None,
    }
    designs.append(entry)

  report = {
      'input': args.pdb or args.coords,
      'chain': args.chain,
      'backbone_length': length,
      'native_sequence': native,
      'model': args.model,
      'temperature': args.temperature,
      'num_steps': num_steps,
      'num_samples': args.num_samples,
      'conditioning': 'all_atom' if args.all_atom else 'backbone_only(N,CA,C,O)',
      'designs': designs,
      'diversity': {
          'mean_pairwise_identity': mean_pairwise,
          'mean_pairwise_diversity': (
              float('nan') if np.isnan(mean_pairwise) else 1.0 - mean_pairwise
          ),
          'unique_designs': len(set(sequences)),
          'identity_matrix': matrix,
      },
  }
  if native:
    recoveries = [d['native_recovery'] for d in designs]
    report['recovery'] = {
        'best': max(recoveries),
        'mean': float(np.mean(recoveries)),
        'worst': min(recoveries),
    }

  eb.write_json(report, args.output)
  _write_fasta(
      _sibling(args.output, '.fasta', args.fasta),
      [
          (
              f'{d["id"]} T={args.temperature}'
              + (
                  f' native_recovery={d["native_recovery"]:.3f}'
                  if d['native_recovery'] is not None
                  else ''
              ),
              d['sequence'],
          )
          for d in designs
      ],
  )

  if native:
    print(
        f'Native-sequence recovery: best {report["recovery"]["best"]:.1%}, '
        f'mean {report["recovery"]["mean"]:.1%}'
    )
  print(
      f'Diversity: {len(set(sequences))}/{len(sequences)} unique, '
      f'mean pairwise identity {mean_pairwise:.1%}'
  )
  print(
      'These designs are UNVALIDATED. Run `self-consistency` before reporting '
      'or ordering any of them.'
  )


def cmd_self_consistency(args) -> None:
  """Designs, re-folds, and ranks by how well each design regenerates the fold."""
  client = BiohubClient()
  coords, native = load_backbone(args)
  length = coords.shape[0]
  num_steps = _resolve_num_steps(args.num_steps, length)

  print(
      f'Designing {args.num_samples} sequence(s) for a {length}-residue '
      f'backbone (T={args.temperature}, num_steps={num_steps})...'
  )
  sequences = design_sequences(
      client,
      coords,
      num_samples=args.num_samples,
      temperature=args.temperature,
      num_steps=num_steps,
      model=args.model,
  )

  print(f'Re-folding {len(sequences)} design(s) with {args.fold_model}...')
  folds = refold(client, sequences, args.fold_model)

  rows = []
  for i, (sequence, fold_result) in enumerate(zip(sequences, folds), start=1):
    sctm, scrmsd = self_consistency_scores(fold_result['coordinates'], coords)
    rows.append({
        'id': f'design_{i}',
        'sequence': sequence,
        'length': len(sequence),
        'sctm': sctm,
        'scrmsd': scrmsd,
        'plddt': fold_result['plddt'],
        'ptm': fold_result['ptm'],
        'native_recovery': identity(sequence, native) if native else None,
        'verdict': _verdict(sctm),
    })

  # Rank by scTM (higher is better), breaking ties on scRMSD (lower is better).
  rows.sort(key=lambda r: (-r['sctm'], r['scrmsd']))
  for rank, row in enumerate(rows, start=1):
    row['rank'] = rank

  _, mean_pairwise = pairwise_identity(sequences)
  best = rows[0]
  report = {
      'input': args.pdb or args.coords,
      'chain': args.chain,
      'backbone_length': length,
      'native_sequence': native,
      'design_model': args.model,
      'fold_model': args.fold_model,
      'temperature': args.temperature,
      'num_steps': num_steps,
      'num_samples': args.num_samples,
      'conditioning': 'all_atom' if args.all_atom else 'backbone_only(N,CA,C,O)',
      'summary': {
          'best_id': best['id'],
          'best_sctm': best['sctm'],
          'best_scrmsd': best['scrmsd'],
          'best_plddt': best['plddt'],
          'best_verdict': best['verdict'],
          'n_passing_sctm_0.8': sum(r['sctm'] > 0.8 for r in rows),
          'n_passing_sctm_0.5': sum(r['sctm'] >= 0.5 for r in rows),
          'mean_pairwise_identity': mean_pairwise,
          'unique_designs': len(set(sequences)),
      },
      'designs': rows,
  }
  if native:
    report['summary']['best_native_recovery'] = max(
        r['native_recovery'] for r in rows
    )

  eb.write_json(report, args.output)

  csv_path = _sibling(args.output, '.csv', args.csv)
  csv_path.parent.mkdir(parents=True, exist_ok=True)
  columns = [
      'rank', 'id', 'sctm', 'scrmsd', 'plddt', 'ptm', 'native_recovery',
      'verdict', 'sequence',
  ]
  with open(csv_path, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
  print(f'Success! Data written to: {csv_path}')

  _write_fasta(
      _sibling(args.output, '.fasta', args.fasta),
      [
          (
              f'{r["id"]} rank={r["rank"]} scTM={r["sctm"]:.3f} '
              f'scRMSD={r["scrmsd"]:.2f} pLDDT={r["plddt"]:.3f}',
              r['sequence'],
          )
          for r in rows
      ],
  )

  print('\nRanked designs (scTM > 0.8 = adopts the target fold):')
  for row in rows:
    print(
        f'  {row["rank"]}. {row["id"]:10s} scTM={row["sctm"]:.3f}  '
        f'scRMSD={row["scrmsd"]:5.2f} A  pLDDT={row["plddt"]:.3f}  '
        f'{row["verdict"]}'
    )


def cmd_recovery(args) -> None:
  """Reports native-sequence recovery per position, with an agreement plot."""
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  from matplotlib.colors import ListedColormap

  client = BiohubClient()
  coords, native = load_backbone(args)
  if not native:
    raise BiohubError(
        'Recovery needs a native sequence to compare against. Pass --pdb (the '
        'sequence is read from its ATOM records) or supply --native.'
    )
  length = coords.shape[0]
  num_steps = _resolve_num_steps(args.num_steps, length)

  print(
      f'Designing {args.num_samples} sequence(s) for a {length}-residue '
      f'backbone (T={args.temperature}, num_steps={num_steps})...'
  )
  sequences = design_sequences(
      client,
      coords,
      num_samples=args.num_samples,
      temperature=args.temperature,
      num_steps=num_steps,
      model=args.model,
  )

  ragged = [len(s) for s in sequences if len(s) != len(native)]
  if ragged:
    raise BiohubError(
        f'Designs must be the same length as the native sequence '
        f'({len(native)}); got {ragged}. Per-position recovery is undefined '
        'without an alignment.'
    )
  matches = np.array(
      [[s[i] == native[i] for i in range(len(native))] for s in sequences],
      dtype=bool,
  )
  per_position = matches.mean(axis=0)
  per_design = matches.mean(axis=1)

  positions = [
      {
          'position': i + 1,
          'native': native[i],
          'agreement': float(per_position[i]),
          'designed': [s[i] for s in sequences],
      }
      for i in range(len(native))
  ]

  report = {
      'input': args.pdb or args.coords,
      'chain': args.chain,
      'backbone_length': length,
      'native_sequence': native,
      'model': args.model,
      'temperature': args.temperature,
      'num_steps': num_steps,
      'num_samples': args.num_samples,
      'overall_recovery': float(matches.mean()),
      'best_design_recovery': float(per_design.max()),
      'worst_design_recovery': float(per_design.min()),
      'per_design_recovery': [float(x) for x in per_design],
      'fully_conserved_positions': int((per_position == 1.0).sum()),
      'never_recovered_positions': [
          i + 1 for i in range(len(native)) if per_position[i] == 0.0
      ],
      'designs': [
          {'id': f'design_{i+1}', 'sequence': s, 'recovery': float(per_design[i])}
          for i, s in enumerate(sequences)
      ],
      'per_position': positions,
  }
  eb.write_json(report, args.output)

  plot_path = _sibling(args.output, '.png', args.plot)
  plot_path.parent.mkdir(parents=True, exist_ok=True)
  height = 2.6 + 0.28 * len(sequences)
  fig, (ax_top, ax_bottom) = plt.subplots(
      2, 1, figsize=(max(9.0, length * 0.13), height),
      sharex=True, height_ratios=[max(1.0, 0.32 * len(sequences)), 1.4],
      constrained_layout=True,
  )

  ax_top.imshow(
      matches, aspect='auto', interpolation='nearest',
      cmap=ListedColormap(['#c0392b', '#27ae60']), vmin=0, vmax=1,
      extent=(0.5, length + 0.5, len(sequences) - 0.5, -0.5),
  )
  ax_top.set_yticks(range(len(sequences)))
  ax_top.set_yticklabels(
      [f'design_{i+1} ({per_design[i]:.0%})' for i in range(len(sequences))],
      fontsize=8,
  )
  ax_top.set_title(
      f'Native-sequence recovery — {args.num_samples} design(s) at T='
      f'{args.temperature}   (green = matches native, red = differs)\n'
      f'overall {matches.mean():.1%}   best design {per_design.max():.1%}',
      fontsize=10,
  )

  ax_bottom.bar(
      np.arange(1, length + 1), per_position, width=1.0,
      color='#2980b9', edgecolor='none',
  )
  ax_bottom.axhline(
      matches.mean(), color='#e67e22', linestyle='--', linewidth=1.2,
      label=f'overall {matches.mean():.1%}',
  )
  ax_bottom.axhline(
      1.0 / len(eb.AA20), color='#7f8c8d', linestyle=':', linewidth=1.2,
      label=f'random baseline {1.0 / len(eb.AA20):.1%}',
  )
  ax_bottom.set_ylim(0, 1.05)
  ax_bottom.set_xlim(0.5, length + 0.5)
  ax_bottom.set_ylabel('fraction of\ndesigns matching', fontsize=9)
  ax_bottom.set_xlabel('residue position', fontsize=9)
  ax_bottom.legend(fontsize=8, loc='lower right', framealpha=0.9)
  fig.savefig(plot_path, dpi=150)
  plt.close(fig)
  print(f'Success! Data written to: {plot_path}')

  print(
      f'Overall recovery {matches.mean():.1%} '
      f'(best design {per_design.max():.1%}, worst {per_design.min():.1%}); '
      f'{report["fully_conserved_positions"]}/{length} positions recovered by '
      'every design.'
  )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _add_backbone_args(parser: argparse.ArgumentParser) -> None:
  source = parser.add_mutually_exclusive_group(required=True)
  source.add_argument('--pdb', help='Input backbone (.pdb / .cif).')
  source.add_argument(
      '--coords',
      help='Input backbone as .npy: (L, 37, 3), (L, 4, 3) N/CA/C/O, or '
      '(L, 3, 3) N/CA/C.',
  )
  parser.add_argument(
      '--chain', help='Chain to read from --pdb. Default: the first chain.'
  )
  parser.add_argument(
      '--num-samples', type=int, required=True,
      help='Number of designs to sample. Required. Diversity comes from '
      'resampling, so ask for several (8-32 is typical).',
  )
  parser.add_argument(
      '--temperature', type=float, default=0.1,
      help='Sampling temperature. 0.1 (default) maximises native recovery; '
      '~1.0 adds diversity; >1.5 destroys the design.',
  )
  parser.add_argument(
      '--num-steps', type=int, default=None,
      help='ESM3 decoding steps. Must be <= backbone length; the API caps it '
      'at 100. Default: min(length, 100).',
  )
  parser.add_argument(
      '--model', default=eb.DEFAULT_ESM3,
      help=f'ESM3 model. Default: {eb.DEFAULT_ESM3}.',
  )
  parser.add_argument(
      '--all-atom', action='store_true',
      help='Condition on every atom in the file instead of stripping to the '
      'N/CA/C/O backbone. Verified to be a no-op (ESM3 reads only the N/CA/C '
      'frame); the default keeps the design provably sidechain-blind.',
  )
  parser.add_argument(
      '--output', required=True, help='Output JSON file path.'
  )


def main() -> None:
  parser = argparse.ArgumentParser(
      description='ESM3 inverse folding: design sequences for a 3D backbone, '
      'and validate them by self-consistency.'
  )
  sub = parser.add_subparsers(dest='command', required=True)

  p_design = sub.add_parser(
      'design', help='Design N sequences for a backbone (FASTA + JSON).'
  )
  _add_backbone_args(p_design)
  p_design.add_argument(
      '--native', help='Native sequence to score recovery against (optional; '
      'read from --pdb automatically).'
  )
  p_design.add_argument(
      '--fasta', help='FASTA output path. Default: --output with .fasta.'
  )
  p_design.set_defaults(func=cmd_design)

  p_sc = sub.add_parser(
      'self-consistency',
      help='Design, re-fold each design with ESMFold2, rank by scTM/scRMSD. '
      'ALWAYS run this before trusting a design.',
  )
  _add_backbone_args(p_sc)
  p_sc.add_argument(
      '--native', help='Native sequence to score recovery against (optional).'
  )
  p_sc.add_argument(
      '--fold-model', default=eb.DEFAULT_ESMFOLD2,
      help=f'ESMFold2 model used to re-fold. Default: {eb.DEFAULT_ESMFOLD2}.',
  )
  p_sc.add_argument(
      '--csv', help='Ranked CSV output path. Default: --output with .csv.'
  )
  p_sc.add_argument(
      '--fasta', help='FASTA output path. Default: --output with .fasta.'
  )
  p_sc.set_defaults(func=cmd_self_consistency)

  p_rec = sub.add_parser(
      'recovery',
      help='Per-position native-sequence recovery + agreement plot.',
  )
  _add_backbone_args(p_rec)
  p_rec.add_argument(
      '--native', help='Native sequence. Required if using --coords.'
  )
  p_rec.add_argument(
      '--plot', help='PNG output path. Default: --output with .png.'
  )
  p_rec.set_defaults(func=cmd_recovery)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
