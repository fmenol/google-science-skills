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

"""Screens candidate binders against a target by folding each complex.

Folds every target:candidate complex with ESMFold2 and ranks the candidates by
interface confidence (iPTM). This is the SELECTION half of the binder-design
protocol: it decides which candidates are worth ordering. It does NOT generate
candidates de novo.
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
import json
import pathlib
import sys
from typing import Any

import numpy as np

from esm_biohub import BiohubClient, BiohubError  # vendored, same directory
import esm_biohub as eb

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

TARGET_CHAIN = 'A'
BINDER_CHAIN = 'B'

# A cross-chain residue pair is "in contact" when its representative atoms (CB;
# CA for glycine, which has none) are within this distance. The 8 A CB-CB rule
# is the conventional definition of a protein-protein interface contact.
CB_CONTACT_CUTOFF = 8.0

# A contact is only believable if the model is also confident about where the
# two chains sit RELATIVE to each other. PAE(i, j) is the expected error in the
# position of residue j when the structure is aligned on residue i, so a
# cross-chain pair with a large PAE is a contact the model does not believe in.
# Measured: for a true complex (barnase:barstar) every geometric contact clears
# this gate; for non-binders essentially none do. See SKILL.md.
CONTACT_PAE_GATE = 15.0

# PAE is reported in Angstrom and saturates around 32 A (the AlphaFold/ESMFold
# convention). Used to fold interface PAE into the composite score.
PAE_SATURATION = 31.75

# The API rejects num_sampling_steps > 100 with HTTP 422 (verified).
MAX_SAMPLING_STEPS = 100

# Composite score weights. iPTM dominates because it is the only metric that
# directly measures interface confidence; binder pLDDT guards against a binder
# that cannot fold at all; interface PAE breaks ties. Contact COUNT is
# deliberately excluded -- it scales with binder size and does not discriminate
# binders from non-binders (see SKILL.md).
SCORE_WEIGHTS = {'iptm': 0.5, 'binder_plddt': 0.3, 'interface_pae': 0.2}

SORT_KEYS = (
    'iptm',
    'selection_score',
    'binder_plddt',
    'interface_pae',
    'confident_interface_contacts',
    'interface_contacts',
    'ptm',
    'isoelectric_point',
)

# Metrics where a SMALLER value is better.
LOWER_IS_BETTER = ('interface_pae', 'isoelectric_point')

CSV_COLUMNS = (
    'rank',
    'candidate',
    'selected',
    'iptm',
    'selection_score',
    'binder_plddt',
    'interface_pae',
    'confident_interface_contacts',
    'interface_contacts',
    'ptm',
    'complex_plddt',
    'target_plddt',
    'n_interface_residues_target',
    'n_interface_residues_binder',
    'isoelectric_point',
    'binder_length',
    'verdict',
    'pdb',
)

# pKa values for isoelectric-point estimation (EMBOSS scale, as used by
# Biopython's ProteinAnalysis, which the reference protocol calls).
PKA_SIDE_CHAIN = {
    'K': 10.8, 'R': 12.5, 'H': 6.5, 'D': 3.9, 'E': 4.1, 'C': 8.5, 'Y': 10.1,
}
PKA_N_TERM = 8.6
PKA_C_TERM = 3.6
POSITIVE_AA = ('K', 'R', 'H')
NEGATIVE_AA = ('D', 'E', 'C', 'Y')


# --------------------------------------------------------------------------
# Interface geometry and metrics
# --------------------------------------------------------------------------


def _chain_letters(complex_data: dict[str, Any]) -> np.ndarray:
  """Maps each token's raw chain id to its author chain letter (e.g. 'A')."""
  lookup = (complex_data.get('metadata') or {}).get('chain_lookup', {})
  return np.array(
      [str(lookup.get(str(c), str(c)))[:1] for c in complex_data['chain_id']]
  )


def representative_coords(complex_data: dict[str, Any]) -> np.ndarray:
  """Returns (L, 3) CB coordinates per token, falling back to CA for glycine.

  Missing atoms stay NaN, which propagates to an infinite distance and is
  therefore never counted as a contact.
  """
  positions = eb.to_array(complex_data['atom_positions'], dtype=np.float64)
  atom_names = complex_data['atom_names']
  token_to_atoms = complex_data['token_to_atoms']

  coords = np.full((len(token_to_atoms), 3), np.nan, dtype=np.float64)
  for token_idx, (start, end) in enumerate(token_to_atoms):
    slots = {atom_names[a]: a for a in range(int(start), int(end))}
    atom = slots.get('CB', slots.get('CA'))
    if atom is not None:
      coords[token_idx] = positions[atom]
  return coords


def compute_interface(
    result: dict[str, Any], binder_sequence: str
) -> dict[str, Any]:
  """Derives every interface metric from one `fold_all_atom` response.

  Args:
    result: The raw API response (must have been requested with include_pae).
    binder_sequence: The candidate's sequence, used for pI and length.

  Returns:
    A dict of scalar metrics plus the per-residue interface detail.
  """
  complex_data = result['complex']
  letters = _chain_letters(complex_data)
  coords = representative_coords(complex_data)
  plddt = eb.to_array(result['plddt'])

  target_idx = np.flatnonzero(letters == TARGET_CHAIN)
  binder_idx = np.flatnonzero(letters == BINDER_CHAIN)
  if target_idx.size == 0 or binder_idx.size == 0:
    raise BiohubError(
        f'Expected chains {TARGET_CHAIN} and {BINDER_CHAIN} in the folded '
        f'complex; found {sorted(set(letters))}.'
    )

  pae = result.get('pae')
  if pae is None:
    raise BiohubError(
        'The response carries no PAE matrix. Fold with include_pae=True.'
    )
  pae = eb.to_array(pae)
  # PAE is asymmetric (row = frame of reference). Symmetrize for the contact
  # gate; average both off-diagonal blocks for the scalar interface PAE.
  pae_sym = (pae + pae.T) / 2.0
  interface_pae = float(
      (
          pae[np.ix_(target_idx, binder_idx)].mean()
          + pae[np.ix_(binder_idx, target_idx)].mean()
      )
      / 2.0
  )

  # Cross-chain CB-CB distances: (n_target, n_binder).
  delta = coords[target_idx][:, None, :] - coords[binder_idx][None, :, :]
  with np.errstate(invalid='ignore'):
    distances = np.linalg.norm(delta, axis=-1)
  distances = np.nan_to_num(distances, nan=np.inf)

  close = distances < CB_CONTACT_CUTOFF
  confident = close & (pae_sym[np.ix_(target_idx, binder_idx)] < CONTACT_PAE_GATE)

  residues = complex_data['sequence']
  target_residues = _interface_residues(
      target_idx, confident, distances, residues, plddt, axis=1, chain=TARGET_CHAIN
  )
  binder_residues = _interface_residues(
      binder_idx, confident, distances, residues, plddt, axis=0, chain=BINDER_CHAIN
  )

  binder_plddt = float(plddt[binder_idx].mean())
  iptm = float(result['interface_ptm'])
  metrics = {
      'iptm': iptm,
      'ptm': float(result['ptm']),
      'binder_plddt': binder_plddt,
      'target_plddt': float(plddt[target_idx].mean()),
      'complex_plddt': float(plddt.mean()),
      'interface_pae': interface_pae,
      'interface_contacts': int(close.sum()),
      'confident_interface_contacts': int(confident.sum()),
      'n_interface_residues_target': len(target_residues),
      'n_interface_residues_binder': len(binder_residues),
      'min_cross_chain_distance': float(
          distances.min() if np.isfinite(distances).any() else np.inf
      ),
      'isoelectric_point': isoelectric_point(binder_sequence),
      'binder_length': len(binder_sequence),
      'target_length': int(target_idx.size),
      'selection_score': selection_score(iptm, binder_plddt, interface_pae),
      'verdict': verdict(iptm),
  }
  return {
      'metrics': metrics,
      'interface_residues_target': target_residues,
      'interface_residues_binder': binder_residues,
  }


def _interface_residues(
    chain_idx: np.ndarray,
    confident: np.ndarray,
    distances: np.ndarray,
    residues: list[str],
    plddt: np.ndarray,
    axis: int,
    chain: str,
) -> list[dict[str, Any]]:
  """Lists the residues of one chain that make a confident cross-chain contact.

  Args:
    chain_idx: Token indices of this chain in the complex.
    confident: (n_target, n_binder) boolean contact matrix, PAE-gated.
    distances: (n_target, n_binder) CB-CB distances.
    residues: Three-letter residue names, one per token.
    plddt: Per-token pLDDT.
    axis: Axis to reduce over -- 1 for the target (rows), 0 for the binder.
    chain: Chain letter to report.

  Returns:
    One record per interface residue, ordered by position in the chain.
  """
  counts = confident.sum(axis=axis)
  masked = np.where(confident, distances, np.inf)
  nearest = masked.min(axis=axis)

  out: list[dict[str, Any]] = []
  for position, token in enumerate(chain_idx):
    if counts[position] == 0:
      continue
    out.append({
        'chain': chain,
        'residue_number': position + 1,  # 1-based within the chain, as in the PDB
        'residue': str(residues[token]),
        'contacts': int(counts[position]),
        'min_distance': round(float(nearest[position]), 2),
        'plddt': round(float(plddt[token]), 3),
    })
  return out


def selection_score(
    iptm: float, binder_plddt: float, interface_pae: float
) -> float:
  """Composite 0-1 score. Documented in SKILL.md; do not recompute by hand.

    selection_score = 0.5 * iPTM
                    + 0.3 * binder_pLDDT
                    + 0.2 * (1 - min(interface_PAE, 31.75) / 31.75)

  Every term is already on a 0-1 scale where higher is better (interface PAE is
  inverted). iPTM carries half the weight because it is the only term that
  measures the INTERFACE; the other two guard against a binder that does not
  fold and an interface the model cannot place. Contact count is excluded on
  purpose: it scales with binder size and does not separate binders from
  non-binders.
  """
  pae_term = 1.0 - min(interface_pae, PAE_SATURATION) / PAE_SATURATION
  return round(
      SCORE_WEIGHTS['iptm'] * iptm
      + SCORE_WEIGHTS['binder_plddt'] * binder_plddt
      + SCORE_WEIGHTS['interface_pae'] * pae_term,
      4,
  )


def verdict(iptm: float) -> str:
  """Turns iPTM into the calibrated band. NOT an affinity prediction."""
  if iptm > 0.8:
    return 'confident interface'
  if iptm >= 0.5:
    return 'possible interface'
  return 'likely no binding'


def isoelectric_point(sequence: str) -> float:
  """Estimates pI by bisection on the Henderson-Hasselbalch net charge.

  The reference protocol filters minibinders to pI < 6 before ordering, as a
  proxy for solubility and expression.
  """
  counts = {aa: sequence.count(aa) for aa in set(PKA_SIDE_CHAIN)}

  def net_charge(ph: float) -> float:
    charge = 1.0 / (1.0 + 10 ** (ph - PKA_N_TERM))
    charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERM - ph))
    for aa in POSITIVE_AA:
      charge += counts.get(aa, 0) / (1.0 + 10 ** (ph - PKA_SIDE_CHAIN[aa]))
    for aa in NEGATIVE_AA:
      charge -= counts.get(aa, 0) / (1.0 + 10 ** (PKA_SIDE_CHAIN[aa] - ph))
    return charge

  low, high = 0.0, 14.0
  for _ in range(100):
    mid = (low + high) / 2.0
    if net_charge(mid) > 0:
      low = mid
    else:
      high = mid
  return round((low + high) / 2.0, 2)


# --------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------


def sort_rows(rows: list[dict[str, Any]], sort_by: str) -> list[dict[str, Any]]:
  """Sorts candidates best-first by `sort_by`."""
  reverse = sort_by not in LOWER_IS_BETTER
  return sorted(rows, key=lambda r: r['metrics'][sort_by], reverse=reverse)


def apply_filters(rows: list[dict[str, Any]], args) -> list[dict[str, Any]]:
  """Drops candidates that fail any threshold the user set."""
  kept = []
  for row in rows:
    m = row['metrics']
    if args.min_iptm is not None and m['iptm'] < args.min_iptm:
      continue
    if (
        args.min_binder_plddt is not None
        and m['binder_plddt'] < args.min_binder_plddt
    ):
      continue
    if (
        args.max_interface_pae is not None
        and m['interface_pae'] > args.max_interface_pae
    ):
      continue
    if (
        args.min_contacts is not None
        and m['confident_interface_contacts'] < args.min_contacts
    ):
      continue
    if (
        args.max_isoelectric_point is not None
        and m['isoelectric_point'] > args.max_isoelectric_point
    ):
      continue
    kept.append(row)
  return kept


def finalize(
    rows: list[dict[str, Any]], sort_by: str, top_n: int
) -> list[dict[str, Any]]:
  """Sorts, assigns ranks, and flags the top-n shortlist."""
  rows = sort_rows(rows, sort_by)
  for i, row in enumerate(rows, start=1):
    row['rank'] = i
    row['selected'] = i <= top_n
  return rows


def write_table(rows: list[dict[str, Any]], path: pathlib.Path) -> None:
  """Writes the ranked table as CSV."""
  path.parent.mkdir(parents=True, exist_ok=True)
  with open(path, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, extrasaction='ignore')
    writer.writeheader()
    for row in rows:
      record = {
          'rank': row['rank'],
          'candidate': row['candidate'],
          'selected': row['selected'],
          'pdb': row.get('pdb', ''),
      }
      for key, value in row['metrics'].items():
        record[key] = round(value, 4) if isinstance(value, float) else value
      writer.writerow(record)
  print(f'Success! Data written to: {path}')


def print_table(rows: list[dict[str, Any]], sort_by: str) -> None:
  """Prints a short ranked summary to stdout."""
  print(f'\n{"rank":>4}  {"candidate":<24} {"iPTM":>6} {"bpLDDT":>7} '
        f'{"iPAE":>6} {"contacts":>9} {"score":>6}  verdict')
  print('-' * 92)
  for row in rows:
    m = row['metrics']
    mark = '*' if row['selected'] else ' '
    print(
        f'{row["rank"]:>3}{mark}  {row["candidate"][:24]:<24} '
        f'{m["iptm"]:>6.3f} {m["binder_plddt"]:>7.3f} '
        f'{m["interface_pae"]:>6.2f} {m["confident_interface_contacts"]:>9d} '
        f'{m["selection_score"]:>6.3f}  {m["verdict"]}'
    )
  print(f'\nRanked by {sort_by}. * = shortlisted (--top-n).')


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------


def _load_target(args) -> tuple[str, str]:
  """Returns (target_name, target_sequence) from --target or --target-fasta."""
  if args.target_fasta:
    records = eb.read_fasta(args.target_fasta)
    if len(records) != 1:
      raise BiohubError(
          f'--target-fasta must hold exactly one record; found {len(records)}.'
      )
    name, sequence = records[0]
    return name, eb.validate_sequence(sequence)
  return 'target', eb.validate_sequence(args.target)


def _check_sampling_steps(steps: int) -> None:
  if not 1 <= steps <= MAX_SAMPLING_STEPS:
    raise BiohubError(
        f'--num-sampling-steps must be between 1 and {MAX_SAMPLING_STEPS}; got '
        f'{steps}. The API rejects anything larger with HTTP 422. (The paper '
        f'uses 200 for its final critic; that is not reachable through the API.)'
    )


def fold_complex(
    client: BiohubClient, target: str, binder: str, args
) -> dict[str, Any]:
  """Folds one target:binder complex. Target is chain A, binder is chain B."""
  return client.fold_all_atom(
      [
          {'type': 'protein', 'id': TARGET_CHAIN, 'sequence': target},
          {'type': 'protein', 'id': BINDER_CHAIN, 'sequence': binder},
      ],
      model=args.model,
      include_pae=True,
      num_loops=args.num_loops,
      num_sampling_steps=args.num_sampling_steps,
  )


def cmd_screen(args) -> None:
  """Folds every target:candidate complex and ranks the candidates."""
  _check_sampling_steps(args.num_sampling_steps)
  target_name, target = _load_target(args)
  candidates = [
      (name, eb.validate_sequence(seq))
      for name, seq in eb.read_fasta(args.candidates)
  ]

  out_dir = pathlib.Path(args.output).expanduser()
  pdb_dir = out_dir / 'structures'
  pdb_dir.mkdir(parents=True, exist_ok=True)

  client = BiohubClient()
  workers = max(1, min(args.max_workers, 4))
  print(
      f'Screening {len(candidates)} candidate(s) against {target_name} '
      f'({len(target)} aa) with {args.model}...',
      file=sys.stderr,
  )

  def screen_one(item: tuple[str, str]) -> dict[str, Any]:
    name, sequence = item
    result = fold_complex(client, target, sequence, args)
    analysis = compute_interface(result, sequence)
    pdb_path = pdb_dir / f'{name}.pdb'
    pdb_path.write_text(eb.complex_to_pdb(result['complex']), encoding='utf-8')
    print(
        f'  folded {name}: iPTM={analysis["metrics"]["iptm"]:.3f}',
        file=sys.stderr,
        flush=True,
    )
    return {
        'candidate': name,
        'binder_sequence': sequence,
        'pdb': str(pdb_path),
        **analysis,
    }

  with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    rows = list(pool.map(screen_one, candidates))

  rows = finalize(rows, args.sort_by, args.top_n)

  payload = {
      'target': {'name': target_name, 'sequence': target, 'length': len(target)},
      'params': {
          'model': args.model,
          'num_loops': args.num_loops,
          'num_sampling_steps': args.num_sampling_steps,
          'sort_by': args.sort_by,
          'top_n': args.top_n,
          'cb_contact_cutoff': CB_CONTACT_CUTOFF,
          'contact_pae_gate': CONTACT_PAE_GATE,
          'score_weights': SCORE_WEIGHTS,
      },
      'ranking_metric': args.sort_by,
      'results': rows,
      'shortlist': [r['candidate'] for r in rows if r['selected']],
  }
  eb.write_json(payload, str(out_dir / 'results.json'))
  write_table(rows, out_dir / 'results.csv')
  print_table(rows, args.sort_by)
  print(f'\nStructures: {pdb_dir}')


def cmd_rank(args) -> None:
  """Re-ranks and filters an existing results.json."""
  with open(args.results, encoding='utf-8') as handle:
    payload = json.load(handle)
  rows = payload['results']
  before = len(rows)
  rows = apply_filters(rows, args)
  rows = finalize(rows, args.sort_by, args.top_n)

  payload['ranking_metric'] = args.sort_by
  payload['results'] = rows
  payload['shortlist'] = [r['candidate'] for r in rows if r['selected']]
  payload['filters'] = {
      'min_iptm': args.min_iptm,
      'min_binder_plddt': args.min_binder_plddt,
      'max_interface_pae': args.max_interface_pae,
      'min_contacts': args.min_contacts,
      'max_isoelectric_point': args.max_isoelectric_point,
  }
  payload['n_candidates_before_filter'] = before
  payload['n_candidates_after_filter'] = len(rows)

  out_path = pathlib.Path(args.output).expanduser()
  eb.write_json(payload, str(out_path))
  write_table(rows, out_path.with_suffix('.csv'))
  print(f'\n{len(rows)}/{before} candidate(s) passed the filters.')
  if rows:
    print_table(rows, args.sort_by)
  else:
    print('No candidate passed. Loosen the thresholds, or accept that none of '
          'these candidates is predicted to bind.')


def cmd_analyze_interface(args) -> None:
  """Deep-dives one target:binder complex."""
  _check_sampling_steps(args.num_sampling_steps)
  target_name, target = _load_target(args)
  if args.binder_fasta:
    records = eb.read_fasta(args.binder_fasta)
    if len(records) != 1:
      raise BiohubError(
          f'--binder-fasta must hold exactly one record; found {len(records)}.'
      )
    binder_name, binder = records[0][0], eb.validate_sequence(records[0][1])
  else:
    binder_name, binder = 'binder', eb.validate_sequence(args.binder)

  out_dir = pathlib.Path(args.output).expanduser()
  out_dir.mkdir(parents=True, exist_ok=True)

  client = BiohubClient()
  result = fold_complex(client, target, binder, args)
  analysis = compute_interface(result, binder)
  metrics = analysis['metrics']

  pdb_path = out_dir / f'{target_name}_{binder_name}.pdb'
  pdb_path.write_text(eb.complex_to_pdb(result['complex']), encoding='utf-8')

  pae = eb.to_array(result['pae'])
  png_path = out_dir / f'{target_name}_{binder_name}_pae.png'
  plot_pae(pae, len(target), target_name, binder_name, metrics, png_path)

  payload = {
      'target': {'name': target_name, 'sequence': target, 'length': len(target)},
      'binder': {'name': binder_name, 'sequence': binder, 'length': len(binder)},
      'params': {
          'model': args.model,
          'num_loops': args.num_loops,
          'num_sampling_steps': args.num_sampling_steps,
          'cb_contact_cutoff': CB_CONTACT_CUTOFF,
          'contact_pae_gate': CONTACT_PAE_GATE,
      },
      'pdb': str(pdb_path),
      'pae_heatmap': str(png_path),
      **analysis,
  }
  eb.write_json(payload, str(out_dir / f'{target_name}_{binder_name}_interface.json'))

  print(f'\nInterface: {target_name} (chain A) : {binder_name} (chain B)')
  print(f'  iPTM                         {metrics["iptm"]:.3f}  '
        f'-> {metrics["verdict"]}')
  print(f'  pTM                          {metrics["ptm"]:.3f}')
  print(f'  binder mean pLDDT            {metrics["binder_plddt"]:.3f}')
  print(f'  interface PAE (A)            {metrics["interface_pae"]:.2f}  '
        f'(lower is better)')
  print(f'  confident contacts           '
        f'{metrics["confident_interface_contacts"]}')
  print(f'  geometric contacts           {metrics["interface_contacts"]}  '
        f'(NOT a binding signal on its own)')
  print(f'  interface residues           '
        f'{metrics["n_interface_residues_target"]} on {target_name}, '
        f'{metrics["n_interface_residues_binder"]} on {binder_name}')
  print(f'  selection_score              {metrics["selection_score"]:.3f}')
  print(f'\nStructure: {pdb_path}\nPAE heatmap: {png_path}')


def plot_pae(
    pae: np.ndarray,
    target_length: int,
    target_name: str,
    binder_name: str,
    metrics: dict[str, Any],
    path: pathlib.Path,
) -> None:
  """Draws the PAE matrix with the chain boundary marked."""
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  fig, ax = plt.subplots(figsize=(7, 6))
  image = ax.imshow(pae, cmap='Greens_r', vmin=0, vmax=PAE_SATURATION,
                    interpolation='nearest')
  ax.axhline(target_length - 0.5, color='crimson', linewidth=1.2)
  ax.axvline(target_length - 0.5, color='crimson', linewidth=1.2)

  total = pae.shape[0]
  ax.set_xticks([target_length / 2, target_length + (total - target_length) / 2])
  ax.set_xticklabels([f'{target_name}\n(chain A)', f'{binder_name}\n(chain B)'])
  ax.set_yticks([target_length / 2, target_length + (total - target_length) / 2])
  ax.set_yticklabels([f'{target_name}', f'{binder_name}'], rotation=90, va='center')
  ax.set_title(
      f'PAE — {target_name}:{binder_name}\n'
      f'iPTM {metrics["iptm"]:.3f} · interface PAE {metrics["interface_pae"]:.1f} Å'
      f' · {metrics["verdict"]}'
  )
  bar = fig.colorbar(image, ax=ax, shrink=0.8)
  bar.set_label('Predicted aligned error (Å) — lower is better')
  fig.text(
      0.5, 0.015,
      'Off-diagonal blocks = cross-chain. Dark = confident relative placement.',
      ha='center', fontsize=8, color='0.35',
  )
  fig.tight_layout(rect=(0, 0.03, 1, 1))
  fig.savefig(path, dpi=150)
  plt.close(fig)
  print(f'Success! Data written to: {path}')


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _add_fold_args(parser: argparse.ArgumentParser) -> None:
  parser.add_argument(
      '--model', default=eb.DEFAULT_ESMFOLD2,
      choices=sorted(eb.ESMFOLD2_MODELS), help='ESMFold2 model.',
  )
  parser.add_argument(
      '--num-loops', type=int, default=3,
      help='Trunk recycling loops. The reference protocol scores with 3.',
  )
  parser.add_argument(
      '--num-sampling-steps', type=int, default=100,
      help=f'Diffusion sampling steps, 1-{MAX_SAMPLING_STEPS} (API cap).',
  )


def _add_target_args(parser: argparse.ArgumentParser) -> None:
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--target', help='Target protein sequence.')
  group.add_argument('--target-fasta', help='FASTA holding one target record.')


def _add_filter_args(parser: argparse.ArgumentParser) -> None:
  parser.add_argument('--min-iptm', type=float,
                      help='Drop candidates below this iPTM (0.8 = confident).')
  parser.add_argument('--min-binder-plddt', type=float,
                      help='Drop candidates whose binder folds below this pLDDT.')
  parser.add_argument('--max-interface-pae', type=float,
                      help='Drop candidates above this interface PAE (Angstrom).')
  parser.add_argument('--min-contacts', type=int,
                      help='Drop candidates with fewer confident contacts.')
  parser.add_argument(
      '--max-isoelectric-point', type=float,
      help='Drop candidates above this pI. The reference protocol uses 6 for '
           'minibinders as a solubility/expression proxy.',
  )


def main() -> None:
  parser = argparse.ArgumentParser(
      description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  sub = parser.add_subparsers(dest='command', required=True)

  screen = sub.add_parser(
      'screen', help='Fold every target:candidate complex and rank them.'
  )
  _add_target_args(screen)
  screen.add_argument('--candidates', required=True,
                      help='FASTA of candidate binder sequences.')
  screen.add_argument('--output', required=True, help='Output directory.')
  screen.add_argument(
      '--top-n', type=int, required=True,
      help='How many candidates to shortlist. Required: choose deliberately.',
  )
  screen.add_argument('--sort-by', default='iptm', choices=SORT_KEYS,
                      help='Ranking metric (default: iptm).')
  screen.add_argument('--max-workers', type=int, default=4,
                      help='Concurrent folds (capped at 4).')
  _add_fold_args(screen)
  screen.set_defaults(func=cmd_screen)

  rank = sub.add_parser(
      'rank', help='Re-rank or filter an existing results.json.'
  )
  rank.add_argument('--results', required=True,
                    help='results.json written by `screen`.')
  rank.add_argument('--output', required=True, help='Output JSON path.')
  rank.add_argument(
      '--top-n', type=int, required=True,
      help='How many candidates to shortlist. Required: choose deliberately.',
  )
  rank.add_argument('--sort-by', default='iptm', choices=SORT_KEYS,
                    help='Ranking metric (default: iptm).')
  _add_filter_args(rank)
  rank.set_defaults(func=cmd_rank)

  interface = sub.add_parser(
      'analyze-interface', help='Deep-dive one target:binder complex.'
  )
  _add_target_args(interface)
  binder_group = interface.add_mutually_exclusive_group(required=True)
  binder_group.add_argument('--binder', help='Binder protein sequence.')
  binder_group.add_argument('--binder-fasta',
                            help='FASTA holding one binder record.')
  interface.add_argument('--output', required=True, help='Output directory.')
  _add_fold_args(interface)
  interface.set_defaults(func=cmd_analyze_interface)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
