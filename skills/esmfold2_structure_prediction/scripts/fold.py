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

"""Predicts all-atom 3D structures from sequence with ESMFold2 (Biohub API).

Subcommands:
  fold          One protein chain            -> PDB + metrics JSON
  fold-complex  A molecular complex          -> PDB + metrics JSON (incl. iPTM)
  analyze       A metrics JSON               -> confidence report + PNG plots

pLDDT from this API is on a 0-1 scale. Everything written here keeps it on the
0-1 scale and says so (`plddt_scale`), except PDB B-factors, which are scaled to
0-100 by convention.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
#   "numpy",
#   "matplotlib",
#   "pyyaml",
# ]
# ///

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any

import matplotlib

matplotlib.use('Agg')  # Headless: no display, write PNGs only.

import matplotlib.patches as mpatches  # pylint: disable=g-import-not-at-top
import matplotlib.pyplot as plt  # pylint: disable=g-import-not-at-top
import numpy as np  # pylint: disable=g-import-not-at-top
import yaml  # pylint: disable=g-import-not-at-top

import esm_biohub as eb  # vendored, same directory
from esm_biohub import BiohubError

# The fast model was not trained with MSAs and IGNORES them silently (verified:
# it returns a normal result and no warning). Only this model actually uses one.
MSA_MODEL = 'esmfold2-2026-05'

# Confidence bands. pLDDT is 0-1 here; the conventional 0-100 cut-offs are
# 90/70/50. Colours are the AlphaFold DB convention.
PLDDT_BANDS = [
    ('very_high', 0.90, 1.01, '#0053D6', 'very high'),
    ('confident', 0.70, 0.90, '#65CBF3', 'confident'),
    ('low', 0.50, 0.70, '#FFDB13', 'low'),
    ('disordered', -0.01, 0.50, '#FF7D45', 'disordered / very low'),
]

# Runs of residues below this are reported as low-confidence regions.
LOW_CONFIDENCE_CUTOFF = 0.70
MIN_REGION_LENGTH = 3

NUCLEIC_ALPHABET = {'dna': set('ACGTN'), 'rna': set('ACGUN')}


# --------------------------------------------------------------------------
# Entity handling
# --------------------------------------------------------------------------


def _parse_id(raw: str) -> str | list[str]:
  """Parses a chain id. 'A' -> 'A'; 'A,B' -> ['A', 'B'] (one homo-oligomer)."""
  ids = [part.strip() for part in raw.split(',') if part.strip()]
  if not ids:
    raise BiohubError(f'Empty chain id in {raw!r}.')
  return ids[0] if len(ids) == 1 else ids


def _split_flag(value: str, flag: str) -> tuple[str, str]:
  """Splits an `ID:VALUE` flag on the first colon."""
  if ':' not in value:
    raise BiohubError(
        f'{flag} expects ID:VALUE (e.g. {flag} A:MKTAYIAK...), got {value!r}.'
    )
  chain, payload = value.split(':', 1)
  if not chain.strip() or not payload.strip():
    raise BiohubError(f'{flag} expects a non-empty ID and VALUE, got {value!r}.')
  return chain.strip(), payload.strip()


def _entity_from_flag(dest: str, value: str) -> dict[str, Any]:
  """Builds one entity dict from a convenience flag."""
  flag = '--' + dest.replace('_', '-')
  chain, payload = _split_flag(value, flag)
  entity_id = _parse_id(chain)

  if dest == 'protein':
    return {
        'type': 'protein',
        'id': entity_id,
        'sequence': eb.validate_sequence(payload),
    }
  if dest in ('dna', 'rna'):
    seq = ''.join(payload.split()).upper()
    bad = sorted(set(seq) - NUCLEIC_ALPHABET[dest])
    if bad:
      raise BiohubError(
          f'{flag} {chain}: not a {dest.upper()} sequence; unexpected '
          f'bases {bad}. Allowed: {"".join(sorted(NUCLEIC_ALPHABET[dest]))}.'
      )
    return {'type': dest, 'id': entity_id, 'sequence': seq}
  if dest == 'ligand_ccd':
    codes = [c.strip().upper() for c in payload.split(',') if c.strip()]
    return {'type': 'ligand', 'id': entity_id, 'ccd': codes}
  if dest == 'ligand_smiles':
    return {'type': 'ligand', 'id': entity_id, 'smiles': payload}
  raise BiohubError(f'Unknown entity flag {dest!r}.')


def _load_spec(path: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
  """Loads entities (and optional covalent bonds) from a JSON or YAML spec.

  Accepts either a bare list of entities, or a mapping with a `sequences` key
  (matching the API's own field name) plus optional `covalent_bonds`.
  """
  text = pathlib.Path(path).expanduser().read_text(encoding='utf-8')
  try:
    # YAML is a superset of JSON, so this parses both.
    data = yaml.safe_load(text)
  except yaml.YAMLError as exc:
    raise BiohubError(f'Could not parse spec {path}: {exc}') from exc

  bonds: list[dict[str, Any]] = []
  if isinstance(data, dict):
    entities = data.get('sequences') or data.get('entities')
    bonds = data.get('covalent_bonds') or []
    if entities is None:
      raise BiohubError(
          f'Spec {path} must have a `sequences` (or `entities`) list.'
      )
  elif isinstance(data, list):
    entities = data
  else:
    raise BiohubError(f'Spec {path} must be a list or a mapping.')
  if not isinstance(entities, list) or not entities:
    raise BiohubError(f'Spec {path} contains no entities.')
  return entities, bonds


def _validate_entities(entities: list[dict[str, Any]]) -> list[dict[str, Any]]:
  """Checks every entity is well formed. Raises BiohubError with the fix."""
  out = []
  for i, entity in enumerate(entities):
    if not isinstance(entity, dict):
      raise BiohubError(f'Entity {i} is not a mapping: {entity!r}')
    kind = str(entity.get('type', '')).lower()
    if kind not in ('protein', 'dna', 'rna', 'ligand'):
      raise BiohubError(
          f'Entity {i}: type must be protein|dna|rna|ligand, got '
          f'{entity.get("type")!r}.'
      )
    if not entity.get('id'):
      raise BiohubError(f'Entity {i}: missing `id` (the chain identifier).')
    if kind == 'ligand':
      if not entity.get('ccd') and not entity.get('smiles'):
        raise BiohubError(
            f'Entity {i} (ligand {entity["id"]}): needs `ccd` (e.g. '
            '["SAH"]) or `smiles`.'
        )
      if entity.get('ccd') and isinstance(entity['ccd'], str):
        entity['ccd'] = [entity['ccd']]  # The API wants a list.
    else:
      if not entity.get('sequence'):
        raise BiohubError(f'Entity {i} ({kind} {entity["id"]}): no `sequence`.')
      if kind == 'protein':
        entity['sequence'] = eb.validate_sequence(entity['sequence'])
      else:
        seq = ''.join(str(entity['sequence']).split()).upper()
        bad = sorted(set(seq) - NUCLEIC_ALPHABET[kind])
        if bad:
          raise BiohubError(
              f'Entity {i} ({kind} {entity["id"]}): unexpected bases {bad}. '
              f'Allowed: {"".join(sorted(NUCLEIC_ALPHABET[kind]))}.'
          )
        entity['sequence'] = seq
    out.append(entity)
  return out


def _entity_chain_ids(entity: dict[str, Any]) -> list[str]:
  """Chain ids an entity occupies. A homo-oligomer declares several."""
  raw = entity['id']
  return [str(c) for c in raw] if isinstance(raw, list) else [str(raw)]


# --------------------------------------------------------------------------
# Result -> metrics
# --------------------------------------------------------------------------


def _token_chain_labels(complex_data: dict[str, Any]) -> list[str]:
  """Chain label for every model *token* (the axis PAE lives on).

  Proteins and nucleic acids contribute one token per residue, but a ligand
  contributes one token per ATOM. Verified: a 76-aa protein plus a 26-atom SAH
  yields a (102, 102) PAE while `complex['sequence']` holds only 77 entries.
  """
  hetero = complex_data.get('atom_hetero') or []
  chain_lookup = (complex_data.get('metadata') or {}).get('chain_lookup', {})
  labels: list[str] = []
  for token_idx, (start, end) in enumerate(complex_data['token_to_atoms']):
    raw = complex_data['chain_id'][token_idx]
    label = str(chain_lookup.get(str(raw), raw))
    is_ligand = bool(hetero[start]) if start < len(hetero) else False
    labels.extend([label] * ((int(end) - int(start)) if is_ligand else 1))
  return labels


def _residue_chain_labels(complex_data: dict[str, Any]) -> list[str]:
  """Chain label for every entry of `complex['plddt']` (the residue axis)."""
  chain_lookup = (complex_data.get('metadata') or {}).get('chain_lookup', {})
  return [
      str(chain_lookup.get(str(raw), raw)) for raw in complex_data['chain_id']
  ]


def _plddt_stats(plddt: np.ndarray) -> dict[str, float]:
  """Summary statistics on the 0-1 scale."""
  return {
      'plddt_mean': float(np.mean(plddt)),
      'plddt_median': float(np.median(plddt)),
      'plddt_min': float(np.min(plddt)),
      'plddt_max': float(np.max(plddt)),
  }


def _save_pae(pae: np.ndarray | None, metrics_path: str) -> str | None:
  """Writes the PAE matrix beside the metrics JSON as .npy. Returns the path."""
  if pae is None:
    return None
  path = pathlib.Path(metrics_path).expanduser()
  out = path.with_name(path.stem + '_pae.npy').resolve()
  out.parent.mkdir(parents=True, exist_ok=True)
  np.save(out, pae.astype(np.float32))
  print(f'Success! PAE matrix {pae.shape} written to: {out}')
  return str(out)


def _write_pdb(text: str, output_pdb: str) -> str:
  """Writes a PDB file. Returns the absolute path, so metrics stay portable."""
  path = pathlib.Path(output_pdb).expanduser()
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(text, encoding='utf-8')
  print(f'Success! Structure written to: {path.resolve()}')
  return str(path.resolve())


# --------------------------------------------------------------------------
# fold
# --------------------------------------------------------------------------


def _read_msa(path: str, query: str) -> dict[str, Any]:
  """Reads an A3M/FASTA MSA, strips insertions, and returns the API payload.

  Lowercase columns in A3M are insertions relative to the query and are removed,
  after which every row must be the same length as the query.
  """
  records = eb.read_fasta(path)
  rows: list[str] = []
  for name, seq in records:
    row = ''.join(ch for ch in seq if not ch.islower() and ch != '.')
    if len(row) != len(query):
      raise BiohubError(
          f'MSA row {name!r} is {len(row)} columns after removing insertions '
          f'but the query is {len(query)}. The alignment must match the query.'
      )
    rows.append(row.upper())
  if rows[0].replace('-', '') != query:
    print(
        '[warning] The first MSA row is not the query sequence. ESMFold2 '
        'expects the query first; results may be degraded.',
        file=sys.stderr,
    )
  print(f'MSA: {len(rows)} sequences x {len(query)} columns from {path}')
  return {'sequences': rows}


def _resolve_msa_model(model: str | None, has_msa: bool) -> str:
  """Picks the model, refusing to let an MSA be silently discarded."""
  if not has_msa:
    return model or eb.DEFAULT_ESMFOLD2
  if model is None:
    print(
        f'[warning] An MSA was supplied, so the model was switched to '
        f'{MSA_MODEL}. The default {eb.DEFAULT_ESMFOLD2} was not trained with '
        'MSAs and would have IGNORED it silently.',
        file=sys.stderr,
    )
    return MSA_MODEL
  if model != MSA_MODEL:
    raise BiohubError(
        f'Model {model!r} does not use MSAs — it would accept your alignment '
        f'and silently ignore it, giving you a single-sequence prediction that '
        f'looks like an MSA one. Use --model {MSA_MODEL}, or drop --msa.'
    )
  return model


def cmd_fold(args) -> None:
  """Folds one protein chain."""
  if bool(args.sequence) == bool(args.fasta):
    raise BiohubError('Provide exactly one of --sequence or --fasta.')

  if args.fasta:
    records = eb.read_fasta(args.fasta)
    if args.record:
      matches = [r for r in records if r[0] == args.record]
      if not matches:
        raise BiohubError(
            f'No record {args.record!r} in {args.fasta}. Found: '
            f'{[r[0] for r in records]}'
        )
      name, sequence = matches[0]
    elif len(records) > 1:
      raise BiohubError(
          f'{args.fasta} holds {len(records)} records: {[r[0] for r in records]}'
          '. `fold` predicts one structure at a time — pass --record ID to pick'
          ' one. To fold many sequences and rank them, use the '
          '`esmfold2-binder-screening` skill.'
      )
    else:
      name, sequence = records[0]
  else:
    name, sequence = 'query', args.sequence

  sequence = eb.validate_sequence(sequence)
  msa = _read_msa(args.msa, sequence) if args.msa else None
  model = _resolve_msa_model(args.model, msa is not None)

  print(
      f'Folding {name} ({len(sequence)} aa) with {model} '
      f'[num_loops={args.num_loops}, num_sampling_steps='
      f'{args.num_sampling_steps}]...'
  )
  client = eb.BiohubClient()
  result = client.fold(
      sequence,
      model=model,
      num_loops=args.num_loops,
      num_sampling_steps=args.num_sampling_steps,
      include_pae=args.include_pae,
      msa=msa,
  )

  coordinates = eb.to_array(result['coordinates'])
  plddt = eb.to_array(result['plddt'])
  pae = eb.to_array(result['pae']) if result.get('pae') is not None else None

  pdb_file = _write_pdb(
      eb.atom37_to_pdb(coordinates, sequence, plddt), args.output_pdb
  )
  pae_file = _save_pae(pae, args.output_metrics)

  metrics: dict[str, Any] = {
      'kind': 'monomer',
      'name': name,
      'model': model,
      'num_loops': args.num_loops,
      'num_sampling_steps': args.num_sampling_steps,
      'used_msa': msa is not None,
      'sequence': sequence,
      'sequence_length': len(sequence),
      'plddt_scale': '0-1',
      'plddt': [float(x) for x in plddt],
      **_plddt_stats(plddt),
      'ptm': float(result['ptm']),
      'interface_ptm': None,  # Meaningless for a single chain.
      'chain_labels': ['A'] * len(plddt),
      'token_chain_labels': ['A'] * (pae.shape[0] if pae is not None else 0),
      'chains': [{
          'id': 'A',
          'type': 'protein',
          'length': len(sequence),
          **_plddt_stats(plddt),
      }],
      'pdb_file': pdb_file,
      'pae_file': pae_file,
  }
  eb.write_json(metrics, args.output_metrics)
  print(
      f'pLDDT {metrics["plddt_mean"]:.3f} (0-1 scale) | '
      f'pTM {metrics["ptm"]:.3f}'
  )
  print(f'Next: analyze --metrics {args.output_metrics}')


# --------------------------------------------------------------------------
# fold-complex
# --------------------------------------------------------------------------


def cmd_fold_complex(args) -> None:
  """Folds a molecular complex (proteins, DNA, RNA, ligands)."""
  flag_entities = [
      _entity_from_flag(dest, value)
      for dest, value in (getattr(args, 'entities', None) or [])
  ]
  bonds: list[dict[str, Any]] = []
  if args.spec:
    if flag_entities:
      raise BiohubError('Use --spec or the entity flags, not both.')
    entities, bonds = _load_spec(args.spec)
  else:
    entities = flag_entities
  if not entities:
    raise BiohubError(
        'No entities. Pass --spec FILE, or flags such as '
        '--protein A:MKT... --ligand-ccd L:SAH'
    )
  entities = _validate_entities(entities)

  if len(entities) == 1 and entities[0]['type'] == 'protein' and (
      len(_entity_chain_ids(entities[0])) == 1
  ):
    print(
        '[note] This is a single protein chain with no partners. `fold` is '
        'the simpler command for that; interface_ptm will be null.',
        file=sys.stderr,
    )

  described = ', '.join(
      f'{e["type"]}:{"/".join(_entity_chain_ids(e))}'
      f'({len(e["sequence"])}nt/aa)' if e['type'] != 'ligand'
      else f'ligand:{"/".join(_entity_chain_ids(e))}'
      f'({e.get("ccd") or "smiles"})'
      for e in entities
  )
  model = args.model or eb.DEFAULT_ESMFOLD2
  print(f'Folding complex [{described}] with {model}...')

  client = eb.BiohubClient()
  result = client.fold_all_atom(
      entities,
      model=model,
      num_loops=args.num_loops,
      num_sampling_steps=args.num_sampling_steps,
      include_pae=args.include_pae,
      covalent_bonds=bonds or None,
  )

  complex_data = result.get('complex')
  if not complex_data:
    raise BiohubError('The API returned no `complex` payload.')

  # `complex['plddt']` is the residue/component axis and lines up with the PDB
  # we write. The top-level `plddt`/`pae` live on the token axis, where a ligand
  # contributes one token per atom.
  plddt = eb.to_array(complex_data['plddt'])
  pae = eb.to_array(result['pae']) if result.get('pae') is not None else None
  chain_labels = _residue_chain_labels(complex_data)
  token_labels = _token_chain_labels(complex_data)

  if pae is not None and len(token_labels) != pae.shape[0]:
    print(
        f'[warning] Token/PAE axis mismatch ({len(token_labels)} vs '
        f'{pae.shape[0]}); chain boundaries will be omitted from the PAE plot.',
        file=sys.stderr,
    )
    token_labels = []

  pdb_file = _write_pdb(eb.complex_to_pdb(complex_data), args.output_pdb)
  pae_file = _save_pae(pae, args.output_metrics)

  chains = []
  for entity in entities:
    for chain_id in _entity_chain_ids(entity):
      mask = np.array([c == chain_id for c in chain_labels])
      entry: dict[str, Any] = {
          'id': chain_id,
          'type': entity['type'],
          'length': int(mask.sum()),
      }
      if entity['type'] != 'ligand':
        entry['sequence'] = entity['sequence']
      else:
        entry['ccd'] = entity.get('ccd')
        entry['smiles'] = entity.get('smiles')
      if mask.any():
        entry.update(_plddt_stats(plddt[mask]))
      chains.append(entry)

  interface_ptm = result.get('interface_ptm')
  metrics: dict[str, Any] = {
      'kind': 'complex',
      'model': model,
      'num_loops': args.num_loops,
      'num_sampling_steps': args.num_sampling_steps,
      'entities': entities,
      'covalent_bonds': bonds,
      'sequence_length': int(len(plddt)),
      'plddt_scale': '0-1',
      'plddt': [float(x) for x in plddt],
      **_plddt_stats(plddt),
      'ptm': float(result['ptm']),
      'interface_ptm': (
          float(interface_ptm) if interface_ptm is not None else None
      ),
      'chain_labels': chain_labels,
      'token_chain_labels': token_labels,
      'chains': chains,
      'pdb_file': pdb_file,
      'pae_file': pae_file,
  }
  eb.write_json(metrics, args.output_metrics)

  iptm = metrics['interface_ptm']
  iptm_text = f'{iptm:.3f}' if iptm is not None else 'n/a'
  print(
      f'pLDDT {metrics["plddt_mean"]:.3f} (0-1 scale) | '
      f'pTM {metrics["ptm"]:.3f} | iPTM {iptm_text}'
  )
  print(f'Next: analyze --metrics {args.output_metrics}')


# --------------------------------------------------------------------------
# analyze
# --------------------------------------------------------------------------


def _band_of(value: float) -> str:
  for name, low, high, _, _ in PLDDT_BANDS:
    if low <= value < high:
      return name
  return 'disordered'


def _band_breakdown(plddt: np.ndarray) -> dict[str, dict[str, float]]:
  total = len(plddt)
  out = {}
  for name, low, high, _, _ in PLDDT_BANDS:
    count = int(np.sum((plddt >= low) & (plddt < high)))
    out[name] = {
        'count': count,
        'fraction': (count / total) if total else 0.0,
    }
  return out


def _low_confidence_regions(
    plddt: np.ndarray, chain_labels: list[str]
) -> list[dict[str, Any]]:
  """Contiguous runs below the cutoff, reported as 1-based ranges per chain."""
  # Where each chain starts, so ranges can be numbered from 1 within the chain.
  chain_start: dict[str, int] = {}
  for i, chain in enumerate(chain_labels):
    chain_start.setdefault(chain, i)

  regions: list[dict[str, Any]] = []
  start = None
  for i in range(len(plddt) + 1):
    # A chain change ends the current run, as does rising above the cutoff.
    at_end = i == len(plddt)
    changed_chain = (
        not at_end
        and i > 0
        and bool(chain_labels)
        and chain_labels[i] != chain_labels[i - 1]
    )
    below = not at_end and plddt[i] < LOW_CONFIDENCE_CUTOFF
    broken = at_end or not below or changed_chain

    if broken:
      if start is not None and (i - start) >= MIN_REGION_LENGTH:
        chain = chain_labels[start] if chain_labels else 'A'
        offset = chain_start.get(chain, 0)
        window = plddt[start:i]
        regions.append({
            'chain': chain,
            'start': int(start - offset + 1),  # 1-based, within the chain.
            'end': int(i - offset),
            'length': int(i - start),
            'mean_plddt': float(window.mean()),
            'band': _band_of(float(window.mean())),
        })
      # A chain change can also immediately begin a new run.
      start = i if below else None
    elif start is None:
      start = i
  return regions


def _verdict(value: float | None, kind: str) -> str:
  """Turns pTM / iPTM into the SPEC's wording. Never leave this to the agent."""
  if value is None:
    return 'not applicable'
  if kind == 'ptm':
    if value > 0.8:
      return 'confident fold'
    if value >= 0.5:
      return 'plausible fold'
    return 'unreliable'
  if value > 0.8:
    return 'confident interface'
  if value >= 0.5:
    return 'possible interface'
  return 'likely no binding'


def _plot_plddt(
    plddt: np.ndarray, chain_labels: list[str], output: str, title: str
) -> None:
  """pLDDT vs residue, with the confidence bands shaded behind it."""
  fig, ax = plt.subplots(figsize=(11, 4))
  x = np.arange(1, len(plddt) + 1)

  for _, low, high, colour, _ in PLDDT_BANDS:
    ax.axhspan(max(low, 0.0), min(high, 1.0), color=colour, alpha=0.18, lw=0)
  ax.plot(x, plddt, color='#222222', lw=1.4, zorder=3)
  ax.axhline(LOW_CONFIDENCE_CUTOFF, color='#666666', ls='--', lw=0.8, zorder=2)

  # Chain separators.
  if chain_labels:
    boundaries = [
        i for i in range(1, len(chain_labels))
        if chain_labels[i] != chain_labels[i - 1]
    ]
    for b in boundaries:
      ax.axvline(b + 0.5, color='black', lw=1.2, zorder=4)
    edges = [0, *boundaries, len(chain_labels)]
    for left, right in zip(edges[:-1], edges[1:]):
      ax.text(
          (left + right) / 2 + 0.5,
          1.02,
          f'chain {chain_labels[left]}',
          ha='center',
          va='bottom',
          fontsize=9,
          transform=ax.get_xaxis_transform(),
      )

  ax.set_xlim(1, max(len(plddt), 2))
  ax.set_ylim(0, 1)
  ax.set_xlabel('Residue (model order)')
  ax.set_ylabel('pLDDT (0-1 scale)')
  # Padded so the per-chain labels drawn at y=1.02 do not collide with it.
  ax.set_title(title, fontweight='bold', pad=18)
  handles = [
      mpatches.Patch(
          facecolor=colour,
          alpha=0.5,
          label=f'{label} ({max(low, 0.0):.2f}-{min(high, 1.0):.2f})',
      )
      for _, low, high, colour, label in PLDDT_BANDS
  ]
  ax.legend(
      handles=handles, loc='lower right', fontsize=8, ncol=2, framealpha=0.9
  )
  fig.tight_layout()
  fig.savefig(output, dpi=150)
  plt.close(fig)
  print(f'Success! pLDDT plot written to: {output}')


def _plot_pae(
    pae: np.ndarray, token_labels: list[str], output: str, title: str
) -> None:
  """PAE heatmap. Dark = confident relative placement."""
  fig, ax = plt.subplots(figsize=(7.2, 6))
  im = ax.imshow(
      pae, cmap='Greens_r', vmin=0, vmax=30, origin='upper', aspect='equal'
  )
  if token_labels and len(token_labels) == pae.shape[0]:
    boundaries = [
        i for i in range(1, len(token_labels))
        if token_labels[i] != token_labels[i - 1]
    ]
    for b in boundaries:
      ax.axhline(b - 0.5, color='black', lw=1.2)
      ax.axvline(b - 0.5, color='black', lw=1.2)
    edges = [0, *boundaries, len(token_labels)]
    ticks, names = [], []
    for left, right in zip(edges[:-1], edges[1:]):
      ticks.append((left + right) / 2)
      names.append(token_labels[left])
    ax.set_xticks(ticks)
    ax.set_xticklabels(names)
    ax.set_yticks(ticks)
    ax.set_yticklabels(names)
  ax.set_xlabel('Scored token')
  ax.set_ylabel('Aligned token')
  ax.set_title(title, fontweight='bold')
  cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
  cbar.set_label('Expected position error (Angstrom)')
  cbar.ax.invert_yaxis()
  fig.tight_layout()
  fig.savefig(output, dpi=150)
  plt.close(fig)
  print(f'Success! PAE heatmap written to: {output}')


def _interchain_pae(
    pae: np.ndarray, token_labels: list[str]
) -> list[dict[str, Any]]:
  """Mean PAE in each off-diagonal (inter-chain) block."""
  if not token_labels or len(token_labels) != pae.shape[0]:
    return []
  labels = np.array(token_labels)
  unique = list(dict.fromkeys(token_labels))
  out = []
  for i, a in enumerate(unique):
    for b in unique[i + 1:]:
      block = pae[np.ix_(labels == a, labels == b)]
      out.append({
          'chains': f'{a}-{b}',
          'mean_pae': float(block.mean()),
          'min_pae': float(block.min()),
      })
  return out


def cmd_analyze(args) -> None:
  """Reads a metrics JSON and writes a heuristic confidence report + plots."""
  metrics = json.loads(
      pathlib.Path(args.metrics).expanduser().read_text(encoding='utf-8')
  )
  if 'plddt' not in metrics:
    raise BiohubError(
        f'{args.metrics} has no `plddt`. Pass a metrics JSON written by '
        '`fold` or `fold-complex`.'
    )
  plddt = np.asarray(metrics['plddt'], dtype=np.float64)
  if plddt.size and plddt.max() > 1.5:
    raise BiohubError(
        'This metrics file has pLDDT above 1.0, so it is not on the 0-1 scale '
        'this API returns. Refusing to analyze a rescaled file.'
    )
  chain_labels = metrics.get('chain_labels') or ['A'] * len(plddt)
  kind = metrics.get('kind', 'monomer')
  name = metrics.get('name') or kind

  bands = _band_breakdown(plddt)
  regions = _low_confidence_regions(plddt, chain_labels)
  stats = _plddt_stats(plddt)
  ptm = metrics.get('ptm')
  iptm = metrics.get('interface_ptm')

  ordered = bands['very_high']['fraction'] + bands['confident']['fraction']
  if stats['plddt_mean'] >= 0.7 and ordered >= 0.7:
    overall = 'well-folded, confidently predicted'
  elif stats['plddt_mean'] >= 0.5:
    overall = 'partially ordered; treat low-confidence regions with caution'
  else:
    overall = (
        'largely disordered or unfoldable — do not use for downstream '
        'structural analysis'
    )

  report: dict[str, Any] = {
      'source_metrics': str(pathlib.Path(args.metrics).expanduser()),
      'kind': kind,
      'model': metrics.get('model'),
      'sequence_length': metrics.get('sequence_length', int(len(plddt))),
      'plddt_scale': '0-1',
      'plddt': {
          **stats,
          'mean_on_0_100_scale': stats['plddt_mean'] * 100.0,
      },
      'confidence_bands': bands,
      'low_confidence_regions': regions,
      'ptm': {'value': ptm, 'verdict': _verdict(ptm, 'ptm')},
      'interface_ptm': {'value': iptm, 'verdict': _verdict(iptm, 'iptm')},
      'chains': metrics.get('chains', []),
      'overall': overall,
  }

  # ---- plots ----
  title = f'{name} — {metrics.get("model", "esmfold2")}'
  _plot_plddt(plddt, chain_labels, args.output_plddt_plot, f'pLDDT — {title}')

  pae_file = metrics.get('pae_file')
  if pae_file and pathlib.Path(pae_file).is_file():
    pae = np.load(pae_file)
    token_labels = metrics.get('token_chain_labels') or []
    _plot_pae(pae, token_labels, args.output_pae_plot, f'PAE — {title}')
    report['pae'] = {
        'file': pae_file,
        'shape': list(pae.shape),
        'mean_pae': float(pae.mean()),
        'interchain': _interchain_pae(pae, token_labels),
    }
  else:
    report['pae'] = None
    print(
        '[warning] No PAE matrix in this metrics file, so no heatmap was '
        'written. Re-run the fold with --include-pae to get one.',
        file=sys.stderr,
    )

  eb.write_json(report, args.output_report)

  # ---- human-readable summary ----
  print()
  print('=' * 70)
  print(f'CONFIDENCE REPORT — {name} ({report["sequence_length"]} residues)')
  print('=' * 70)
  print(
      f'  mean pLDDT     {stats["plddt_mean"]:.3f}  (0-1 scale; '
      f'{stats["plddt_mean"] * 100:.1f} on the 0-100 scale)'
  )
  print(f'  median pLDDT   {stats["plddt_median"]:.3f}')
  print('  bands:')
  for band_name, low, high, _, label in PLDDT_BANDS:
    entry = bands[band_name]
    print(
        f'    {label:22s} {low:.2f}-{min(high, 1.0):.2f}  '
        f'{entry["count"]:5d} residues  {entry["fraction"] * 100:5.1f}%'
    )
  if ptm is not None:
    print(f'  pTM            {ptm:.3f}  -> {report["ptm"]["verdict"]}')
  if iptm is not None:
    print(
        f'  iPTM           {iptm:.3f}  -> '
        f'{report["interface_ptm"]["verdict"]}'
    )
  if report.get('pae') and report['pae']['interchain']:
    for block in report['pae']['interchain']:
      print(
          f'  PAE {block["chains"]:>7s}   mean {block["mean_pae"]:.1f} A  '
          f'(min {block["min_pae"]:.1f} A)'
      )
  if regions:
    print(f'  low-confidence regions (pLDDT < {LOW_CONFIDENCE_CUTOFF:.2f}):')
    for region in regions:
      print(
          f'    chain {region["chain"]}  {region["start"]:>4d}-'
          f'{region["end"]:<4d} ({region["length"]:>3d} aa)  mean '
          f'{region["mean_plddt"]:.3f}  [{region["band"]}]'
      )
  else:
    print('  low-confidence regions: none')
  print(f'  OVERALL: {overall}')
  print('=' * 70)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


class _EntityAction(argparse.Action):
  """Collects entity flags into one list, preserving command-line order."""

  def __call__(self, parser, namespace, values, option_string=None):
    items = list(getattr(namespace, 'entities', None) or [])
    items.append((self.dest, values))
    setattr(namespace, 'entities', items)


def _add_fold_flags(parser: argparse.ArgumentParser) -> None:
  parser.add_argument(
      '--num-loops',
      type=int,
      default=20,
      help='Refinement loops. 20 = API default (accurate); 10 is ~2x faster.',
  )
  parser.add_argument(
      '--num-sampling-steps',
      type=int,
      default=100,
      help='Diffusion steps. 100 = API default; 50 is faster and usually fine.',
  )
  parser.add_argument(
      '--include-pae',
      action='store_true',
      help='Also return the PAE matrix (needed for the `analyze` heatmap).',
  )
  parser.add_argument(
      '--model',
      default=None,
      help=(
          f'ESMFold2 model. Default {eb.DEFAULT_ESMFOLD2} (fast, no MSA). '
          f'{MSA_MODEL} is slower and is the only one that uses an MSA.'
      ),
  )
  parser.add_argument('--output-pdb', required=True, help='Output .pdb path.')
  parser.add_argument(
      '--output-metrics', required=True, help='Output metrics .json path.'
  )


def main() -> None:
  parser = argparse.ArgumentParser(
      description='Predict 3D structures from sequence with ESMFold2.'
  )
  sub = parser.add_subparsers(dest='command', required=True)

  p_fold = sub.add_parser('fold', help='Fold a single protein chain.')
  p_fold.add_argument('--sequence', help='Protein sequence (one letter).')
  p_fold.add_argument('--fasta', help='FASTA file holding the sequence.')
  p_fold.add_argument('--record', help='Which FASTA record to fold, by id.')
  p_fold.add_argument(
      '--msa',
      help=(
          'A3M/FASTA alignment. Forces --model '
          f'{MSA_MODEL}; the fast model ignores MSAs silently.'
      ),
  )
  _add_fold_flags(p_fold)
  p_fold.set_defaults(func=cmd_fold)

  p_cx = sub.add_parser(
      'fold-complex',
      help='Fold a complex: proteins, DNA, RNA and/or small-molecule ligands.',
  )
  p_cx.add_argument('--spec', help='JSON/YAML file listing the entities.')
  p_cx.add_argument(
      '--protein',
      action=_EntityAction,
      metavar='ID:SEQ',
      help='Protein chain, e.g. A:MKTAYIAK. A homodimer is one entity: A,B:SEQ',
  )
  p_cx.add_argument(
      '--dna', action=_EntityAction, metavar='ID:SEQ', help='DNA chain, B:GATC'
  )
  p_cx.add_argument(
      '--rna', action=_EntityAction, metavar='ID:SEQ', help='RNA chain, C:GAUC'
  )
  p_cx.add_argument(
      '--ligand-ccd',
      action=_EntityAction,
      metavar='ID:CCD',
      help='Ligand by PDB CCD code, e.g. L:SAH',
  )
  p_cx.add_argument(
      '--ligand-smiles',
      action=_EntityAction,
      metavar='ID:SMILES',
      help='Ligand by SMILES, e.g. L:CC(=O)Oc1ccccc1C(=O)O',
  )
  _add_fold_flags(p_cx)
  p_cx.set_defaults(func=cmd_fold_complex)

  p_an = sub.add_parser(
      'analyze', help='Interpret a metrics JSON: report + pLDDT/PAE plots.'
  )
  p_an.add_argument(
      '--metrics', required=True, help='Metrics JSON from fold/fold-complex.'
  )
  p_an.add_argument(
      '--output-report', required=True, help='Output report .json path.'
  )
  p_an.add_argument(
      '--output-plddt-plot', required=True, help='Output pLDDT .png path.'
  )
  p_an.add_argument(
      '--output-pae-plot',
      required=True,
      help='Output PAE heatmap .png path (needs a fold run with --include-pae).',
  )
  p_an.set_defaults(func=cmd_analyze)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
