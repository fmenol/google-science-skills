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

"""Multimodal, promptable protein design with ESM3.

Generates protein sequences conditioned on any combination of a partial
sequence, a 3D structural motif, an SS8 secondary-structure string and a
per-residue SASA profile. Every design is folded with ESMFold2 and ranked by
pTM before it is written out: a sequence without its quality control is not a
design, it is a guess.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
#   "numpy",
# ]
# ///

from __future__ import annotations

import argparse
import concurrent.futures
import json
import pathlib
import sys
from typing import Any, Callable, Sequence

import numpy as np

from esm_biohub import BiohubClient, BiohubError  # vendored, same directory
import esm_biohub as eb

# The ONLY ESM3 model this API key can reach. `esm3-medium-*` / `esm3-large-*`,
# which the upstream cookbook tutorials use, return HTTP 403.
ESM3_MODEL = eb.DEFAULT_ESM3  # 'esm3-open-2024-03'
FOLD_MODEL = eb.DEFAULT_ESMFOLD2  # 'esmfold2-fast-2026-05'

MAX_WORKERS = 8
CA = eb.ATOM37.index('CA')

# Residue tokens that are not one of the 20 canonical amino acids. Banning them
# during sampling guarantees a synthesizable design (no 'X', no gaps).
NON_CANONICAL_IDS = [eb.VOCAB[t] for t in ('X', 'B', 'U', 'Z', 'O', '.', '-', '|')]


# --------------------------------------------------------------------------
# Core generation + QC
# --------------------------------------------------------------------------


def num_steps_for(n_masked: int) -> int:
  """Decoding steps for `n_masked` positions.

  One step per two masked positions, which the ESM3 cookbook uses throughout.
  The API rejects num_steps > sequence length and caps it at 100.
  """
  return max(1, min(n_masked // 2, 100))


def generate_sequence(
    client: BiohubClient,
    prompt: str,
    *,
    coordinates: np.ndarray | None = None,
    secondary_structure: str | None = None,
    sasa: list[float | None] | None = None,
    temperature: float,
) -> str:
  """Decodes the sequence track from a (possibly multimodal) prompt."""
  n_masked = prompt.count(eb.MASK_CHAR)
  if n_masked == 0:
    raise BiohubError(
        'The sequence prompt has no "_" positions, so there is nothing to '
        'design. Mask the positions you want ESM3 to fill in.'
    )
  out = client.generate(
      'sequence',
      ESM3_MODEL,
      sequence=prompt,
      coordinates=coordinates,
      secondary_structure=secondary_structure,
      sasa=sasa,
      num_steps=num_steps_for(n_masked),
      temperature=temperature,
      invalid_ids=NON_CANONICAL_IDS,
  )
  design = out['outputs']['sequence']
  if design is None:
    raise BiohubError('ESM3 returned no sequence for this prompt.')
  return design


def fold_and_score(client: BiohubClient, sequence: str) -> dict[str, Any]:
  """Folds one design with ESMFold2.

  Note that `fold` returns its payload at the TOP level of the response, not
  under an `outputs` key like `generate` does.
  """
  out = client.fold(sequence, FOLD_MODEL)
  plddt = eb.to_array(out['plddt'])
  return {
      'ptm': float(out['ptm']),
      'plddt_mean': float(np.nanmean(plddt)),
      'plddt': plddt,
      'coordinates': eb.to_array(out['coordinates']),
  }


def sample_designs(build: Callable[[], str], num_samples: int) -> list[str]:
  """Draws `num_samples` independent designs, in parallel.

  Design is stochastic and the spread is wide, so the caller is required to say
  how many samples it wants rather than silently taking one.
  """
  if num_samples < 1:
    raise BiohubError('--num-samples must be at least 1.')
  workers = min(MAX_WORKERS, num_samples)
  with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    return list(pool.map(lambda _: build(), range(num_samples)))


def qc_and_rank(
    client: BiohubClient,
    designs: Sequence[str],
    annotate: Callable[[dict[str, Any], str], None] | None = None,
    extras: Sequence[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
  """Folds every design, ranks by pTM (descending) and returns records.

  This is the mandatory QC step. `annotate` may attach extra per-design metrics
  (for example a motif RMSD) given the record and its folded coordinates;
  `extras[i]` is merged into the record for `designs[i]`.
  """
  workers = min(MAX_WORKERS, len(designs))
  with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    folded = list(pool.map(lambda s: fold_and_score(client, s), designs))

  records = []
  for i, (sequence, metrics) in enumerate(zip(designs, folded)):
    record: dict[str, Any] = {
        'sequence': sequence,
        'length': len(sequence),
        'ptm': metrics['ptm'],
        'plddt_mean': metrics['plddt_mean'],
        'plddt_min': float(np.nanmin(metrics['plddt'])),
        '_coordinates': metrics['coordinates'],
        '_plddt': metrics['plddt'],
    }
    if extras is not None:
      record.update(extras[i])
    if annotate is not None:
      annotate(record, sequence)
    records.append(record)

  records.sort(key=lambda r: r['ptm'], reverse=True)
  for rank, record in enumerate(records, start=1):
    record['rank'] = rank
    record['id'] = f'design_{rank}'
  return records


def write_outputs(
    records: list[dict[str, Any]],
    prefix: str,
    metadata: dict[str, Any],
) -> None:
  """Writes ranked designs as FASTA + metrics JSON + one PDB per design."""
  base = pathlib.Path(prefix).expanduser()
  base.parent.mkdir(parents=True, exist_ok=True)

  # Append, never with_suffix(): a prefix like "out/run.v2" must not become
  # "out/run.fasta".
  fasta_path = base.parent / f'{base.name}.fasta'
  json_path = base.parent / f'{base.name}.json'

  extra_keys = [
      k for k in ('motif_rmsd_ca', 'motif_verbatim')
      if any(k in r for r in records)
  ]

  with open(fasta_path, 'w', encoding='utf-8') as handle:
    for record in records:
      tags = [
          f'rank={record["rank"]}',
          f'len={record["length"]}',
          f'ptm={record["ptm"]:.3f}',
          f'plddt={record["plddt_mean"]:.3f}',
      ]
      for key in extra_keys:
        if key in record:
          value = record[key]
          tags.append(
              f'{key}={value:.3f}' if isinstance(value, float)
              else f'{key}={value}'
          )
      handle.write(f'>{record["id"]} {" ".join(tags)}\n')
      sequence = record['sequence']
      for i in range(0, len(sequence), 60):
        handle.write(sequence[i : i + 60] + '\n')

  clean: list[dict[str, Any]] = []
  for record in records:
    pdb_path = base.parent / f'{base.name}_{record["id"]}.pdb'
    pdb_path.write_text(
        eb.atom37_to_pdb(
            record['_coordinates'], record['sequence'], record['_plddt']
        ),
        encoding='utf-8',
    )
    entry = {k: v for k, v in record.items() if not k.startswith('_')}
    entry['pdb_file'] = str(pdb_path)
    clean.append(entry)

  payload = dict(metadata)
  payload['esm3_model'] = ESM3_MODEL
  payload['fold_model'] = FOLD_MODEL
  payload['num_designs'] = len(clean)
  payload['designs'] = clean

  with open(json_path, 'w', encoding='utf-8') as handle:
    json.dump(payload, handle, indent=2, default=_json_default)

  best = clean[0]
  print(f'Success! {len(clean)} design(s) written.')
  print(f'  FASTA:   {fasta_path}')
  print(f'  Metrics: {json_path}')
  print(f'  Best:    {best["id"]}  pTM={best["ptm"]:.3f}  '
        f'pLDDT={best["plddt_mean"]:.3f}'
        + (f'  motif RMSD={best["motif_rmsd_ca"]:.2f} A'
           if 'motif_rmsd_ca' in best else ''))


def _json_default(obj: Any) -> Any:
  if isinstance(obj, np.ndarray):
    return obj.tolist()
  if isinstance(obj, (np.floating, np.integer)):
    return obj.item()
  if isinstance(obj, np.bool_):
    return bool(obj)
  raise TypeError(f'Not JSON serializable: {type(obj)}')


def build_prompt(sequence: str | None, length: int | None) -> str:
  """Returns a sequence prompt: an explicit one, or `length` masked positions."""
  if sequence:
    return eb.validate_sequence(sequence, allow_mask=True)
  if not length or length < 1:
    raise BiohubError('Provide either --sequence (with "_" masks) or --length.')
  return eb.MASK_CHAR * length


# --------------------------------------------------------------------------
# Subcommands
# --------------------------------------------------------------------------


def cmd_generate(args) -> None:
  """De novo design, optionally inpainting a partial sequence."""
  client = BiohubClient()
  prompt = build_prompt(args.sequence, args.length)
  print(f'Designing {args.num_samples} sequence(s), length {len(prompt)}, '
        f'{prompt.count(eb.MASK_CHAR)} masked position(s), T={args.temperature}')

  designs = sample_designs(
      lambda: generate_sequence(client, prompt, temperature=args.temperature),
      args.num_samples,
  )
  records = qc_and_rank(client, designs)
  write_outputs(records, args.output_prefix, {
      'command': 'generate',
      'prompt': prompt,
      'length': len(prompt),
      'temperature': args.temperature,
      'num_samples': args.num_samples,
  })


def cmd_scaffold_motif(args) -> None:
  """Grafts a 3D structural motif from a PDB into a larger designed scaffold."""
  client = BiohubClient()

  try:
    lo, hi = (int(x) for x in args.motif_range.split('-'))
  except ValueError as exc:
    raise BiohubError(
        f'--motif-range must look like "124-146", got {args.motif_range!r}.'
    ) from exc

  sequence, coords, res_ids = eb.parse_pdb_atom37(args.motif_pdb, args.chain)
  picked = [i for i, r in enumerate(res_ids) if lo <= r <= hi]
  if not picked:
    raise BiohubError(
        f'No residues numbered {lo}-{hi} in {args.motif_pdb}'
        f'{f" chain {args.chain}" if args.chain else ""}. '
        f'This file spans {res_ids[0]}-{res_ids[-1]}.'
    )
  motif_seq = ''.join(sequence[i] for i in picked)
  motif_coords = coords[picked]
  motif_len = len(motif_seq)

  # --motif-start is 1-indexed, like --motif-range.
  start = args.motif_start - 1
  if start < 0 or start + motif_len > args.scaffold_length:
    raise BiohubError(
        f'A {motif_len}-residue motif starting at position {args.motif_start} '
        f'does not fit inside a {args.scaffold_length}-residue scaffold.'
    )

  prompt_chars = [eb.MASK_CHAR] * args.scaffold_length
  prompt_chars[start : start + motif_len] = list(motif_seq)
  prompt = ''.join(prompt_chars)

  coord_prompt = np.full((args.scaffold_length, 37, 3), np.nan)
  coord_prompt[start : start + motif_len] = motif_coords

  print(f'Motif: {motif_len} residues ({lo}-{hi}) = {motif_seq}')
  print(f'Scaffold: {args.scaffold_length} residues, motif at '
        f'{args.motif_start}-{args.motif_start + motif_len - 1} (1-indexed)')
  print(f'Designing {args.num_samples} scaffold(s), T={args.temperature}')

  designs = sample_designs(
      lambda: generate_sequence(
          client, prompt, coordinates=coord_prompt,
          temperature=args.temperature,
      ),
      args.num_samples,
  )

  reference_ca = motif_coords[:, CA]

  def annotate(record: dict[str, Any], sequence: str) -> None:
    grafted = sequence[start : start + motif_len]
    record['motif_verbatim'] = grafted == motif_seq
    record['motif_sequence'] = grafted
    design_ca = record['_coordinates'][start : start + motif_len, CA]
    rmsd, _ = eb.kabsch_rmsd(design_ca, reference_ca)
    record['motif_rmsd_ca'] = rmsd

  records = qc_and_rank(client, designs, annotate)
  write_outputs(records, args.output_prefix, {
      'command': 'scaffold_motif',
      'motif_pdb': args.motif_pdb,
      'motif_range': args.motif_range,
      'motif_sequence': motif_seq,
      'motif_length': motif_len,
      'motif_start': args.motif_start,
      'scaffold_length': args.scaffold_length,
      'temperature': args.temperature,
      'num_samples': args.num_samples,
  })


def cmd_design_with_ss(args) -> None:
  """Designs a sequence that adopts a requested SS8 secondary structure."""
  client = BiohubClient()
  ss8 = args.ss8.strip().upper()
  bad = sorted(set(ss8) - set(eb.SS8_VOCAB))
  if bad:
    raise BiohubError(
        f'--ss8 contains characters outside the SS8 alphabet '
        f'{eb.SS8_VOCAB}: {bad}'
    )

  prompt = build_prompt(args.sequence, len(ss8))
  if len(prompt) != len(ss8):
    raise BiohubError(
        f'--sequence has length {len(prompt)} but --ss8 has length {len(ss8)}. '
        'They must match: the design is exactly as long as the SS8 string.'
    )

  print(f'Designing {args.num_samples} sequence(s) for a {len(ss8)}-residue '
        f'SS8 target, T={args.temperature}')
  print(f'  SS8: {ss8}')

  designs = sample_designs(
      lambda: generate_sequence(
          client, prompt, secondary_structure=ss8,
          temperature=args.temperature,
      ),
      args.num_samples,
  )
  records = qc_and_rank(client, designs)
  write_outputs(records, args.output_prefix, {
      'command': 'design_with_ss',
      'ss8': ss8,
      'length': len(ss8),
      'temperature': args.temperature,
      'num_samples': args.num_samples,
  })


def parse_sasa(spec: str | None, path: str | None, length: int) -> list:
  """Builds a length-L SASA prompt; None marks an unconstrained position.

  `spec` is "start-end:value,..." with 1-indexed inclusive residue ranges, in
  square Angstroms. `path` is a JSON list of length L allowing nulls.
  """
  if path:
    with open(path, encoding='utf-8') as handle:
      values = json.load(handle)
    if not isinstance(values, list) or len(values) != length:
      raise BiohubError(
          f'--sasa-json must hold a list of exactly {length} entries '
          '(numbers, or null for unconstrained).'
      )
    return [None if v is None else float(v) for v in values]

  sasa: list[float | None] = [None] * length
  if not spec:
    raise BiohubError('Provide --sasa or --sasa-json.')
  for chunk in spec.split(','):
    chunk = chunk.strip()
    if not chunk:
      continue
    try:
      span, raw = chunk.split(':')
      value = float(raw)
      lo, hi = (int(x) for x in span.split('-')) if '-' in span else (
          int(span), int(span)
      )
    except ValueError as exc:
      raise BiohubError(
          f'Cannot parse --sasa entry {chunk!r}. Use "20-30:5.0,45:80.0" '
          '(1-indexed inclusive residue ranges : SASA in A^2).'
      ) from exc
    if lo < 1 or hi > length or lo > hi:
      raise BiohubError(
          f'--sasa range {lo}-{hi} is outside 1-{length}.'
      )
    for i in range(lo - 1, hi):
      sasa[i] = value
  return sasa


def cmd_design_with_sasa(args) -> None:
  """Designs a sequence with a requested solvent-accessibility profile."""
  client = BiohubClient()
  prompt = build_prompt(args.sequence, args.length)
  sasa = parse_sasa(args.sasa, args.sasa_json, len(prompt))
  n_constrained = sum(v is not None for v in sasa)
  if not n_constrained:
    raise BiohubError('Every SASA entry is unconstrained; nothing to condition on.')

  print(f'Designing {args.num_samples} sequence(s), length {len(prompt)}, '
        f'{n_constrained} SASA-constrained position(s), T={args.temperature}')

  designs = sample_designs(
      lambda: generate_sequence(
          client, prompt, sasa=sasa, temperature=args.temperature
      ),
      args.num_samples,
  )
  records = qc_and_rank(client, designs)
  write_outputs(records, args.output_prefix, {
      'command': 'design_with_sasa',
      'length': len(prompt),
      'sasa': sasa,
      'num_constrained': n_constrained,
      'temperature': args.temperature,
      'num_samples': args.num_samples,
  })


def cmd_chain_of_thought(args) -> None:
  """Decodes tracks one after another, each conditioned on all the previous.

  The default order (secondary_structure -> structure -> sequence) lets ESM3
  commit to a fold plan before it writes any residues. A single spot check at
  length 60 gave pTM 0.65 against a mean of 0.32 for direct sequence sampling
  (n=8), so this is promising but not established -- draw several samples and
  trust the measured pTM, not the anecdote.
  """
  client = BiohubClient()
  tracks = [t.strip() for t in args.tracks.split(',') if t.strip()]
  if tracks[-1] != 'sequence':
    raise BiohubError(
        f'--tracks must end with "sequence" (got {tracks[-1]!r}); the point of '
        'the chain is to arrive at a designed sequence.'
    )
  known = {'sequence', 'structure', 'secondary_structure', 'sasa'}
  bad = [t for t in tracks if t not in known]
  if bad:
    raise BiohubError(f'Unknown track(s): {bad}. Known: {sorted(known)}')

  prompt = build_prompt(args.sequence, args.length)
  length = len(prompt)
  print(f'Designing {args.num_samples} sequence(s) of length {length} via '
        f'{" -> ".join(tracks)}, T={args.temperature}')

  def one_chain() -> tuple[str, dict[str, Any]]:
    state: dict[str, Any] = {
        'sequence': prompt,
        'coordinates': None,
        'secondary_structure': None,
        'sasa': None,
    }
    trace: dict[str, Any] = {}
    for track in tracks:
      # Intermediate tracks get a full decoding pass (one token per step);
      # the sequence track keeps the standard n_masked // 2 schedule.
      steps = (
          num_steps_for(prompt.count(eb.MASK_CHAR))
          if track == 'sequence'
          else max(1, min(length, 100))
      )
      out = client.generate(
          track,
          ESM3_MODEL,
          sequence=state['sequence'],
          coordinates=state['coordinates'],
          secondary_structure=state['secondary_structure'],
          sasa=state['sasa'],
          num_steps=steps,
          temperature=args.temperature,
          invalid_ids=NON_CANONICAL_IDS if track == 'sequence' else (),
      )['outputs']

      if track == 'structure':
        state['coordinates'] = eb.to_array(out['coordinates'])
        trace['structure_ptm'] = out.get('ptm')
      elif track == 'sequence':
        state['sequence'] = out['sequence']
      else:
        state[track] = out[track]
        trace[track] = out[track]

    return state['sequence'], trace

  workers = min(MAX_WORKERS, args.num_samples)
  with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    chains = list(pool.map(lambda _: one_chain(), range(args.num_samples)))

  designs = [seq for seq, _ in chains]
  traces = [trace for _, trace in chains]
  records = qc_and_rank(client, designs, extras=traces)

  write_outputs(records, args.output_prefix, {
      'command': 'chain_of_thought',
      'tracks': tracks,
      'length': length,
      'temperature': args.temperature,
      'num_samples': args.num_samples,
  })


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main() -> None:
  parser = argparse.ArgumentParser(
      description=__doc__,
      formatter_class=argparse.RawDescriptionHelpFormatter,
  )
  sub = parser.add_subparsers(dest='command', required=True)

  def common(p, *, temperature=0.5):
    p.add_argument(
        '--num-samples', type=int, required=True,
        help='Number of independent designs to draw. Design is stochastic: '
             'draw several and keep the best by pTM.',
    )
    p.add_argument(
        '--temperature', type=float, default=temperature,
        help=f'Sampling temperature (default {temperature}). Higher is more '
             'diverse and less foldable; > 0.7 degrades quality sharply.',
    )
    p.add_argument(
        '--output-prefix', required=True,
        help='Output path prefix; writes <prefix>.fasta, <prefix>.json and '
             'one PDB per design.',
    )

  p = sub.add_parser(
      'generate', help='De novo design, or inpaint a partial sequence.')
  p.add_argument('--length', type=int, help='Length of a fully de novo design.')
  p.add_argument(
      '--sequence',
      help='Partial sequence prompt; "_" marks positions to design. Use '
           'instead of --length to inpaint a scaffold you already have.',
  )
  common(p)
  p.set_defaults(func=cmd_generate)

  p = sub.add_parser(
      'scaffold-motif',
      help='Graft a 3D structural motif from a PDB into a designed scaffold.')
  p.add_argument('--motif-pdb', required=True, help='Path to a .pdb/.cif file.')
  p.add_argument(
      '--motif-range', required=True,
      help='Author residue numbers of the motif, 1-indexed inclusive, e.g. '
           '"124-146" (the numbers you read off in PyMOL).',
  )
  p.add_argument('--chain', default=None, help='Chain id (default: first).')
  p.add_argument(
      '--scaffold-length', type=int, required=True,
      help='Total length of the designed protein.',
  )
  p.add_argument(
      '--motif-start', type=int, required=True,
      help='Position in the scaffold where the motif begins, 1-indexed.',
  )
  common(p)
  p.set_defaults(func=cmd_scaffold_motif)

  p = sub.add_parser(
      'design-with-ss',
      help='Design a sequence that folds to a requested SS8 string.')
  p.add_argument(
      '--ss8', required=True,
      help=f'SS8 secondary structure over the alphabet {eb.SS8_VOCAB} '
           '(H=alpha helix, E=beta strand, C=coil, T=turn, G/I=3-10/pi helix, '
           'B=bridge, S=bend). The design is exactly this long.',
  )
  p.add_argument(
      '--sequence',
      help='Optional partial sequence of the same length; "_" marks positions '
           'to design.',
  )
  common(p)
  p.set_defaults(func=cmd_design_with_ss)

  p = sub.add_parser(
      'design-with-sasa',
      help='Design a sequence with a requested burial / exposure profile.')
  p.add_argument('--length', type=int, help='Length of the design.')
  p.add_argument('--sequence', help='Partial sequence prompt; "_" masks.')
  p.add_argument(
      '--sasa',
      help='Per-residue SASA in A^2 as "start-end:value,...", 1-indexed '
           'inclusive; unlisted positions are unconstrained. Buried < 20, '
           'exposed > 50. Example: "20-30:5,45-50:90".',
  )
  p.add_argument(
      '--sasa-json',
      help='Path to a JSON list of exactly --length entries; null means '
           'unconstrained.',
  )
  common(p)
  p.set_defaults(func=cmd_design_with_sasa)

  p = sub.add_parser(
      'chain-of-thought',
      help='Decode tracks one at a time, each conditioned on the last.')
  p.add_argument('--length', type=int, help='Length of the design.')
  p.add_argument('--sequence', help='Partial sequence prompt; "_" masks.')
  p.add_argument(
      '--tracks', default='secondary_structure,structure,sequence',
      help='Comma-separated decoding order. Must end with "sequence". '
           'Default: secondary_structure,structure,sequence.',
  )
  common(p)
  p.set_defaults(func=cmd_chain_of_thought)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
