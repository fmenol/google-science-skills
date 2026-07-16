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

"""Predicts functional domain annotations from sequence alone with ESM3.

Uses the ESM3 `function` track, which decodes InterPro-style annotations
directly from the residue sequence with no homology search and no database hit.
Every annotation is a MODEL PREDICTION, not a curated fact.

Annotation ranges are 1-indexed and inclusive: length = end - start + 1.
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
import re
import sys

from esm_biohub import BiohubClient, BiohubError
import esm_biohub as eb

# InterPro accessions are 'IPR' + exactly 6 digits, e.g. IPR000626. The function
# track emits them inside the label, as 'Ubiquitin-like domain (IPR000626)'.
INTERPRO_RE = re.compile(r'\bIPR\d{6}\b')
INTERPRO_SUFFIX_RE = re.compile(r'\s*\(\s*IPR\d{6}\s*\)\s*')

MAX_WORKERS = 8  # esm_common/SPEC.md: keep batch concurrency <= 8.
MAX_NUM_STEPS = 100  # Server-side cap; num_steps must also be <= len(sequence).
DEFAULT_NUM_STEPS = 8
DEFAULT_TEMPERATURE = 1.0

CAVEAT = (
    'These are ESM3 MODEL PREDICTIONS decoded from sequence alone, not curated '
    'annotations. No homology search was performed and no database was '
    'consulted. Always present them as predictions, and resolve any predicted '
    'InterPro accession against the authoritative entry using the '
    '`interpro_database` skill before relying on it.'
)


# --------------------------------------------------------------------------
# Annotation parsing
# --------------------------------------------------------------------------


def parse_annotation(raw: list) -> dict:
  """Parses one raw `[label, start, end]` triple into a structured annotation.

  The function track returns a mix of free-text InterPro keywords ('ubiquitin
  domain', 'core', 'site') and full InterPro entries that carry an accession
  ('Ubiquitin-like domain (IPR000626)'). We pull the accession out so the agent
  can hand it to the `interpro_database` skill.

  Args:
    raw: A `[label, start, end]` triple. `start`/`end` are 1-indexed inclusive.

  Returns:
    A dict with label, name, start, end, length, interpro_id and kind.
  """
  label = str(raw[0])
  start, end = int(raw[1]), int(raw[2])
  match = INTERPRO_RE.search(label)
  interpro_id = match.group(0) if match else None
  name = INTERPRO_SUFFIX_RE.sub(' ', label).strip() if interpro_id else label
  return {
      'label': label,
      'name': name,
      'start': start,
      'end': end,
      'length': end - start + 1,  # 1-indexed, inclusive.
      'interpro_id': interpro_id,
      'kind': 'interpro_entry' if interpro_id else 'keyword',
  }


def merge_spans(annotations: list[dict], length: int) -> list[dict]:
  """Merges annotation intervals into the non-overlapping domain architecture.

  Overlapping *and* directly adjacent intervals are merged, so each output span
  is a maximal run of consecutively annotated residues.

  Args:
    annotations: Parsed annotations.
    length: Sequence length, used to clamp any out-of-range range.

  Returns:
    Non-overlapping spans sorted by position, each listing the labels that
    contributed to it.
  """
  if not annotations or length <= 0:
    return []

  intervals = []
  for ann in annotations:
    start = max(1, min(int(ann['start']), length))
    end = max(1, min(int(ann['end']), length))
    if start > end:
      start, end = end, start
    intervals.append((start, end))
  intervals.sort()

  merged: list[list[int]] = []
  cur_start, cur_end = intervals[0]
  for start, end in intervals[1:]:
    if start <= cur_end + 1:  # Overlapping or directly adjacent.
      cur_end = max(cur_end, end)
    else:
      merged.append([cur_start, cur_end])
      cur_start, cur_end = start, end
  merged.append([cur_start, cur_end])

  spans = []
  for start, end in merged:
    inside = [
        a for a in annotations if a['start'] <= end and a['end'] >= start
    ]
    inside.sort(key=lambda a: (-a['length'], a['label']))
    spans.append({
        'start': start,
        'end': end,
        'length': end - start + 1,
        'n_annotations': len(inside),
        'interpro_ids': sorted(
            {a['interpro_id'] for a in inside if a['interpro_id']}
        ),
        'labels': [a['label'] for a in inside],
    })
  return spans


def effective_num_steps(num_steps: int, length: int, seq_id: str) -> int:
  """Clamps num_steps into the range the API accepts: 1 .. min(length, 100).

  The server rejects anything outside this with HTTP 422. Clamping (rather than
  erroring) keeps one short record from killing an entire batch.
  """
  ceiling = max(1, min(length, MAX_NUM_STEPS))
  steps = max(1, min(int(num_steps), ceiling))
  if steps != int(num_steps):
    print(
        f'[note] {seq_id}: num_steps {num_steps} -> {steps} '
        f'(must be 1..min(len={length}, {MAX_NUM_STEPS})).',
        file=sys.stderr,
    )
  return steps


def check_temperature(temperature: float) -> None:
  """Guards the two degenerate ends of the temperature range.

  Verified against the live API:
    * temperature = 0 returns ZERO annotations even for real ubiquitin. The
      function track collapses onto the null token. Silently reads as 'this
      protein has no function', so we refuse it outright.
    * temperature > 1.0 confabulates. At 1.5 a *scrambled* ubiquitin produced 36
      annotations including an invented 'NAD(P)-binding domain superfamily
      (IPR036291)'.
  """
  if temperature <= 0:
    raise BiohubError(
        'temperature must be > 0. The ESM3 function track collapses to ZERO '
        'annotations at temperature=0 (verified on real ubiquitin), which is '
        'indistinguishable from "no function found". Use 0.1-1.0; '
        f'the default is {DEFAULT_TEMPERATURE}.'
    )
  if temperature > 1.0:
    print(
        f'[warning] temperature={temperature} > 1.0. The function track '
        'confabulates above 1.0: a scrambled sequence yielded 36 annotations '
        'and an invented InterPro accession at temperature=1.5. Treat every '
        'annotation from this run as unreliable.',
        file=sys.stderr,
    )


# --------------------------------------------------------------------------
# Prediction
# --------------------------------------------------------------------------


def predict_one(
    client: BiohubClient,
    seq_id: str,
    sequence: str,
    model: str,
    num_steps: int,
    temperature: float,
) -> dict:
  """Predicts function annotations for one sequence."""
  sequence = eb.validate_sequence(sequence)
  length = len(sequence)
  steps = effective_num_steps(num_steps, length, seq_id)

  result = client.generate(
      'function',
      model=model,
      sequence=sequence,
      num_steps=steps,
      temperature=temperature,
  )
  raw = (result.get('outputs') or {}).get('function') or []

  annotations = [parse_annotation(a) for a in raw]
  annotations.sort(key=lambda a: (a['start'], -a['length'], a['label']))
  spans = merge_spans(annotations, length)
  covered = sum(s['length'] for s in spans)

  return {
      'id': seq_id,
      'sequence': sequence,
      'length': length,
      'model': model,
      'num_steps': steps,
      'temperature': temperature,
      'n_annotations': len(annotations),
      'annotations': annotations,
      'interpro_ids': sorted(
          {a['interpro_id'] for a in annotations if a['interpro_id']}
      ),
      'domain_architecture': spans,
      'covered_residues': covered,
      'coverage': round(covered / length, 4) if length else 0.0,
      'caveat': CAVEAT,
  }


def _resolve_input(args) -> list[tuple[str, str]]:
  """Returns [(id, sequence)] from --sequence or --fasta."""
  if bool(args.sequence) == bool(args.fasta):
    raise BiohubError('Provide exactly one of --sequence or --fasta.')
  if args.sequence:
    return [(args.id or 'query', args.sequence)]
  records = eb.read_fasta(args.fasta)
  if len(records) > 1:
    print(
        f'[note] {args.fasta} holds {len(records)} records; `predict` uses the '
        'first. Use the `batch` subcommand for all of them.',
        file=sys.stderr,
    )
  return records[:1]


def cmd_predict(args) -> None:
  """Predicts annotations for a single sequence."""
  check_temperature(args.temperature)
  (seq_id, sequence), = _resolve_input(args)

  client = BiohubClient()
  record = predict_one(
      client, seq_id, sequence, args.model, args.num_steps, args.temperature
  )
  eb.write_json(record, args.output)

  arch = ', '.join(
      f'{s["start"]}-{s["end"]}' for s in record['domain_architecture']
  ) or 'none'
  print(
      f'{record["id"]}: {record["n_annotations"]} predicted annotation(s), '
      f'coverage {record["coverage"]:.1%} of {record["length"]} aa'
  )
  print(f'  domain architecture: {arch}')
  print(f'  InterPro accessions: {record["interpro_ids"] or "none"}')
  if not record['n_annotations']:
    print(
        '  No annotation returned. ESM3 does not recognise this sequence: it '
        'may be a genuine novel/dark-matter protein, or not a real protein.'
    )
  print(f'  NOTE: {CAVEAT}')


def cmd_batch(args) -> None:
  """Predicts annotations for every record in a FASTA file, concurrently."""
  check_temperature(args.temperature)
  records = eb.read_fasta(args.fasta)
  workers = max(1, min(args.max_workers, MAX_WORKERS))
  client = BiohubClient()

  # Keyed by index, not id: a FASTA may legally repeat an id, and keying by it
  # would silently drop records.
  results: dict[int, dict] = {}
  with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    futures = {
        pool.submit(
            predict_one, client, seq_id, seq, args.model, args.num_steps,
            args.temperature,
        ): index
        for index, (seq_id, seq) in enumerate(records)
    }
    for future in concurrent.futures.as_completed(futures):
      index = futures[future]
      seq_id = records[index][0]
      try:
        results[index] = future.result()
      except Exception as exc:  # noqa: BLE001  One bad record must not sink the batch.
        print(f'[error] {seq_id}: {exc}', file=sys.stderr)
        results[index] = {'id': seq_id, 'error': str(exc)}

  ordered = [results[index] for index in range(len(records))]
  payload = {
      'model': args.model,
      'num_steps': args.num_steps,
      'temperature': args.temperature,
      'n_records': len(ordered),
      'n_failed': sum(1 for r in ordered if 'error' in r),
      'records': ordered,
      'caveat': CAVEAT,
  }
  eb.write_json(payload, args.output)

  if args.summary:
    write_summary_csv(ordered, args.summary)

  ok = [r for r in ordered if 'error' not in r]
  print(f'Predicted function for {len(ok)}/{len(ordered)} record(s).')
  for record in ordered:
    if 'error' in record:
      print(f'  {record["id"]:<20s} FAILED: {record["error"][:60]}')
      continue
    print(
        f'  {record["id"]:<20s} {record["length"]:>5d} aa  '
        f'{record["n_annotations"]:>3d} annotation(s)  '
        f'coverage {record["coverage"]:>6.1%}  '
        f'IPR: {",".join(record["interpro_ids"]) or "none"}'
    )
  print(f'NOTE: {CAVEAT}')


def write_summary_csv(records: list[dict], output_file: str) -> None:
  """Writes the one-row-per-record summary table."""
  path = pathlib.Path(output_file).expanduser()
  path.parent.mkdir(parents=True, exist_ok=True)
  columns = [
      'id', 'length', 'n_annotations', 'covered_residues', 'coverage',
      'n_interpro', 'interpro_ids', 'domain_architecture', 'top_labels',
      'error',
  ]
  with open(path, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    for record in records:
      if 'error' in record:
        writer.writerow({'id': record['id'], 'error': record['error']})
        continue
      top = sorted(
          record['annotations'], key=lambda a: -a['length']
      )[:3]
      writer.writerow({
          'id': record['id'],
          'length': record['length'],
          'n_annotations': record['n_annotations'],
          'covered_residues': record['covered_residues'],
          'coverage': f'{record["coverage"]:.4f}',
          'n_interpro': len(record['interpro_ids']),
          'interpro_ids': ';'.join(record['interpro_ids']),
          'domain_architecture': ';'.join(
              f'{s["start"]}-{s["end"]}' for s in record['domain_architecture']
          ),
          'top_labels': ';'.join(a['label'] for a in top),
          'error': '',
      })
  print(f'Success! Data written to: {path}')


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------


def cmd_plot(args) -> None:
  """Draws a domain-architecture diagram from a `predict` or `batch` JSON."""
  import matplotlib
  matplotlib.use('Agg')  # Headless; must precede the pyplot import.
  import matplotlib.patches as mpatches
  import matplotlib.pyplot as plt

  with open(args.input, encoding='utf-8') as handle:
    data = json.load(handle)

  if 'records' in data:  # A `batch` payload.
    candidates = [r for r in data['records'] if 'error' not in r]
    if not candidates:
      raise BiohubError(f'{args.input} holds no successful records to plot.')
    if args.record:
      matches = [r for r in candidates if r['id'] == args.record]
      if not matches:
        raise BiohubError(
            f'No record {args.record!r} in {args.input}. Available: '
            f'{[r["id"] for r in candidates]}'
        )
      record = matches[0]
    else:
      record = candidates[0]
  else:
    record = data

  annotations = record.get('annotations') or []
  length = int(record['length'])

  # Longest bars at the top: read like a domain diagram, general -> specific.
  ordered = sorted(annotations, key=lambda a: (-a['length'], a['start']))
  n_rows = max(len(ordered), 1)

  fig_height = max(2.4, 0.32 * n_rows + 1.6)
  fig, axis = plt.subplots(figsize=(12, fig_height))

  interpro_color = '#1f77b4'
  keyword_color = '#a6c8e6'

  # The sequence backbone.
  axis.add_patch(
      mpatches.Rectangle(
          (1, n_rows + 0.15), length, 0.35,
          facecolor='#d9d9d9', edgecolor='#8c8c8c', linewidth=0.6,
      )
  )
  axis.text(
      length / 2, n_rows + 0.33, f'{record.get("id", "sequence")} ({length} aa)',
      ha='center', va='center', fontsize=8, color='#333333',
  )

  for row, ann in enumerate(ordered):
    y = n_rows - 1 - row
    is_ipr = ann['interpro_id'] is not None
    axis.add_patch(
        mpatches.Rectangle(
            (ann['start'], y + 0.15), ann['length'], 0.7,
            facecolor=interpro_color if is_ipr else keyword_color,
            edgecolor='#404040', linewidth=0.5, alpha=0.95 if is_ipr else 0.8,
        )
    )
    axis.text(
        ann['start'] + ann['length'] / 2, y + 0.5,
        f'{ann["label"]}  [{ann["start"]}-{ann["end"]}]',
        ha='center', va='center', fontsize=7,
        color='white' if is_ipr else '#1a1a1a',
        clip_on=True,
    )

  axis.set_xlim(0, length + 1)
  axis.set_ylim(-0.2, n_rows + 0.8)
  axis.set_yticks([])
  axis.set_xlabel('Residue position (1-indexed, inclusive)')
  coverage = record.get('coverage', 0.0)
  axis.set_title(
      f'ESM3 PREDICTED domain architecture — {record.get("id", "sequence")}\n'
      f'{len(annotations)} predicted annotation(s), '
      f'{coverage:.1%} residue coverage — predictions, NOT curated annotations',
      fontsize=10,
  )
  # Below the axes: an in-axes legend would sit on top of the backbone bar.
  axis.legend(
      handles=[
          mpatches.Patch(
              facecolor=interpro_color,
              label='InterPro entry (accession -> verify with interpro_database)',
          ),
          mpatches.Patch(facecolor=keyword_color, label='Free-text keyword'),
      ],
      loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=2,
      fontsize=7, frameon=False,
  )
  for spine in ('top', 'right', 'left'):
    axis.spines[spine].set_visible(False)

  fig.tight_layout()
  path = pathlib.Path(args.output).expanduser()
  path.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
  plt.close(fig)
  print(f'Success! Data written to: {path}')


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          'Predict InterPro-style function annotations from sequence alone '
          'with the ESM3 function track. Output ranges are 1-indexed '
          'inclusive. Predictions, not curated annotations.'
      )
  )
  sub = parser.add_subparsers(dest='command', required=True)

  def add_model_flags(sub_parser) -> None:
    sub_parser.add_argument(
        '--model', default=eb.DEFAULT_ESM3,
        help=f'ESM3 model (default: {eb.DEFAULT_ESM3}; the only one reachable).',
    )
    sub_parser.add_argument(
        '--num-steps', type=int, default=DEFAULT_NUM_STEPS,
        help=(
            f'Decoding steps, 1..min(len, {MAX_NUM_STEPS}) '
            f'(default: {DEFAULT_NUM_STEPS}). More steps = more annotations.'
        ),
    )
    sub_parser.add_argument(
        '--temperature', type=float, default=DEFAULT_TEMPERATURE,
        help=(
            f'Sampling temperature (default: {DEFAULT_TEMPERATURE}). Use '
            '0.1-1.0. 0 returns NO annotations; >1.0 confabulates.'
        ),
    )

  p = sub.add_parser(
      'predict', help='Predict annotations for one sequence.'
  )
  p.add_argument('--sequence', help='Protein sequence (one-letter).')
  p.add_argument('--fasta', help='FASTA file; the first record is used.')
  p.add_argument('--id', help='Name for the sequence (with --sequence).')
  p.add_argument('--output', required=True, help='Output JSON file path.')
  add_model_flags(p)
  p.set_defaults(func=cmd_predict)

  p = sub.add_parser(
      'batch', help='Predict annotations for every record in a FASTA.'
  )
  p.add_argument('--fasta', required=True, help='Input FASTA file.')
  p.add_argument('--output', required=True, help='Output JSON file path.')
  p.add_argument('--summary', help='Optional summary table (CSV) path.')
  p.add_argument(
      '--max-workers', type=int, default=MAX_WORKERS,
      help=f'Concurrent requests, capped at {MAX_WORKERS}.',
  )
  add_model_flags(p)
  p.set_defaults(func=cmd_batch)

  p = sub.add_parser(
      'plot', help='Draw a domain-architecture diagram from a predict/batch JSON.'
  )
  p.add_argument(
      '--input', required=True, help='JSON written by `predict` or `batch`.'
  )
  p.add_argument('--output', required=True, help='Output PNG file path.')
  p.add_argument(
      '--record', help='With a `batch` JSON: which record id to draw.'
  )
  p.add_argument('--dpi', type=int, default=150, help='Figure DPI.')
  p.set_defaults(func=cmd_plot)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
