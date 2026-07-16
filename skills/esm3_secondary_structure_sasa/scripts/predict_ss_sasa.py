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

"""Predict per-residue secondary structure (SS8) and solvent accessibility
(SASA) directly from a protein sequence with ESM3 — no folding, no DSSP."""

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
import sys

import numpy as np

from esm_biohub import BiohubClient, BiohubError
import esm_biohub as eb

# DSSP 8-state secondary-structure alphabet, as ESM3 emits it.
SS8_MEANING = {
    'H': 'alpha helix',
    'G': '3-10 helix',
    'I': 'pi helix',
    'E': 'beta strand',
    'B': 'beta bridge',
    'T': 'turn',
    'S': 'bend',
    'C': 'coil / loop',
}
# Coarse 3-state collapse (helix / strand / coil).
SS8_TO_SS3 = {
    'H': 'H', 'G': 'H', 'I': 'H',
    'E': 'E', 'B': 'E',
    'T': 'C', 'S': 'C', 'C': 'C',
}

# Burial bands on absolute SASA (Angstrom^2). Absolute thresholds are coarse —
# a large residue exposes more area than a small one at the same relative
# burial — so these are reported as a convenience, not a substitute for
# relative accessibility.
BURIED_MAX = 20.0
EXPOSED_MIN = 50.0

MAX_WORKERS = 8


def _num_steps(length: int) -> int:
  """Iterative-decoding steps. The API caps num_steps at min(length, 100)."""
  return max(1, min(length // 2, 100))


def predict_one(
    client: BiohubClient,
    sequence: str,
    model: str,
    temperature: float,
) -> dict:
  """Predicts SS8 and SASA for one sequence.

  Returns:
    A dict with the sequence, the SS8 string, the SS3 collapse, per-residue
    SASA, burial labels, and summary statistics.
  """
  sequence = eb.validate_sequence(sequence)
  length = len(sequence)
  steps = _num_steps(length)

  ss8 = client.generate(
      'secondary_structure', model, sequence=sequence,
      num_steps=steps, temperature=temperature,
  )['outputs']['secondary_structure']
  sasa = client.generate(
      'sasa', model, sequence=sequence,
      num_steps=steps, temperature=temperature,
  )['outputs']['sasa']

  if ss8 is None or len(ss8) != length:
    raise BiohubError(
        f'ESM3 returned a secondary-structure string of length '
        f'{0 if ss8 is None else len(ss8)} for a {length}-residue sequence.'
    )
  if sasa is None or len(sasa) != length:
    raise BiohubError(
        f'ESM3 returned {0 if sasa is None else len(sasa)} SASA values for a '
        f'{length}-residue sequence.'
    )

  sasa_arr = np.array([float(x) for x in sasa], dtype=np.float64)
  ss3 = ''.join(SS8_TO_SS3.get(c, 'C') for c in ss8)
  burial = [
      'buried' if v <= BURIED_MAX else 'exposed' if v >= EXPOSED_MIN
      else 'intermediate'
      for v in sasa_arr
  ]

  def _frac(alphabet: str, string: str) -> dict:
    return {c: round(string.count(c) / len(string), 4) for c in alphabet
            if c in string}

  return {
      'sequence': sequence,
      'length': length,
      'secondary_structure_ss8': ss8,
      'secondary_structure_ss3': ss3,
      'sasa': [round(float(v), 3) for v in sasa_arr],
      'burial': burial,
      'summary': {
          'ss8_composition': _frac('HGIEBTSC', ss8),
          'ss3_composition': _frac('HEC', ss3),
          'percent_helix': round(100 * ss3.count('H') / length, 1),
          'percent_strand': round(100 * ss3.count('E') / length, 1),
          'percent_coil': round(100 * ss3.count('C') / length, 1),
          'sasa_mean': round(float(sasa_arr.mean()), 2),
          'sasa_min': round(float(sasa_arr.min()), 2),
          'sasa_max': round(float(sasa_arr.max()), 2),
          'n_buried': int(sum(b == 'buried' for b in burial)),
          'n_exposed': int(sum(b == 'exposed' for b in burial)),
      },
  }


def cmd_predict(args) -> None:
  client = BiohubClient()
  if bool(args.sequence) == bool(args.fasta):
    raise BiohubError('Provide exactly one of --sequence or --fasta.')

  if args.sequence:
    records = [('query', args.sequence)]
  else:
    records = eb.read_fasta(args.fasta)

  results = []
  for name, seq in records:
    result = predict_one(client, seq, args.model, args.temperature)
    result['id'] = name
    results.append(result)
    s = result['summary']
    print(
        f'{name}: len={result["length"]} '
        f'helix={s["percent_helix"]}% strand={s["percent_strand"]}% '
        f'coil={s["percent_coil"]}% | SASA mean={s["sasa_mean"]} '
        f'buried={s["n_buried"]} exposed={s["n_exposed"]}',
        file=sys.stderr,
    )

  eb.write_json(results if len(results) > 1 else results[0], args.output)


def cmd_batch(args) -> None:
  client = BiohubClient()
  records = eb.read_fasta(args.fasta)

  def _run(item):
    name, seq = item
    result = predict_one(client, seq, args.model, args.temperature)
    result['id'] = name
    return result

  results: list[dict] = [None] * len(records)  # type: ignore[list-item]
  with concurrent.futures.ThreadPoolExecutor(
      max_workers=min(MAX_WORKERS, len(records))
  ) as pool:
    futures = {pool.submit(_run, r): i for i, r in enumerate(records)}
    for future in concurrent.futures.as_completed(futures):
      results[futures[future]] = future.result()

  eb.write_json(results, args.output)

  if args.summary_csv:
    lines = ['id,length,percent_helix,percent_strand,percent_coil,'
             'sasa_mean,n_buried,n_exposed']
    for r in results:
      s = r['summary']
      lines.append(
          f'{r["id"]},{r["length"]},{s["percent_helix"]},'
          f'{s["percent_strand"]},{s["percent_coil"]},{s["sasa_mean"]},'
          f'{s["n_buried"]},{s["n_exposed"]}'
      )
    with open(args.summary_csv, 'w', encoding='utf-8') as handle:
      handle.write('\n'.join(lines) + '\n')
    print(f'Summary table written to: {args.summary_csv}')


def cmd_plot(args) -> None:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  from matplotlib.patches import Patch

  data = _load(args.input)
  ss8 = data['secondary_structure_ss8']
  sasa = np.array(data['sasa'], dtype=np.float64)
  length = data['length']

  colors = {'H': '#d1495b', 'E': '#edae49', 'C': '#66a182'}
  ss3 = data['secondary_structure_ss3']

  fig, (ax_ss, ax_sasa) = plt.subplots(
      2, 1, figsize=(max(8, length / 12), 3.4), height_ratios=[1, 3],
      sharex=True,
  )
  for i, state in enumerate(ss3):
    ax_ss.axvspan(i - 0.5, i + 0.5, color=colors.get(state, '#cccccc'))
  ax_ss.set_yticks([])
  ax_ss.set_ylabel('SS', rotation=0, ha='right', va='center')
  ax_ss.set_title(f'{data.get("id", "protein")} — ESM3 predicted SS8 & SASA')
  ax_ss.legend(
      handles=[Patch(color=c, label=l) for l, c in
               (('helix', colors['H']), ('strand', colors['E']),
                ('coil', colors['C']))],
      ncol=3, loc='upper right', fontsize=7, frameon=False,
  )

  ax_sasa.fill_between(range(length), sasa, color='#2e4057', alpha=0.7)
  ax_sasa.axhline(BURIED_MAX, color='#888', ls=':', lw=0.8)
  ax_sasa.axhline(EXPOSED_MIN, color='#888', ls=':', lw=0.8)
  ax_sasa.set_ylabel('SASA (Å²)')
  ax_sasa.set_xlabel('residue')
  ax_sasa.set_xlim(-0.5, length - 0.5)

  fig.tight_layout()
  fig.savefig(args.output, dpi=130)
  print(f'Plot written to: {args.output}')


def _load(path: str) -> dict:
  import json
  with open(path, encoding='utf-8') as handle:
    data = json.load(handle)
  if isinstance(data, list):
    if not data:
      raise BiohubError(f'{path} holds an empty list.')
    return data[0]
  return data


def main() -> None:
  parser = argparse.ArgumentParser(
      description='Predict secondary structure (SS8) and solvent accessibility '
                  '(SASA) from sequence with ESM3.'
  )
  sub = parser.add_subparsers(dest='command', required=True)

  def add_common(p):
    p.add_argument('--model', default=eb.DEFAULT_ESM3,
                   help=f'ESM3 model (default {eb.DEFAULT_ESM3}).')
    p.add_argument('--temperature', type=float, default=0.7,
                   help='Sampling temperature. Lower = more conservative.')

  p_pred = sub.add_parser('predict', help='Predict for one sequence or FASTA.')
  p_pred.add_argument('--sequence', help='A single amino-acid sequence.')
  p_pred.add_argument('--fasta', help='A FASTA file (records handled in order).')
  p_pred.add_argument('--output', required=True, help='Output JSON path.')
  add_common(p_pred)
  p_pred.set_defaults(func=cmd_predict)

  p_batch = sub.add_parser('batch', help='Predict for every record in a FASTA.')
  p_batch.add_argument('--fasta', required=True, help='Input FASTA path.')
  p_batch.add_argument('--output', required=True, help='Output JSON path.')
  p_batch.add_argument('--summary-csv', help='Optional per-record summary CSV.')
  add_common(p_batch)
  p_batch.set_defaults(func=cmd_batch)

  p_plot = sub.add_parser('plot', help='Plot an SS8 ribbon over a SASA profile.')
  p_plot.add_argument('--input', required=True,
                      help='A prediction JSON (single record, or the first of '
                           'a list).')
  p_plot.add_argument('--output', required=True, help='Output PNG path.')
  p_plot.set_defaults(func=cmd_plot)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
