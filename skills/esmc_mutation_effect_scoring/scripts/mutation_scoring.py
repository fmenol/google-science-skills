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

"""Zero-shot mutation-effect scoring with ESM C (leave-one-out masked scan)."""

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
import math
import pathlib
import re
import sys
from typing import Any, Sequence

import matplotlib

matplotlib.use('Agg')  # Headless: never try to open a window.

import matplotlib.pyplot as plt  # pylint: disable=g-import-not-at-top
import numpy as np  # pylint: disable=g-import-not-at-top

from esm_biohub import BiohubClient, BiohubError  # pylint: disable=g-import-not-at-top
import esm_biohub as eb  # pylint: disable=g-import-not-at-top

# The rate limiter is cross-process, so concurrency is already bounded; keeping
# the pool small avoids piling up retries behind a 429.
MAX_WORKERS = 8

# Standard biology variant notation: wild-type residue, 1-indexed position,
# mutant residue. e.g. 'K48R' = lysine 48 to arginine.
VARIANT_RE = re.compile(r'^([A-Za-z])(\d+)([A-Za-z])$')

# Isoelectric points, used only to order the heatmap rows so that chemically
# similar residues sit together (as in the ESM mutation-scoring tutorial).
AA_ISOELECTRIC_POINT = {
    'A': 6.11, 'R': 10.76, 'N': 5.43, 'D': 2.98, 'C': 5.15, 'E': 3.08,
    'Q': 5.65, 'G': 6.06, 'H': 7.64, 'I': 6.04, 'L': 6.04, 'K': 9.47,
    'M': 5.71, 'F': 5.76, 'P': 6.30, 'S': 5.70, 'T': 5.60, 'W': 5.88,
    'Y': 5.63, 'V': 6.02,
}


# --------------------------------------------------------------------------
# Core numerics
# --------------------------------------------------------------------------


def leave_one_out_logprobs(
    client: BiohubClient,
    sequence: str,
    model: str,
    positions: Sequence[int] | None = None,
    max_workers: int = MAX_WORKERS,
) -> np.ndarray:
  """Masks each requested position in turn and returns its log-probabilities.

  For position `i` the model is shown `sequence[:i] + '_' + sequence[i+1:]`, so
  its prediction at `i` is made without ever seeing the wild-type residue there.
  That is the leave-one-out protocol, and it costs one API request per position.

  The sequence is encoded once and the mask token is then swapped in locally.
  This was verified against the live API to be byte-identical to encoding the
  '_'-masked string server-side, and it halves the request count (L+1 instead
  of 2L).

  Args:
    client: A BiohubClient.
    sequence: Wild-type sequence of length L.
    model: An ESMC model name.
    positions: 0-indexed positions to mask. None means every position.
    max_workers: Thread-pool size, capped at MAX_WORKERS.

  Returns:
    (L, 64) float64 log-softmax rows. Row `i` is the distribution predicted at
    residue `i`, read from token index `i + 1` because BOS occupies index 0.
    Rows for positions that were not requested are NaN.

  Raises:
    BiohubError: If the API returns tensors of an unexpected shape.
  """
  length = len(sequence)
  wanted = list(range(length)) if positions is None else sorted(set(positions))

  base_tokens = list(client.encode(sequence, model))
  if len(base_tokens) != length + 2:
    raise BiohubError(
        f'Expected {length + 2} tokens (BOS + {length} residues + EOS) but the '
        f'API returned {len(base_tokens)}.'
    )

  out = np.full((length, eb.LOGIT_DIM), np.nan, dtype=np.float64)

  def score_one(i: int) -> tuple[int, np.ndarray]:
    tokens = list(base_tokens)
    tokens[i + 1] = eb.SEQUENCE_MASK_TOKEN  # +1: BOS occupies token index 0.
    response = client.logits(tokens, model, sequence=True)
    logits = eb.to_array(response['logits']['sequence'])  # (L+2, 64)
    if logits.shape != (length + 2, eb.LOGIT_DIM):
      raise BiohubError(
          f'Expected logits of shape ({length + 2}, {eb.LOGIT_DIM}), got '
          f'{logits.shape}.'
      )
    return i, eb.log_softmax(logits[i + 1])  # +1: BOS offset. Never use [i].

  workers = max(1, min(int(max_workers), MAX_WORKERS, len(wanted)))
  with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    for i, row in pool.map(score_one, wanted):
      out[i] = row
  return out


def wt_log_probs(logp: np.ndarray, sequence: str) -> np.ndarray:
  """(L,) log-probability the model assigns to the actual wild-type residue."""
  index = np.array([eb.VOCAB[aa] for aa in sequence])
  return logp[np.arange(len(sequence)), index]


def llr_matrix(logp: np.ndarray, sequence: str) -> np.ndarray:
  """(L, 20) log-likelihood ratios over the canonical amino acids.

  LLR(i, aa) = log p(aa | masked context) - log p(wt_i | masked context), both
  read from the same leave-one-out distribution at position i. The wild-type
  column is therefore exactly 0 by construction.

  Columns follow `eb.AA20` (alphabetical: ACDEFGHIKLMNPQRSTVWY).
  """
  return logp[:, eb.AA20_IDX] - wt_log_probs(logp, sequence)[:, None]


def entropy_bits(logp: np.ndarray) -> np.ndarray:
  """(L,) Shannon entropy in bits of each leave-one-out distribution.

  Computed over the full 64-wide output distribution (the vocabulary is
  zero-padded to 64), matching the ESM tutorial. In practice the padding tokens
  carry ~1e-15 of the probability mass, so this is indistinguishable from the
  entropy over the real tokens.
  """
  return eb.shannon_entropy_bits(np.exp(logp), axis=-1)


def pseudo_perplexity(logp: np.ndarray, sequence: str) -> float:
  """exp(-mean_i log p(seq[i] | leave-one-out context)). Lower = more protein-like."""
  return float(math.exp(-float(np.mean(wt_log_probs(logp, sequence)))))


def wt_recovery(logp: np.ndarray, sequence: str) -> float:
  """Fraction of positions whose most-likely residue is the wild type."""
  predicted = np.argmax(logp, axis=1)
  actual = np.array([eb.VOCAB[aa] for aa in sequence])
  return float(np.mean(predicted == actual))


def fraction_deleterious(llr: np.ndarray, sequence: str) -> np.ndarray:
  """(L,) fraction of the 19 non-wild-type substitutions with LLR < 0.

  The wild-type column is excluded (it is exactly 0 and would otherwise cap the
  value at 19/20), so 1.0 means *every* substitution is predicted deleterious.
  """
  negative = (llr < 0).sum(axis=1).astype(np.float64)
  denominator = np.full(len(sequence), 19.0)
  # A non-canonical wild type (X, B, U, Z, O) has no zero column, so all 20
  # canonical substitutions are genuine substitutions there.
  for i, aa in enumerate(sequence):
    if aa not in eb.AA20:
      denominator[i] = 20.0
  return negative / denominator


# --------------------------------------------------------------------------
# Input handling
# --------------------------------------------------------------------------


def resolve_sequences(args) -> list[tuple[str, str]]:
  """Returns [(id, sequence)] from --sequence or --fasta."""
  if getattr(args, 'sequence', None) and getattr(args, 'fasta', None):
    raise BiohubError('Pass either --sequence or --fasta, not both.')
  if getattr(args, 'fasta', None):
    return [(name, eb.validate_sequence(seq)) for name, seq in
            eb.read_fasta(args.fasta)]
  if getattr(args, 'sequence', None):
    return [('query', eb.validate_sequence(args.sequence))]
  raise BiohubError('One of --sequence or --fasta is required.')


def resolve_one_sequence(args) -> tuple[str, str]:
  """Returns a single (id, sequence); errors if a FASTA holds several."""
  records = resolve_sequences(args)
  if len(records) > 1:
    raise BiohubError(
        f'{args.fasta} contains {len(records)} records but this subcommand '
        'scores one sequence. Split the FASTA.'
    )
  return records[0]


def check_model(model: str) -> str:
  """Validates the model name against the ESMC models this API key can reach."""
  if model not in eb.ESMC_MODELS:
    raise BiohubError(
        f'Unknown ESMC model {model!r}. Available: '
        f'{", ".join(sorted(eb.ESMC_MODELS))}'
    )
  return model


def parse_variants(specs: Sequence[str], sequence: str) -> list[dict[str, Any]]:
  """Parses and validates 1-indexed variant strings such as 'K48R'.

  Every problem is collected and reported together, because a mismatched
  wild-type residue almost always means the caller used 0-indexing or is
  holding the wrong isoform — and they will want to see all of them at once.

  Args:
    specs: Variant strings, e.g. ['A123G', 'K48R'].
    sequence: The wild-type sequence the positions refer to.

  Returns:
    [{'variant', 'wt', 'position' (1-indexed), 'mutant', 'index' (0-indexed)}]

  Raises:
    BiohubError: If any variant is unparseable, out of range, names a
      non-canonical mutant residue, or states a wild-type residue that does not
      match the sequence.
  """
  parsed: list[dict[str, Any]] = []
  problems: list[str] = []

  for raw in specs:
    spec = raw.strip()
    if not spec:
      continue
    match = VARIANT_RE.match(spec)
    if not match:
      problems.append(
          f'{spec!r}: cannot parse. Expected <WT><POS><MUT>, e.g. "K48R" '
          '(1-indexed).'
      )
      continue

    wt, position, mutant = (
        match.group(1).upper(), int(match.group(2)), match.group(3).upper()
    )

    if not 1 <= position <= len(sequence):
      problems.append(
          f'{spec!r}: position {position} is outside the sequence '
          f'(1..{len(sequence)}). Variants are 1-INDEXED.'
      )
      continue
    if mutant not in eb.AA20:
      problems.append(
          f'{spec!r}: mutant residue {mutant!r} is not one of the 20 canonical '
          'amino acids.'
      )
      continue

    actual = sequence[position - 1]
    if actual != wt:
      zero_indexed_hint = ''
      if position < len(sequence) and sequence[position] == wt:
        zero_indexed_hint = (
            f' NOTE: the sequence does have {wt} at 1-indexed position '
            f'{position + 1}, so this looks like a 0-indexing error.'
        )
      problems.append(
          f'{spec!r}: states wild-type {wt} at position {position}, but the '
          f'sequence has {actual} there (1-indexed).{zero_indexed_hint}'
      )
      continue

    parsed.append({
        'variant': f'{wt}{position}{mutant}',
        'wt': wt,
        'position': position,
        'mutant': mutant,
        'index': position - 1,
    })

  if problems:
    raise BiohubError(
        'Variant validation failed — refusing to score anything.\n  '
        + '\n  '.join(problems)
        + '\n\nVariants are 1-INDEXED and the wild-type residue must match the '
        'sequence you passed. A mismatch usually means 0-indexing was used, or '
        'that the position numbers came from a different isoform / a construct '
        'with a tag or cleaved signal peptide.'
    )
  if not parsed:
    raise BiohubError('No variants supplied.')
  return parsed


def read_variant_file(path: str) -> list[str]:
  """Reads variant specs, one per line. '#' comments and blanks are skipped."""
  specs: list[str] = []
  with open(path, encoding='utf-8') as handle:
    for line in handle:
      line = line.split('#', 1)[0].strip()
      if line:
        specs.extend(part for part in line.replace(',', ' ').split() if part)
  if not specs:
    raise BiohubError(f'No variants found in {path}')
  return specs


# --------------------------------------------------------------------------
# Scan persistence
# --------------------------------------------------------------------------


def write_scan(
    out_dir: pathlib.Path,
    sequence: str,
    model: str,
    logp: np.ndarray,
    top: int,
) -> dict[str, Any]:
  """Writes llr.npy, llr.csv, positions.csv and summary.json. Returns summary."""
  out_dir.mkdir(parents=True, exist_ok=True)
  length = len(sequence)

  llr = llr_matrix(logp, sequence)
  entropy = entropy_bits(logp)
  deleterious = fraction_deleterious(llr, sequence)
  mean_llr = llr.mean(axis=1)

  np.save(out_dir / 'llr.npy', llr.astype(np.float32))

  with open(out_dir / 'llr.csv', 'w', encoding='utf-8', newline='') as handle:
    writer = csv.writer(handle)
    writer.writerow(['position', 'wt', 'mutant', 'variant', 'llr'])
    for i in range(length):
      for j, aa in enumerate(eb.AA20):
        writer.writerow([
            i + 1, sequence[i], aa, f'{sequence[i]}{i + 1}{aa}',
            f'{llr[i, j]:.6f}',
        ])

  with open(
      out_dir / 'positions.csv', 'w', encoding='utf-8', newline=''
  ) as handle:
    writer = csv.writer(handle)
    writer.writerow([
        'position', 'wt', 'entropy_bits', 'fraction_deleterious', 'mean_llr',
        'best_substitution', 'best_llr', 'worst_substitution', 'worst_llr',
    ])
    for i in range(length):
      order = np.argsort(llr[i])
      # Best and worst *substitutions*: the wild-type column is 0 by
      # construction and is not a substitution, so it is skipped in both.
      worst_aa = next(
          eb.AA20[int(j)] for j in order if eb.AA20[int(j)] != sequence[i]
      )
      best_aa = next(
          eb.AA20[int(j)] for j in order[::-1] if eb.AA20[int(j)] != sequence[i]
      )
      writer.writerow([
          i + 1, sequence[i], f'{entropy[i]:.6f}', f'{deleterious[i]:.4f}',
          f'{mean_llr[i]:.6f}',
          f'{sequence[i]}{i + 1}{best_aa}',
          f'{llr[i, eb.AA20.index(best_aa)]:.6f}',
          f'{sequence[i]}{i + 1}{worst_aa}',
          f'{llr[i, eb.AA20.index(worst_aa)]:.6f}',
      ])

  by_entropy = np.argsort(entropy)

  def position_record(i: int) -> dict[str, Any]:
    return {
        'position': i + 1,
        'wt': sequence[i],
        'entropy_bits': round(float(entropy[i]), 4),
        'fraction_deleterious': round(float(deleterious[i]), 4),
        'mean_llr': round(float(mean_llr[i]), 4),
    }

  # Every single-residue substitution, ranked. WT columns are dropped (LLR 0).
  substitutions = [
      {
          'variant': f'{sequence[i]}{i + 1}{eb.AA20[j]}',
          'position': i + 1,
          'wt': sequence[i],
          'mutant': eb.AA20[j],
          'llr': round(float(llr[i, j]), 4),
      }
      for i in range(length)
      for j in range(len(eb.AA20))
      if eb.AA20[j] != sequence[i]
  ]
  substitutions.sort(key=lambda record: record['llr'])

  summary: dict[str, Any] = {
      'model': model,
      'sequence': sequence,
      'length': length,
      'method': 'leave-one-out masked marginal (one API request per position)',
      'amino_acid_columns': eb.AA20,
      'files': {
          'llr_matrix_npy': str((out_dir / 'llr.npy').resolve()),
          'llr_long_csv': str((out_dir / 'llr.csv').resolve()),
          'positions_csv': str((out_dir / 'positions.csv').resolve()),
      },
      'pseudo_perplexity': round(pseudo_perplexity(logp, sequence), 4),
      'mean_wt_log_prob': round(float(np.mean(wt_log_probs(logp, sequence))), 4),
      'wt_recovery_fraction': round(wt_recovery(logp, sequence), 4),
      'entropy_bits': [round(float(x), 4) for x in entropy],
      'fraction_deleterious': [round(float(x), 4) for x in deleterious],
      'mean_llr_per_position': [round(float(x), 4) for x in mean_llr],
      'entropy_summary': {
          'mean': round(float(entropy.mean()), 4),
          'std': round(float(entropy.std()), 4),
          'min': round(float(entropy.min()), 4),
          'max': round(float(entropy.max()), 4),
      },
      'most_constrained_positions': [
          position_record(int(i)) for i in by_entropy[:top]
      ],
      'most_tolerant_positions': [
          position_record(int(i)) for i in by_entropy[::-1][:top]
      ],
      'most_deleterious_substitutions': substitutions[:top],
      'best_tolerated_substitutions': substitutions[::-1][:top],
  }
  eb.write_json(summary, str(out_dir / 'summary.json'))
  return summary


def load_scan(scan_dir: str) -> tuple[str, np.ndarray, np.ndarray, dict]:
  """Loads a previous scan. Returns (sequence, llr (L,20), entropy (L,), summary)."""
  path = pathlib.Path(scan_dir).expanduser()
  summary_path, matrix_path = path / 'summary.json', path / 'llr.npy'
  for required in (summary_path, matrix_path):
    if not required.is_file():
      raise BiohubError(
          f'{required} not found. Run `scan --output-dir {scan_dir}` first.'
      )
  with open(summary_path, encoding='utf-8') as handle:
    summary = json.load(handle)
  llr = np.load(matrix_path)
  entropy = np.asarray(summary['entropy_bits'], dtype=np.float64)
  sequence = summary['sequence']
  if llr.shape != (len(sequence), len(eb.AA20)):
    raise BiohubError(
        f'{matrix_path} has shape {llr.shape}, expected '
        f'({len(sequence)}, {len(eb.AA20)}).'
    )
  return sequence, llr, entropy, summary


def scan_or_load(args) -> tuple[str, np.ndarray, np.ndarray, dict]:
  """Reuses --scan-dir if given, otherwise runs a fresh scan from --sequence."""
  if getattr(args, 'scan_dir', None):
    return load_scan(args.scan_dir)
  _, sequence = resolve_one_sequence(args)
  model = check_model(args.model)
  client = BiohubClient()
  announce(sequence, model)
  logp = leave_one_out_logprobs(
      client, sequence, model, max_workers=args.max_workers
  )
  llr = llr_matrix(logp, sequence)
  return sequence, llr, entropy_bits(logp), {'model': model}


def announce(sequence: str, model: str) -> None:
  """Tells the user what the scan is about to cost, before it costs it."""
  print(
      f'Leave-one-out scan: {len(sequence)} positions -> '
      f'{len(sequence) + 1} API requests to {model} ...',
      file=sys.stderr,
  )


# --------------------------------------------------------------------------
# Subcommands
# --------------------------------------------------------------------------


def cmd_scan(args) -> None:
  """Full deep-mutational scan: every residue mutated to every amino acid."""
  _, sequence = resolve_one_sequence(args)
  model = check_model(args.model)
  client = BiohubClient()
  announce(sequence, model)

  logp = leave_one_out_logprobs(
      client, sequence, model, max_workers=args.max_workers
  )
  summary = write_scan(
      pathlib.Path(args.output_dir).expanduser(), sequence, model, logp,
      args.top,
  )

  print(
      f'Scanned {summary["length"]} positions x 20 amino acids. '
      f'pseudo-perplexity={summary["pseudo_perplexity"]:.3f}, '
      f'WT recovery={summary["wt_recovery_fraction"]:.1%}, '
      f'entropy {summary["entropy_summary"]["min"]:.2f}-'
      f'{summary["entropy_summary"]["max"]:.2f} bits.'
  )
  print('Read summary.json for the ranked positions and substitutions.')


def cmd_score_variants(args) -> None:
  """Scores specific 1-indexed variants such as K48R."""
  _, sequence = resolve_one_sequence(args)
  model = check_model(args.model)

  specs: list[str] = []
  if args.variants:
    specs.extend(
        part for part in args.variants.replace(',', ' ').split() if part
    )
  if args.variants_file:
    specs.extend(read_variant_file(args.variants_file))
  if not specs:
    raise BiohubError('Pass --variants and/or --variants-file.')

  # Validate BEFORE spending a single API request.
  variants = parse_variants(specs, sequence)

  positions = sorted({v['index'] for v in variants})
  client = BiohubClient()
  print(
      f'Scoring {len(variants)} variant(s) at {len(positions)} position(s) -> '
      f'{len(positions) + 1} API requests to {model} ...',
      file=sys.stderr,
  )
  logp = leave_one_out_logprobs(
      client, sequence, model, positions=positions,
      max_workers=args.max_workers,
  )

  results = []
  for variant in variants:
    i = variant['index']
    row = logp[i]
    wt_lp = float(row[eb.VOCAB[variant['wt']]])
    mut_lp = float(row[eb.VOCAB[variant['mutant']]])
    llr_row = row[eb.AA20_IDX] - wt_lp  # (20,), wild-type column is 0.

    # Where does this substitution sit among the 19 alternatives available at
    # this position? Rank 1 = the least deleterious substitution here. This is
    # the only calibration available without a full scan: raw LLR magnitudes
    # are not comparable across positions, let alone across proteins.
    alternatives = sorted(
        (
            (float(llr_row[j]), aa)
            for j, aa in enumerate(eb.AA20)
            if aa != variant['wt']
        ),
        reverse=True,
    )
    rank = 1 + [aa for _, aa in alternatives].index(variant['mutant'])

    results.append({
        'variant': variant['variant'],
        'position': variant['position'],
        'wt': variant['wt'],
        'mutant': variant['mutant'],
        'llr': round(mut_lp - wt_lp, 4),
        'wt_log_prob': round(wt_lp, 4),
        'mutant_log_prob': round(mut_lp, 4),
        'position_entropy_bits': round(
            float(eb.shannon_entropy_bits(np.exp(row))), 4
        ),
        'position_fraction_deleterious': round(
            float((llr_row < 0).sum() / 19.0), 4
        ),
        'rank_at_position': rank,
        'rank_out_of': len(alternatives),
        'best_substitution_at_position': (
            f'{variant["wt"]}{variant["position"]}{alternatives[0][1]}'
        ),
        'best_llr_at_position': round(alternatives[0][0], 4),
        'predicted': 'deleterious' if mut_lp < wt_lp else 'tolerated',
    })

  results.sort(key=lambda record: record['llr'])
  payload = {
      'model': model,
      'sequence_length': len(sequence),
      'indexing': '1-indexed (standard biology convention)',
      'method': 'leave-one-out masked marginal',
      'n_variants': len(results),
      'variants_ranked_most_to_least_deleterious': results,
  }
  eb.write_json(payload, args.output)

  print(
      f'\n{"variant":<10s} {"LLR":>9s}  {"entropy":>8s}  {"rank@pos":>9s}  '
      'prediction'
  )
  for record in results:
    print(
        f'{record["variant"]:<10s} {record["llr"]:>9.3f}  '
        f'{record["position_entropy_bits"]:>8.2f}  '
        f'{record["rank_at_position"]:>4d}/{record["rank_out_of"]:<4d}  '
        f'{record["predicted"]}'
    )
  print(
      '\nLLR is unitless and only comparable WITHIN this protein. '
      'Rank the variants; do not read the magnitudes as kcal/mol.'
  )


def cmd_pseudo_perplexity(args) -> None:
  """Sequence fitness: how protein-like is this sequence to ESM C?"""
  records = resolve_sequences(args)
  model = check_model(args.model)
  client = BiohubClient()

  results = []
  for name, sequence in records:
    announce(sequence, model)
    logp = leave_one_out_logprobs(
        client, sequence, model, max_workers=args.max_workers
    )
    results.append({
        'id': name,
        'length': len(sequence),
        'pseudo_perplexity': round(pseudo_perplexity(logp, sequence), 4),
        'mean_wt_log_prob': round(
            float(np.mean(wt_log_probs(logp, sequence))), 4
        ),
        'wt_recovery_fraction': round(wt_recovery(logp, sequence), 4),
        'mean_entropy_bits': round(float(np.mean(entropy_bits(logp))), 4),
    })

  results.sort(key=lambda record: record['pseudo_perplexity'])
  eb.write_json(
      {
          'model': model,
          'method': (
              'exp(-mean log p(residue | leave-one-out masked context)); '
              'lower is more protein-like'
          ),
          'results_ranked_best_to_worst': results,
      },
      args.output,
  )

  print(f'\n{"id":<24s} {"len":>5s} {"pseudo-perplexity":>18s} {"WT recovery":>12s}')
  for record in results:
    print(
        f'{record["id"]:<24s} {record["length"]:>5d} '
        f'{record["pseudo_perplexity"]:>18.4f} '
        f'{record["wt_recovery_fraction"]:>11.1%}'
    )


def cmd_heatmap(args) -> None:
  """Renders the (L, 20) LLR matrix as a diverging heatmap centred at 0."""
  sequence, llr, _, _ = scan_or_load(args)
  length = len(sequence)

  # Rows ordered by descending isoelectric point so that basic residues sit at
  # the top and acidic at the bottom (as in the ESM tutorial).
  row_order = sorted(eb.AA20, key=lambda aa: -AA_ISOELECTRIC_POINT[aa])
  image = np.stack([llr[:, eb.AA20.index(aa)] for aa in row_order])  # (20, L)

  vmax = (
      float(args.vmax) if args.vmax is not None
      else float(np.nanmax(np.abs(image)))
  )
  if not np.isfinite(vmax) or vmax <= 0:
    vmax = 1e-6

  width = max(8.0, min(40.0, 0.16 * length + 4.0))
  fig, axes = plt.subplots(figsize=(width, 5.0))
  # 'bwr_r': negative (deleterious) -> red, 0 -> white, positive -> blue.
  mesh = axes.imshow(
      image, aspect='auto', cmap='bwr_r', vmin=-vmax, vmax=vmax,
      interpolation='nearest',
  )

  axes.set_xlabel('Sequence position (1-indexed)')
  axes.set_ylabel('Substituted amino acid (descending pI)')
  axes.set_yticks(np.arange(len(row_order)))
  axes.set_yticklabels(row_order, fontsize=8)

  ticks = [0] + [i - 1 for i in range(10, length, 10)] + [length - 1]
  axes.set_xticks(ticks)
  axes.set_xticklabels([str(t + 1) for t in ticks], fontsize=8)

  # Mark the wild-type residue at each position.
  for i in range(length):
    if sequence[i] in row_order:
      axes.text(
          i, row_order.index(sequence[i]), '.', ha='center', va='center',
          color='black', fontsize=7, fontweight='bold',
      )

  fig.colorbar(mesh, ax=axes, label='Log-likelihood ratio (LLR)', pad=0.01)
  axes.set_title(
      f'ESM C zero-shot mutation effects — L={length} '
      '(red = deleterious, blue = tolerated, dot = wild type)'
  )
  fig.tight_layout()

  out_path = pathlib.Path(args.output).expanduser()
  out_path.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(out_path, dpi=args.dpi)
  plt.close(fig)
  print(f'Success! Heatmap written to: {out_path}')
  print(
      f'Colour scale is symmetric about 0 at +/-{vmax:.2f} LLR. '
      'Red columns are positions where nearly every substitution is bad.'
  )


def cmd_entropy(args) -> None:
  """Per-position constraint profile: low entropy = evolutionarily constrained."""
  sequence, llr, entropy, _ = scan_or_load(args)
  length = len(sequence)
  deleterious = fraction_deleterious(llr, sequence)

  def uniform_entropy(n: int) -> float:
    """Entropy of a uniform distribution over n amino acids, in bits."""
    return math.log2(n)

  fig, axes = plt.subplots(
      figsize=(max(8.0, min(40.0, 0.16 * length + 4.0)), 4.5)
  )
  axes.bar(np.arange(1, length + 1), entropy, color='steelblue', width=0.85)
  for n, alpha in ((2, 0.25), (4, 0.4), (8, 0.55), (16, 0.7)):
    axes.axhline(
        uniform_entropy(n), color='grey', linestyle='--', linewidth=1.2,
        alpha=alpha, label=f'{n} AA uniform ({uniform_entropy(n):.0f} bits)',
    )
  axes.set_xlabel('Position (1-indexed)')
  axes.set_ylabel('Entropy (bits)')
  axes.set_title(
      'ESM C per-position constraint — low entropy = evolutionarily constrained'
  )
  axes.set_xlim(0.5, length + 0.5)
  axes.legend(fontsize=7, loc='upper right')

  # Annotate the most constrained positions with their wild-type residue.
  for i in np.argsort(entropy)[: min(args.annotate, length)]:
    axes.text(
        i + 1, entropy[i] + 0.02, sequence[i], ha='center', va='bottom',
        fontsize=6,
    )
  fig.tight_layout()

  out_path = pathlib.Path(args.output).expanduser()
  out_path.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(out_path, dpi=args.dpi)
  plt.close(fig)

  csv_path = out_path.with_suffix('.csv')
  with open(csv_path, 'w', encoding='utf-8', newline='') as handle:
    writer = csv.writer(handle)
    writer.writerow(
        ['position', 'wt', 'entropy_bits', 'fraction_deleterious']
    )
    for i in range(length):
      writer.writerow(
          [i + 1, sequence[i], f'{entropy[i]:.6f}', f'{deleterious[i]:.4f}']
      )

  order = np.argsort(entropy)
  print(f'Success! Entropy plot written to: {out_path}')
  print(f'Success! Data written to: {csv_path}')
  print(
      f'Entropy {entropy.min():.2f}-{entropy.max():.2f} bits '
      f'(mean {entropy.mean():.2f}, std {entropy.std():.2f}).'
  )
  print(
      '  Most constrained: '
      + ', '.join(
          f'{sequence[i]}{i + 1} ({entropy[i]:.2f})' for i in order[:5]
      )
  )
  print(
      '  Most tolerant:    '
      + ', '.join(
          f'{sequence[i]}{i + 1} ({entropy[i]:.2f})' for i in order[::-1][:5]
      )
  )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def add_sequence_args(parser: argparse.ArgumentParser) -> None:
  parser.add_argument('--sequence', help='Protein sequence (one-letter codes)')
  parser.add_argument('--fasta', help='FASTA file holding the sequence')


def add_common_args(parser: argparse.ArgumentParser) -> None:
  parser.add_argument(
      '--model', default=eb.DEFAULT_ESMC,
      help=f'ESMC model (default: {eb.DEFAULT_ESMC})',
  )
  parser.add_argument(
      '--max-workers', type=int, default=MAX_WORKERS,
      help=f'Parallel API requests, capped at {MAX_WORKERS}',
  )


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          'Zero-shot mutation-effect scoring with ESM C. No MSA, no structure, '
          'no labels. Variant positions are 1-INDEXED.'
      )
  )
  sub = parser.add_subparsers(dest='command', required=True)

  p_scan = sub.add_parser(
      'scan', help='Deep mutational scan: (L, 20) LLR matrix + entropy profile'
  )
  add_sequence_args(p_scan)
  add_common_args(p_scan)
  p_scan.add_argument(
      '--output-dir', required=True,
      help='Directory for llr.npy, llr.csv, positions.csv, summary.json',
  )
  p_scan.add_argument(
      '--top', type=int, required=True,
      help='How many ranked positions/substitutions to list in summary.json',
  )
  p_scan.set_defaults(func=cmd_scan)

  p_variants = sub.add_parser(
      'score-variants', help='Score specific 1-indexed variants, e.g. K48R'
  )
  add_sequence_args(p_variants)
  add_common_args(p_variants)
  p_variants.add_argument(
      '--variants', help='Comma-separated 1-indexed variants, e.g. A123G,K48R'
  )
  p_variants.add_argument(
      '--variants-file', help='File of variants, one per line'
  )
  p_variants.add_argument('--output', required=True, help='Output JSON path')
  p_variants.set_defaults(func=cmd_score_variants)

  p_pppl = sub.add_parser(
      'pseudo-perplexity',
      help='Sequence fitness: how protein-like is this sequence?',
  )
  add_sequence_args(p_pppl)
  add_common_args(p_pppl)
  p_pppl.add_argument('--output', required=True, help='Output JSON path')
  p_pppl.set_defaults(func=cmd_pseudo_perplexity)

  p_heatmap = sub.add_parser(
      'heatmap', help='(L, 20) LLR heatmap PNG, diverging colormap centred at 0'
  )
  add_sequence_args(p_heatmap)
  add_common_args(p_heatmap)
  p_heatmap.add_argument(
      '--scan-dir', help='Reuse a previous `scan` --output-dir (no API calls)'
  )
  p_heatmap.add_argument('--output', required=True, help='Output PNG path')
  p_heatmap.add_argument(
      '--vmax', type=float,
      help='Symmetric colour limit. Default: max |LLR| in the matrix.',
  )
  p_heatmap.add_argument('--dpi', type=int, default=150, help='PNG resolution')
  p_heatmap.set_defaults(func=cmd_heatmap)

  p_entropy = sub.add_parser(
      'entropy', help='Per-position constraint profile + plot'
  )
  add_sequence_args(p_entropy)
  add_common_args(p_entropy)
  p_entropy.add_argument(
      '--scan-dir', help='Reuse a previous `scan` --output-dir (no API calls)'
  )
  p_entropy.add_argument('--output', required=True, help='Output PNG path')
  p_entropy.add_argument(
      '--annotate', type=int, default=25,
      help='Label this many of the most constrained positions',
  )
  p_entropy.add_argument('--dpi', type=int, default=150, help='PNG resolution')
  p_entropy.set_defaults(func=cmd_entropy)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
