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

"""Objective-guided ESM3 protein generation (Soft Value-Based Decoding).

Reimplements the SVDD guided-decoding loop of Li et al. 2024 (arXiv:2408.08252)
directly against the Biohub REST API. The `esm` SDK's `ESM3GuidedDecoding` is
deliberately NOT used: it requires a HuggingFace login to download gated
tokenizer assets and targets `esm3-medium-2024-08`, which this API key cannot
reach (HTTP 403).
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
import dataclasses
import math
import pathlib
import random
import re
import sys
import time
from typing import Any, Callable

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # pylint: disable=g-import-not-at-top

from esm_biohub import BiohubClient, BiohubError  # vendored, same directory
import esm_biohub as eb  # pylint: disable=g-import-not-at-top

# --------------------------------------------------------------------------
# Cheap sequence-level biophysics (numpy only; no biopython)
# --------------------------------------------------------------------------

# Kyte & Doolittle (1982) hydropathy index. GRAVY = mean over residues.
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

# Bjellqvist pKa set -- the scale used by ExPASy Compute pI/Mw and Biopython's
# `IsoelectricPoint`. Verified in evals/eval_esm3_guided_generation.py to match
# Biopython to < 0.01 pH on real proteins. Do NOT swap in the EMBOSS values:
# they shift pI by up to 1.3 pH units and silently break `--constraint pi<=5`.
PKA_POSITIVE = {'Nterm': 7.5, 'K': 10.0, 'R': 12.0, 'H': 5.98}
PKA_NEGATIVE = {'Cterm': 3.55, 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0}
# The terminal pKa depends on which residue occupies the terminus.
PKA_NTERMINAL = {
    'A': 7.59, 'M': 7.00, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44, 'E': 7.70,
}
PKA_CTERMINAL = {'D': 4.55, 'E': 4.75}


def gravy(sequence: str) -> float:
  """Kyte-Doolittle GRAVY: mean hydropathy. Positive = hydrophobic."""
  if not sequence:
    return 0.0
  values = [KYTE_DOOLITTLE[a] for a in sequence if a in KYTE_DOOLITTLE]
  return float(np.mean(values)) if values else 0.0


def _terminal_pkas(sequence: str) -> tuple[dict[str, float], dict[str, float]]:
  """Bjellqvist tables with the N-/C-terminal pKa set by the terminal residue."""
  positive = dict(PKA_POSITIVE)
  negative = dict(PKA_NEGATIVE)
  positive['Nterm'] = PKA_NTERMINAL.get(sequence[0], PKA_POSITIVE['Nterm'])
  negative['Cterm'] = PKA_CTERMINAL.get(sequence[-1], PKA_NEGATIVE['Cterm'])
  return positive, negative


def net_charge(sequence: str, ph: float) -> float:
  """Net charge of a sequence at a given pH (Henderson-Hasselbalch)."""
  if not sequence:
    return 0.0
  positive, negative = _terminal_pkas(sequence)
  counts = {a: sequence.count(a) for a in set(sequence)}
  charge = 0.0
  # Positively ionisable groups: the protonated fraction carries +1.
  for group, pka in positive.items():
    n = 1 if group == 'Nterm' else counts.get(group, 0)
    if n:
      charge += n / (10.0 ** (ph - pka) + 1.0)
  # Negatively ionisable groups: the deprotonated fraction carries -1.
  for group, pka in negative.items():
    n = 1 if group == 'Cterm' else counts.get(group, 0)
    if n:
      charge -= n / (10.0 ** (pka - ph) + 1.0)
  return float(charge)


def isoelectric_point(sequence: str) -> float:
  """pI: the pH at which net charge is zero.

  Bisection over the full pH 0-14 range. (Biopython brackets its search to
  [4.05, 12] and therefore reports the bracket edge for extremely acidic or
  basic sequences; we do not.)
  """
  if not sequence:
    return 7.0
  low, high = 0.0, 14.0
  for _ in range(100):
    mid = 0.5 * (low + high)
    charge = net_charge(sequence, mid)
    if abs(charge) < 1e-8:
      return float(mid)
    if charge > 0:
      low = mid  # still net-positive: the pI is at a higher pH
    else:
      high = mid
  return float(0.5 * (low + high))


def radius_of_gyration(ca_coords: np.ndarray) -> float:
  """Rg of a CA trace in Angstrom. Lower = more compact/globular."""
  ca = np.asarray(ca_coords, dtype=np.float64)
  ca = ca[np.all(np.isfinite(ca), axis=1)]
  if len(ca) < 2:
    return float('nan')
  centred = ca - ca.mean(axis=0)
  return float(np.sqrt(np.mean(np.sum(centred**2, axis=1))))


# --------------------------------------------------------------------------
# Metrics for one complete sequence (folds lazily, at most once)
# --------------------------------------------------------------------------


@dataclasses.dataclass
class Metrics:
  """Everything an objective or constraint can read off a finished sequence."""

  sequence: str
  n_cysteine: int
  gravy: float
  pi: float
  ptm: float | None = None
  plddt: float | None = None
  rg: float | None = None
  coordinates: np.ndarray | None = None
  plddt_per_residue: np.ndarray | None = None

  def as_dict(self) -> dict[str, Any]:
    out = {
        'sequence': self.sequence,
        'length': len(self.sequence),
        'n_cysteine': self.n_cysteine,
        'gravy': round(self.gravy, 4),
        'isoelectric_point': round(self.pi, 3),
    }
    if self.ptm is not None:
      out['ptm'] = round(float(self.ptm), 4)
    if self.plddt is not None:
      out['mean_plddt'] = round(float(self.plddt), 4)
    if self.rg is not None:
      out['radius_of_gyration'] = round(float(self.rg), 3)
    return out


class Evaluator:
  """Computes `Metrics`. Folds only when an objective/constraint needs it."""

  def __init__(
      self,
      client: BiohubClient,
      *,
      needs_fold: bool,
      fold_model: str = eb.DEFAULT_ESMFOLD2,
      num_loops: int = 4,
      num_sampling_steps: int = 20,
  ):
    self.client = client
    self.needs_fold = needs_fold
    self.fold_model = fold_model
    self.num_loops = num_loops
    self.num_sampling_steps = num_sampling_steps
    self.n_folds = 0

  def evaluate(self, sequence: str, force_fold: bool = False) -> Metrics:
    metrics = Metrics(
        sequence=sequence,
        n_cysteine=sequence.count('C'),
        gravy=gravy(sequence),
        pi=isoelectric_point(sequence),
    )
    if not (self.needs_fold or force_fold):
      return metrics

    out = self.client.fold(
        sequence,
        self.fold_model,
        num_loops=self.num_loops,
        num_sampling_steps=self.num_sampling_steps,
    )
    self.n_folds += 1
    coords = eb.to_array(out['coordinates'])
    plddt = eb.to_array(out['plddt'])
    metrics.ptm = float(out['ptm'])
    metrics.plddt = float(np.nanmean(plddt))
    metrics.plddt_per_residue = plddt
    metrics.coordinates = coords
    metrics.rg = radius_of_gyration(coords[:, eb.ATOM37.index('CA'), :])
    return metrics


# --------------------------------------------------------------------------
# Objective registry. Every objective is MAXIMISED.
# --------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Objective:
  name: str
  description: str
  needs_fold: bool
  score: Callable[[Metrics], float]
  # Human-readable value reported alongside the (signed) score.
  report: Callable[[Metrics], float]
  units: str


OBJECTIVES: dict[str, Objective] = {
    'none': Objective(
        name='none',
        description=(
            'Unguided baseline. Runs the identical decoding loop with selection '
            'switched off (one candidate per step, kept unconditionally). Use '
            'to measure what guidance actually bought you.'
        ),
        needs_fold=False,
        score=lambda m: 0.0,
        report=lambda m: 0.0,
        units='',
    ),
    'no-cysteine': Objective(
        name='no-cysteine',
        description=(
            'Minimise the cysteine count (score = -count). CHEAP and '
            'deterministically verifiable: a converged run has zero C.'
        ),
        needs_fold=False,
        score=lambda m: -float(m.n_cysteine),
        report=lambda m: float(m.n_cysteine),
        units='cysteines',
    ),
    'hydrophobicity': Objective(
        name='hydrophobicity',
        description=(
            'Kyte-Doolittle GRAVY (mean hydropathy). CHEAP. maximize = '
            'hydrophobic/membrane-like; minimize = polar/soluble.'
        ),
        needs_fold=False,
        score=lambda m: m.gravy,
        report=lambda m: m.gravy,
        units='GRAVY',
    ),
    'isoelectric-point': Objective(
        name='isoelectric-point',
        description=(
            'Isoelectric point (EMBOSS pKa set). CHEAP. maximize = basic '
            'protein; minimize = acidic protein.'
        ),
        needs_fold=False,
        score=lambda m: m.pi,
        report=lambda m: m.pi,
        units='pH',
    ),
    'ptm': Objective(
        name='ptm',
        description=(
            'Maximise predicted TM-score of the folded candidate (foldability / '
            'generation quality). EXPENSIVE: one ESMFold2 fold per candidate.'
        ),
        needs_fold=True,
        score=lambda m: float(m.ptm),
        report=lambda m: float(m.ptm),
        units='pTM',
    ),
    'radius-of-gyration': Objective(
        name='radius-of-gyration',
        description=(
            'Maximise globularity by minimising Rg of the folded CA trace '
            '(score = -Rg). EXPENSIVE: one ESMFold2 fold per candidate. '
            'Pair with --constraint ptm>=0.7 or you will get compact garbage.'
        ),
        needs_fold=True,
        score=lambda m: -float(m.rg),
        report=lambda m: float(m.rg),
        units='Angstrom',
    ),
}

# Metric names usable in --constraint, and whether reading them needs a fold.
CONSTRAINT_METRICS: dict[str, tuple[Callable[[Metrics], float], bool]] = {
    'ptm': (lambda m: float(m.ptm), True),
    'plddt': (lambda m: float(m.plddt), True),
    'rg': (lambda m: float(m.rg), True),
    'gravy': (lambda m: m.gravy, False),
    'pi': (lambda m: m.pi, False),
    'n_cysteine': (lambda m: float(m.n_cysteine), False),
}

_CONSTRAINT_RE = re.compile(
    r'^\s*([A-Za-z_]+)\s*(>=|<=|==|=|>|<)\s*(-?\d+(?:\.\d+)?)\s*$'
)


@dataclasses.dataclass(frozen=True)
class Constraint:
  """A hard constraint, enforced by feasibility-first candidate ranking."""

  metric: str
  op: str
  value: float
  needs_fold: bool

  def violation(self, metrics: Metrics) -> float:
    """0.0 when satisfied; otherwise the (positive) distance to feasibility."""
    getter, _ = CONSTRAINT_METRICS[self.metric]
    x = getter(metrics)
    if not math.isfinite(x):
      return float('inf')
    if self.op in ('>=', '>'):
      return max(0.0, self.value - x)
    if self.op in ('<=', '<'):
      return max(0.0, x - self.value)
    return abs(x - self.value)

  def __str__(self) -> str:
    return f'{self.metric}{self.op}{self.value:g}'


def parse_constraint(text: str) -> Constraint:
  """Parses `ptm>=0.7`, `rg<=15`, `n_cysteine<=0`, ..."""
  match = _CONSTRAINT_RE.match(text)
  if not match:
    raise BiohubError(
        f'Cannot parse --constraint {text!r}. Expected METRIC OP VALUE, e.g. '
        "'ptm>=0.7'. Metrics: " + ', '.join(sorted(CONSTRAINT_METRICS))
    )
  metric, op, value = match.group(1).lower(), match.group(2), match.group(3)
  if metric not in CONSTRAINT_METRICS:
    raise BiohubError(
        f'Unknown constraint metric {metric!r}. Available: '
        + ', '.join(sorted(CONSTRAINT_METRICS))
    )
  op = '==' if op == '=' else op
  return Constraint(
      metric=metric,
      op=op,
      value=float(value),
      needs_fold=CONSTRAINT_METRICS[metric][1],
  )


# --------------------------------------------------------------------------
# The SVDD guided-decoding loop
# --------------------------------------------------------------------------


@dataclasses.dataclass
class Candidate:
  prompt: str          # the candidate's partially-decoded state (may hold '_')
  denoised: str        # a COMPLETE completion of `prompt`; what gets scored
  metrics: Metrics
  score: float         # signed objective (always maximised)
  violation: float     # 0.0 if feasible

  def rank_key(self) -> tuple[int, float]:
    """Feasible candidates always outrank infeasible ones."""
    if self.violation <= 0.0:
      return (1, self.score)
    return (0, -self.violation)


def _n_masked(prompt: str) -> int:
  return prompt.count(eb.MASK_CHAR)


def guided_generate(
    client: BiohubClient,
    *,
    length: int,
    objective: Objective,
    direction: str,
    num_decoding_steps: int,
    num_samples_per_step: int,
    temperature: float,
    completion_steps: int,
    constraint: Constraint | None,
    evaluator: Evaluator,
    model: str,
    seed: int,
    max_workers: int = 8,
    verbose: bool = True,
) -> dict[str, Any]:
  """Runs Soft Value-Based Decoding (SVDD) on the ESM3 sequence track.

  Each step, for each of `num_samples_per_step` candidates:
    1. PROPOSE  — a full iterative ESM3 completion of the current prompt, at
       `temperature`. (1 API call.)
    2. COMMIT   — lock in a random slice of the newly-filled positions. This is
       the candidate's partially-decoded state; everything else returns to '_'.
    3. DENOISE  — complete the candidate again, so the value estimate reflects
       what the candidate ACTUALLY locked in. (1 API call; skipped on the final
       step, where nothing is left masked.)
    4. SCORE    — evaluate the objective on the denoised prediction (folding it
       first if the objective or the constraint reads a structural metric).
  The best candidate's prompt carries into the next step. The final step commits
  everything, so the returned sequence is scored on the TRUE objective.

  Two corrections to the SDK's `ESM3GuidedDecoding`, both load-bearing
  ---------------------------------------------------------------------
  * `predict_denoised` uses temperature 0. That is a BIASED estimator of the
    soft value V(x_t) = E[r(x_0) | x_t]: greedy ESM3 completions of a
    mostly-masked prompt collapse toward poly-leucine and look nothing like the
    temperature-1 samples that will actually be committed. Ranking by the greedy
    completion optimises the wrong quantity — measured on this API it drove
    GRAVY *down* when asked to maximise it. We sample the completion at the
    generation temperature, making r(denoised) an unbiased single-sample
    Monte-Carlo estimate of V.
  * The proposal is a full iterative completion, not a single forward pass, so
    the committed tokens come from a coherent protein rather than from the
    per-position marginal. Same one API call either way.

  Step 3 is what gives the loop its teeth: because the denoised prediction is
  conditioned on the candidate's committed slice, a candidate that locked in a
  cysteine keeps it and is penalised, while a candidate that did not can escape.
  Scoring the *proposal* instead (one call cheaper) cannot tell those apart, and
  in simulation that drops the zero-cysteine success rate from ~100% to ~55%.

  `objective='none'` is the unguided ablation: draw a single completion per step,
  commit a slice, never score. Drawing one and keeping it is statistically
  identical to drawing N and picking one uniformly at random, so it ablates
  *selection* exactly, at a fraction of the cost.
  """
  sign = -1.0 if direction == 'minimize' else 1.0
  unguided = objective.name == 'none'
  rng = random.Random(seed)

  prompt = eb.MASK_CHAR * length
  trajectory: list[dict[str, Any]] = []
  started = time.time()

  for step in range(num_decoding_steps):
    masked = [i for i, ch in enumerate(prompt) if ch == eb.MASK_CHAR]
    if not masked:
      break
    remaining_steps = num_decoding_steps - step
    n_unmask = math.ceil(len(masked) / remaining_steps)
    n_draw = 1 if unguided else num_samples_per_step

    # Draw every candidate's committed slice up front, on this thread:
    # random.Random is not thread-safe and we want `--seed` to be meaningful.
    slices = [
        rng.sample(masked, min(n_unmask, len(masked))) for _ in range(n_draw)
    ]

    # PROPOSE: draw `n_draw` completions of the SAME prompt. They are i.i.d.
    # samples, so the assignment of a completion to an (indexed) committed slice
    # is statistically arbitrary. We nevertheless fix it deterministically by
    # sorting the completions before pairing. This makes the selected winner —
    # and therefore the whole trajectory — reproducible under any thread
    # completion order, which is what lets a recorded run replay exactly.
    def _propose(_idx: int) -> str:
      return client.generate(
          'sequence', model, sequence=prompt,
          num_steps=min(completion_steps, length), temperature=temperature,
      )['outputs']['sequence']

    if n_draw == 1:
      completions = [_propose(0)]
    else:
      with concurrent.futures.ThreadPoolExecutor(
          max_workers=min(max_workers, n_draw)
      ) as pool:
        completions = sorted(pool.map(_propose, range(n_draw)))

    def _sample_and_score(idx: int) -> Candidate:
      """COMMIT -> DENOISE -> SCORE for one candidate's proposed completion."""
      completion = completions[idx]

      chars = list(prompt)
      for i in slices[idx]:
        chars[i] = completion[i]
      cand_prompt = ''.join(chars)

      if _n_masked(cand_prompt) == 0:
        # Final step: the candidate is already complete. No second call, and the
        # score is the TRUE objective rather than an estimate of it.
        denoised = cand_prompt
      elif unguided:
        denoised = cand_prompt  # never scored; skip the call entirely
      else:
        denoised = client.generate(
            'sequence', model, sequence=cand_prompt,
            num_steps=min(completion_steps, length), temperature=temperature,
        )['outputs']['sequence']

      if unguided:
        # Cheap metrics only: the unguided arm is never scored, and folding it
        # would burn a credit per candidate for nothing.
        cheap = Metrics(
            denoised, denoised.count('C'), gravy(denoised),
            isoelectric_point(denoised),
        )
        return Candidate(cand_prompt, denoised, cheap, 0.0, 0.0)
      metrics = evaluator.evaluate(denoised)
      score = sign * objective.score(metrics)
      violation = constraint.violation(metrics) if constraint else 0.0
      return Candidate(cand_prompt, denoised, metrics, score, violation)

    if n_draw == 1:
      candidates = [_sample_and_score(0)]
    else:
      with concurrent.futures.ThreadPoolExecutor(
          max_workers=min(max_workers, n_draw)
      ) as pool:
        candidates = list(pool.map(_sample_and_score, range(n_draw)))

    best = max(candidates, key=Candidate.rank_key)
    prompt = best.prompt

    scores = [c.score for c in candidates]
    trajectory.append({
        'step': step,
        'positions_committed': n_unmask,
        'positions_remaining': _n_masked(prompt),
        'best_score': round(best.score, 5),
        'best_value': round(objective.report(best.metrics), 5),
        'candidate_scores': [round(s, 5) for s in scores],
        'mean_candidate_score': round(float(np.mean(scores)), 5),
        'feasible': bool(best.violation <= 0.0),
        'violation': round(best.violation, 5),
        'denoised_prediction': best.denoised,
    })
    if verbose:
      tag = '' if constraint is None else (
          f'  {constraint}: '
          + ('OK' if best.violation <= 0 else f'VIOLATED by {best.violation:.3g}')
      )
      print(
          f'  step {step + 1:>2}/{num_decoding_steps}  '
          f'committed {n_unmask:>3} ({_n_masked(prompt):>3} left)  '
          f'{objective.name}={objective.report(best.metrics):+.4g} '
          f'{objective.units}{tag}',
          file=sys.stderr, flush=True,
      )

  # The final step commits every remaining position, so `prompt` is complete.
  if _n_masked(prompt) > 0:  # defensive: only reachable if the loop broke early
    prompt = client.generate(
        'sequence', model, sequence=prompt,
        num_steps=min(completion_steps, length), temperature=temperature,
    )['outputs']['sequence']

  final_sequence = eb.validate_sequence(prompt)
  return {
      'sequence': final_sequence,
      'trajectory': trajectory,
      'elapsed_s': round(time.time() - started, 1),
  }


def plain_generate(
    client: BiohubClient, *, length: int, num_steps: int, temperature: float,
    model: str,
) -> str:
  """The no-skill baseline: one plain iterative ESM3 generation."""
  out = client.generate(
      'sequence', model, sequence=eb.MASK_CHAR * length,
      num_steps=min(num_steps, length), temperature=temperature,
  )
  return eb.validate_sequence(out['outputs']['sequence'])


# --------------------------------------------------------------------------
# Cost accounting
# --------------------------------------------------------------------------


def needs_fold(objective: Objective, constraint: Constraint | None) -> bool:
  """True when scoring a candidate requires a structure."""
  return objective.needs_fold or (
      constraint is not None and constraint.needs_fold
  )


def estimate_calls(
    *, objective: Objective, constraint: Constraint | None,
    num_decoding_steps: int, num_samples_per_step: int,
    final_fold: bool = True,
) -> int:
  """API calls for one run. This is the cost formula quoted in SKILL.md.

    guided:    steps x samples x (2 + F)  -  samples  [+ 1 final fold]
    unguided:  steps                                  [+ 1 final fold]

  Per candidate: 1 PROPOSE call + 1 DENOISE call + F fold calls, where F = 1 if
  the objective or the constraint reads a structural metric (pTM / pLDDT / Rg)
  and 0 otherwise. The DENOISE call is skipped on the final step (the candidate
  is already complete), which is the `- samples` term. The unguided arm draws
  one candidate per step and never scores it, so it costs one call per step.
  """
  if objective.name == 'none':
    return num_decoding_steps + (1 if final_fold else 0)
  fold = 1 if needs_fold(objective, constraint) else 0
  total = num_decoding_steps * num_samples_per_step * (2 + fold)
  total -= num_samples_per_step  # final step needs no DENOISE call
  return total + (1 if final_fold else 0)


# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------


def plot_trajectory(
    trajectory: list[dict[str, Any]],
    objective: Objective,
    direction: str,
    constraint: Constraint | None,
    output: str,
    baselines: list[tuple[str, float]] | None = None,
) -> None:
  """The guided objective per decoding step, with the candidate spread.

  `baselines` draws horizontal reference lines (the unguided and plain-generate
  finals). We plot those as lines rather than as trajectories because the
  unguided arm is never scored per-step — inventing a per-step curve for it
  would be fiction.
  """
  fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)

  if trajectory:
    steps = [t['step'] + 1 for t in trajectory]
    values = [t['best_value'] for t in trajectory]
    ax.plot(steps, values, 'o-', color='#1f77b4', label='guided (best/step)',
            lw=2, ms=6, zorder=3)
    # Candidate spread: every candidate this step, mapped back onto the
    # human-readable value axis. `_signed` is its own inverse (it only ever
    # multiplies by +/-1), so it converts score -> value as well as value ->
    # score.
    for i, t in enumerate(trajectory):
      vals = [_signed(objective, direction, s) for s in t['candidate_scores']]
      ax.plot(
          [t['step'] + 1] * len(vals), vals, '.', color='#1f77b4', alpha=0.35,
          ms=8, zorder=2,
          label='candidates scored' if i == 0 else None,
      )

  for colour, (label, value) in zip(
      ('#d62728', '#7f7f7f'), baselines or []
  ):
    if value is not None and math.isfinite(value):
      ax.axhline(value, ls='--', lw=1.6, color=colour,
                 label=f'{label} = {value:+.3g}', zorder=1)

  if constraint is not None and constraint.metric == _objective_metric(objective):
    ax.axhline(
        constraint.value, ls=':', color='k', lw=1.4,
        label=f'constraint {constraint}', zorder=1,
    )

  ax.set_xlabel('decoding step')
  ax.set_ylabel(f'{objective.name} ({objective.units})')
  ax.set_title(
      f'SVDD guided generation — {objective.name} ({direction})', fontsize=11
  )
  ax.grid(alpha=0.3)
  ax.legend(fontsize=8)
  fig.savefig(output, dpi=150)
  plt.close(fig)
  print(f'Success! Plot written to: {output}')


def rescore_unguided_trajectory(
    trajectory: list[dict[str, Any]], objective: Objective, direction: str
) -> None:
  """Rewrites an unguided trajectory to carry the REAL objective's values.

  The unguided arm runs under the `none` objective, whose reported value is a
  constant 0.0. Leaving that in a `compare` report whose objective is, say,
  hydrophobicity would be actively misleading. The per-step completions are
  stored, so for a sequence-only objective we can recover the true values for
  free. A structural objective would need a fold per step, which we refuse to
  pay for, so we null the values out rather than invent them.
  """
  sign = -1.0 if direction == 'minimize' else 1.0
  for entry in trajectory:
    if objective.needs_fold:
      entry['best_value'] = None
      entry['best_score'] = None
      entry['candidate_scores'] = []
      entry['note'] = 'unguided arm not scored per-step (would need a fold)'
      continue
    seq = entry['denoised_prediction']
    metrics = Metrics(
        seq, seq.count('C'), gravy(seq), isoelectric_point(seq)
    )
    entry['best_value'] = round(objective.report(metrics), 5)
    entry['best_score'] = round(sign * objective.score(metrics), 5)
    entry['candidate_scores'] = [entry['best_score']]


def _objective_metric(objective: Objective) -> str:
  return {
      'ptm': 'ptm',
      'radius-of-gyration': 'rg',
      'hydrophobicity': 'gravy',
      'isoelectric-point': 'pi',
      'no-cysteine': 'n_cysteine',
  }.get(objective.name, '')


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------


def cmd_list_objectives(args) -> None:
  del args
  print('Guided-generation objectives (all are MAXIMISED; use --direction '
        'minimize to flip):\n')
  rows = []
  for obj in OBJECTIVES.values():
    cost = 'EXPENSIVE (folds every candidate)' if obj.needs_fold else 'cheap (sequence-only)'
    rows.append((obj.name, cost, obj.description))
  width = max(len(r[0]) for r in rows)
  for name, cost, desc in rows:
    print(f'  {name:<{width}}  [{cost}]')
    for line in _wrap(desc, 72):
      print(f'  {"":<{width}}  {line}')
    print()
  print('Constraint metrics for --constraint (e.g. --constraint "ptm>=0.7"):')
  for metric, (_, needs_fold) in sorted(CONSTRAINT_METRICS.items()):
    tag = 'EXPENSIVE (folds every candidate)' if needs_fold else 'cheap'
    print(f'  {metric:<12} [{tag}]')


def _wrap(text: str, width: int) -> list[str]:
  words, lines, cur = text.split(), [], ''
  for word in words:
    if len(cur) + len(word) + 1 > width:
      lines.append(cur)
      cur = word
    else:
      cur = f'{cur} {word}'.strip()
  if cur:
    lines.append(cur)
  return lines


def _setup(args) -> tuple[Objective, Constraint | None, Evaluator, BiohubClient]:
  """Shared validation + client/evaluator construction."""
  if args.objective not in OBJECTIVES:
    raise BiohubError(
        f'Unknown objective {args.objective!r}. Available: '
        + ', '.join(OBJECTIVES)
    )
  objective = OBJECTIVES[args.objective]
  constraint = parse_constraint(args.constraint) if args.constraint else None

  if not 1 <= args.num_decoding_steps <= args.length:
    raise BiohubError(
        f'--num-decoding-steps must be in 1..{args.length} (the sequence '
        f'length); got {args.num_decoding_steps}.'
    )
  if args.num_samples_per_step < 1:
    raise BiohubError('--num-samples-per-step must be >= 1.')

  client = BiohubClient()
  evaluator = Evaluator(
      client,
      needs_fold=needs_fold(objective, constraint),
      fold_model=args.fold_model,
      num_loops=args.fold_loops,
      num_sampling_steps=args.fold_sampling_steps,
  )
  return objective, constraint, evaluator, client


# Biohub bills one credit per API call. Free-tier keys are capped at 100 credits
# PER DAY -- a single tutorial-scale guided run (L=256, 32 steps x 10 samples)
# would need ~320 credits, i.e. three days of quota. Refuse to silently burn it.
DAILY_CREDIT_LIMIT = 100
COST_LIMIT = 60


def _cost_gate(calls: int, args, label: str = '') -> None:
  """Warns loudly, and refuses runs that would eat the daily quota."""
  suffix = f' [{label}]' if label else ''
  print(
      f'Estimated API calls{suffix}: ~{calls} '
      f'({args.num_decoding_steps} steps x {args.num_samples_per_step} '
      f'samples/step). Biohub bills 1 credit per call; free-tier keys get '
      f'{DAILY_CREDIT_LIMIT} credits PER DAY.',
      file=sys.stderr,
  )
  if calls > COST_LIMIT and not args.yes:
    raise BiohubError(
        f'This run needs ~{calls} API calls, above the {COST_LIMIT}-call safety '
        f'limit (free-tier keys get only {DAILY_CREDIT_LIMIT} credits per DAY, '
        'and exceeding it returns HTTP 429 until the quota resets). Guided '
        'decoding costs steps x samples x (2 + fold) - samples calls. Reduce '
        '--num-decoding-steps / --num-samples-per-step, choose a cheap '
        'objective (see `list-objectives`), or pass --yes to proceed anyway. '
        'Note that --length does NOT affect cost.'
    )


def cmd_generate(args) -> None:
  objective, constraint, evaluator, client = _setup(args)
  calls = estimate_calls(
      objective=objective, constraint=constraint,
      num_decoding_steps=args.num_decoding_steps,
      num_samples_per_step=args.num_samples_per_step,
  )
  _cost_gate(calls, args)

  print(
      f'Guided generation: length={args.length} objective={objective.name} '
      f'({args.direction})'
      + (f' constraint={constraint}' if constraint else ''),
      file=sys.stderr,
  )
  run = guided_generate(
      client,
      length=args.length,
      objective=objective,
      direction=args.direction,
      num_decoding_steps=args.num_decoding_steps,
      num_samples_per_step=args.num_samples_per_step,
      temperature=args.temperature,
      completion_steps=args.completion_steps,
      constraint=constraint,
      evaluator=evaluator,
      model=args.model,
      seed=args.seed,
  )

  # Always fold the final design so the report carries pTM / pLDDT / Rg.
  final = evaluator.evaluate(run['sequence'], force_fold=True)
  report = {
      'command': 'generate',
      'model': args.model,
      'fold_model': args.fold_model,
      'objective': objective.name,
      'objective_units': objective.units,
      'direction': args.direction,
      'constraint': str(constraint) if constraint else None,
      'constraint_satisfied': (
          bool(constraint.violation(final) <= 0.0) if constraint else None
      ),
      'length': args.length,
      'num_decoding_steps': args.num_decoding_steps,
      'num_samples_per_step': args.num_samples_per_step,
      'temperature': args.temperature,
      'completion_steps': args.completion_steps,
      'seed': args.seed,
      'sequence': run['sequence'],
      'final_metrics': final.as_dict(),
      'final_objective_value': round(objective.report(final), 5),
      'trajectory': run['trajectory'],
      'objective_trajectory': [t['best_value'] for t in run['trajectory']],
      'api_calls_estimated': calls,
      'tokens_used': client.tokens_used,
      'elapsed_s': run['elapsed_s'],
  }
  eb.write_json(report, args.output)

  out_path = pathlib.Path(args.output)
  plot_path = args.plot or str(out_path.with_suffix('.png'))
  plot_trajectory(
      run['trajectory'], objective, args.direction, constraint, plot_path,
  )

  if final.coordinates is not None:
    pdb_path = out_path.with_suffix('.pdb')
    pdb_path.write_text(
        eb.atom37_to_pdb(
            final.coordinates, run['sequence'], final.plddt_per_residue
        ),
        encoding='utf-8',
    )
    print(f'Success! Structure written to: {pdb_path}')

  print(
      f'\nFinal sequence ({args.length} aa): {run["sequence"]}\n'
      f'{objective.name}: {objective.report(final):+.4g} {objective.units}   '
      f'pTM: {final.ptm:.3f}   mean pLDDT: {final.plddt:.3f}   '
      f'Rg: {final.rg:.1f} A   Cys: {final.n_cysteine}'
  )


def cmd_compare(args) -> None:
  """Guided vs unguided with identical settings — an exact selection ablation."""
  objective, constraint, evaluator, client = _setup(args)
  if objective.name == 'none':
    raise BiohubError(
        "`compare` needs a real objective; 'none' IS the unguided arm."
    )
  # `compare` folds the three final sequences only when the objective actually
  # needs a structure. For a cheap objective that would be 3 wasted credits.
  fold_finals = needs_fold(objective, constraint)
  guided_calls = estimate_calls(
      objective=objective, constraint=constraint,
      num_decoding_steps=args.num_decoding_steps,
      num_samples_per_step=args.num_samples_per_step,
      final_fold=fold_finals,
  )
  none_obj = OBJECTIVES['none']
  unguided_calls = estimate_calls(
      objective=none_obj, constraint=None,
      num_decoding_steps=args.num_decoding_steps,
      num_samples_per_step=args.num_samples_per_step,
      final_fold=fold_finals,
  )
  plain_calls = 1 + (1 if fold_finals else 0)
  _cost_gate(
      guided_calls + unguided_calls + plain_calls, args, label='all three arms'
  )

  arms: dict[str, dict[str, Any]] = {}
  for label, obj in (('guided', objective), ('unguided', none_obj)):
    print(f'\n[{label}] objective={obj.name}', file=sys.stderr)
    run = guided_generate(
        client,
        length=args.length,
        objective=obj,
        direction=args.direction,
        num_decoding_steps=args.num_decoding_steps,
        num_samples_per_step=args.num_samples_per_step,
        temperature=args.temperature,
        completion_steps=args.completion_steps,
        constraint=constraint if label == 'guided' else None,
        evaluator=evaluator,
        model=args.model,
        seed=args.seed,
    )
    if label == 'unguided':
      # The unguided arm ran under `none`, whose reported value is constant 0.
      # Restate its trajectory in terms of the objective under test.
      rescore_unguided_trajectory(
          run['trajectory'], objective, args.direction
      )
    # Score BOTH arms with the SAME objective.
    metrics = evaluator.evaluate(run['sequence'], force_fold=fold_finals)
    arms[label] = {
        'sequence': run['sequence'],
        'objective_value': round(objective.report(metrics), 5),
        'metrics': metrics.as_dict(),
        'trajectory': run['trajectory'],
        'elapsed_s': run['elapsed_s'],
    }

  # A third reference point: what a plain, non-SVDD generation gives you.
  plain_seq = plain_generate(
      client, length=args.length, num_steps=args.num_decoding_steps,
      temperature=args.temperature, model=args.model,
  )
  plain_metrics = evaluator.evaluate(plain_seq, force_fold=fold_finals)

  g_val = arms['guided']['objective_value']
  u_val = arms['unguided']['objective_value']
  # "Better" is direction-aware and objective-aware: compare signed scores.
  better = _signed(objective, args.direction, g_val) > _signed(
      objective, args.direction, u_val
  )

  report = {
      'command': 'compare',
      'model': args.model,
      'fold_model': args.fold_model,
      'objective': objective.name,
      'objective_units': objective.units,
      'direction': args.direction,
      'constraint': str(constraint) if constraint else None,
      'length': args.length,
      'num_decoding_steps': args.num_decoding_steps,
      'num_samples_per_step': args.num_samples_per_step,
      'temperature': args.temperature,
      'seed': args.seed,
      'guided': arms['guided'],
      'unguided': arms['unguided'],
      'plain_generate_baseline': {
          'sequence': plain_seq,
          'objective_value': round(objective.report(plain_metrics), 5),
          'metrics': plain_metrics.as_dict(),
          'note': (
              'A single plain ESM3 iterative generation (no SVDD loop at all) '
              '— what the esm3-protein-design skill would give you.'
          ),
      },
      'guided_better_than_unguided': bool(better),
      'delta_guided_minus_unguided': round(g_val - u_val, 5),
      'api_calls_estimated': guided_calls + unguided_calls + plain_calls,
      'tokens_used': client.tokens_used,
  }
  eb.write_json(report, args.output)

  plot_path = args.plot or str(pathlib.Path(args.output).with_suffix('.png'))
  plot_trajectory(
      arms['guided']['trajectory'], objective, args.direction, constraint,
      plot_path,
      baselines=[
          ('unguided (selection off)', u_val),
          ('plain generate', objective.report(plain_metrics)),
      ],
  )

  print(
      f'\n{objective.name} ({args.direction}, {objective.units}):\n'
      f'  guided            {g_val:+.4g}\n'
      f'  unguided          {u_val:+.4g}\n'
      f'  plain generate    {objective.report(plain_metrics):+.4g}\n'
      f'  guided better than unguided: {better}'
  )


def _signed(objective: Objective, direction: str, value: float) -> float:
  """Maps a human-readable objective value onto the maximised score axis."""
  sign = -1.0 if direction == 'minimize' else 1.0
  if objective.name in ('no-cysteine', 'radius-of-gyration'):
    return -value * sign
  return value * sign


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _add_common(parser: argparse.ArgumentParser) -> None:
  parser.add_argument('--length', type=int, required=True,
                      help='Length of the protein to design, in residues.')
  parser.add_argument('--objective', required=True,
                      help='Objective name. See `list-objectives`.')
  parser.add_argument('--num-decoding-steps', type=int, required=True,
                      help='REQUIRED. Decoding steps (1..length). Cost is '
                           'linear in this.')
  parser.add_argument('--num-samples-per-step', type=int, required=True,
                      help='REQUIRED. Candidates proposed and scored per step. '
                           'Cost is linear in this.')
  parser.add_argument('--direction', choices=('maximize', 'minimize'),
                      default='maximize',
                      help='Optimise the objective up or down (default: '
                           'maximize).')
  parser.add_argument('--temperature', type=float, default=1.0,
                      help='Sampling temperature for the candidate completions '
                           '(default: 1.0). Do NOT set 0: greedy ESM3 '
                           'completions of a mostly-masked prompt collapse to '
                           'poly-leucine and the value estimate becomes '
                           'meaningless.')
  parser.add_argument('--completion-steps', type=int, default=8,
                      help='ESM3 iterative decoding steps used to build each '
                           'candidate completion (default: 8). Server-side '
                           'work only -- still one API call.')
  parser.add_argument('--constraint', default=None,
                      help="Hard constraint, e.g. 'ptm>=0.7'. Feasible "
                           'candidates always outrank infeasible ones.')
  parser.add_argument('--seed', type=int, default=0,
                      help='Seed for LOCAL position selection only. ESM3 '
                           'sampling happens server-side and is NOT seedable, '
                           'so runs are not bit-reproducible.')
  parser.add_argument('--model', default=eb.DEFAULT_ESM3,
                      help=f'ESM3 model (default: {eb.DEFAULT_ESM3}).')
  parser.add_argument('--fold-model', default=eb.DEFAULT_ESMFOLD2,
                      help=f'ESMFold2 model (default: {eb.DEFAULT_ESMFOLD2}).')
  parser.add_argument('--fold-loops', type=int, default=4,
                      help='ESMFold2 num_loops when scoring (default: 4).')
  parser.add_argument('--fold-sampling-steps', type=int, default=20,
                      help='ESMFold2 num_sampling_steps (default: 20).')
  parser.add_argument('--plot', default=None,
                      help='Trajectory plot path (default: <output>.png).')
  parser.add_argument('--output', required=True, help='Output JSON file path.')
  parser.add_argument('--yes', action='store_true',
                      help=f'Proceed past the {COST_LIMIT}-API-call safety '
                           'limit.')


def main() -> None:
  parser = argparse.ArgumentParser(
      description='Objective-guided ESM3 generation (SVDD, Li et al. 2024).'
  )
  sub = parser.add_subparsers(dest='command', required=True)

  p_list = sub.add_parser(
      'list-objectives',
      help='Print the objective registry with descriptions and API cost.',
  )
  p_list.set_defaults(func=cmd_list_objectives)

  p_gen = sub.add_parser(
      'generate', help='Design a protein that optimises an objective.'
  )
  _add_common(p_gen)
  p_gen.set_defaults(func=cmd_generate)

  p_cmp = sub.add_parser(
      'compare',
      help='Guided vs unguided with identical settings (selection ablation).',
  )
  _add_common(p_cmp)
  p_cmp.set_defaults(func=cmd_compare)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
