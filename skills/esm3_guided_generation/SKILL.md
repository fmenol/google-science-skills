---
name: esm3-guided-generation
description: >
  Designs proteins that OPTIMISE a user-specified objective using ESM3 with
  Soft Value-Based Decoding (SVDD, Li et al. 2024). Use when the user wants to
  maximise or minimise a measurable property of a designed protein — "design a
  protein with no cysteines", "make it more hydrophobic / more soluble", "design
  an acidic protein", "maximise pTM / foldability", "make it globular / compact",
  "steer generation toward X", "optimise this objective", "guided generation",
  "SVDD" — or to impose a hard constraint such as "keep pTM above 0.7". Do NOT
  use when no objective is being optimised: plain unconditional generation,
  motif scaffolding, or generation conditioned on structure / secondary
  structure / SASA belongs to the `esm3-protein-design` skill, which does the
  same job in ONE API call instead of hundreds.
---

# Guided Protein Generation with ESM3 (SVDD)

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions.
2.  **User Notification**: If `.licenses/esm3_guided_generation_LICENSE.txt` does
    not exist in the workspace root then (1) prominently notify the user to check
    the terms at https://biohub.org/acceptable-use-policy/ and https://biohub.ai/,
    then (2) create the file recording the notification text and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for it.

> **Compute options.** This skill uses the hosted Biohub API (`BIOHUB_API_KEY` in `~/.env`). To run equivalent capabilities on the open
> weights instead — on [Modal](https://modal.com), or on a local GPU if one is verified available — see
> [compute-options.md](../esm_common/references/compute-options.md).


## READ THIS FIRST — the cost warning

Guided decoding turns generation into optimisation by *sampling many candidates
and throwing most of them away*. It is the most API-expensive skill in the ESM
suite by one to two orders of magnitude.

```
guided run:    calls = steps x samples x (2 + F)  -  samples   [+ 1 final fold]
unguided run:  calls = steps                                   [+ 1 final fold]

F = 1 if the objective OR the constraint reads a structural metric
    (ptm / plddt / rg), because every candidate must then be folded.
F = 0 for the sequence-only objectives (no-cysteine, hydrophobicity,
    isoelectric-point).
```

**Biohub bills one credit per API call, and free-tier keys get 100 credits PER
DAY.** Exceeding that returns `HTTP 429: You have exceeded your daily credit
limit of 100 credits` and *everything* stops until the quota resets — including
every other ESM skill sharing the key.

| Run | Calls | Fits in a 100-credit day? |
|---|---|---|
| `no-cysteine`, 4 steps x 3 samples | 22 | yes |
| `hydrophobicity`, 3 steps x 4 samples | 21 | yes |
| `ptm`, 3 steps x 2 samples | 17 | yes |
| `ptm`, 8 steps x 5 samples | 116 | **no** |
| The tutorial's settings (L=256, 32 steps x 10 samples, pTM) | **951** | **no — 10 days of quota** |

The script refuses any run above **60 calls** unless you pass `--yes`. Do not
reflexively add `--yes`. **Tell the user the estimated call count and get their
agreement before launching a long run.**

Cost does **not** depend on `--length`. A 300-residue design costs the same
number of credits as a 48-residue one — only latency grows. Scale up length
before you scale up steps or samples.

## Overview

`scripts/guided_design.py` implements Soft Value-Based Decoding on ESM3's
sequence track. Starting from a fully-masked chain, each decoding step:

1.  **PROPOSE** — draw `--num-samples-per-step` full ESM3 completions of the
    current partially-decoded prompt (1 call each).
2.  **COMMIT** — lock in a random slice of one candidate's newly-filled
    positions; the rest return to `_`.
3.  **DENOISE** — complete each candidate again so the value estimate reflects
    what that candidate *actually locked in* (1 call each; free on the last
    step).
4.  **SCORE** — evaluate the objective on the denoised prediction, folding it
    first if the objective needs a structure. Keep the best candidate.

The final step commits everything, so the returned sequence is scored on the
**true** objective rather than an estimate of it.

**What this skill does NOT do.** It does not do plain or conditioned generation
(use `esm3-protein-design`), inverse folding (`esm3-inverse-folding`), or
structure prediction (`esmfold2-structure-prediction`). It optimises the
*sequence* track only.

### Why the `esm` SDK's `ESM3GuidedDecoding` is not used

The vendor SDK ships `esm.sdk.experimental.ESM3GuidedDecoding`. This skill
**deliberately reimplements it against the raw API** for three reasons:

*   It calls `get_esm3_model_tokenizers()`, which **downloads gated tokenizer
    assets and therefore requires a HuggingFace login.** No ESM skill may
    download model assets.
*   The tutorial points it at `esm3-medium-2024-08`, which this API key **cannot
    reach (HTTP 403)**. Only `esm3-open-2024-03` is available.
*   Its `predict_denoised()` scores candidates at **temperature 0**, which is a
    *biased* estimator of the soft value function and silently breaks
    composition objectives — see below.

## Core Rules

-   **ALWAYS `uv run --no-project`.** Without `--no-project`, `uv` walks up the
    directory tree, finds an unrelated `pyproject.toml`, and tries to build that
    project instead.
-   **ALWAYS quote the estimated call count to the user before a long run.** The
    daily quota is 100 credits and it is shared with every other ESM skill.
-   **NEVER set `--temperature 0`.** Greedy ESM3 completions of a mostly-masked
    prompt collapse to poly-leucine (`MLLLLLLL...`, measured). The value estimate
    becomes meaningless and guidance optimises the wrong thing.
-   **NEVER interpret `--seed` as making a run reproducible.** It seeds only the
    *local* choice of which positions to commit. ESM3 samples server-side and is
    not seedable, so two runs with the same seed give different sequences. Do not
    promise the user reproducibility.
-   **Do not compute GRAVY, pI, cysteine counts, Rg, pTM or "is it better than
    unguided" yourself; always use the script's output.** The JSON report carries
    every number you need.
-   **ALWAYS run `compare` before claiming guidance worked.** A single guided run
    proves nothing on its own — the objective could have landed there by chance.
-   **Before interpreting a run, read `docs/interpretation-guide.md` and review
    the worked example in `docs/examples/objective_guided/`** — including its
    honest guided-vs-unguided result, so you know what a real (and a merely lucky)
    improvement looks like.
-   **Pair `radius-of-gyration` with `--constraint ptm>=0.7`.** Minimising Rg
    without a foldability constraint produces compact garbage.
-   **Notification**: If this skill is used, ensure this is mentioned in the
    output.

## Utility Scripts

### `list-objectives` — the registry and its cost (0 API calls)

```bash
uv run --no-project scripts/guided_design.py list-objectives
```

| Objective | Cost | Meaning (all are MAXIMISED; `--direction minimize` flips) |
|---|---|---|
| `none` | cheap | Unguided baseline: the same loop with selection switched off. |
| `no-cysteine` | cheap | score = −(cysteine count). Converged run has **zero** C. |
| `hydrophobicity` | cheap | Kyte-Doolittle GRAVY. maximize = membrane-like; minimize = soluble. |
| `isoelectric-point` | cheap | pI (Bjellqvist). maximize = basic; minimize = acidic. |
| `ptm` | **EXPENSIVE** | pTM of the folded candidate. One ESMFold2 fold per candidate. |
| `radius-of-gyration` | **EXPENSIVE** | score = −Rg of the folded CA trace. Globularity. |

### `generate` — design a protein that optimises an objective

```bash
# Cheap, sequence-only objective. 22 calls.
uv run --no-project scripts/guided_design.py generate \
  --length 48 --objective no-cysteine \
  --num-decoding-steps 4 --num-samples-per-step 3 \
  --seed 0 --output design.json
```

```bash
# Expensive: folds every candidate. 17 calls.
uv run --no-project scripts/guided_design.py generate \
  --length 100 --objective ptm \
  --num-decoding-steps 3 --num-samples-per-step 2 \
  --output ptm_design.json
```

```bash
# Globular design with a hard foldability constraint. 17 calls.
uv run --no-project scripts/guided_design.py generate \
  --length 100 --objective radius-of-gyration --constraint 'ptm>=0.7' \
  --num-decoding-steps 3 --num-samples-per-step 2 \
  --output globular.json
```

`--num-decoding-steps` and `--num-samples-per-step` are **required** — there is
no safe default for something that multiplies your bill. Writes
`design.json` (sequence, per-step trajectory, final pTM/pLDDT/Rg),
`design.png` (trajectory plot) and `design.pdb` (the folded design).

### `compare` — guided vs unguided, same settings (the honest test)

```bash
uv run --no-project scripts/guided_design.py compare \
  --length 48 --objective hydrophobicity \
  --num-decoding-steps 3 --num-samples-per-step 4 \
  --seed 0 --output compare.json
```

Runs three arms and reports the objective for each:

*   **guided** — full SVDD.
*   **unguided** — *the identical decoding loop with selection switched off*
    (one candidate per step, kept unconditionally). This is an exact ablation of
    the selection step, not a different algorithm, so any difference is
    attributable to guidance alone.
*   **plain generate** — one ordinary ESM3 generation, i.e. what
    `esm3-protein-design` would give you for 1 credit.

### Hard constraints

`--constraint 'METRIC OP VALUE'` where METRIC is `ptm`, `plddt`, `rg`, `gravy`,
`pi` or `n_cysteine`, e.g. `'ptm>=0.7'`, `'rg<=15'`, `'n_cysteine<=0'`.

Enforced by **feasibility-first ranking**: feasible candidates always outrank
infeasible ones; among feasible candidates the best objective wins; if *none* is
feasible the least-violating one is kept so the run can still make progress. This
replaces the SDK's MDMM dual-variable ascent, which needs a learning rate and a
damping term to tune and can diverge. `constraint_satisfied` in the JSON tells
you whether the final design actually met it — **check it; a run can finish
infeasible.**

## Interpreting the Output

| Quantity | Reading |
|---|---|
| pTM | > 0.8 confident fold · 0.5–0.8 plausible · < 0.5 unreliable |
| pLDDT (0–1) | > 0.9 very high · 0.7–0.9 confident · 0.5–0.7 low · < 0.5 disordered |
| GRAVY | > 0 hydrophobic/membrane-like · ~ −0.4 typical soluble protein · < −1 very polar |
| pI | < 5 acidic · ~7 neutral · > 9 basic |
| Rg | a compact globular domain runs ~ 2.2·L^0.38 Å (≈11 Å at L=48). Much larger = extended/unfolded. |

**Reference points measured on this API.** Unguided ESM3 48-mers carry a mean of
**1.67 cysteines** (83% have at least one), GRAVY **−0.275 ± 0.317**, and fold to
pTM **0.09–0.26** with Rg ≈ 36 Å — i.e. short unconditional designs are
essentially extended chains. Judge a guided design against *those* numbers, not
against a natural protein.

**Read the trajectory, not just the final number.** `objective_trajectory` in the
JSON is the best candidate's objective at each step. A healthy guided run climbs
and then plateaus. A flat trajectory at the optimum (e.g. `[0,0,0,0]` for
`no-cysteine`) means the objective was easy at those settings. A trajectory that
*declines* means the value estimate is fighting you — raise `--num-samples-per-step`.

**Do not over-claim on a small budget.** With 2–3 steps and 2 samples the
selection sees only 4–6 candidates. That is enough to demonstrate the machinery
and to hit a hard, cheap target like zero-cysteine, but it is **not** enough to
reliably beat unguided generation on a noisy structural metric like pTM. Say so
rather than reporting a lucky draw as a result.

## Common Mistakes

1.  **Using this skill when nothing is being optimised.** Plain generation,
    motif scaffolding, or conditioning on structure/SS/SASA is
    `esm3-protein-design` — **one** API call. Reaching for guided decoding there
    burns hundreds of credits for an identical result. Guided decoding is only
    worth it when you have a scoring function you genuinely want to maximise.
2.  **Scaling steps x samples without doing the arithmetic.** `--num-decoding-steps
    16 --num-samples-per-step 10 --objective ptm` is 471 calls — nearly five days
    of quota — and it will die mid-run with HTTP 429, wasting everything spent up
    to that point. Compute the cost first; the script prints it.
3.  **Setting `--temperature 0` "to make it deterministic".** It does the
    opposite of what you want: greedy completions of a mostly-masked prompt
    collapse to poly-leucine, so every candidate scores the same and guidance
    degenerates to a random walk.
4.  **Believing a guided run without its control.** ESM3's unguided output is
    already variable; one guided sequence with a good score proves nothing. Run
    `compare`, and for a stochastic objective average over ≥ 2 seeds.
5.  **Minimising `radius-of-gyration` with no `--constraint ptm>=…`.** The most
    compact sequence is not a protein. You will get a dense, unfoldable blob.

## Dependencies

*   `uv` — environment setup.
*   `credentials` — safe handling of `BIOHUB_API_KEY`.
*   `esm3-protein-design` — **use it instead** for plain/conditioned generation.
*   `esmfold2-structure-prediction` — for folding a design in more depth than the
    single fold this skill reports.

## References

*   [Worked example: objective-guided generation](docs/examples/objective_guided/report.md)
    — `no-cysteine` hits exactly 0 cysteines; guided GRAVY beats the unguided
    ablation on the 2-seed mean (+0.132 vs +0.005) but only ties a single plain
    generation. Review it before running your own, especially the honest
    guided-vs-unguided story.
*   [Interpretation guide](docs/interpretation-guide.md) — reading the objective
    trajectory, choosing cheap sequence-level vs expensive fold-based objectives,
    the cost warning, and why the SDK's `ESM3GuidedDecoding` is reimplemented on
    the API.
*   [ESM/Biohub API reference](references/esm-biohub-api.md) — shared endpoints,
    reachable models, the daily credit model, and the record/replay cassettes.
*   `references/citation.bib` — cite **both** ESM3 (Hayes et al. 2025) and
    **SVDD (Li et al. 2024, arXiv:2408.08252)** if you use this skill.
