# Interpretation Guide — Guided Generation (SVDD)

How to read what `guided_design.py` gives you, how to choose an objective without
blowing the daily credit budget, and how to stay honest about a stochastic
result. Read this before interpreting a run. The worked example lives in
[examples/objective_guided/report.md](examples/objective_guided/report.md).

--------------------------------------------------------------------------------

## 1. Read the objective trajectory, not just the final number

The JSON's `objective_trajectory` (and the blue line in the PNG) is the
best-ranked candidate's objective value at each decoding step. What the shape
tells you:

| Trajectory shape | Reading |
|---|---|
| Climbs, then plateaus | Healthy guided run. Selection is finding and keeping gains. |
| Flat at the optimum (e.g. `[0,0,0,0]` for no-cysteine) | The target was **easy at these settings** — a handful of candidates already sufficed. Not a failure; judge against the *unguided baseline*, not the flat line. |
| Declines across steps | The value estimate is fighting you. Raise `--num-samples-per-step` so selection has more candidates to choose from. |
| Non-monotonic wobble (down then flat) | Normal on a noisy objective and a small budget. Each per-step value is a single Monte-Carlo estimate on a *partial* sequence; only the final step scores the true, complete sequence. Trust the gap to the baseline, not the wiggle. |

The faint candidate dots at each step are every candidate that was scored. A wide
spread with the best point well above the pack is selection doing real work
(see step 1 of the no-cysteine example: candidates at 0, 1, 2 cysteines, best
kept at 0). A tight cluster means the candidates barely differ and guidance has
little to grip.

**The per-step value is an estimate; the final value is real.** The last step
commits every position, so `final_objective_value` is the true objective on the
returned sequence. When they disagree, the final number wins.

## 2. Reading the metrics table

`final_metrics` always carries the cheap sequence-level numbers; pTM/pLDDT/Rg
appear only when the run folded the design. Thresholds:

| Quantity | Reading |
|---|---|
| pTM | > 0.8 confident fold · 0.5–0.8 plausible · < 0.5 unreliable |
| pLDDT (0–1 scale) | > 0.9 very high · 0.7–0.9 confident · 0.5–0.7 low · < 0.5 disordered |
| GRAVY | > 0 hydrophobic/membrane-like · ~ −0.4 typical soluble protein · < −1 very polar |
| pI | < 5 acidic · ~7 neutral · > 9 basic |
| Rg | a compact globular domain runs ~ 2.2·L^0.38 Å (≈ 11 Å at L=48); much larger = extended/unfolded |

**Anchor every judgement to unguided ESM3, measured on this API:** short
unconditional 48-mers carry a mean of **1.67 cysteines** (83% have ≥ 1), GRAVY
**−0.275 ± 0.317**, and fold to **pTM 0.09–0.26** with **Rg ≈ 36 Å**. In other
words, a short unconditional design is essentially an extended chain. A guided
48-mer that folds to pTM 0.31 has not "failed to fold well" relative to that
baseline — a length-48 ESM3 design was never going to be a crisp domain. Judge a
guided design against *these* numbers, never against a natural protein.

## 3. Choosing an objective: cheap vs expensive

The single most important cost decision is whether your objective reads a
**structure**. That sets `F` in the cost formula and can multiply your bill.

| Objective | Cost | Reads a fold? |
|---|---|---|
| `no-cysteine`, `hydrophobicity`, `isoelectric-point` | **cheap** | no (F = 0) |
| `ptm`, `radius-of-gyration` | **EXPENSIVE** | yes (F = 1: one ESMFold2 fold per candidate) |

```
guided run:    calls = steps x samples x (2 + F)  -  samples   [+ 1 final fold]
unguided run:  calls = steps                                   [+ 1 final fold]
```

A structural objective (or a structural `--constraint`, e.g. `ptm>=0.7`) forces
`F = 1`, folding **every** candidate. The same 2 x 2 run costs 11 credits with
`ptm` but only ~7 with a cheap objective. The script refuses any run above **60
calls** unless you pass `--yes`; do not reflexively add it. Quote the estimated
call count to the user and get agreement before a long run.

Practical guidance:

-   Prefer a **cheap, sequence-level** objective whenever it captures what you
    want. Composition, charge, and hydropathy targets need no fold.
-   Reach for `ptm`/`radius-of-gyration` only when foldability or globularity is
    genuinely the goal, and budget for `steps x samples` folds.
-   Always pair `radius-of-gyration` with `--constraint ptm>=0.7`. The most
    compact sequence is a dense, unfoldable blob, not a protein.
-   `--length` does **not** affect cost — only latency. Scale length freely;
    scale `steps` and `samples` only after doing the arithmetic.

## 4. The cost warning, restated

Biohub bills **one credit per API call**, and free-tier keys get **100 credits
per day**, shared with every other ESM skill. Exceeding it returns
`HTTP 429: You have exceeded your daily credit limit` and stops *everything*
until the quota resets at 00:00 UTC. A run that dies mid-way wastes every credit
it already spent. This is the most credit-hungry skill in the ESM suite by one
to two orders of magnitude — the tutorial's own settings (L=256, 32 x 10, pTM)
would cost **951 credits**, ten days of quota. Compute the cost first; the script
prints it.

## 5. Negative results & scientific integrity

A single guided run proves nothing on its own. Server-side ESM3 sampling is
stochastic and `--seed` seeds only the *local* choice of committed positions, so
two runs with the same seed give different sequences — never promise the user
reproducibility. Before you claim guidance worked:

-   **Run `compare`.** It runs the guided arm against the unguided ablation — the
    identical loop with candidate selection switched off — so any difference is
    attributable to selection alone, not to two different algorithms.
-   **Average a stochastic objective over ≥ 2 seeds.** In the worked example the
    two seeds *disagree*: guided won by +0.89 GRAVY at seed 0 and *lost* by −0.64
    at seed 1. Only the 2-seed mean (+0.132 guided vs +0.005 unguided) supports
    "guidance helps." Reporting seed 1 alone would have "proven" the opposite.
-   **Distinguish beating the ablation from beating plain generation.** On a
    small budget guidance reliably beats the unguided ablation but may only *tie*
    a single plain generation (guided −0.190 vs plain −0.183 at seed 0). Say
    which one you beat.
-   **Do not compute GRAVY, pI, cysteine counts, Rg or pTM yourself.** The JSON
    report carries every number; recomputing invites subtle scale errors (pLDDT
    is 0–1 here, not 0–100; the pI uses the Bjellqvist pKa set).
-   **Check `constraint_satisfied`.** A run can finish infeasible — the loop keeps
    the least-violating candidate so it can still make progress. A `false` there
    means the final design did *not* meet your `--constraint`.

## 6. Why this skill reimplements the SDK's guided decoder

The `esm` SDK ships `esm.sdk.experimental.ESM3GuidedDecoding`. This skill
deliberately does **not** use it, and reimplements SVDD against the raw Biohub
API, for three verified reasons:

-   It calls `get_esm3_model_tokenizers()`, which **downloads gated tokenizer
    assets and therefore requires a HuggingFace login.** No ESM skill may
    download model assets.
-   The vendor tutorial points it at `esm3-medium-2024-08`, which this API key
    **cannot reach (HTTP 403)**. Only `esm3-open-2024-03` is available.
-   Its `predict_denoised()` scores candidates at **temperature 0**, a *biased*
    estimator of the soft value function: greedy ESM3 completions of a
    mostly-masked prompt collapse toward poly-leucine and look nothing like the
    temperature-1 samples that actually get committed. Measured on this API, that
    bias drove GRAVY *down* when asked to maximise it. This skill scores at the
    generation temperature, making the denoised score an unbiased single-sample
    Monte-Carlo estimate — which is why **you must never set `--temperature 0`**.

See the shared [ESM/Biohub API reference](../references/esm-biohub-api.md) for the
reachable models, the credit model, and the record/replay cassettes that let you
reproduce the worked example for free.
