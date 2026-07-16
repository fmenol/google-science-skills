---
name: esmfold2-binder-screening
description: >
  Screens candidate protein binders against a target by folding every
  target:binder complex with ESMFold2 and ranking the candidates by interface
  confidence (iPTM, interface PAE, interface contacts). Use when the user has a
  target protein and a set of candidate binders, minibinders, nanobodies or
  scFvs and asks which of them actually binds, which to shortlist, which to
  order, or wants to triage / rank / filter / score binder designs, or to
  inspect a predicted protein-protein interface. Do not use when the user wants
  to GENERATE binder sequences from scratch (use `esm3-protein-design`, or
  `esm3-inverse-folding` to design a sequence onto a backbone) or to fold a
  single chain with no binding partner (use `esmfold2-structure-prediction`).
---

# Binder Screening with ESMFold2

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If `.licenses/esmfold2_binder_screening_LICENSE.txt`
    does not already exist in the workspace root directory then (1) prominently
    notify the user to check the terms at
    https://biohub.org/acceptable-use-policy/ and https://biohub.ai/, then
    (2) create the file recording the notification text and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for and request
    this key.

## Overview

Given one **target** protein and a FASTA of **candidate binders**, this skill
folds every `target:candidate` complex with ESMFold2 and ranks the candidates by
how confident the model is in the *interface*. It answers exactly one question:

> **Of the candidates I already have, which ones actually bind — and which
> should I order?**

This is the **selection / triage** stage of the binder-design protocol from
*Language Modeling Materializes a World Model of Protein Biology* — the stage
that decides which designs are worth spending money on.

### Scope — read this before you start

**This skill SCREENS and RANKS candidate binders. It does NOT generate them.**

The published protocol has two halves. Only one is reachable here:

| Half | What it does | Available? |
|---|---|---|
| **Design** (Alg. 11–14) | Backpropagates structural losses *through the folding trunk* to invent a binder sequence | **No.** Needs local H100 weights, 27–51 GB VRAM, Modal, and a `torch.autograd` path into the model. There is no remote-gradient endpoint. Do not attempt it. |
| **Selection** (the "critics") | Folds target+binder, ranks candidates by interface confidence | **Yes — this skill.** |

So: **produce** candidates with `esm3-protein-design` (generate new sequences) or
`esm3-inverse-folding` (design a sequence onto a backbone), or bring your own
(homologs, mutant panels, antibody variants, literature binders) — then **come
back here to rank them.**

Two further consequences of the API surface, so nobody wastes time:

*   `include_distogram=True` returns **HTTP 422 — not implemented server-side.**
    The distogram-based iPTM proxy from the paper (Alg. 15) is therefore
    **unreachable**. This skill uses the **real `interface_ptm`** the API
    returns, which is a better metric anyway. Do not try to reconstruct the
    proxy.
*   `num_sampling_steps` is **capped at 100** by the API (HTTP 422 above that).
    The paper's final critic uses 200; that exact setting is not reachable. It
    does not matter — the controls separate cleanly at 50–100.

## Core Rules

-   **ALWAYS `uv run --no-project`.** Without `--no-project`, `uv` walks up the
    directory tree, finds an unrelated `pyproject.toml`, and tries to build
    *that* project instead. Never `pip install`, never bare `python`/`python3`.
-   **NEVER download model weights.** No `torch`, `transformers`,
    `huggingface_hub`, `esm` or `modal`. Everything runs remotely on the Biohub
    API.
-   **ALWAYS rank on `iptm`** unless the user asks otherwise. It is the only
    metric that directly measures the *interface*, and it is by far the most
    stable one (measured: barnase:barstar returns iPTM 0.963–0.964 across every
    `num_loops`/`num_sampling_steps` setting tried).
-   **NEVER judge binding by contact count.** See *Common Mistakes* — a
    non-binder routinely shows **more** raw contacts than a true binder.
-   **Do not compute any of these metrics yourself; always use the script's
    output.** Do not eyeball the PDB, re-derive iPTM, or estimate an interface
    by hand.
-   **iPTM is a structural-confidence score, NOT an affinity.** It does not give
    you a Kd. See *Interpreting the Output*.
-   **Before interpreting results, read
    [`docs/interpretation-guide.md`](docs/interpretation-guide.md).** It defines
    the iPTM bands, the interface-PAE reading, the contact-count trap, and the
    "not an affinity" caveat.
-   **Review a worked example in [`docs/examples/`](docs/examples/) before
    writing a report** — including the negative
    [`decoy_controls`](docs/examples/decoy_controls/report.md), so you know what
    "no binding" looks like (a well-folded non-binder, a composition-matched
    decoy).
-   **Write the report using
    [`docs/report-templates.md`](docs/report-templates.md).**
-   **This skill RANKS candidates; it does NOT design them.** To generate binder
    sequences use `esm3-protein-design` or `esm3-inverse-folding`, then return
    here to rank. See *Scope*.
-   If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All commands run from the skill directory.

### 1. `screen` — fold every complex and rank the candidates

Folds each `target:candidate` complex (target = chain A, candidate = chain B),
derives the interface metrics, and writes a ranked table plus one PDB per
candidate.

```bash
uv run --no-project scripts/screen_binders.py screen \
  --target-fasta target.fasta \
  --candidates candidates.fasta \
  --output results/ \
  --top-n 10
```

`--target SEQ` works instead of `--target-fasta`. **`--top-n` is required** — you
must decide explicitly how many candidates you would actually order; there is no
silent default. Other flags: `--sort-by` (default `iptm`), `--num-loops`
(default 3), `--num-sampling-steps` (default 100, max 100), `--max-workers`
(default 4, capped at 4), `--model`.

Writes into `--output`: `results.json` (everything), `results.csv` (the ranked
table), and `structures/<candidate>.pdb` (two chains: A = target, B = binder;
B-factor column carries pLDDT × 100).

### 2. `rank` — re-rank or filter an existing screen

Cheap and offline — reuses `results.json`, folds nothing.

```bash
uv run --no-project scripts/screen_binders.py rank \
  --results results/results.json \
  --output shortlist.json \
  --min-iptm 0.8 \
  --top-n 5
```

Filters: `--min-iptm`, `--min-binder-plddt`, `--max-interface-pae`,
`--min-contacts` (confident contacts), `--max-isoelectric-point`. Re-sort with
`--sort-by`. The reference protocol filters minibinders to **pI < 6** before
ordering, as a solubility/expression proxy — pass `--max-isoelectric-point 6`
(do not apply it to antibodies).

### 3. `analyze-interface` — deep-dive one complex

```bash
uv run --no-project scripts/screen_binders.py analyze-interface \
  --target-fasta target.fasta \
  --binder-fasta best_candidate.fasta \
  --output interface/
```

Lists the interface residues on **both** chains (with per-residue contact counts,
minimum distance and pLDDT), reports the contact and PAE numbers, writes the
complex PDB, and renders a **PAE heatmap PNG with the chain boundary marked**.

## Interpreting the Output

### iPTM — the primary metric

| iPTM | Reading |
|---|---|
| **> 0.8** | Confident interface. A real hypothesis worth testing. |
| **0.5 – 0.8** | Possible interface. Ambiguous — do not over-claim. |
| **< 0.5** | Likely no binding. |

> [!WARNING]
> **iPTM is a structural-confidence proxy, NOT an affinity prediction.** It tells
> you how confident ESMFold2 is that it has placed the two chains correctly
> relative to each other. It does **not** give you a Kd, a ΔG, an IC50, or a
> rank-order of affinities. A high iPTM is a **hypothesis to test
> experimentally**, not a measurement. Never report an iPTM as if it were a
> binding constant, and never tell a user that iPTM 0.9 "binds more tightly"
> than iPTM 0.85.

### The other metrics

| Metric | Reading |
|---|---|
| `binder_plddt` (0–1) | Does the **binder itself** fold? > 0.7 confident. **A well-folded binder is not a binding binder** — see the negative example. |
| `interface_pae` (Å) | Mean cross-chain predicted aligned error. **Lower is better.** < 5 Å ≈ a confidently placed interface; > 20 Å means the model has no idea where the chains sit relative to each other. |
| `confident_interface_contacts` | Cross-chain residue pairs with CB–CB (CA for Gly) < 8 Å **and** PAE < 15 Å. Use this one. |
| `interface_contacts` | The same count **without** the PAE gate. Reported for transparency only — **it does not discriminate binders.** |
| `selection_score` | Composite, 0–1. See below. |

`selection_score = 0.5·iPTM + 0.3·binder_pLDDT + 0.2·(1 − min(interface_PAE, 31.75)/31.75)`

Every term is on a 0–1 scale where higher is better. iPTM carries half the weight
because it is the only term that measures the *interface*; binder pLDDT guards
against a binder that cannot fold; interface PAE breaks ties. **Contact count is
deliberately excluded** — it scales with binder size and does not separate
binders from non-binders. Use `selection_score` as a summary, but **rank on
`iptm`**: a well-folded non-binder can score a middling composite on the
pLDDT term alone (lysozyme scores 0.375 below).

### Worked example — a TRUE binder (barnase : barstar)

Barstar is barnase's natural inhibitor, one of the tightest known
protein–protein complexes (Kd ~ 10⁻¹⁴ M).

```
iPTM                       0.964   -> confident interface
pTM                        0.967
binder mean pLDDT          0.940
interface PAE              2.59 Å
confident contacts         52
geometric contacts         52       (all 52 survive the PAE gate)
interface residues         18 on barnase, 19 on barstar
selection_score            0.948
```

The interface it finds is the real one: barnase's active-site hotspots
**Lys27, Arg59, Arg83, His102** contacting barstar's **Asp35 / Asp39** binding
loop — Asp39 is the residue that plugs barnase's active site and mimics the
substrate phosphate. In the PAE heatmap **the entire matrix is dark**, including
the off-diagonal cross-chain blocks.

### Worked example — a NON-binder (barnase : lysozyme)  ← learn this one

Hen lysozyme is a well-folded, completely unrelated protein. It does not bind
barnase. This is what "no" looks like:

```
iPTM                       0.131   -> likely no binding
pTM                        0.572
binder mean pLDDT          0.887   <-- HIGH! it folds beautifully. It still does not bind.
interface PAE             24.86 Å  <-- the chains are not placed relative to each other
confident contacts         0
geometric contacts        53       <-- MORE than the true binder's 52 (!)
interface residues         0 on barnase, 0 on lysozyme
selection_score            0.375
```

Three lessons, all of them traps:

1.  **`binder_plddt` was 0.887.** The binder folds *better* than many real
    designs. Folding is not binding. Never conclude "it binds" from pLDDT.
2.  **It had 53 raw geometric contacts — more than barstar's 52.** ESMFold2 must
    output *some* coordinates, so it packs the two chains against each other
    regardless. A big non-binder racks up incidental proximity. **Raw contact
    count is not evidence of binding.** Once gated on PAE: **zero**.
3.  **iPTM (0.131) and interface PAE (24.9 Å) both called it correctly.** Trust
    these two.

In the PAE heatmap the two *diagonal* blocks are dark (each chain folds fine on
its own) while the *off-diagonal* cross-chain blocks are washed out and pale.
**That signature — confident chains, uncertain relative placement — is what a
non-binder looks like.** Learn to recognise it.

For reference, the third control (barstar with its sequence scrambled — same
amino-acid composition, fold destroyed) gives iPTM 0.139, binder pLDDT 0.315,
interface PAE 23.2 Å, 0 confident contacts. It fails on *both* folding and
binding.

## Common Mistakes

1.  **Ranking by contact count.** The single biggest trap, and it is measured:
    across folding settings, the scrambled decoy produced up to **72** raw
    cross-chain contacts and lysozyme up to **58**, against the true binder's
    stable **51–53**. Rank by raw contacts and you will order the non-binder
    over the nanomolar binder. Always use `confident_interface_contacts`, and
    rank on `iptm`.
2.  **Reading a high `binder_plddt` as binding.** It only says the binder folds.
    See the lysozyme example above.
3.  **Quoting iPTM as an affinity.** It is a confidence score. There is no Kd
    here. Say "predicted confident interface — a hypothesis to test", never
    "binds with high affinity".
4.  **Trying to design binders with this skill,** or trying to reach the paper's
    gradient-descent loop / distogram iPTM proxy through the API. Both are
    impossible remotely (see *Scope*). Generate candidates with
    `esm3-protein-design` / `esm3-inverse-folding` first.
5.  **Screening a huge FASTA without thinking.** Every candidate is one fold
    call. Concurrency is capped at 4 workers on purpose. Triage a long list
    first, or expect to wait.

## Dependencies

*   **`esm3-protein-design`** and **`esm3-inverse-folding`** — to *generate* the
    candidate binders that this skill *ranks*.
*   **`esmfold2-structure-prediction`** — to fold a single chain with no partner.
*   **`credentials`** — the safe protocol for `BIOHUB_API_KEY`.
*   **`uv`** — environment setup.

## References

*   [Interpretation guide](docs/interpretation-guide.md) — how to read iPTM,
    interface PAE and contacts; the "iPTM is not an affinity" caveat; the
    Negative Results & Scientific Integrity section; and the pre-report
    checklist. **Read this before interpreting a screen.**
*   [Report template](docs/report-templates.md) — the scaffold for the written
    report (ranked table + per-candidate verdicts).
*   [Worked example — a TRUE binder (barnase : barstar)](docs/examples/barnase_barstar/report.md)
    — iPTM 0.964, interface PAE 2.65 Å, 52 confident contacts; the positive
    control to read every screen against.
*   [Worked example — decoy controls](docs/examples/decoy_controls/report.md) —
    what "no binding" looks like: a well-folded non-binder (lysozyme, iPTM 0.230)
    and a composition-matched scramble (iPTM 0.144) ranked below the true binder
    and filtered out.
*   [ESM/Biohub API reference](references/esm-biohub-api.md) — endpoints,
    reachable models, the daily-credit model, and verified gotchas shared across
    every ESM skill.
*   `references/citation.bib` — cite the ESMFold2/ESMC paper if you use this.
*   Source protocol: `esm/cookbook/tutorials/binder_design.ipynb` (the selection
    stage, "Pick the designs to order").
