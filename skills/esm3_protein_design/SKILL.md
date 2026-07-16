---
name: esm3-protein-design
description: >
  Designs new protein sequences with ESM3, conditioned on any combination of a
  partial sequence, a 3D structural motif, a secondary-structure (SS8) plan and
  a per-residue solvent-accessibility profile. Use when the user wants to
  design, generate, invent, or create a protein; scaffold or graft a functional
  motif (an active site, a binding epitope) into a new protein; build a de novo
  binder backbone; inpaint or redesign a loop or a region of an existing
  protein; or make a sequence that adopts a requested fold (e.g. "design a
  44-residue helix-loop-helix"). Every design is folded and QC'd before it is
  reported. Do not use when the user wants the structure of a sequence they
  already have (use `esmfold2-structure-prediction`); when they have a backbone
  and want a sequence for it, with no new scaffold to invent (use
  `esm3-inverse-folding`); or when they want to score or rank mutations of an
  existing protein (use `esmc-mutation-effect-scoring`).
---

# Protein Design using ESM3

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If `.licenses/esm3_protein_design_LICENSE.txt` does
    not already exist in the workspace root directory then (1) prominently
    notify the user to check the terms at
    https://biohub.org/acceptable-use-policy/ and https://biohub.ai/, then
    (2) create the file recording the notification text and timestamp.
3.  **`BIOHUB_API_KEY`**: This skill requires an API key to function. Register
    at https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for and request
    this key if this skill looks relevant to the user's request.

## Overview

ESM3 is a generative model over three tracks at once — sequence, structure, and
function. You give it a partly-specified protein and it fills in the blanks.
This skill exposes that as promptable design: mask what you want invented, pin
down what you already know.

What you can condition on, in any combination:

| Prompt | How it is expressed | Masked with |
|---|---|---|
| Partial sequence | `MKT___________AVL` | `_` |
| 3D structural motif | atom37 coordinates `(L, 37, 3)` | `NaN` |
| Secondary structure | SS8 string over `GHITEBSC` | — |
| Solvent accessibility | per-residue SASA in Å² | `None` |

What this skill does **not** do:

*   It does **not** predict the structure of a sequence you already have. That
    is `esmfold2-structure-prediction`.
*   It does **not** do plain structure → sequence recovery on a backbone you
    already have. That is `esm3-inverse-folding`.
*   It does **not** score or rank point mutations. That is
    `esmc-mutation-effect-scoring`.
*   It does **not** predict function or GO terms. That is
    `esm3-function-prediction`.

## Core Rules

-   **ALWAYS fold and QC a design before reporting it. NEVER present a generated
    sequence as a "design" without its pTM and pLDDT.** A raw ESM3 sample is a
    guess, not a design. `design.py` folds every sample with ESMFold2 and ranks
    by pTM automatically — report those numbers. A sequence quoted without them
    is a scientific error.
-   **Before interpreting any result, read
    [`docs/interpretation-guide.md`](docs/interpretation-guide.md).** It defines
    how to read pTM/pLDDT/motif RMSD, the designed-vs-random floor logic, the
    temperature and low-complexity checks, and the integrity rules.
-   **Review at least one worked example in [`docs/examples/`](docs/examples/)
    before writing a report — including the negative one
    ([`random_baseline`](docs/examples/random_baseline/report.md)) so you know
    what the foldability floor looks like** and never oversell a design that
    barely clears it.
-   **Write the report using
    [`docs/report-templates.md`](docs/report-templates.md).** Report the whole
    ranked set with QC metrics, compare de novo designs to the random floor, and
    report the motif CA RMSD for scaffolds.
-   **The only ESM3 model this API key can reach is `esm3-open-2024-03`**
    (aliases `esm3-sm-open-v1`, `esm3-open`). The upstream cookbook tutorials
    (`esm3_generate.ipynb`, `gfp_design.ipynb`) call `esm3-medium-2024-08` and
    `esm3-medium-2024-03`; **those return HTTP 403 — no access.** Everything
    here is re-pointed at `esm3-open-2024-03`. Do not "fix" a 403 by retrying a
    medium/large model.
-   **NEVER run `python3` or `python3 -c` directly, and NEVER `pip install`.**
    ALWAYS `uv run --no-project`. The `--no-project` flag is mandatory: without
    it `uv` walks up the directory tree, finds an unrelated `pyproject.toml`,
    and tries to build that project instead.
-   **NEVER download model weights.** No `torch`, no `transformers`, no
    `huggingface_hub`, no `esm` PyPI package. All inference is remote.
-   **ALWAYS draw more than one sample.** Design is stochastic and the spread is
    large. `--num-samples` is required for exactly this reason. Draw 4–8, keep
    the best by pTM, and say how many you drew.
-   **Do not compute pTM, pLDDT or RMSD yourself; always use the script's
    output.** The metrics JSON is the source of truth.
-   **Notification**: If this skill is used, ensure this is mentioned in the
    output.

## Utility Scripts

All commands write `<prefix>.fasta` (designs ranked by pTM), `<prefix>.json`
(full metrics) and one PDB per design.

### `generate` — de novo design, or inpaint a partial sequence

```bash
uv run --no-project scripts/design.py generate \
  --length 60 --num-samples 8 --output-prefix out/denovo
```

Inpainting: keep what you know, mask what you want invented.

```bash
uv run --no-project scripts/design.py generate \
  --sequence "MKTAYIAKQR____________________QITLGGRLSHQ" \
  --num-samples 8 --output-prefix out/loop
```

### `scaffold-motif` — graft a 3D motif into a bigger protein

This is the flagship workflow: hold an active site or epitope rigid in 3D and
let ESM3 invent a whole protein around it.

```bash
uv run --no-project scripts/design.py scaffold-motif \
  --motif-pdb 1ITU.pdb --chain A --motif-range 124-146 \
  --scaffold-length 200 --motif-start 73 \
  --num-samples 8 --output-prefix out/scaffold
```

`--motif-range` is **author residue numbering, 1-indexed inclusive** — the
numbers you read off in PyMOL. `--motif-start` is **1-indexed** too: the
position in the *new* protein where the motif begins. The motif's residues are
pinned in the sequence prompt and its atom37 coordinates are pinned in the
structure prompt; everything else is invented.

### `design-with-ss` — design a sequence for a requested fold

```bash
uv run --no-project scripts/design.py design-with-ss \
  --ss8 CCHHHHHHHHHHHHHHHHCCCCCCHHHHHHHHHHHHHHHHHHCC \
  --num-samples 8 --output-prefix out/hlh
```

SS8 alphabet is `GHITEBSC`: **H** α-helix, **E** β-strand, **C** coil, **T**
turn, **G** 3₁₀-helix, **I** π-helix, **B** β-bridge, **S** bend. The design is
exactly as long as the SS8 string.

### `design-with-sasa` — design a burial / exposure profile

```bash
uv run --no-project scripts/design.py design-with-sasa \
  --length 100 --sasa "30-40:5,60-70:90" \
  --num-samples 8 --output-prefix out/sasa
```

SASA is per-residue, in Å²: **buried < 20**, intermediate 20–50, **exposed >
50**. Unlisted positions are unconstrained. Use `--sasa-json` for a full
length-L list where `null` means unconstrained.

### `chain-of-thought` — decode tracks one at a time

```bash
uv run --no-project scripts/design.py chain-of-thought \
  --length 60 --tracks secondary_structure,structure,sequence \
  --num-samples 8 --output-prefix out/cot
```

Each track conditions on all the previous ones, so ESM3 commits to a fold plan
before it writes a single residue. In a spot check at length 60 (**n = 1**, so
treat it as indicative, not established) the three-track chain reached pTM 0.65
/ pLDDT 0.78, against a mean of pTM 0.32 for direct sequence sampling (n = 8).
It is worth trying for unconditioned de novo design of small proteins — but draw
several samples and judge by the pTM you actually measure, not by this number.
The chain must end with `sequence`.

## Interpreting the Output

Read these off the metrics JSON. Do not recompute them.

| Quantity | Reading |
|---|---|
| pLDDT (0–1) | > 0.9 very high · 0.7–0.9 confident · 0.5–0.7 low · < 0.5 disordered |
| pTM | > 0.8 confident fold · 0.5–0.8 plausible · < 0.5 unreliable |
| Motif CA RMSD | < 1.0 Å excellent · 1.0–2.0 Å good · > 3.0 Å the graft failed |

**Calibration, measured against this API — quote these when judging a design.**
Sample sizes are given because they are small; treat them as reference points,
not constants.

| Setting (length 60 unless noted) | pTM | pLDDT | n |
|---|---|---|---|
| Ubiquitin, a real protein (the ceiling) | 0.78 | 0.82 | 1 |
| **Uniformly-random 60-mer (the floor)** | **0.15** | **0.41** | 2 |
| `generate`, unconditioned, T=0.5 | 0.32 (range 0.19–0.49) | 0.62 | 8 |
| `generate`, unconditioned, T=0.7 | 0.28 (range 0.12–0.41) | 0.57 | 8 |
| `chain-of-thought`, 3 tracks | 0.65 | 0.78 | **1** |
| `scaffold-motif` into a 200-mer | 0.89–0.98 | 0.87–0.97 | 4 |

Conditioning is what makes ESM3 good: the more you pin down, the better the
design. Motif CA RMSD for those scaffolds was 0.22–0.31 Å.

**Therefore: a low pTM on a short unconditioned design is expected, not a bug.**
Say so plainly instead of overselling the design. If the user needs a
confidently-folded protein, condition it (motif / SS8 / SASA) or use
`chain-of-thought`, and design longer — pTM is harshly penalised below ~80
residues.

A design is only worth reporting if it clears **all** of these:

1.  pTM and pLDDT are quoted, and pTM is above the random floor for its length.
2.  For `scaffold-motif`: the motif is verbatim in the sequence **and** its CA
    RMSD is low. A design that keeps the residues but loses the geometry has
    failed.
3.  It is not low-complexity junk (see below).

## Workflow Checklist

```
Protein Design Progress:
- [ ] Step 1: Create output folder; confirm BIOHUB_API_KEY via the credentials skill
- [ ] Step 2: Decide what is FIXED and what is DESIGNED. Pick the subcommand from
              the conditioning you actually have (motif -> scaffold-motif;
              fold sketch -> design-with-ss; nothing -> chain-of-thought)
- [ ] Step 3: Draw MULTIPLE samples (--num-samples 4-8). Never design once.
- [ ] Step 4: QC (automatic): every design folded with ESMFold2, ranked by pTM
- [ ] Step 5: Read the metrics JSON. Compare pTM/pLDDT against the random floor
              and the calibration table above (MANDATORY before any claim)
- [ ] Step 6: Sanity-check the top design for low-complexity / membrane artefacts
- [ ] Step 7: Report the top design WITH its pTM, pLDDT (and motif RMSD), the
              number of samples drawn, and the model name (esm3-open-2024-03)
- [ ] Step 8: State the limitations. These are computational designs, not
              validated proteins.
```

## Common Mistakes

-   **Reporting a sequence with no pTM/pLDDT.** The single most common failure.
    An unfolded, unranked ESM3 sample is not a design. Always QC.
-   **Designing once and believing it.** Sampling variance is large: at length
    60, individual de novo designs range from pTM 0.12 to 0.49. One sample tells
    you nothing. Draw 4–8 and rank.
-   **Expecting good unconditional designs at short lengths.** ESM3's
    unconditioned 60-mers are frequently low-complexity (poly-Ala runs like
    `LALAALAATHLAARLLGLLAAAALAAALG...`) or look like signal peptides and
    transmembrane helices. This is a real property of the model, not a bug in
    the script. Inspect the top design for repeats before reporting it; if it is
    junk, condition the design or use `chain-of-thought`.
-   **Cranking up the temperature for "more creative" designs.** It does the
    opposite. Measured at length 60: T=0.5 → mean pTM 0.32 (n=8); T=0.7 → 0.28
    (n=8); a single T=1.0 draw gave 0.19, barely above the random floor. The
    default (0.5) is both the tutorials' value and the best measured here. Do
    not exceed ~0.7.
-   **Retrying a 403 against `esm3-medium-*` or `esm3-large-*`.** The key has no
    access to them and never will. Use `esm3-open-2024-03`.
-   **Off-by-one on the motif.** `--motif-range` uses PDB *author* residue
    numbers and `--motif-start` is 1-indexed. Verify against
    `motif_sequence` in the output JSON before trusting the graft.
-   **Burning the daily credit budget.** The API enforces a **daily credit
    limit (100 on this key)** and returns HTTP 429 with
    `"You have exceeded your daily credit limit"` once it is spent. Each
    `generate` and each `fold` costs a credit, so `--num-samples 8` costs 16.
    `encode` is *not* billed. Note that `esm_biohub` reports this 429 with the
    hint "Rate limited even after backoff. Lower BIOHUB_QPS and retry" — that
    hint is **wrong** for a quota 429: lowering QPS will not help and retrying
    just wastes time. Read the server message, not the hint. Budget your samples
    before you start; the quota is shared across every ESM skill using the key.

## Dependencies

*   `credentials` — for `BIOHUB_API_KEY`.
*   `uv` — for running the scripts.
*   Sibling skills, cross-referenced rather than reimplemented:
    `esmfold2-structure-prediction` (fold a known sequence),
    `esm3-inverse-folding` (backbone → sequence),
    `esmc-mutation-effect-scoring` (score mutations),
    `esm3-function-prediction` (predict function).

## References

*   [Interpretation guide](docs/interpretation-guide.md) — how to QC a design and
    read pTM/pLDDT/motif RMSD; the designed-vs-random floor logic; the
    "Negative Results & Scientific Integrity" section; the pre-report checklist.
*   [Report template](docs/report-templates.md) — the scaffold to fill for a de
    novo campaign and for a motif-scaffolding campaign.
*   [Worked example: de novo design](docs/examples/denovo_design/report.md)
    (POSITIVE) — five 60-mers folded; designed pTM 0.259 vs random 0.196.
*   [Worked example: motif scaffolding](docs/examples/motif_scaffold/report.md)
    (POSITIVE) — the 1ITU motif grafted into a 200-mer at pTM 0.981, motif CA
    RMSD 0.26 Å, motif verbatim.
*   [Worked example: random baseline](docs/examples/random_baseline/report.md)
    (NEGATIVE) — random strings fold to pTM 0.196; the floor a real design must
    clear.
*   [ESM/Biohub API reference](references/esm-biohub-api.md) — endpoints,
    reachable models, credit model, and verified gotchas for the platform.
*   [`scripts/design.py`](scripts/design.py) — the CLI. Five subcommands:
    `generate`, `scaffold-motif`, `design-with-ss`, `design-with-sasa`,
    `chain-of-thought`.
*   [`references/citation.bib`](references/citation.bib) — cite
    `hayes2024simulating` (ESM3) for any design produced with this skill.
