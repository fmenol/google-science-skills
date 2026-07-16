---
name: esm3-inverse-folding
description: >
  Design amino-acid sequences that fold into a given 3D backbone (inverse
  folding / fixed-backbone sequence design), then validate them by re-folding.
  Use when the user has a structure file (PDB/mmCIF) or a backbone coordinate
  array and wants sequences for it — "design a sequence for this backbone",
  "what sequence would fold into this structure", "redesign this protein",
  "sequence recovery", "inverse fold this PDB", "self-consistency / scTM".
  Do not use when the user has a sequence and wants its structure — that is
  sequence-to-structure prediction, use the `esmfold2-structure-prediction`
  skill. Do not use when the user wants a novel protein invented from scratch
  with no input backbone — use the `esm3-protein-design` skill.
---

# ESM3 Inverse Folding

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If `.licenses/esm3_inverse_folding_LICENSE.txt` does
    not exist in the workspace root then (1) prominently notify the user to
    check the terms at https://biohub.org/acceptable-use-policy/ and
    https://biohub.ai/, then (2) create the file recording the notification text
    and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for it.

## Overview

Structure in, sequence out. Given a 3D backbone this skill samples amino-acid
sequences predicted to fold into it (ESM3), then **proves** they do by re-folding
each design with an independent model (ESMFold2) and measuring how well the
refolded design reproduces the target backbone — the self-consistency loop
(scTM / scRMSD). This is the standard way the field validates an inverse-folding
design, and it is the only thing that makes a design trustworthy.

**Do NOT use when:**

-   The user has a **sequence** and wants its **structure**. That is the
    opposite direction — use **`esmfold2-structure-prediction`**.
-   The user wants a **novel protein designed from scratch** with no input
    backbone (unconditional or prompt-only generation) — use
    **`esm3-protein-design`**.
-   The user wants to score point mutations in an existing sequence — use
    **`esmc-mutation-effect-scoring`**.

## Core Rules

-   **The `/inverse_fold` endpoint does not exist for this API key.** Every
    reachable model (`esm3-open-2024-03`, `esmfold2-*`, `esmc-*`) answers
    HTTP 422: *"Model '…' does not support the 'inverse_fold' endpoint."* The ESM
    SDK's `inverse_fold()` / `InverseFoldingConfig` path is therefore
    **unreachable** — do not reach for it, and do not "fix" it. The supported
    route, which this skill's script uses, is the ESM3 generate endpoint with
    coordinates supplied and the sequence track left masked:

    ```python
    client.generate('sequence', model='esm3-open-2024-03',
                    coordinates=coords_atom37, num_steps=L, temperature=0.1)
    # -> result['outputs']['sequence']
    ```

-   **NEVER report, rank, or hand over a design without its self-consistency
    scores.** An unvalidated inverse-folding design is worthless: the model will
    always return a plausible-looking sequence, including for a backbone it has
    completely failed to solve. Run `self-consistency`, not just `design`, before
    you tell the user anything about a sequence.
-   **Interpret scTM strictly**:
    -   **scTM > 0.8** — the design very likely adopts the target fold.
    -   **scTM 0.5–0.8** — uncertain; the fold is only partially recovered. Do
        not present these as successes.
    -   **scTM < 0.5** — a **failed** design. Say so plainly.
-   **Use the Wrapper**: ALWAYS run `scripts/inverse_fold.py`. Do not compute
    recovery, scTM, scRMSD or diversity yourself — always use the script's
    output. Do not hand-roll HTTP calls to the API.
-   **Never download model weights.** No `torch`, no `transformers`, no `esm`
    PyPI package. Everything runs remotely on the Biohub Platform.
-   **A CA-only trace silently produces garbage.** ESM3 builds its structure
    representation from each residue's **N–CA–C frame**. Hand it a CA-only model
    and the API does *not* error — it returns a fluent sequence with ~5% identity
    to the native, i.e. pure noise (measured on ubiquitin: 0.046). The script
    refuses such input; **do not work around that check**, fix the input.
-   **Get diversity by resampling, not by raising the temperature.** Each call is
    an independent stochastic decode, so ask for more samples. Raising the
    temperature past ~1.5 destroys the design rather than diversifying it.
-   **Before interpreting or reporting, read
    [`docs/interpretation-guide.md`](docs/interpretation-guide.md) and review the
    worked example in
    [`docs/examples/self_consistency/`](docs/examples/self_consistency/report.md)** —
    it shows the full design → self-consistency → recovery loop with real numbers
    and the decoy control that makes an scTM of 1.000 meaningful.
-   If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All three subcommands take the backbone as either `--pdb` (with optional
`--chain`) or `--coords` (a `.npy` of shape `(L, 37, 3)` atom37, `(L, 4, 3)`
N/CA/C/O, or `(L, 3, 3)` N/CA/C — the layouts backbone generators emit).
`--num-samples` is **required**: you must choose it explicitly.

**1. `self-consistency` — design and validate. This is the one you want.**

Designs N sequences, re-folds each with ESMFold2, scores scTM / scRMSD / pLDDT
against the input backbone, and ranks them best-first. Writes JSON + a ranked CSV
+ FASTA.

```bash
uv run --no-project scripts/inverse_fold.py self-consistency \
  --pdb ./target_backbone.pdb --chain A \
  --num-samples 8 --temperature 0.1 \
  --output ./results/designs.json
```

**2. `design` — sequences only (fast, UNVALIDATED).**

Use only when you will validate separately, or when you genuinely just need
candidates. Reports per-design native-sequence recovery (if the input PDB carries
a sequence) and pairwise diversity. Writes JSON + FASTA.

```bash
uv run --no-project scripts/inverse_fold.py design \
  --pdb ./target_backbone.pdb \
  --num-samples 16 --temperature 1.0 \
  --output ./results/candidates.json
```

**3. `recovery` — how well do the designs reproduce the native sequence?**

Per-position and overall native-sequence recovery, plus a per-position agreement
plot (PNG). Needs a native sequence: read from `--pdb`, or given via `--native`.
This is a benchmark/diagnostic, only meaningful when the backbone came from a
real protein.

```bash
uv run --no-project scripts/inverse_fold.py recovery \
  --pdb ./native_protein.pdb \
  --num-samples 8 --temperature 0.1 \
  --output ./results/recovery.json
```

Useful options: `--temperature` (default `0.1`), `--num-steps` (default
`min(L, 100)`; must be ≤ L, and the API caps it at 100), `--model` (default
`esm3-open-2024-03`), `--fold-model` (default `esmfold2-fast-2026-05`),
`--all-atom` (condition on sidechains too — verified to be a no-op, see below).

## Interpreting the Output

Read the script's numbers; do not recompute them.

| Quantity | Reading |
|---|---|
| **scTM** | **> 0.8** design very likely adopts the target fold · **0.5–0.8** uncertain · **< 0.5** failed design |
| **scRMSD** (CA, Å) | < 2 Å is the practical success bar used in the design field · > 5 Å means the fold was not reproduced |
| **pLDDT** (0–1) | Confidence of the *re-fold*. > 0.9 very high · 0.7–0.9 confident · < 0.5 disordered. A design with high scTM but pLDDT < 0.7 is a weak result |
| **native recovery** | Fraction of positions matching the native sequence. Chance is ~5%. ProteinMPNN/ESM-IF report ~50% on native backbones; on an ESM-predicted backbone at T=0.1, ESM3 returns the native sequence essentially exactly |
| **mean pairwise identity** | Across designs. High (> 0.95) = the backbone strongly determines the sequence, or the temperature is too low for the diversity you asked for |

Judge a design on **scTM first**, then scRMSD, then pLDDT. Recovery is *not* a
quality score — for a de novo backbone there is no native sequence to recover,
and a design can be excellent at 30% recovery.

**Temperature ladder** (measured on ubiquitin's own predicted backbone):

| T | Native recovery | Use |
|---|---|---|
| 0.1 | ~100% | Maximum fidelity. The default. |
| 1.0 | ~98% | Adds real but modest diversity. |
| 1.5 | ~90% | Aggressive. |
| 2.0 | ~40% | **Broken** — barely above chance. Do not use. |

## Common Mistakes

1.  **Reaching for `inverse_fold()`.** Everyone who reads the ESM SDK tries the
    dedicated inverse-folding client first, and every model rejects it with
    HTTP 422. Use `generate(track='sequence', coordinates=...)`. See Core Rules.
2.  **Reporting a design without re-folding it.** `design` alone tells you
    nothing about whether the sequence works. ESM3 returns a confident-looking
    sequence no matter how badly it has failed. Always run `self-consistency`.
3.  **Feeding a CA-only trace or a broken backbone.** The API accepts it and
    returns noise with no warning. The script blocks this; heed it.
4.  **Cranking the temperature for diversity.** Diversity comes from drawing more
    samples (`--num-samples`), not from a hotter decode. Past T≈1.5 you are
    destroying the design, not exploring it.
5.  **Treating high native recovery as proof of a good design.** It only means
    the design matches the *native* sequence. For a de novo backbone recovery is
    undefined and irrelevant — scTM is the score that matters.

## Notes on the conditioning

ESM3 reads only the **N, CA, C** frame of each residue. Verified against the live
API on ubiquitin: conditioning on all 37 atoms, on N/CA/C/O, or on N/CA/C alone
returns byte-identical designs at 100% recovery. The script therefore strips the
input to the backbone by default, which is a numeric no-op but *guarantees* that
no sidechain identity can leak into the design — so a high recovery number is
real inverse folding, not the model reading the answer off the sidechains. Pass
`--all-atom` to send everything anyway.

## Dependencies

-   **`esmfold2-structure-prediction`** — the reverse direction (sequence →
    structure), and the model this skill uses internally to validate designs.
-   **`esm3-protein-design`** — de novo generation when there is no input
    backbone.
-   **`credentials`** — for handling `BIOHUB_API_KEY`.
-   **`uv`** — for running the scripts.

## References

-   [Interpretation guide](docs/interpretation-guide.md) — how to read scTM /
    scRMSD / pLDDT, recovery vs chance, diversity via resampling, the
    `/inverse_fold`-unsupported fact, and a pre-report checklist.
-   [Worked example: ubiquitin self-consistency](docs/examples/self_consistency/README.md)
    — the full `design` → `self-consistency` → `recovery` loop with real numbers
    (native recovery 100%, scTM 1.000) and the decoy control (scTM 0.078); the
    analysis is in its [report.md](docs/examples/self_consistency/report.md).
-   [ESM/Biohub API reference](references/esm-biohub-api.md) — shared endpoints,
    reachable models, request/response shapes, the credit model, and the verified
    gotchas (including `/inverse_fold`) that apply to every ESM skill.
-   [`references/citation.bib`](references/citation.bib) — ESM3, ESMFold, the scTM
    protocol, and the inverse-folding baselines to compare against.
