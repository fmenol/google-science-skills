---
name: esm3-secondary-structure-sasa
description: >
  Predicts per-residue secondary structure (8-state DSSP: helix H/G/I, strand
  E/B, turn/bend/coil T/S/C) and solvent-accessible surface area (SASA, in Å²)
  directly from a protein sequence with ESM3 — no structure and no folding
  required. Use when the user wants secondary-structure content, a helix/strand
  topology, buried-vs-exposed residues, or a solvent-accessibility profile from
  sequence alone, or asks "is this region helical/a strand?" or "which residues
  are buried?". Do NOT use to fold a protein or get 3D coordinates (that is
  esmfold2-structure-prediction), and do NOT use to extract SS/SASA from a PDB
  you already have (that is esm-protein-tracks).
---

# ESM3 Secondary Structure & SASA Prediction

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If
    .licenses/esm3_secondary_structure_sasa_LICENSE.txt does not already exist
    in the workspace root directory then (1) prominently notify the user to
    check the terms at https://biohub.org/acceptable-use-policy/ and
    https://biohub.ai/, then (2) create the file recording the notification text
    and timestamp.
3.  **`BIOHUB_API_KEY`**: This skill requires an API key. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for it before
    running anything.

## Overview

ESM3 was trained on secondary structure and solvent accessibility as native
prediction tracks. This skill queries those tracks directly, so it annotates a
sequence with:

-   **8-state secondary structure (SS8)** — the full DSSP alphabet: `H` α-helix,
    `G` 3₁₀-helix, `I` π-helix, `E` β-strand, `B` β-bridge, `T` turn, `S` bend,
    `C` coil. A 3-state (helix/strand/coil) collapse is also reported.
-   **Solvent-accessible surface area (SASA)** — per-residue, in Å², with a
    coarse buried / intermediate / exposed label.

Because it does not fold the protein, it is fast and works on sequences of any
length, including ones with no known or predictable structure.

**Do NOT use when:**

-   The user wants a 3D structure or coordinates → **esmfold2-structure-prediction**.
-   The user already has a PDB/mmCIF and wants SS/SASA measured from it →
    **esm-protein-tracks** (which computes them geometrically from coordinates).
-   The user wants domain/function annotation → **esm3-function-prediction**.

## Core Rules

-   **NEVER run `python`/`python3` directly.** Always `uv run --no-project`.
    Without `--no-project`, `uv` finds an unrelated parent `pyproject.toml` and
    fails.
-   **These are model PREDICTIONS, not measurements.** ESM3 predicts SS8 and
    SASA from sequence; it does not compute them from a structure. Present them
    as predictions and, when a structure is available, prefer measuring them
    from it (via `esm-protein-tracks`).
-   **SS8 here is genuine 8-state**, unlike `esm-protein-tracks`, whose biotite
    fallback only yields a 3-state (H/E/C) approximation. If a user needs true
    SS8 without a structure, this is the skill.
-   **Absolute SASA depends on residue size.** A large residue exposes more area
    than a small one at the same relative burial, so treat the buried/exposed
    labels as a convenience, not ground truth. Report the numeric Å² values.
-   Do not eyeball or hand-compute compositions; use the script's summary.
-   If this skill is used, mention it in the output.

## Utility Scripts

All commands write results to a file and print a one-line status to stderr.

**1. Predict for one sequence (or a small FASTA)**

```bash
uv run --no-project scripts/predict_ss_sasa.py predict \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --output ubiquitin_ss.json
```

The JSON holds `secondary_structure_ss8`, `secondary_structure_ss3`, the
per-residue `sasa` and `burial` arrays, and a `summary` with helix/strand/coil
percentages and buried/exposed counts.

**2. Batch a FASTA**

```bash
uv run --no-project scripts/predict_ss_sasa.py batch \
  --fasta proteins.fasta --output preds.json --summary-csv preds.csv
```

**3. Plot the SS8 ribbon over the SASA profile**

```bash
uv run --no-project scripts/predict_ss_sasa.py plot \
  --input ubiquitin_ss.json --output ubiquitin_ss.png
```

Optional flags: `--model` (default `esm3-open-2024-03`, the only reachable ESM3
model), `--temperature` (default 0.7; lower is more conservative).

## Interpreting the Output

Before interpreting, read the [interpretation
guide](docs/interpretation-guide.md) and skim the worked example in
[`docs/examples/fold_topology/`](docs/examples/fold_topology/report.md) so you
know what a real β-grasp and an all-α fold look like.

-   **SS8 → topology.** Runs of `H` are α-helices; runs of `E` are β-strands.
    `G`/`I` are rarer helix types; `T`/`S`/`B` are short local features. Collapse
    to SS3 (helix/strand/coil) for a coarse fold class — e.g. mostly `H` with no
    `E` is an all-α fold; alternating `E` and `C` is a β-sheet protein.
-   **SASA → burial.** Low SASA (≤ 20 Å²) marks a buried core residue; high SASA
    (≥ 50 Å²) marks an exposed surface residue. A protein with a real hydrophobic
    core will show many buried residues; an all-exposed profile suggests a
    peptide, an extended/disordered region, or low confidence.
-   Hand off a predicted SS8/SASA to **esm-protein-tracks** `build-prompt` to
    condition an **esm3-protein-design** run on a desired topology.

## Dependencies

-   `esm-protein-tracks` — the complementary skill that *measures* SS8/SASA from
    a known structure and builds masked design prompts.
-   `esmfold2-structure-prediction` — fold the sequence if 3D coordinates (and
    structure-derived SS/SASA) are wanted instead of a sequence-only prediction.

## Common Mistakes

-   Using this to "get the structure" — it returns per-residue labels, not
    coordinates. Fold with `esmfold2-structure-prediction` for 3D.
-   Reporting the buried/exposed labels as if they were relative accessibility —
    they are thresholds on absolute Å² and vary with residue size.
-   Treating SS8 as measured DSSP — it is ESM3's prediction from sequence.

## References

-   [Interpretation guide](docs/interpretation-guide.md) — how to read SS8 and
    SASA, the burial bands, and the prediction-vs-measurement caveat.
-   [Worked example: fold topology (β-grasp ubiquitin vs all-α myoglobin)](docs/examples/fold_topology/report.md)
    — contrasts a mixed α/β and an all-α protein with real numbers.
-   [ESM/Biohub API reference](references/esm-biohub-api.md) — endpoints, the
    reachable ESM3 model, the `generate` track shapes, and the credit model.
-   `references/citation.bib` — cite ESM3 (Hayes et al. 2025) and DSSP (Kabsch &
    Sander 1983) if you use this skill.
