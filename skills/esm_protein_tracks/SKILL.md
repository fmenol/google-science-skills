---
name: esm-protein-tracks
description: >
  Extract ESM3's conditioning tracks — sequence, 3D coordinates, secondary
  structure (SS8), solvent accessibility (SASA) and function annotations — from
  an experimental structure, inspect them, and build a masked PROMPT for ESM3.
  Use when the user has a PDB ID or a local .pdb/.cif and wants to "keep the
  active site and mask the rest", scaffold a motif, pin a binding site, see
  which residues are buried or exposed, or otherwise prepare conditioning
  inputs for a design run. This is the INPUT-PREPARATION step and nothing else.
  Do not use it to predict a structure (use esmfold2-structure-prediction), to
  design or generate a sequence (use esm3-protein-design), to recover a
  sequence from a backbone (use esm3-inverse-folding), or to analyse protein
  function (use esm3-function-prediction) — this skill only PREPARES the inputs
  those skills consume.
---

# ESM Protein Tracks

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If `.licenses/esm_protein_tracks_LICENSE.txt` does not
    exist in the workspace root then (1) prominently notify the user to check the
    terms at https://biohub.org/acceptable-use-policy/ and https://biohub.ai/,
    then (2) create the file recording the notification text and timestamp.
3.  **`BIOHUB_API_KEY`**: **NOT required for this skill.** Extracting, inspecting
    and masking tracks is pure local computation — no model, no network, no
    weights. A key is needed *only* for the optional `--predict-function` flag on
    `extract`. If you do use that flag, you **MUST** use the safe credentials
    protocol in the `credentials` skill to check for the key. Register at
    https://biohub.ai/developer-console/api-keys.

## Overview

ESM3 is promptable on five parallel, residue-aligned tracks. You can hold any
subset fixed and let the model fill in the rest:

| Track | Type | Meaning |
|---|---|---|
| `sequence` | string, length L | The amino acids. |
| `coordinates` | `(L, 37, 3)` | All-atom (atom37) 3D structure. |
| `secondary_structure` | string, length L | SS8 alphabet `GHITEBSC`. |
| `sasa` | list, length L | Solvent-accessible area per residue, Å². |
| `function` | `[label, start, end]` | Function/domain annotations. |

This skill turns a **real experimental structure** into those tracks, lets you
look at them, and masks them into a **prompt** — the exact input that
`esm3-protein-design` and `esm3-inverse-folding` consume. It is the natural
first step of any conditioned-design workflow: *"take this PDB, keep the active
site, mask everything else, give me the prompt."*

**This skill does NOT:**

-   Predict, fold or generate any structure → use **`esmfold2-structure-prediction`**.
-   Design, generate, mutate or optimise any sequence → use **`esm3-protein-design`**.
-   Recover a sequence from a backbone → use **`esm3-inverse-folding`**.
-   Analyse or interpret protein function → use **`esm3-function-prediction`**.
-   Compute embeddings → use **`esmc-protein-embeddings`**.

It only **prepares inputs**. If the user's goal is a designed protein, this skill
produces the prompt and a sibling skill does the design.

## Dependencies

| Sibling skill | Relationship |
|---|---|
| **`esm3-protein-design`** | **Consumer.** Feeds a masked prompt from `build-prompt` to ESM3 to generate sequence and/or structure. |
| **`esm3-inverse-folding`** | **Consumer.** Takes a coordinates-only prompt (`--condition-on coordinates`) and recovers a sequence. |
| **`esm3-function-prediction`** | **Owns the function track.** Use it for any real function analysis. `extract --predict-function` here is only a convenience that fills the `function_annotations` field; it is not a substitute. |

## Core Rules

-   **Use the script.** ALWAYS run `scripts/tracks.py`. Never parse the PDB,
    compute SASA, assign secondary structure, or apply a mask by hand.
-   **Start from the worked example.** For an end-to-end run and the masking /
    SS8 conventions, see [`docs/examples/ubiquitin_tracks/`](docs/examples/ubiquitin_tracks/report.md)
    and [`docs/interpretation-guide.md`](docs/interpretation-guide.md).
-   **Do not compute the masking yourself; always use the script's output.** The
    script reports exactly how many positions it masked on each track. Report
    those numbers. Do not recount them by eye.
-   **Ranges are 1-INDEXED and INCLUSIVE.** `--keep 10-20` keeps **11** residues:
    the 10th through the 20th. This is the biology convention, and it is the
    single easiest thing to get wrong. The script prints the kept positions and
    the kept author residue IDs — check them.
-   **Know which numbering you are using.** By default ranges are **track
    positions**, `1..L` along the extracted sequence. Many PDB files do not start
    at residue 1 (calmodulin 1CM4 starts at 4). To address the file's own author
    numbering, pass `--index residue-id`. If the user quotes residue numbers
    *from a paper or the PDB entry*, they almost certainly mean `--index
    residue-id`.
-   **The mask values are per-track and non-negotiable:**

    | Track | Masked value |
    |---|---|
    | `sequence` | `'_'` |
    | `coordinates` | `NaN` (serialised as JSON `null`) |
    | `secondary_structure` | `'_'` |
    | `sasa` | `None` / JSON `null` |

-   **SS8 here is an SS3 APPROXIMATION, not DSSP.** It comes from biotite's
    `annotate_sse` (the P-SEA algorithm), which only distinguishes helix / strand
    / coil. Those are written into the SS8 alphabet as `H` / `E` / `C`. The
    states `G`, `I`, `T`, `B` and `S` are **never** produced. Say so when you
    report secondary structure. If the user needs true SS8, they need DSSP, and
    this skill does not provide it.
-   **SASA is absolute, in Å².** That is what ESM3 was trained on, and it is what
    goes in the track. The buried/exposed *classification* in `inspect` is a
    separate relative-accessibility calculation (see below) and is never written
    into the track.
-   **No API key is needed** except for `--predict-function`. If a key is absent
    or invalid, `extract` still succeeds and simply leaves the function track
    empty.
-   If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All commands MUST be run with `uv run --no-project`.

### 1. `extract` — pull every track out of a structure

Fetches a structure by PDB ID (cached on disk, so re-runs are offline) or reads a
local `.pdb` / `.cif`. Writes a `.json` (sequence, SS8, SASA, residue IDs,
per-track summary) plus a companion `.npz` (coordinates, SASA) alongside it.
**Keep the two files together** — `inspect` and `build-prompt` need both.

```bash
# By PDB ID, from RCSB
uv run --no-project scripts/tracks.py extract \
  --pdb-id 1UBQ --chain A --output ./out/1ubq_tracks.json

# From a local file
uv run --no-project scripts/tracks.py extract \
  --pdb ./data/mystructure.cif --chain B --output ./out/mine_tracks.json

# Also fill the function track (this one call needs BIOHUB_API_KEY)
uv run --no-project scripts/tracks.py extract \
  --pdb-id 1UBQ --chain A --predict-function --output ./out/1ubq_tracks.json
```

If `--chain` is omitted the **first chain in the file** is used. For a multi-chain
entry, always name the chain explicitly.

### 2. `inspect` — summarise and plot the tracks

Reports length, SS8 composition, the SASA distribution (buried / intermediate /
exposed), and residues with missing coordinates. Writes a summary JSON and a
two-panel PNG: the SS8 ribbon over the SASA profile.

```bash
uv run --no-project scripts/tracks.py inspect \
  --tracks ./out/1ubq_tracks.json \
  --output ./out/1ubq_inspect.json \
  --plot ./out/1ubq_tracks.png
```

### 3. `build-prompt` — mask the tracks into an ESM3 prompt

`--keep` names what SURVIVES; `--mask` names what is DESTROYED. Pass exactly one.
Both are **1-indexed, inclusive**, comma-separated (`10-20,45-50`), and a bare
number selects one residue.

```bash
# Scaffolding: keep a motif, mask everything else on every track
uv run --no-project scripts/tracks.py build-prompt \
  --tracks ./out/1ubq_tracks.json --keep 10-20,45-50 \
  --output ./out/1ubq_prompt.json

# Redesign a loop: mask just that loop, keep the rest of the protein
uv run --no-project scripts/tracks.py build-prompt \
  --tracks ./out/1ubq_tracks.json --mask 30-40 \
  --output ./out/loop_prompt.json

# Inverse folding: keep ALL coordinates, condition on nothing else
uv run --no-project scripts/tracks.py build-prompt \
  --tracks ./out/1ubq_tracks.json --keep 1-76 --condition-on coordinates \
  --output ./out/if_prompt.json

# Address the PDB file's own residue numbering instead of track positions
uv run --no-project scripts/tracks.py build-prompt \
  --tracks ./out/1cm4_tracks.json --keep 10-20 --index residue-id \
  --output ./out/1cm4_prompt.json
```

`--condition-on` selects which tracks the region applies to (default: all four).
**Any track you do not name is emitted FULLY masked** — same length, zero
information — which is how you tell ESM3 to condition on nothing there.

The prompt JSON's `sequence`, `secondary_structure`, `sasa` and `coordinates`
keys map 1:1 onto the `BiohubClient.generate(...)` keyword arguments, so a
sibling skill can hand them straight to the API.

## Interpreting the Output

**Report the numbers the script prints. Do not re-derive them.**

-   **Masked counts.** `build-prompt` prints `masked per track: sequence 65/76,
    coordinates 65/76, ...`. State these. A track showing `0/L` masked is fully
    conditioned; `L/L` means it conditions on nothing.
-   **SS8 composition.** Reported as % helix / strand / coil. Always add the
    caveat that this is a P-SEA SS3 approximation, not DSSP. A globin should come
    out ~65-70% helix and ~0% strand; a β-barrel the opposite. If a well-folded
    protein reports 100% coil, the SSE annotation failed — say so rather than
    reporting it as fact.
-   **SASA.** Absolute Å² per residue. The buried/exposed call uses **relative**
    accessibility (RSA = SASA ÷ the residue's theoretical maximum, Tien et al.
    2013): **buried** RSA < 0.25, **exposed** RSA > 0.50, intermediate between.
    A real folded domain has a substantial buried set — that is its hydrophobic
    core. If nothing is buried, the structure is probably a fragment, a single
    extended chain, or a peptide.
-   **Missing coordinates.** `residues_with_no_atoms` lists residues present in
    the sequence but unresolved in the structure. These are *already* masked on
    the coordinate and SASA tracks (NaN / null) and cannot be conditioned on.
    Mention them; they are usually disordered loops or termini.
-   **Function annotations** (only with `--predict-function`) are
    `[label, start, end]` with start/end **1-indexed and inclusive**, the same
    convention as the ranges.

## Common Mistakes

1.  **Off-by-one on the ranges.** `--keep 10-20` is **11** residues, not 10, and
    not residues 11-21. The conversion is 1-indexed inclusive → `sequence[9:20]`.
    Trust the script's `kept_positions` list; never hand-count.
2.  **Confusing track position with PDB residue number.** They coincide only when
    the file happens to start at residue 1. 1UBQ does; 1CM4 does not (it starts at
    4, so track position 10 is author residue 13). Numbers taken from a paper or
    the PDB entry are author numbers → use `--index residue-id`.
3.  **Reporting the SS8 track as if it were DSSP.** It is helix/strand/coil only.
    Never claim a 3-10 helix (`G`), a π-helix (`I`), a turn (`T`), a β-bridge
    (`B`) or a bend (`S`) from this track — those states are never emitted.
4.  **Losing the `.npz`.** The coordinates live in the companion array file next
    to the JSON. Move or copy them together, or `inspect` and `build-prompt` will
    fail with "companion array file not found".
5.  **Using this skill to do the design.** It stops at the prompt. Hand the prompt
    to `esm3-protein-design` or `esm3-inverse-folding`; do not try to generate
    anything here.

## References

-   [Worked example: ubiquitin tracks (1UBQ)](docs/examples/ubiquitin_tracks/report.md)
    — `extract`, `inspect` and `build-prompt --keep 10-20` end to end, with the
    real track numbers (76 residues, 14.5% H / 34.2% E / 51.3% C, 28 buried, 11
    kept / 65 masked) and the masking arithmetic.
-   [Interpretation guide](docs/interpretation-guide.md) — the masking
    conventions (`'_'` / `NaN` / `None`), the SS8-is-an-SS3-approximation caveat,
    the 1-indexed-inclusive ranges, and how the prompt feeds `esm3-protein-design`
    and `esm3-inverse-folding`.
-   [ESM/Biohub API reference](references/esm-biohub-api.md) — endpoints,
    reachable models, the daily credit model, and the shared gotchas.
