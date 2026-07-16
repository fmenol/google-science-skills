# Interpretation Guide — SS8 & SASA from sequence

Read this before interpreting a prediction. This skill returns two ESM3 tracks
for a sequence: an 8-state secondary-structure string and a per-residue SASA
array. Both are **predictions from sequence** — ESM3 never saw a structure. The
worked example ([`docs/examples/fold_topology/report.md`](examples/fold_topology/report.md))
shows the numbers below on two real proteins; skim it first.

--------------------------------------------------------------------------------

## SS8 → topology

ESM3 emits the full DSSP 8-state alphabet. Group the states, then read the runs.

| State | Meaning | 3-state | Reads as |
|---|---|---|---|
| `H` | α-helix | helix | a run of `H` is an α-helix |
| `G` | 3₁₀-helix | helix | short helix, often capping an α-helix |
| `I` | π-helix | helix | rare; wide helix, usually a few residues |
| `E` | β-strand | strand | a run of `E` is a β-strand |
| `B` | β-bridge | strand | a single-residue strand pairing |
| `T` | turn | coil | hydrogen-bonded turn |
| `S` | bend | coil | high-curvature loop |
| `C` | coil / loop | coil | no regular structure |

**Read runs, not single residues.** A block of `H` is one helix; a block of `E`
is one strand. The script's `secondary_structure_ss3` collapses the eight states
to helix / strand / coil (`GHI→H`, `EB→E`, `TSC→C`) — use it for the coarse fold
class, and go back to SS8 when helix *type* or short features matter.

**Fold-class heuristics** (from the SS3 composition in `summary`):

-   **All-α**: high `%helix`, `%strand ≈ 0`. Globins, four-helix bundles. In the
    worked example myoglobin is 76.5% helix / 0% strand across eight helices.
-   **All-β**: high `%strand`, `%helix ≈ 0`. Immunoglobulin folds, β-barrels.
-   **Mixed α/β** (α+β or α/β): both present. Ubiquitin's β-grasp is 23.7% helix
    / 34.2% strand, with one α-helix at residues 23–34 gripped by a β-sheet.
-   **Mostly coil**: little `H` or `E`. Genuine for a peptide or a disordered
    region — but also what a low-confidence prediction looks like (see below).

**Locate a feature.** To answer "is residue *i* helical / in a strand?", index
the SS8 string at *i−1* (it is 0-indexed, aligned 1:1 with the sequence). To list
a helix's extent, find the contiguous `H`/`G`/`I` run around it.

--------------------------------------------------------------------------------

## SASA → burial

SASA is per-residue solvent-accessible area in Å². The script bins each residue:

| Band | Absolute SASA | Interpretation |
|---|---|---|
| **buried** | ≤ 20 Å² | core residue, packed away from solvent |
| **intermediate** | 20–50 Å² | partially exposed |
| **exposed** | ≥ 50 Å² | surface residue, solvent-facing |

`summary` reports `sasa_mean/min/max` and the buried/exposed counts; the `plot`
draws the dotted 20 and 50 Å² guides over the profile.

**Absolute Å² varies with residue size.** A large residue (Trp, Arg) exposes
more area than a small one (Gly, Ala) at the *same* relative burial, so the
band labels are a convenience, not ground truth — **report the numeric Å²**. If
you need true relative accessibility, measure it from a structure (below).

**What a real buried core looks like.** A folded globular protein has some
near-zero residues *and* some highly exposed ones — a wide range. In the worked
example both proteins span 0.4 → 227.1 Å², and the chemistry splits as expected:

-   Hydrophobic residues (A I L M F W V C) average **~36 Å²** (ubiquitin 35.7,
    myoglobin 26.4) — pushed into the core.
-   Charged residues (D E K R) average **~100 Å²** (ubiquitin 100.8, myoglobin
    98.9) — held at the surface.

If instead the whole profile is uniformly high with no buried residues, suspect a
peptide, an extended / disordered region, or a low-confidence prediction — not a
folded core.

--------------------------------------------------------------------------------

## This is a PREDICTION, not a measurement

ESM3 infers SS8 and SASA from the amino-acid sequence. It did not run DSSP and it
did not compute areas from coordinates. Consequences:

-   **Prefer measuring when you have a structure.** If you already hold a
    PDB/mmCIF — or you fold the sequence with `esmfold2-structure-prediction` —
    compute SS/SASA from the coordinates with **`esm-protein-tracks`** instead.
    Geometry beats a sequence-only guess whenever geometry is available.
-   **Say "ESM3 predicts…".** Never present these labels as observed DSSP or
    crystallographic accessibility.
-   **Sanity-check against known biology.** These are the bundle's canonical
    anchors (see the shared API reference §7): ubiquitin is a β-grasp with one
    helix at ~23–34; myoglobin is an all-α globin. A prediction that flips those
    is wrong regardless of how confident it looks. Lower `--temperature` (e.g.
    0.1, as in the example) for the most conservative call.

--------------------------------------------------------------------------------

## How it complements `esm-protein-tracks`

The two skills are mirror images — pick by what you have:

| | `esm3_secondary_structure_sasa` (this) | `esm-protein-tracks` |
|---|---|---|
| Input | sequence only | a structure (PDB/mmCIF) |
| SS/SASA are | **predicted** by ESM3 | **measured** from coordinates |
| SS resolution | genuine 8-state DSSP | 8-state via DSSP, else 3-state (H/E/C) biotite fallback |
| Use when | you have no structure | you have (or folded) a structure |

Because this skill predicts genuine SS8 even with no structure, it is the right
tool when a user needs true 8-state labels for a sequence they cannot fold. When
a structure exists, `esm-protein-tracks` measures the ground truth and can also
build masked SS/SASA prompts to condition an `esm3-protein-design` run on a
target topology — hand a predicted SS8/SASA from here into that `build-prompt`
step to design toward a desired fold.

--------------------------------------------------------------------------------

## Quick reading checklist

1.  State up front that these are **ESM3 predictions from sequence**.
2.  Read `summary`: `%helix / %strand / %coil` → fold class; `sasa_mean`, buried
    and exposed counts → is there a real core?
3.  Walk the SS8 runs to name and place the helices and strands (1-indexed).
4.  Quote numeric SASA (Å²), not just the buried/exposed labels; call out the
    most buried and most exposed residues.
5.  If a structure is or could be available, note that `esm-protein-tracks` would
    *measure* these instead — and prefer it.
