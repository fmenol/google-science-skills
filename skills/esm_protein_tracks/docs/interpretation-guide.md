# Interpretation Guide: Reading Tracks and Prompts

This skill prepares inputs; it never interprets biology. But you still have to
read its output correctly before handing a prompt to a design run. This guide
covers the masking conventions, the SS8 caveat, and how the prompt reaches
`esm3-protein-design` and `esm3-inverse-folding`. Numbers are quoted from the
worked example, [`examples/ubiquitin_tracks/report.md`](examples/ubiquitin_tracks/report.md),
and from the shared [ESM/Biohub API reference](../references/esm-biohub-api.md).

## The five tracks

ESM3 is promptable on five parallel, residue-aligned tracks. `extract` fills
four of them locally; only the fifth needs the API.

| Track | Type | Filled by |
|---|---|---|
| `sequence` | string, length L | PDB parse (local) |
| `coordinates` | `(L, 37, 3)` atom37 | PDB parse (local) |
| `secondary_structure` | SS8 string, length L | biotite P-SEA (local) |
| `sasa` | list, length L, Å² | biotite SASA (local) |
| `function` | `[label, start, end]` | ESM3 via `--predict-function` (billed) |

## Masking conventions (non-negotiable)

A prompt is the tracks with most positions *masked* — blanked so ESM3 fills
them in. Each track has one mask value, and it is fixed:

| Track | Masked value | On the wire |
|---|---|---|
| `sequence` | `'_'` | the underscore character |
| `secondary_structure` | `'_'` | the underscore character |
| `coordinates` | `NaN` | serialised to JSON `null` |
| `sasa` | `None` | JSON `null` |

`build-prompt` reports exactly how many positions it masked on each track
(`masked per track: sequence 65/76, …`). **Report those counts; never recount by
eye.** Two readings:

-   `0/L` masked → the track is fully conditioned (ESM3 must respect all of it).
-   `L/L` masked → the track conditions on nothing. A track you leave out of
    `--condition-on` is emitted this way on purpose, so every prompt keeps the
    same shape and the masking counts stay self-documenting.

## Ranges are 1-indexed and inclusive

`--keep 10-20` keeps **11** residues, the 10th through the 20th — the biology
convention. Internally that is `sequence[9:20]`, not `[10:20]` and not `[10:21]`.
Off-by-one here is the most likely mistake in the whole skill. The script prints
`kept_positions` and `kept_residue_ids`; check them instead of hand-counting.

**Position vs. residue-id.** By default ranges count *track positions*, 1..L
along the extracted sequence. They coincide with the file's own numbering only
when the structure starts at residue 1. Ubiquitin 1UBQ does (both lists are
[10…20]); calmodulin 1CM4 does not — it starts at residue 4, so track position
10 is author residue 13. A number quoted from a paper or a PDB entry is an author
number → pass `--index residue-id` so the range lands where you mean it.

## Reading the `inspect` summary

**SS8 composition.** Reported as % helix / strand / coil. Ubiquitin gives 14.5%
H / 34.2% E / 51.3% C — a mixed alpha/beta fold with a real strand fraction.
Sanity anchors: a globin should come out ~65–70% helix and ~0% strand; a
β-barrel the opposite. If a well-folded protein reports **100% coil**, the SSE
annotation failed — say so rather than reporting it as fact.

**SS8 is an SS3 approximation, not DSSP.** The track comes from biotite's
`annotate_sse`, the P-SEA algorithm, which only distinguishes helix / strand /
coil. Those three states are written into the SS8 alphabet as `H` / `E` / `C`.
The states `G` (3-10 helix), `I` (π-helix), `T` (turn), `B` (β-bridge) and `S`
(bend) are **never** emitted here. Always add that caveat when you report
secondary structure, and never claim one of those five states from this track.
If the user needs true SS8, they need DSSP, which this skill does not provide.
(By contrast, when ESM3 itself *predicts* `secondary_structure`, it returns a
genuine 8-state DSSP string — that is a different skill, not this one.)

**SASA.** The track is *absolute* accessibility in Å² — what ESM3 was trained
on. The buried / exposed *classification* is a separate, relative calculation:
RSA = SASA ÷ the residue's theoretical maximum (Tien et al. 2013), with
**buried RSA < 0.25**, **exposed RSA > 0.50**, intermediate between. Ubiquitin:
28 buried / 31 intermediate / 17 exposed, mean 63.19 Å². A real folded domain
has a substantial buried set — its hydrophobic core. If nothing is buried, the
structure is probably a fragment, a single extended chain, or a peptide.

**Missing coordinates.** `residues_with_no_atoms` lists residues present in the
sequence but unresolved in the structure. They are *already* masked on the
coordinate and SASA tracks (`NaN` / `null`) and cannot be conditioned on. 1UBQ
is fully resolved, so this list is empty; disordered loops and flexible termini
are the usual occupants elsewhere. Mention them when present.

## How the prompt feeds the consumers

The prompt JSON's `sequence`, `secondary_structure`, `sasa` and `coordinates`
keys map 1:1 onto `BiohubClient.generate(...)` keyword arguments, so a sibling
skill hands them straight to ESM3:

-   **`esm3-protein-design`** consumes the full masked prompt and fills the
    masked positions — scaffolding, motif grafting, partial redesign. The kept
    region is the constraint; the masked region is what ESM3 invents.
-   **`esm3-inverse-folding`** consumes a coordinates-only prompt. Build it with
    `--keep 1-76 --condition-on coordinates`: that keeps all backbones and masks
    the sequence entirely (`sequence` masked = L, `coordinates` masked = 0). On
    ubiquitin's own backbone the consumer recovers the native sequence at **100%
    identity** — the reference anchor for a correct coordinates-only prompt. (For
    context, ESMC recovers 75/76 of ubiquitin's residues under plain masking; the
    structure-conditioned path is the stronger one.)

This skill stops at the prompt. Do not fold, design, generate, or interpret
function here — those belong to the consumer skills above and to
`esmfold2-structure-prediction` / `esm3-function-prediction`.
