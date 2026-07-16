# Interpretation Guide

This guide covers how to read an ESMFold2 prediction: the confidence metrics
(pLDDT, pTM, iPTM), the pLDDT-vs-residue plot, the PAE heatmap, when a fold is
unreliable, and the pre-report reasoning checklist. **Read this before writing
any structure-prediction report.** Do not eyeball the raw arrays — run
`analyze`, and take every threshold, band, and verdict from its report.

--------------------------------------------------------------------------------

## 0. pLDDT FROM THIS API IS ON A 0-1 SCALE, NOT 0-100

> [!CRITICAL] This is the single most common misread of ESMFold2 output, and it
> inverts the conclusion. A mean pLDDT of **0.82 is a good, confident structure**
> — the equivalent of 82 in AlphaFold's 0-100 convention — **not** a
> catastrophic one.

-   The metrics JSON carries `"plddt_scale": "0-1"`. Read it. Every pLDDT value
    in the JSON (`plddt_mean`, the per-residue `plddt` array, per-chain stats) is
    between 0 and 1.
-   The **one** place the value is on the 0-100 scale is the **B-factor column of
    the PDB file**, which is rescaled by convention so you can colour by
    confidence in PyMOL/ChimeraX. `analyze` also prints the mean "on the 0-100
    scale" as a convenience, and the report JSON echoes it as
    `plddt.mean_on_0_100_scale` — but the authoritative arrays stay 0-1.
-   If you ever write "pLDDT 0.82, very low confidence", you have made the worst
    error this skill can produce. When in doubt, state the scale explicitly in
    the report: "0.82 on the 0-1 scale, i.e. 82/100."

The bands below are stated on the 0-1 scale throughout. The conventional 0-100
cut-offs (90 / 70 / 50) are simply these ×100.

--------------------------------------------------------------------------------

## 1. Confidence Metrics Reference

Use these tables to turn the numbers `analyze` reports into a verdict. **Take the
numbers from the report; do not recompute them.** The bands and the verdict
strings are exactly the ones `fold.py` emits, so your prose and the machine
report never disagree.

### pLDDT — per-residue local confidence (0-1 scale)

Confidence that each residue is placed correctly relative to its **local**
neighbourhood. Reported per residue and as a mean/median. Colours are the
AlphaFold DB convention used in the pLDDT plot.

| Band                     | pLDDT (0-1) | 0-100 | Plot colour | Reading                                        |
| :----------------------- | :---------- | :---- | :---------- | :--------------------------------------------- |
| **Very high**            | > 0.90      | > 90  | dark blue   | Backbone and often side chains reliable.       |
| **Confident**            | 0.70 - 0.90 | 70-90 | cyan        | Backbone reliable; good for most analysis.     |
| **Low**                  | 0.50 - 0.70 | 50-70 | yellow      | Treat with caution; local geometry uncertain.  |
| **Disordered / very low**| < 0.50      | < 50  | orange      | Likely disordered, flexible, or unfoldable.    |

`analyze` reports the **fraction of residues in each band** and flags any run of
≥ 3 consecutive residues below **0.70** as a low-confidence region, numbered as a
1-based residue range within its chain. A protein is called *well-folded,
confidently predicted* when mean pLDDT ≥ 0.70 **and** ≥ 70 % of residues are in
the confident-or-better bands.

### pTM — global fold confidence (single number, 0-1)

Predicted TM-score: confidence in the **overall fold / global topology**, not
just local geometry. One number for the whole prediction.

| pTM         | Verdict          | Reading                                             |
| :---------- | :--------------- | :-------------------------------------------------- |
| **> 0.8**   | confident fold   | The model is confident in the global topology.      |
| **0.5 - 0.8** | plausible fold | Topology is plausible but not locked in.            |
| **< 0.5**   | unreliable       | No confident global fold. Do not trust the topology.|

### iPTM — interface confidence, complexes only (0-1)

Interface predicted TM-score: for a **complex**, confidence in how the chains sit
relative to each other. This is the number that answers "are these predicted to
form a complex?". `null`/`n/a` for a single chain.

| iPTM        | Verdict             | Reading                                          |
| :---------- | :------------------ | :----------------------------------------------- |
| **> 0.8**   | confident interface | The model is confident these chains form *this* interface. |
| **0.5 - 0.8** | possible interface| An interface is possible but uncertain.          |
| **< 0.5**   | likely no binding   | No confident interface. Probably do not associate.|

> [!WARNING] **iPTM is confidence, not affinity.** A high iPTM means "I am
> confident these chains form this interface", **not** "this binds tightly", and
> **not** "this happens in a cell". It will happily report a confident interface
> for two proteins that never meet physiologically. For binding strength you need
> experiment or a dedicated affinity predictor — see §7.

### PAE — Predicted Aligned Error (per residue-pair, Ångström)

Expected position error of residue *j* when the structure is aligned on residue
*i*. An `L × L` matrix (tokens, so a ligand contributes one entry per atom). Low
(dark in the heatmap) = the two residues are confidently placed **relative to
each other**. This is the only metric that reads *relative* placement, which is
what tells domains and rigid bodies apart (see §3). `analyze` reports the mean
PAE and, for a complex, the mean/min PAE of each inter-chain block.

--------------------------------------------------------------------------------

## 2. Reference Measurements (calibration anchors)

Measured on this API at `--num-loops 10 --num-sampling-steps 50`. Use them as
sanity anchors; a well-behaved small protein should land near the ubiquitin row,
and a genuine tight complex near the barnase+barstar row.

| System                             | mean pLDDT | pTM  | iPTM | Reading                          |
| :--------------------------------- | :--------- | :--- | :--- | :------------------------------- |
| Ubiquitin (76 aa, real)            | **0.82**   | **0.78** | —    | a real, well-folded domain       |
| Ubiquitin, **sequence shuffled**   | **0.48**   | **0.24** | —    | same composition, no fold        |
| Barnase + barstar (real complex)   | **0.95**   | 0.97 | **0.96** | a genuine, tight interface       |

The scramble row is the load-bearing calibration: a sequence with **identical
amino-acid composition** to ubiquitin, shuffled, collapses from pLDDT 0.82 to
0.48 and pTM 0.78 to 0.24. Both numbers fall, together, and both fall below the
"confident" and "plausible" cut-offs. **If a designed or unfamiliar sequence
scores like that scramble, the model is telling you it has no confident fold.**
Do not proceed to docking, Foldseek, or MD with it. The two worked examples in
`docs/examples/` are exactly this real-vs-scramble pair.

--------------------------------------------------------------------------------

## 3. Reading the pLDDT-vs-Residue Plot

`analyze` writes pLDDT vs residue with the four confidence bands shaded behind
the trace, a dashed line at the 0.70 low-confidence cut-off, and black vertical
lines at chain boundaries (for a complex). Read it like this:

-   **Where does the trace sit?** A trace riding in the cyan/blue bands
    (≥ 0.70) across the whole length is a confidently folded chain. A trace that
    lives in the yellow/orange bands (< 0.70) throughout has no confident fold.
-   **Dips below the dashed line are the interesting part.** A localised dip
    marks a flexible loop, a mobile terminus, or a disordered linker. `analyze`
    lists these as residue ranges — quote them, and quote *which* structural
    element they correspond to when you can identify it.
-   **Termini often dip and that is usually real.** Chain ends are genuinely more
    mobile. A short low-confidence tail at the N- or C-terminus of an otherwise
    confident domain is a feature of the molecule, not a failure of the
    prediction. (In the ubiquitin example, residues **72-76** — the C-terminal
    `RLRGG`, whose `LRGG` is the conjugation motif that is mobile in solution —
    dip to mean pLDDT **0.56** while the globular body 1-71 stays ≥ 0.70. The
    plot is *correctly* reporting a known flexible tail.)
-   **A complex:** compare the per-chain levels. One confident chain next to a
    low-confidence partner means the model folded one but not the other — treat
    them separately.

## 4. Reading the PAE Heatmap for Domains and Flexibility

The PAE heatmap is how you find **domains, rigid bodies, and inter-chain
placement** — things pLDDT (which is purely local) cannot see. Dark = low error =
confidently placed relative to each other.

-   **One solid dark block = one rigid unit.** If the whole matrix is dark, every
    residue is confidently placed relative to every other residue: a single rigid
    domain. (Ubiquitin's heatmap is one dark block — its globular fold is a single
    rigid domain — with a pale fringe along the last few rows/columns for the
    mobile C-terminal tail, whose position relative to the body is uncertain.)
-   **Two (or more) dark blocks on the diagonal with a bright off-diagonal =
    multiple domains.** Each block is internally rigid, but the bright
    off-diagonal says the model is unsure how the domains are oriented *relative
    to each other* — i.e. a flexible inter-domain hinge. Report those as separate
    rigid domains and warn that their relative orientation is not determined.
-   **Bright rows/columns = flexibility.** A residue that is bright against
    everything else is not confidently placed relative to the rest of the
    structure — a flexible loop or tail. This matches the pLDDT dip for the same
    residues.
-   **For a complex, read the off-diagonal chain blocks.** The block relating
    chain A to chain B is the interface confidence. Dark off-diagonal blocks (low
    inter-chain PAE) mean the model is sure how the chains dock. `analyze` prints
    the mean/min PAE per inter-chain block; quote it alongside iPTM. (Barnase +
    barstar: inter-chain A-B mean PAE **2.4 Å**, min **0.3 Å** — a tightly, dark
    off-diagonal block consistent with iPTM 0.96.)

> [!NOTE] The PAE axis is in **tokens**, and a ligand contributes one token per
> atom, so a protein+ligand PAE is larger than the residue count. `analyze`
> labels the chain boundaries for you; do not try to line the axis up with the
> residue index by hand.

--------------------------------------------------------------------------------

## 5. When a Fold Is Unreliable

Stop and treat the prediction as unreliable — and say so in the report — when any
of these hold:

-   **pTM < 0.5** (`unreliable` verdict). The global topology is not to be
    trusted, even if a few local stretches have decent pLDDT.
-   **Mean pLDDT < 0.5**, or the disordered band dominates. `analyze`'s overall
    verdict will read *largely disordered or unfoldable — do not use for
    downstream structural analysis*. Take it at its word.
-   **For a complex, iPTM < 0.5** (`likely no binding`). The chains are not
    predicted to form the interface; do not describe them as "a complex".
-   **The whole chain is one low-confidence region.** When `analyze` reports a
    single low-confidence range covering the entire chain (as it does for the
    scramble control, 1-76), there is no confident core to salvage.

In every one of these cases the correct output is a clearly-stated negative
result, not a hedged structural story. See §6.

--------------------------------------------------------------------------------

## 6. Negative Results & Scientific Integrity

> [!CRITICAL] **A low-confidence fold is a real answer, not a failure.** "ESMFold2
> has no confident structure for this sequence" is a legitimate, useful, and
> often correct scientific finding. Report it plainly.

-   **Do not stretch a low-pLDDT blob into a structural story.** If pTM is 0.24
    and the whole chain is disordered, do not describe secondary-structure
    elements, a "core", or a binding pocket. There is no reliable structure to
    describe.
-   **The value of the negative control.** A scrambled sequence with the same
    composition collapsing to pLDDT 0.48 / pTM 0.24 is *how you know the pipeline
    is reading the sequence at all*. Reproduce it when a designed or unfamiliar
    sequence matters: if your candidate scores like the scramble, the model is
    telling you it has no fold. This is the whole point of the `scramble_control`
    worked example.
-   **Confidence is not validity.** A **high** pLDDT does **not** mean the
    sequence is a real, expressible, or functional protein. ESMFold2 will
    confidently fold sequences that never express, and will confidently fold a
    single mutant that abolishes function. A confident prediction of a mutant is
    not proof of a functional protein — sanity-check designs against a negative
    control, and against experiment where it matters.
-   **iPTM is not evidence of a biological complex.** A confident interface (high
    iPTM) says the model can place the chains together, not that they meet, bind,
    or function together in a cell. Do not upgrade "confident interface" to "these
    proteins interact" without orthogonal evidence.
-   **Do not average pLDDT by hand or invent thresholds.** Use `analyze`. It does
    the banding, region-finding, and pTM/iPTM verdicts numerically, so two reports
    on the same numbers read the same way.

--------------------------------------------------------------------------------

## 7. What ESMFold2 Does Not Tell You (Model Limitations)

A confident structure answers "what does this sequence fold to?" and nothing
else. It does **not** report:

-   **Binding affinity.** iPTM is interface *confidence*, not a Kd. A confident
    interface can be a nanomolar binder or a crystallographic artefact; the number
    is the same. Use experiment or a dedicated affinity predictor.
-   **Whether the protein is real, expressible, or stable in a cell.** See §6:
    confidence ≠ validity.
-   **Function, catalysis, or mechanism.** A confident fold of an enzyme does not
    tell you it is catalytically active, nor which residues are catalytic.
-   **Conformational ensembles / dynamics.** You get one static model. Flexible
    regions show up as low pLDDT / high PAE, but the model does not enumerate
    alternative conformations, allosteric states, or folding pathways.
-   **The effect of ligands, ions, PTMs, or partners you did not include.** The
    prediction is of exactly the entities you gave it. A pocket may look
    different with its cofactor bound.
-   **Anything an MSA would add, unless you asked for it.** The default
    `esmfold2-fast-2026-05` is single-sequence and **silently ignores** an MSA. If
    evolutionary information matters, use `esmfold2-2026-05` with `--msa`
    (`fold.py` enforces this switch for you).

When the user's real question is any of the above, say so, and point them at the
right tool: existing solved/predicted structures →
`alphafold-database-fetch-and-analyze` (UniProt) or `pdb_database` (PDB ID);
structural homologs → `foldseek-structural-search`; ranking many candidate
binders → `esmfold2-binder-screening`.

--------------------------------------------------------------------------------

## 8. Pre-Report Reasoning Checklist

Complete this BEFORE writing the report. Ground every statement in a number from
the `analyze` report or a feature you can point to in a plot.

### 1. Scale and provenance

-   [ ] **Confirm the pLDDT scale.** State it explicitly ("0.82 on the 0-1 scale,
    i.e. 82/100"). Never report a 0-1 value as if it were 0-100.
-   [ ] **Record the model and settings** actually used (`esmfold2-fast-2026-05`
    vs `esmfold2-2026-05`, `--num-loops`, `--num-sampling-steps`, whether an MSA
    was used). A fast-model single-sequence prediction is a triage result.

### 2. Global confidence verdict

-   [ ] **Quote mean pLDDT and pTM** with their verdict strings from the report.
-   [ ] **State the overall verdict** (`well-folded` / `partially ordered` /
    `largely disordered`) and let it, not vibes, decide whether a structural
    story is warranted at all (§5).
-   [ ] **For a complex, quote iPTM and its verdict**, and the inter-chain PAE.
    That is the number that answers "do these bind?" — with the affinity caveat
    (§7) stated.

### 3. Local confidence and topology

-   [ ] **List every low-confidence region as a residue range**, and identify the
    structural element where you can (mobile terminus, loop, linker, disordered
    segment).
-   [ ] **Read the PAE heatmap for domain structure**: one rigid body or several?
    Any flexible hinges (bright off-diagonal blocks)? Report domains as separate
    rigid ranges when the model shows them.
-   [ ] **Reconcile the two plots**: pLDDT dips and PAE-bright regions should
    agree. Note it if they do not.

### 4. Integrity checks

-   [ ] **Is this a negative result?** If pTM < 0.5 or the chain is largely
    disordered, report the negative clearly — do not manufacture structure (§6).
-   [ ] **Compare to a control where it matters.** For a designed or unfamiliar
    sequence, is it clearly better than its own scramble on **both** pLDDT and
    pTM? If not, say the model has no confident fold.
-   [ ] **Confidence ≠ validity/affinity/function.** Do not claim the protein is
    real/functional, or that a complex binds, on the strength of the confidence
    numbers alone.

### 5. Report readiness

-   [ ] **Embed both plots** (pLDDT-vs-residue and, if `--include-pae` was set,
    the PAE heatmap), each with a caption that states what it shows.
-   [ ] **Verify every figure referenced in the report exists** in the report
    folder.
-   [ ] **Note downstream suitability**: which residue ranges are safe for
    docking / Foldseek / MD, and which (low-confidence, disordered) must be
    excluded.
-   [ ] **Remind the user** that per-residue pLDDT is in the PDB B-factor column
    (0-100), so they can colour by confidence in PyMOL/ChimeraX.
-   [ ] **State that this skill (ESMFold2 on the Biohub Platform) was used.**

--------------------------------------------------------------------------------
