# Interpretation Guide

How to read a binder screen: the metric bands, the PAE heatmap signatures,
the traps that make a non-binder look like a hit, and the pre-report reasoning
checklist. **Read this before you interpret a screen or write a report.**

This skill answers exactly one question:

> **Of the candidates I already have, which ones does ESMFold2 predict actually
> bind the target — and which should I order?**

It is the *selection* stage. It ranks candidates you already hold; it does not
design them, and it does not measure affinity. Everything below serves that one
question.

--------------------------------------------------------------------------------

## The one caveat that governs everything

> [!WARNING]
> **iPTM is a structural-confidence proxy, NOT an affinity prediction.** It
> tells you how confident ESMFold2 is that it has placed the two chains
> correctly relative to each other. It does **not** give you a Kd, a ΔG, an
> IC50, or a calibrated rank-order of affinities. A high iPTM is a **hypothesis
> to test experimentally**, not a measurement.
>
> Never report an iPTM as if it were a binding constant. Never tell a user that
> iPTM 0.90 "binds more tightly" than iPTM 0.85 — the model does not resolve
> that difference, and neither number is a Kd. The correct verb is **"predicts
> a confident interface"**, never "binds with affinity X".

Everything this skill emits is a *confidence*, on a 0–1 or Ångström scale. Keep
that framing in every sentence of every report.

--------------------------------------------------------------------------------

## iPTM — the primary metric

**Interface predicted TM-score.** The single number to rank on. It is the only
metric that directly scores the *interface* (as opposed to how well each chain
folds on its own), and it is by far the most stable across sampling settings.

| iPTM        | Verdict                | What to do                             |
| :---------- | :--------------------- | :------------------------------------- |
| **> 0.8**   | Confident interface    | A real hypothesis worth ordering.      |
| **0.5–0.8** | Possible interface     | Ambiguous. Do not over-claim. Deeper look, or deprioritise. |
| **< 0.5**   | Likely no binding      | The model places no confident interface. |

These bands are the same ones the script's `verdict` column applies. They are
calibrated against controls, not invented — see the reference panel below.

> [!NOTE]
> **Anchor value.** Barnase + barstar — one of the tightest protein–protein
> complexes known (Kd ~ 10⁻¹⁴ M) — folds to **iPTM 0.96** with **52 confident
> interface contacts**. That is what an unambiguous, real, high-affinity
> interface looks like through this skill. It is the calibration point in the
> shared API reference §7 and the top of every worked example here.

--------------------------------------------------------------------------------

## Interface PAE — the confirming metric

**Mean cross-chain predicted aligned error (Å). Lower is better.** PAE(i,
j) is the expected error in the position of residue *j* when the structure is
aligned on residue *i*; the interface PAE averages the two off-diagonal
(cross-chain) blocks. It answers: *does the model know where chain B sits
relative to chain A?*

| Interface PAE | Reading                                                        |
| :------------ | :------------------------------------------------------------- |
| **< 5 Å**     | Confidently placed interface. Consistent with a real complex. |
| **5–15 Å**    | Uncertain. The chains are not locked relative to each other.  |
| **> 20 Å**    | The model has no idea where the two chains sit. No interface.  |

iPTM and interface PAE should **agree**. When they do, trust the call. Barnase +
barstar: iPTM 0.96, interface PAE **2.7 Å** — both say the same thing. A true
non-binder sits at iPTM < 0.3 **and** interface PAE > 20 Å (see the decoys). If
they ever disagree, do not force a verdict; report the ambiguity.

--------------------------------------------------------------------------------

## Interface contacts — read the *gated* count, never the raw one

The script reports two contact counts. They are not interchangeable, and
confusing them is the single most expensive mistake in binder screening.

| Field                            | Definition                                          | Use it? |
| :------------------------------- | :-------------------------------------------------- | :------ |
| `confident_interface_contacts`   | Cross-chain residue pairs with CB–CB < 8 Å **and** PAE < 15 Å | **Yes** |
| `interface_contacts`             | The same, **without** the PAE gate (raw geometry)   | No — transparency only |

> [!CRITICAL]
> **Raw contact count does not discriminate binders from non-binders — it is a
> trap.** ESMFold2 must emit *some* coordinates for two chains, so it packs them
> against each other whether or not they bind. A large, unrelated protein racks
> up incidental proximity. In the controls here, lysozyme shows **38** raw
> geometric contacts and the scrambled decoy **39** — both a large fraction of
> the true binder's 52 — yet once gated on PAE they collapse to **0** and **1**.
> Rank on raw contacts and you will order the non-binder over the nanomolar
> binder.

Always rank on `iptm`. Use `confident_interface_contacts` as the corroborating
count — the true binder's 52 gated contacts survive the PAE gate intact (all 52
of its geometric contacts are confident), while the decoys keep essentially
none.

--------------------------------------------------------------------------------

## Binder pLDDT — folding is not binding

`binder_plddt` (0–1) is the mean pLDDT of the **binder chain alone**: does the
candidate fold into *a* structure? Above ~0.7 it does.

> [!WARNING]
> **A well-folded binder is not a binding binder.** In the controls, lysozyme
> folds beautifully — `binder_plddt` **0.897**, better than many real designs —
> and still does not bind barnase (iPTM 0.230, 0 confident contacts). pLDDT says
> the chain is confident *about its own shape*; it says nothing about the
> interface. Never conclude "it binds" from a high binder pLDDT.

pLDDT earns its place as a guard on the *other* side: a candidate that cannot
even fold (the scrambled decoy, `binder_plddt` 0.321) is not orderable
regardless of any interface signal. Use it to reject, not to select.

--------------------------------------------------------------------------------

## selection_score — the tie-breaking composite

A single 0–1 summary, higher is better:

```
selection_score = 0.5·iPTM
                + 0.3·binder_pLDDT
                + 0.2·(1 − min(interface_PAE, 31.75)/31.75)
```

iPTM carries half the weight because it is the only interface term; binder pLDDT
guards against a candidate that cannot fold; the inverted interface PAE breaks
ties. **Contact count is deliberately excluded** — it scales with binder size
and does not separate binders from non-binders.

Use `selection_score` as a convenience summary, but **rank on `iptm`**. The
composite can be inflated by the folding term alone: lysozyme, a confident fold
that does not bind, still scores **0.441** on its pLDDT contribution. The
composite is a summary, not the verdict.

--------------------------------------------------------------------------------

## Reference measurements — the calibrated control panel

Every number below is a real value produced by this skill's `screen` /
`analyze-interface` commands on barnase (target, chain A) at the settings the
eval uses (`--num-loops 3 --num-sampling-steps 50`). Treat these as the sanity
anchors your own screens should straddle.

| Candidate            | iPTM      | interface PAE | binder pLDDT | confident contacts | geometric contacts | verdict             |
| :------------------- | :-------- | :------------ | :----------- | :----------------- | :----------------- | :------------------ |
| **barstar** (true)   | **0.964** | **2.65 Å**    | 0.937        | **52**             | 52                 | confident interface |
| lysozyme (decoy)     | 0.230     | 22.76 Å       | 0.897        | 0                  | 38                 | likely no binding   |
| scrambled barstar    | 0.144     | 23.02 Å       | 0.321        | 1                  | 39                 | likely no binding   |

The signal is not subtle: the true binder sits **~0.73 iPTM above** the better
decoy, and the two metrics that matter (iPTM, interface PAE) call all three
correctly. The two metrics that *don't* discriminate (binder pLDDT, raw
geometric contacts) would each mislead you on their own — lysozyme's pLDDT beats
many real designs, and its 38 raw contacts approach barstar's 52.

> [!NOTE]
> Exact decoy values drift a little with sampling depth (the shared reference
> and SKILL.md quote lysozyme in the 0.13–0.23 range across settings); the true
> binder is rock-stable at **iPTM ~0.96**. What never changes is the *ordering*
> and the *bands*: barstar » both decoys, and both decoys land firmly below 0.5.
> Report the ordering and the band, and quote the exact number you measured.

--------------------------------------------------------------------------------

## Reading the PAE heatmap

`analyze-interface` renders the PAE matrix with the chain boundary marked in
red. The diagonal blocks are each chain against itself; the **off-diagonal
blocks are the cross-chain interface**. There are three signatures, and learning
to tell them apart is most of the skill.

1.  **Binder — the whole matrix is dark.** Both diagonal blocks *and* both
    off-diagonal blocks are dark (low PAE). The model is confident about each
    chain and about their relative placement. This is barnase:barstar.
2.  **Folds-but-does-not-bind — dark diagonals, pale off-diagonals.** Each chain
    folds confidently (dark diagonal blocks) but the cross-chain blocks are
    washed out and pale. The model knows each shape and has no idea how they fit
    together. This is barnase:lysozyme — the canonical non-binder signature.
    **Learn to recognise it.**
3.  **Neither folds nor binds — one diagonal also pale.** The target block stays
    dark, but the candidate's own diagonal block is pale (it cannot fold) *and*
    the cross-chain blocks are pale. This is barnase:scrambled-barstar: a fold
    destroyed by shuffling, failing on both counts.

--------------------------------------------------------------------------------

## Negative Results & Scientific Integrity

> [!CRITICAL] **Most candidate binders do not bind.** A screen whose honest
> answer is "none of these is predicted to bind" is a valuable, publishable
> result. Do not stretch to find a binder where the model shows none.

-   **The negative case is the point.** This skill exists to separate binders
    from non-binders. Reporting "ESMFold2 predicts no confident interface for
    any candidate" — and shortlisting nothing — is a correct, complete answer.
    An empty shortlist is a finding, not a failure.
-   **Do not rescue a non-binder with the wrong metric.** A high `binder_plddt`,
    a large raw `interface_contacts`, a decent `selection_score` from the
    folding term — none of these is binding. If iPTM is below the band and
    interface PAE is high, the candidate does not bind, however nicely it folds.
-   **Strict anti-speculation.** Do not invent an interface, a hotspot, or a
    binding mode the model did not place. If iPTM is 0.5–0.8, say "possible,
    ambiguous", not "binds". Never upgrade a confidence into an affinity.
-   **iPTM is not a Kd — restate it in every report.** A confident interface is
    a hypothesis for the wet lab, not a binding constant. The final decision to
    order still rests on experiment.
-   **Model limitations — what a confident interface does *not* tell you:**
    -   **Affinity / Kd / kinetics.** No quantitative binding strength, on/off
        rates, or ΔG.
    -   **Specificity.** iPTM scores this target:candidate pair. It does not say
        the candidate won't also stick to ten other proteins.
    -   **Expression, solubility, aggregation, stability.** Predicted only
        indirectly (`isoelectric_point` is a coarse solubility proxy; the
        reference protocol filters minibinders to pI < 6). A confident interface
        on a candidate that will not express is not orderable.
    -   **Conformational change, allostery, induced fit, post-translational
        modification, cofactor- or membrane-dependent binding.** ESMFold2 folds
        one static complex from sequence; dynamic or context-dependent binding
        is out of scope.
    -   **The design half of the protocol.** This skill cannot invent binders
        (see *Scope* in SKILL.md); it can only rank the ones you bring.

--------------------------------------------------------------------------------

## Pre-Report Reasoning Checklist

Complete this BEFORE writing the report. Ground every claim in a real number
from the script's output — never eyeball the PDB or re-derive a metric by hand.

### 1. Rank and bands

-   [ ] **Ranked on iPTM?** Confirm the table is sorted by `iptm` (or state and
    justify the alternative). Rank is the deliverable.
-   [ ] **Assigned each candidate a band** (> 0.8 / 0.5–0.8 / < 0.5) and used
    the matching verb — "confident interface" / "possible" / "likely no
    binding" — never an affinity word.
-   [ ] **Checked iPTM and interface PAE agree** for every candidate you call a
    hit. If they disagree, reported the ambiguity instead of forcing a verdict.

### 2. Traps cleared

-   [ ] **Did not rank or argue from raw `interface_contacts`.** Used
    `confident_interface_contacts` only, and noted where the raw count would
    have misled.
-   [ ] **Did not read a high `binder_plddt` as binding.** Confirmed folding and
    binding are reported as separate facts.
-   [ ] **Did not quote any iPTM as a Kd, ΔG, or affinity.**

### 3. Evidence

-   [ ] **Embedded the PAE heatmap** for each candidate discussed, and described
    its signature (whole-matrix dark / dark-diagonals-pale-off / one-diagonal
    pale). Every figure referenced exists in the folder.
-   [ ] **Named the interface residues** from the script's per-residue output
    when claiming a specific binding mode — do not assert a hotspot the output
    does not list.

### 4. The shortlist

-   [ ] **Chose `--top-n` deliberately** — how many would you actually order?
    Justified the cutoff against the iPTM band, not a round number.
-   [ ] **Stated the filters applied** (`--min-iptm`, `--max-isoelectric-point`,
    …) and how many candidates each dropped.
-   [ ] **Reported an empty shortlist honestly** if nothing cleared the band.

### 5. Framing

-   [ ] **Every verdict is a prediction to test**, not a result. The report says
    so explicitly.
-   [ ] **Stated that this skill was used** to produce the screen.
-   [ ] **Directly answered the user's question** — which to order, which to
    drop, and why.
