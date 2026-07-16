# Interpretation Guide

This guide covers how to quality-control an ESM3 design, how to read its fold
metrics, how to judge a scaffolded motif, and the pre-report reasoning checklist.
Read it before writing any design report. Every number quoted here was measured
against this API in the skill's eval; the worked examples in `docs/examples/`
were produced from the same runs.

--------------------------------------------------------------------------------

## The QC Contract (read this first)

> [!CRITICAL] A raw ESM3 sample is a **guess**, not a design. Never present a
> generated sequence as a "design" without folding it and quoting its **pTM** and
> **pLDDT**. A sequence reported without its QC metrics is a scientific error, not
> a shortcut.

ESM3 is a generative model: it will happily emit a plausible-looking amino-acid
string for any prompt, including prompts it has no idea how to satisfy. The only
evidence that the string is a *protein* — something that would fold to a defined
structure — is an independent structure prediction. `design.py` folds every
sample with ESMFold2 and ranks by pTM automatically, so the QC is already done;
your job is to **read the metrics JSON and report it honestly**, not to recompute
it and not to skip it.

Three questions decide whether a design is worth reporting at all:

1.  **Does it fold?** Is pTM above the random floor for its length (see below),
    and is pLDDT in the confident range?
2.  **For a scaffold: did the motif survive?** Is the motif verbatim in the
    sequence *and* is its CA RMSD to the reference low?
3.  **Is it real protein, not junk?** Is it free of low-complexity runs and
    membrane/signal-peptide artifacts?

If the answer to any of these is no, the design is a *rejected* design. Say so.

--------------------------------------------------------------------------------

## Score Interpretation

Read these off the metrics JSON `design.py` writes. Do not recompute them.

| Quantity | Reading |
| :------- | :------ |
| **pLDDT** (0-1) | > 0.9 very high · 0.7-0.9 confident · 0.5-0.7 low · < 0.5 disordered |
| **pTM** | > 0.8 confident global fold · 0.5-0.8 plausible · < 0.5 unreliable |
| **Motif CA RMSD** | < 1.0 Å excellent · 1.0-2.0 Å good · > 3.0 Å the graft failed |

> [!NOTE] **pLDDT is on a 0-1 scale here**, not 0-100. A design at pLDDT 0.71 is
> at the conventional "70" confidence, not "0.71". Multiply by 100 only when
> writing B-factors. The skill already reports the 0-1 value.

### Calibration — quote these when judging a design

These are the reference points the skill is measured against. Sample sizes are
small and are stated because they are small; treat them as anchors, not
constants. All folded with `esmfold2-fast-2026-05`; all designs from
`esm3-open-2024-03`.

| Setting (length 60 unless noted) | pTM | pLDDT | n |
| :------------------------------- | :-- | :---- | :- |
| Ubiquitin, a real 76-mer (the ceiling) | 0.78 | 0.82 | 1 |
| **Uniformly-random 60-mer (the floor)** | **0.20** | **0.43** | 5 |
| `generate`, unconditioned, T=0.5 | 0.26 (range 0.18-0.38) | 0.59 | 5 |
| `chain-of-thought`, 3 tracks | 0.35 | 0.68 | 1 |
| `scaffold-motif` into a 200-mer | 0.98 | 0.96 | 2 |

The floor, the `generate` row, and the `scaffold-motif` row are read directly
from this skill's shipped worked examples (`docs/examples/`), so you can verify
them against the JSON. The ceiling (ubiquitin 0.78 / 0.82) is the shared API
reference's §7 measurement for a real protein of comparable length. The
`chain-of-thought` row is a single stochastic sample and is only indicative — a
separate single draw recorded in `SKILL.md` gave pTM 0.65 where this cassette's
gave 0.35; a 0.30 spread on n=1 is exactly why you must draw several and judge by
the pTM you measure. `SKILL.md`'s "Interpreting the Output" table carries the
broader calibration (a temperature sweep at n=8, scaffolds at n=4); the numbers
here agree with it within the small sampling spread.

**The single most important lesson in this table: conditioning is what makes
ESM3 good.** An unconditioned 60-mer barely clears the random floor (mean pTM
0.26 vs 0.20); the *same* model, given a rigid motif to build around, folds the
whole 200-residue protein to pTM 0.98. If a user needs a confidently-folded
protein, the answer is almost never "sample more unconditioned designs" — it is
"condition the design" (motif / SS8 / SASA) or use `chain-of-thought`, and design
longer. pTM is harshly penalised below ~80 residues.

--------------------------------------------------------------------------------

## The Designed-vs-Random Foldability Logic

This is the central scientific test of a de novo campaign, and the one most
often skipped. A pTM of 0.26 sounds bad in isolation. It is only interpretable
against the **floor** — what a *random* amino-acid string of the same length and
composition folds to. That is why the negative control is not optional:

-   **Random 60-mers fold to pTM ~0.20, pLDDT ~0.43** (measured, n=5). This is
    not zero: ESMFold2 assigns some confidence to any string, so "the design has
    pTM 0.3, therefore it folds" is a false inference.
-   **A de novo design is only doing something if it beats this floor.** In the
    worked example the designs average pTM 0.26 vs 0.20 and pLDDT 0.59 vs 0.43 —
    a real but *modest* margin, exactly what a short unconditioned prompt should
    give.
-   **The distributions overlap.** In the same worked example the *best* random
    control (pTM 0.261) edged above the *mean* design (0.259). At length 60 the
    signal is weak enough that a single lucky design or a single lucky control
    tells you nothing. This is why the skill requires multiple samples and
    compares *means*, and why the eval draws five, not two.

> [!IMPORTANT] When you report a de novo design, report it **relative to the
> floor**, not in absolute terms. "Mean pTM 0.26 vs a random-control floor of
> 0.20 (n=5 each)" is honest. "pTM 0.38, a well-folded design" — quoting only the
> best sample, with no floor — is not.

--------------------------------------------------------------------------------

## Reading Motif RMSD (scaffolding)

`scaffold-motif` grafts a rigid 3D motif (an active site, an epitope) into a
larger designed protein. Two things must both be true for the graft to have
worked, and they are independent:

1.  **The motif sequence is verbatim.** The exact motif residues must appear at
    the requested offset in the design. The JSON records `motif_verbatim: true`
    and the grafted substring; check it. `design.py` pins the motif in the
    sequence prompt, so this is normally true, but verify it — an off-by-one in
    `--motif-range` or `--motif-start` silently grafts the wrong residues.
2.  **The motif geometry is preserved.** After folding the whole design, the
    motif's CA trace must superpose on the *original crystal* motif to low RMSD.
    This is the `motif_rmsd_ca` field, a Kabsch superposition the skill computes
    for you.

> [!CRITICAL] A design that keeps the motif *residues* but loses the motif
> *geometry* has **failed**, even at high pTM. The whole point of scaffolding is
> to hold the functional site in its active conformation. Sequence identity in
> the motif is necessary but not sufficient; the RMSD is what certifies the
> graft. Report both.

In the worked example (1ITU renal-dipeptidase motif, 23 residues, into a 200-mer)
both designs kept the motif verbatim and superposed at **CA RMSD 0.26 Å and
0.28 Å** — sub-Ångström, i.e. the active-site geometry is essentially identical
to the crystal. Anything under ~1.0 Å is excellent; 1-2 Å is usable; above ~3 Å
the graft has failed regardless of the design's overall pTM.

--------------------------------------------------------------------------------

## Temperature Guidance

Sampling temperature trades foldability for diversity, and the trade is not in
your favour for quality. Raising the temperature to get "more creative" designs
does the opposite:

-   **T=0.5 is the default and the best measured setting.** It is the value the
    ESM3 cookbook uses and the one behind every calibration number above.
-   **Higher temperature lowers pTM.** Measured at length 60: T=0.5 gave mean pTM
    0.32 (n=8, a separate calibration run); T=0.7 gave 0.28; a single T=1.0 draw
    gave 0.19 — indistinguishable from the random floor.
-   **Do not exceed ~0.7.** If designs look repetitive or low-complexity, the fix
    is more conditioning or more samples, not more temperature.

--------------------------------------------------------------------------------

## Low-Complexity and Artifact Sanity Check

A high-ranked design can still be junk. ESM3's unconditioned short designs are
frequently **low-complexity** (poly-Ala or poly-Lys runs) or resemble **signal
peptides** and **transmembrane helices**. This is a real, documented property of
the model at short lengths, not a bug in the script — but it means pTM alone does
not certify a *useful* design.

Inspect the top design's sequence before reporting it. Real examples from the
worked de novo run (all five are genuine top-5 designs):

-   `design_1`, the **best by pTM (0.376)**, ends in
    `...KKEEDKKKKEKEKKKKKEKEKKKKKQDKKK` — a long, low-complexity poly-Lys tail.
    The highest pTM in the set is still not a clean sequence.
-   `design_3` (`...PLLPALAP...ALLAAALAAAP...`) is Ala/Leu-repetitive.
-   The `chain-of-thought` design (pTM 0.35) ends in `...KAAAAAAAAAA`, a
    ten-residue poly-Ala run.

> [!CAUTION] If the top design is low-complexity, do not launder it through its
> pTM. Report the artifact, and either draw more samples, add conditioning
> (SS8 / SASA / motif), or design longer. A poly-Ala helix with pTM 0.9 is a
> well-predicted poly-Ala helix, not a useful protein.

--------------------------------------------------------------------------------

## Negative Results & Scientific Integrity

> [!CRITICAL] **A design that folds poorly is a rejected design, and reporting it
> as rejected is a real result.** Do not stretch a low-pTM sample into a
> "candidate" it is not.

-   **Value of negative results.** "ESM3's unconditioned 60-mer designs beat a
    random-string floor only modestly (mean pTM 0.26 vs 0.20) and are frequently
    low-complexity; conditioning is required for a confident fold" is a correct,
    useful finding. It tells the user what to do next. A fabricated "we designed a
    novel mini-protein (pTM 0.38)" does not.
-   **Quote the floor, always.** A de novo pTM is meaningless without the random
    control it is measured against. If you did not run the control, say the pTM
    is uncalibrated.
-   **Best-sample cherry-picking is the classic trap.** Design is stochastic with
    a wide spread; the best of five samples is an upward-biased estimate. Report
    the number of samples drawn and the mean, then the best. Never report only
    the best as if it were typical.
-   **Occupancy ≠ success.** A motif being present in the sequence does not mean
    the graft worked (check the RMSD). A high pTM does not mean the sequence is
    useful (check for low-complexity). A design existing does not mean it is a
    protein (check it against the floor).
-   **These are computational designs, not validated proteins.** Every report
    must say so. A confident pTM/pLDDT and a sub-Ångström motif RMSD are
    *in-silico* predictions of foldability, not evidence of expression, stability,
    solubility, or function. Wet-lab validation is a separate step this skill does
    not perform.

### What ESM3 design does *not* establish

-   **Function.** A scaffolded active-site motif is held in the right *geometry*;
    whether the designed protein is catalytically active is untested. Predicting
    function is `esm3-function-prediction`.
-   **Binding / affinity.** A de novo binder backbone is a backbone, not a
    validated binder. Interface quality (iPTM) needs a complex fold
    (`esmfold2-binder-screening`); the reference iPTM for a true high-affinity
    complex (barnase-barstar) is 0.96, far above anything a raw backbone implies.
-   **Expressibility / stability / solubility.** Not modeled at all.

--------------------------------------------------------------------------------

## Model & API Limitations

-   **Only `esm3-open-2024-03` is reachable.** The cookbook's `esm3-medium-*` /
    `esm3-large-*` return HTTP 403. A 403 is a permissions wall, not a transient
    error; do not retry a larger model.
-   **pTM is length-sensitive.** Short designs (< ~80 aa) are penalised; do not
    compare a 60-mer's pTM to a 200-mer's as if the scale were the same.
-   **The daily credit budget is a hard wall (100/key).** Each `generate` and each
    `fold` costs one credit, so `--num-samples 8` costs 16. Budget before you
    start; the quota is shared across every ESM skill. A quota-exhaustion 429 is
    reported by the client with a *misleading* "lower BIOHUB_QPS and retry" hint —
    that hint is wrong for a quota 429; read the server body.

--------------------------------------------------------------------------------

## Pre-Report Reasoning Checklist

Complete this BEFORE writing the report. Ground every claim in a number from the
metrics JSON.

### 1. QC completeness

-   [ ] **Every reported design is folded.** pTM and pLDDT are quoted for the top
    design (and ideally the whole ranked set). No raw, unfolded sequence is
    presented as a design.
-   [ ] **The model is named.** `esm3-open-2024-03` for the design,
    `esmfold2-fast-2026-05` for the QC fold.
-   [ ] **The sample count is stated.** How many designs were drawn, ranked how.

### 2. Calibration against the floor

-   [ ] **De novo designs are compared to the random floor** (pTM ~0.20 /
    pLDDT ~0.43 at length 60), not reported in absolute terms.
-   [ ] **Means, not just the best sample**, are reported for a campaign.
-   [ ] **The pTM is read against the length.** A modest pTM on a short
    unconditioned design is expected and is stated as such.

### 3. Scaffold-specific (if applicable)

-   [ ] **Motif verbatim** confirmed from the JSON (`motif_verbatim`, grafted
    substring), and the offset checked against `--motif-range` / `--motif-start`.
-   [ ] **Motif CA RMSD quoted** and judged (< 1 Å excellent, > 3 Å failed). Both
    the sequence and the geometry checks are reported.

### 4. Sanity

-   [ ] **Top design inspected for low-complexity / TM / signal-peptide
    artifacts.** If junk, flagged, not laundered through pTM.
-   [ ] **Temperature** stated; if raised above 0.5, justified.

### 5. Integrity

-   [ ] **A poorly-folding design is reported as rejected**, with a concrete next
    step (condition it / more samples / longer), not oversold.
-   [ ] **The "computational design, not a validated protein" limitation is
    stated.**
-   [ ] **The report answers the user's actual request** (design what they asked
    for), and cross-references the right sibling skill if their real need is
    folding, inverse folding, mutation scoring, or function prediction.

--------------------------------------------------------------------------------
