# Interpretation Guide

This guide covers how to read SAE feature output: the two rankings, how to read a
feature description, the Jaccard feature-set comparison, the calibrated
thresholds, and a **"Model Limitations & Scientific Integrity"** section you must
internalise before writing anything. It closes with a pre-report reasoning
checklist. Read this before interpreting any run.

> [!CRITICAL] **SAE feature descriptions are auto-generated HYPOTHESES, not
> curated annotations.** They were produced by an automated agent that read each
> feature's activation pattern across large protein databases. A high activation
> means *the model associates this concept with this region* — it is **NOT**
> evidence that the protein truly has that function. Never state a feature
> description as established fact. This rule governs every sentence of every
> report; the [Model Limitations](#model-limitations--scientific-integrity)
> section below makes it operational.

--------------------------------------------------------------------------------

## What the SAE gives you

ESMC turns each residue into a dense embedding whose individual dimensions are
*polysemantic* — one axis mixes many unrelated concepts. The Sparse Autoencoder
re-expresses that vector in a **16,384-feature basis where exactly k=64 features
are active at every residue** (a top-k SAE read off layer 60 of ESM C 6B). Each
feature *tends* to correspond to one recognisable concept, so the features that
fire on your protein are candidate answers to "what does the model see here?".

The script (`sae_features.py`) computes every statistic directly from the sparse
`(L, 64)` arrays. **Never densify** the `(L, 16384)` matrix, never re-derive
activations by hand, and never re-trim BOS/EOS — the arrays are already trimmed,
so **row `i` is residue `i`** (0-based); residue `i` is reported at 1-based
position `i+1`.

--------------------------------------------------------------------------------

## The two rankings: read BOTH

`top-features` and `describe` emit two rankings. They answer different questions,
and reporting only one is a misreading.

| Ranking | What it surfaces | Track shape |
| :--- | :--- | :--- |
| **Max activation** | Motif-like / **local** features: candidate catalytic residues, binding motifs, short conserved patterns | A sharp spike over a few residues, ~zero elsewhere |
| **Prevalence** (residues above threshold) | Domain- / **family-level** features: fold class, domain identity, protein-family or compartment signal | Broadly on across much of the chain |

**A feature high on _both_ lists is usually the protein's dominant domain.** In
the worked lysozyme example, feature `2773` ("Cell wall glycan hydrolase cores")
is #1 by max activation (1.097) and #9 by prevalence (82/129) — the dominant
catalytic domain. Feature `8581` ("GlcNAc glycan-degrading domains") is #2 on
both lists.

**Why you cannot skip the prevalence list.** Prevalence surfaces features that
are quietly on across the whole chain but never spike — precisely the ones
max-activation ranking buries. Lysozyme's prevalence list is topped by feature
`4916` ("Secretory luminal ectodomains", 128/129 residues, peak only 0.397): a
compartment-level signal invisible on the max-activation list. That broad, flat
plateau is often the single most informative statement about what the protein
*is* (here: a secreted, disulfide-rich ectodomain), even though it never spikes.

> [!NOTE] The two lists overlap on the dominant domain and diverge on everything
> else. If you report only the sharp motif features, you miss the domain-level
> feature that tells you the protein's identity; if you report only the broad
> features, you miss the candidate catalytic/binding motifs. Read both, and say
> which list each cited feature came from.

--------------------------------------------------------------------------------

## Reading a feature description

`describe` fetches, for each ranked feature, a record from the unauthenticated
feature API and joins it onto the ranking. The fields, in decreasing order of how
much weight you should put on them:

| Field | What it is | How to use it |
| :--- | :--- | :--- |
| `top_swissprot` | The real UniProt entries that activate this feature most strongly, with activation values | **The single most useful field.** Check whether those proteins actually resemble yours. If they do, the hypothesis is credible; if they are unrelated, distrust the label. |
| `threshold` | The per-feature activation above which this feature is "meaningfully on" | Calibrated per feature (observed range ~0.14–0.71). More reliable than any global cutoff for deciding whether a feature is genuinely active. |
| `activation_pattern` | Where the feature tends to fire (domain-level vs residue-level, which residue types) | Sanity-check against your own track (`plot`): does the observed peak match the described pattern? |
| `exemplar_protein_families` | Protein families the automating agent associated with the feature | Context, not proof. |
| `category` | One-word bucket (e.g. `Catalytic function`, `Domain`, `Membrane-associated`) | Coarse; do not over-read. `category` can disagree with the biology (see below). |
| `summary` / `description` | The prose hypothesis | The lead to verify. Quote it with a hedge, never as fact. |
| `label` | The short human-readable name | Convenient handle; the weakest evidence. |

**`top_swissprot` is your trust check.** For lysozyme's #1 feature `2773`, the
strongest SwissProt activation is **P00698 at 18.985** — and P00698 *is* hen
egg-white lysozyme, the exact entry our sequence corresponds to. When the top
exemplar is the very protein you fed in, the hypothesis has earned a lot of
credibility. When it is not, downgrade the label to a guess.

> [!CAUTION] **`category` is not a curated ontology.** Feature `1924`
> ("Peptidoglycan hydrolase catalytic helix"), whose description is squarely
> about muramidase catalytic helices, is filed under `category:
> Post-translational modification` — a mismatch. Read the `summary` and
> `top_swissprot`, not the bucket.

--------------------------------------------------------------------------------

## Activation magnitude

With the default TF-IDF normalization on (`normalize_features=True`), activations
run **~0–1.2**. As a global rule of thumb:

-   **≥ ~1.0** — a strong association (the model is confident this concept is
    present in this region).
-   **~0.2–1.0** — a moderate association; lean on the per-feature `threshold`.
-   **< ~0.2** — weak; usually below the feature's own `threshold`.

Prefer the **per-feature `threshold`** from `describe` over any global cutoff: it
is calibrated per feature. Observed thresholds in the worked examples span
`8581` = 0.144 up to `3995` = 0.707. A feature whose peak on your protein is
below its own threshold is not meaningfully on, however high it ranks.

> [!NOTE] Normalization changes the rankings. `normalize_features=True` is
> **rejected by the API for ESMC 300M SAE models** — pass `--no-normalize`
> there. Never mix normalized and unnormalized runs in one comparison.

--------------------------------------------------------------------------------

## Jaccard feature-set comparison

`compare` takes the top-N features (by peak activation) of two sequences and
reports the **Jaccard overlap** of the two sets, plus the shared and distinctive
features. `--top-n 50` is the calibrated default.

| Jaccard (top-50) | Reading |
| :--- | :--- |
| **> 0.7** | Near-identical machinery (point mutants, close homologs) |
| **0.3 – 0.7** | Substantial: shared domains, fold, or family |
| **0.1 – 0.3** | Modest: some shared machinery, largely distinct proteins |
| **< 0.1** | The model represents these as unrelated proteins |

**Reference points measured against this API** (reproduced in replay for the
worked examples):

-   Lysozyme vs. its own **E35Q point mutant**: Jaccard = **1.00** (50/50 shared).
-   Lysozyme vs. **ubiquitin** (unrelated): Jaccard = **0.02** (2/98 shared).

Unrelated proteins genuinely share almost nothing; the two features lysozyme and
ubiquitin *do* share are generic compartment/terminus signals ("Mature secretory
lumenal domains", "C-terminal helix–tail motif"), not biology-specific ones.

> [!IMPORTANT] **What Jaccard does and does not tell you.** The feature *set* is
> a coarse fingerprint of **identity / fold / family**. It is deliberately robust
> to a single substitution: the classic activity-abolishing lysozyme mutation
> E35Q leaves the top-50 set *completely unchanged* (Jaccard 1.00), even though
> the per-residue activations do move — and move most at residue 35 (the
> mutated site). A high Jaccard therefore means "the model sees the same kind of
> protein", **not** "the mutation is harmless". To score a mutation's functional
> effect, use **`esmc-mutation-effect-scoring`**, not this skill.

--------------------------------------------------------------------------------

## Model Limitations & Scientific Integrity

> [!CRITICAL] **The descriptions are hypotheses. Treat every one as a lead to
> verify, never as an established fact.** This is the single most important rule
> of the skill.

**Why the descriptions are hypotheses, not annotations.** They were generated
by an automated agent that inspected each feature's activation pattern across
protein databases and wrote a natural-language guess. They are not from UniProt,
InterPro, or any curated source. They can be incomplete, mis-bucketed, or simply
wrong — most often for rare biology, where the feature fired on few training
proteins.

**Language rules — apply to every sentence that cites a feature:**

-   **Right:** "ESMC's SAE most strongly activates a feature *hypothesised* to
    detect peptidoglycan hydrolase cores; its strongest database exemplar is
    lysozyme itself (P00698)."
-   **Wrong:** "This protein is a peptidoglycan hydrolase." (States a hypothesis
    as fact and substitutes for a real annotation.)
-   Prefer verbs of association: *"the model associates…"*, *"a feature
    hypothesised to…"*, *"activates a feature whose exemplars are…"*.
-   Never let a feature label stand in for a database annotation. If the user
    needs ground truth (GO / EC / InterPro / catalytic activity), direct them to
    UniProt / InterPro, or to **`esm3-function-prediction`**.

**What the SAE does NOT tell you:**

-   **Ground-truth function.** A feature firing is an association, not an assay.
-   **The effect of a mutation.** The feature set is near-invariant to single
    substitutions (see Jaccard, above). Use `esmc-mutation-effect-scoring`.
-   **A 3D structure or contacts.** Use `esmfold2-structure-prediction`.
-   **Novel biology absent from the codebook.** Only 16,384 concepts exist; a
    genuinely novel motif may have no dedicated feature and will be split across
    generic ones.

**The negative control that makes descriptions trustworthy.** Lysozyme's terms
(peptidoglycan / muramidase / lysozyme / glycoside / hydrolase / cell wall) are
diagnostic *only because* they do not fire on an unrelated protein: **9/10** of
lysozyme's top features hit them, **0/10** of ubiquitin's do. If you suspect a
label is generic boilerplate, run an unrelated protein and check whether the same
label appears. If it does, it is not diagnostic. (This is exactly why the eval
rejected "glycan" and "disulfide" as keywords: they leak into ubiquitin's
ERAD/N-glycanase feature `7865`.)

--------------------------------------------------------------------------------

## Sanity checks before you interpret anything

-   **The #1 feature should make sense for the protein.** Ubiquitin's is
    "Ubiquitin-like domain detector"; lysozyme's are cell-wall / peptidoglycan
    hydrolase features. If the top features look unrelated to a protein you know,
    suspect the **input** (wrong sequence, truncation, a scramble), not the
    biology.
-   **Rows == residues.** BOS/EOS are already stripped. Do not apply `[1:-1]`
    again — that offset is only for raw `logits`/`embeddings` payloads.
-   **Both rankings read.** Confirm you looked at max-activation *and* prevalence.
-   **`top_swissprot` checked** for every feature you plan to quote.
-   **Peak matches the description.** Overlay the track (`plot`) against the
    feature's `activation_pattern`; a mismatch is a reason to distrust the label.

--------------------------------------------------------------------------------

## Pre-Report Reasoning Checklist

Complete this BEFORE writing the report. Ground every claim in a real number
from the run and a real field from the feature record.

### 1. Extraction sanity

-   [ ] **Rows == residues**: the summary reports `sequence_length` == your input
    length and 64 active features/residue. No re-trimming applied.
-   [ ] **Model recorded**: `sae_model` and `esmc_model` noted; normalization
    state noted (and `--no-normalize` used for any 300M SAE).

### 2. Both rankings

-   [ ] **Max-activation list** read: which features are sharp/local, and where
    do they peak (1-based residue)?
-   [ ] **Prevalence list** read: which features are broadly on, and which appear
    *only* here (not on the max list)?
-   [ ] **Dominant domain** identified: which feature is high on both lists?

### 3. Descriptions as hypotheses

-   [ ] **`top_swissprot` checked** for each quoted feature — do the exemplars
    resemble the protein?
-   [ ] **Per-feature `threshold`** compared against the observed peak.
-   [ ] **Hedged language** used everywhere ("the model associates…"). No feature
    stated as fact; no feature substituted for a database annotation.
-   [ ] **The one miss noted**: which top feature is generic/compartment-level
    rather than protein-specific? (For lysozyme, `12315` "Mature secretory
    lumenal domains" — true but not diagnostic.)

### 4. Comparison (if `compare` was run)

-   [ ] **Jaccard cited** with the top-N and the shared/union counts.
-   [ ] **Reading matches the calibration table**, and the "identity, not
    mutation effect" caveat is stated if a point mutant was involved.

### 5. Report readiness

-   [ ] **Figures exist and are referenced** with relative paths; every embedded
    PNG is present in the folder.
-   [ ] **The hypotheses-not-facts disclaimer** appears prominently in the
    report.
-   [ ] **The user's question is answered directly**, and the right sibling skill
    is named if they actually need function, mutation effect, or structure.

--------------------------------------------------------------------------------
