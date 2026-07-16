---
name: esmc-sae-feature-interpretation
description: >
  Mechanistic interpretability for proteins. Decomposes ESMC's internal
  representation of a protein into ~16,384 sparse, human-interpretable features
  via a Sparse Autoencoder, finds which features fire on your sequence and
  where, and fetches a natural-language description of each. Use when the user
  asks "what does the model see in this protein", "which motifs/domains does
  ESM detect", "why does the model think this is an enzyme", "what features are
  active at residue N", or wants to compare two proteins by shared internal
  machinery. Do not use when the user wants dense embedding vectors for
  similarity search or clustering (use `esmc-protein-embeddings`), a real
  function annotation such as GO/InterPro terms (use `esm3-function-prediction`
  — SAE descriptions are auto-generated hypotheses, not curated annotations),
  or a mutation's fitness effect (use `esmc-mutation-effect-scoring`).
---

# ESMC SAE Feature Interpretation

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If
    `.licenses/esmc_sae_feature_interpretation_LICENSE.txt` does not already
    exist in the workspace root directory then (1) prominently notify the user
    to check the terms at https://biohub.org/acceptable-use-policy/ and
    https://biohub.ai/, then (2) create the file recording the notification text
    and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for and request
    this key.

## Overview

When ESMC reads a protein, each residue becomes a dense embedding vector. Those
vectors are information-rich but uninterpretable: any single dimension mixes many
unrelated concepts (*polysemanticity*). A **Sparse Autoencoder (SAE)** is an
interpretability model trained to re-express that dense vector in a much larger
but sparse basis, where individual coordinates tend to correspond to one
recognisable biological concept.

The released codebook is a **top-k SAE**: 16,384 learned features, of which
**exactly k=64 are active at every residue**. This skill extracts those
activations, ranks the features that fire on your protein, and joins them to
natural-language descriptions.

It answers: *which learned concepts does ESMC detect in this protein, where
along the chain, and do two proteins share internal machinery?*

**Do NOT use when:**

-   The user wants dense embedding vectors for similarity search, clustering, or
    as ML features — use **`esmc-protein-embeddings`**.
-   The user wants a genuine functional annotation (GO terms, EC numbers,
    InterPro domains) — use **`esm3-function-prediction`**. SAE descriptions are
    hypotheses, not annotations, and must never be substituted for a real one.
-   The user wants to score the effect of a mutation — use
    **`esmc-mutation-effect-scoring`**.
-   The user wants a 3D structure — use **`esmfold2-structure-prediction`**.

## Core Rules

-   **SAE feature descriptions are auto-generated HYPOTHESES, not curated
    annotations.** They were produced by an automated agent that inspected each
    feature's activation pattern across large protein databases. A high
    activation means *the model associates this concept with this region* — it
    is **NOT** evidence that the protein truly has that function. **ALWAYS**
    present them with that caveat ("the model associates…", "a hypothesised
    feature…"). **NEVER** state a feature description as established fact, and
    never let one stand in for a real database annotation. If the user needs
    ground truth, direct them to UniProt/InterPro.
-   **Before interpreting results, read
    [`docs/interpretation-guide.md`](docs/interpretation-guide.md).** It covers
    the two rankings, how to read a feature description (and use `top_swissprot`
    as the trust check), the Jaccard calibration, and the pre-report checklist.
-   **Review at least one worked example in
    [`docs/examples/`](docs/examples/) before writing a report** — including the
    negative
    [`specificity_control`](docs/examples/specificity_control/report.md) case, so
    you know what a non-match looks like and why the descriptions are hypotheses,
    not facts (ubiquitin scores 0/10 on lysozyme terms; a point mutant and an
    unrelated protein are separated only by Jaccard, 1.00 vs 0.02).
-   **Write the report using
    [`docs/report-templates.md`](docs/report-templates.md).** Keep the
    hypotheses-not-facts disclaimer near the top and hedge every functional
    claim.
-   **BOS/EOS are already trimmed.** `eb.sae_features()` (and therefore every
    array this skill writes) has the BOS and EOS tokens removed: **row `i` is
    residue `i`** (0-based), so residue `i` is reported at 1-based position
    `i+1`. **NEVER** re-trim, and never apply the `[1:-1]` offset that raw
    `logits` payloads need.
-   **ALWAYS use `uv run --no-project`.** Without `--no-project`, `uv` walks up
    the directory tree, finds an unrelated `pyproject.toml`, and tries to build
    that project instead.
-   **Do not compute rankings, prevalence, or Jaccard overlaps yourself; always
    use the script's output.** Do not attempt to densify the feature matrix or
    re-derive activations by hand — the sparse arithmetic is easy to get wrong.
-   **Read BOTH rankings.** Peak activation and prevalence answer different
    questions (see *Interpreting the Output*). Reporting only one is a
    misreading.
-   **NEVER download model weights.** No `torch`, `transformers`, or the `esm`
    PyPI package. Everything runs remotely against the Biohub API.
-   `normalize_features=True` (the default TF-IDF normalization) is **rejected by
    the API for ESMC 300M SAE models**. Pass `--no-normalize` for those.
-   If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All commands run from the skill directory. The default SAE is
`esmc-6b-2024-12-sae-layer60-k64-codebook16384`, which runs on `esmc-6b-2024-12`.
SAE codebooks are named `{esmc_model}-sae-layer{L}-k{k}-codebook{C}`.

**1. `extract` — sequence → sparse SAE features**

Writes a compact sparse `.npz` (the `(L, 64)` indices + values, never a dense
`(L, 16384)` matrix) plus a JSON summary. Every other subcommand consumes this
`.npz`.

```bash
uv run --no-project scripts/sae_features.py extract \
  --sequence KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL \
  --output results/lysozyme.npz
```

Accepts `--fasta` instead of `--sequence` (first record is used).

**2. `top-features` — rank the features that fire**

Emits **two** rankings. `--limit` is required.

```bash
uv run --no-project scripts/sae_features.py top-features \
  --features results/lysozyme.npz \
  --limit 15 \
  --output results/lysozyme_top.json
```

`--threshold` (default `0.01`) sets the activation above which a feature counts
as "on" at a residue, for the prevalence ranking.

**3. `describe` — fetch the natural-language descriptions**

Fetches each feature's description from the (unauthenticated) feature API,
caches it on disk permanently — the descriptions are static — and joins it onto
the rankings. `--report` additionally writes a human-readable markdown report
carrying the hypothesis disclaimer.

```bash
uv run --no-project scripts/sae_features.py describe \
  --features results/lysozyme.npz \
  --limit 10 \
  --output results/lysozyme_described.json \
  --report results/lysozyme_report.md
```

To look up specific features directly, with no `.npz`:

```bash
uv run --no-project scripts/sae_features.py describe \
  --feature-ids 2773,8581,1924 \
  --output results/features.json
```

**4. `plot` — per-residue activation track**

```bash
uv run --no-project scripts/sae_features.py plot \
  --features results/lysozyme.npz \
  --feature-ids 2773,8581,1924 \
  --output results/lysozyme_tracks.png
```

One panel per feature, activation vs. residue position, peak residue annotated.

**5. `compare` — do two proteins share machinery?**

Jaccard overlap of the two top-N feature sets (ranked by peak activation), plus
the shared and distinctive features with their labels. `--top-n` is required.

```bash
uv run --no-project scripts/sae_features.py compare \
  --sequence-a <SEQ_A> --sequence-b <SEQ_B> \
  --top-n 50 \
  --output results/comparison.json
```

## Interpreting the Output

**The two rankings mean different things. Always read both.**

| Ranking | What it surfaces | Shape of the track |
|---|---|---|
| **Max activation** | Motif-like / **local** features: candidate catalytic residues, binding motifs, short conserved patterns | A sharp spike over a few residues, ~zero elsewhere |
| **Prevalence** (residues above threshold) | Domain- / **family-level** features: fold class, domain identity, protein-family or taxonomic signal | Broadly on across much of the chain |

A feature high on *both* lists is usually the protein's dominant domain.

**Activation magnitude.** With TF-IDF normalization on (the default), activations
are ~0–1.2. Peak activation near or above ~1.0 is a strong association; below
~0.2 is weak. Use the per-feature `threshold` returned by `describe` as the
feature's own "meaningfully on" level (typically 0.14–0.5) — it is calibrated
per feature and is more reliable than any global cutoff.

**Jaccard overlap** from `compare`, on top-50 sets — calibrated on this SAE:

| Jaccard | Reading |
|---|---|
| > 0.7 | Near-identical machinery (point mutants, close homologs) |
| 0.3 – 0.7 | Substantial: shared domains, fold, or family |
| 0.1 – 0.3 | Modest: some shared machinery, largely distinct proteins |
| < 0.1 | The model represents these as unrelated proteins |

Reference points measured against this API: a single-point mutant of lysozyme
(E35Q) scores **1.00** against wild-type lysozyme; **lysozyme vs. ubiquitin
scores 0.02**. Unrelated proteins genuinely do share almost nothing.

**Sanity check before you interpret anything:** the top feature should make sense
for the protein. Ubiquitin's #1 feature is *"Ubiquitin-like domain detector"*;
lysozyme's are *"Cell wall glycan hydrolase cores"*, *"GlcNAc glycan-degrading
domains"*, *"Lysozyme-like peptidoglycan hydrolases"*. If the top features look
unrelated to the protein, suspect the input, not the biology.

**`top_swissprot` in the `describe` output** lists the real UniProt entries that
activate a feature most strongly. This is the single most useful field for
judging whether a hypothesised description is trustworthy — check whether those
proteins actually resemble yours.

## Common Mistakes

1.  **Reporting a feature description as fact.** *"This protein is a
    peptidoglycan hydrolase"* is wrong. *"ESMC's SAE most strongly activates a
    feature hypothesised to detect peptidoglycan hydrolase cores"* is right. The
    descriptions are auto-generated and can be incomplete or simply incorrect,
    especially for rare biology.
2.  **Off-by-one on residue numbering.** BOS/EOS are already stripped. Row `i`
    of the arrays is residue `i` (0-based) = position `i+1` (1-based). Do **not**
    apply `[1:-1]` again — that is only needed for raw `logits`/`embeddings`
    payloads, not for `sae_features`.
3.  **Ranking only by max activation.** You will see the sharp motif features and
    completely miss the domain-level feature that is quietly on across the whole
    chain — often the most informative one about what the protein *is*.
4.  **Passing `normalize_features=True` to a 300M SAE.** The API rejects it. Use
    `--no-normalize` with any `esmc-300m-*` codebook. Note that normalization
    changes the rankings, so never mix normalized and unnormalized runs in one
    comparison.
5.  **Densifying the feature matrix.** `(L, 16384)` floats is 65 MB for a
    1000-residue protein and is never necessary — the scripts compute every
    statistic directly from the sparse `(L, 64)` arrays.

## Dependencies

-   **`uv`** — environment management (required).
-   **`credentials`** — safe handling of `BIOHUB_API_KEY` (required).
-   Sibling skills, cross-referenced rather than reimplemented:
    **`esmc-protein-embeddings`** (dense vectors),
    **`esmc-embedding-layer-sweep`** (which layer carries the signal),
    **`esmc-mutation-effect-scoring`** (variant effects),
    **`esm3-function-prediction`** (real function annotation),
    **`esmfold2-structure-prediction`** (map features onto a 3D structure).

## References

-   [Interpretation guide](docs/interpretation-guide.md) — how to read the two
    rankings, a feature description, and the Jaccard comparison; the calibrated
    thresholds; the **Model Limitations & Scientific Integrity** section; and the
    pre-report reasoning checklist. **Read this before interpreting any run.**
-   [Report template](docs/report-templates.md) — the scaffold for the written
    report (single-protein and two-protein comparison).
-   [Worked example: lysozyme (positive)](docs/examples/lysozyme_features/report.md)
    — top features are cell-wall / peptidoglycan hydrolases; the #1 feature's
    strongest exemplar is lysozyme itself (9/10 top features hit real biology).
-   [Worked example: specificity control (negative / limitation)](docs/examples/specificity_control/report.md)
    — ubiquitin scores 0/10 on lysozyme terms, and Jaccard separates a point
    mutant (1.00) from an unrelated protein (0.02); shows what a non-match looks
    like and why a feature set is not a mutation-effect score.
-   [ESM/Biohub API reference](references/esm-biohub-api.md) — endpoints, the SAE
    output format, the feature-description API, credits, and gotchas.
-   [`references/citation.bib`](references/citation.bib) — ESM C/ESM3, top-k
    sparse autoencoders, and protein-LM interpretability.
