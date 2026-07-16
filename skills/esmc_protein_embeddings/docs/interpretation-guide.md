# Interpretation Guide: Reading ESMC Embeddings

Read this before you report a similarity, a cluster assignment, or a Rand index.
Embeddings are easy to compute and easy to over-read; the two traps below — the
cosine floor and layer choice — cause most wrong conclusions. Every number here
is measured on `esmc-600m-2024-12`, mean-pooled, and is reproduced by the worked
example in [examples/clustering/](examples/clustering/).

--------------------------------------------------------------------------------

## Trap 1 — The cosine floor (~0.5 is not 0)

Cosine similarity over ESMC embeddings has a **high floor**. Two proteins with
no shared fold, function, or homology do **not** score near 0 — they score
around **0.5**. Use these measured reference points as your scale, not
intuition:

| Pair | Cosine | Reading |
|---|---|---|
| ubiquitin vs itself | `1.000` | ESMC is deterministic: identical input -> identical vector |
| ubiquitin vs ubiquitin `R42F` | `0.984` | a single substitution barely moves the mean vector |
| ubiquitin vs lysozyme | `0.507` | **two unrelated real proteins — the floor** |
| ubiquitin vs scrambled ubiquitin | `-0.103` | composition preserved, biology destroyed |

**Rules for reading cosine:**

-   **Never read cosine as percent identity.** A cosine of `0.6` does not mean
    "60% similar" and does not mean the proteins are related. The floor for
    arbitrary real proteins is ~0.5, so anything from ~0.5 up to ~0.98 is the
    *entire* dynamic range you have to work with.
-   **Judge a pair only against the spread of its own matrix.** The `similarity`
    JSON reports `mean_cosine`, `min_cosine` and `max_cosine` for exactly this.
    The meaningful statement is "this pair is the highest of the N pairs here,"
    not "this pair scored 0.71."
-   **Below ~0 means not protein-like.** A scramble (same amino-acid
    composition, destroyed sequence) is the only thing that falls to ~0 or
    negative. A cosine near or below zero flags input that is not a real protein
    — a useful FASTA sanity check, not evidence of biological "anti-similarity."
-   **A point mutation is nearly invisible to the mean vector.** Mean-pooling
    averages one changed residue over the whole length (~0.016 cosine for a
    76-residue protein). Do **not** use `embed` or `similarity` to score a
    mutation; use `--per-residue` to localise a difference, or
    `esmc-mutation-effect-scoring` for a proper position-wise score.

--------------------------------------------------------------------------------

## Trap 2 — Choosing a layer

The layer you pool from is the **single most consequential choice** in this
skill, and the default (final output embedding) is often not the best answer.

Hidden states have `n_layers + 1` rows: row `0` is the raw embedding layer
(essentially no context) and row `n_layers` is the final block.

| Model | Blocks | Hidden rows | `--layer -2` resolves to | D |
|---|---|---|---|---|
| `esmc-300m-2024-12` | 30 | 0..30 | 29 | 960 |
| `esmc-600m-2024-12` | 36 | 0..36 | 35 | 1152 |
| `esmc-6b-2024-12` | 80 | 0..80 | 79 | 2560 |

-   **The last layer is tuned for masked-token prediction, not for your task.**
    Omitting `--layer` uses that final output embedding; it is a reasonable
    default only when classes are far apart (as in the worked example, where two
    different folds separate at any layer).
-   **Intermediate layers often beat the last layer.** In the ESM tutorial this
    skill is built from, layer 12 of the 30-block 300M model separates adenylate
    kinase structural sub-classes in the top principal components far better
    than layer 30 does. When your classes are subtle — sub-classes of one
    enzyme, not two unrelated proteins — expect an intermediate layer to win.
-   **`--layer -2` is the recommended blind default** when you have no labels.
    The next best blind guess is roughly two-thirds of the way through the
    network (e.g. `--layer 24` on the 600M model).
-   **With labels, don't guess — sweep.** The `esmc-embedding-layer-sweep` skill
    evaluates every layer against your labelled task in a *single* request and
    picks the winner. Prefer it whenever labels exist.
-   **Negative `--layer` is resolved locally**, so `--layer -2` is safe on every
    model including 6B. Passing a negative index straight to the API is an HTTP
    500 — always go through `scripts/embed.py`.
-   **One space per comparison.** Vectors from different models or different
    layers live in different spaces; their cosine is meaningless. Re-embed
    *everything* you compare with the same `--model` and `--layer`, and check
    the `model` and `layer` fields in each sidecar before combining files.

--------------------------------------------------------------------------------

## Reading the clustering metrics

From `clusters.json` (thresholds as in SKILL.md; measured values from the worked
example in parentheses):

-   **`adjusted_rand_index`** (only with `--labels`): `> 0.8` recovers the known
    groups, `0.5-0.8` partial recovery, `< 0.2` no better than chance. It is 0
    for random labelling and 1 for a perfect match. *(Worked example: 1.00.)*
-   **`silhouette`**: `> 0.5` well-separated, `0.25-0.5` weak structure,
    `< 0.25` the clusters are not really there. Reported even without labels, so
    it is your only quality signal in the unsupervised case. *(Worked example:
    +0.863.)*
-   **`explained_variance_ratio`**: if PC1+PC2 capture only a small fraction of
    the variance, the 2-D scatter is hiding structure — raise `--pca-components`
    before concluding the sequences do not separate. *(Worked example: PC1
    92.9%, PC2 4.4%, so 97.2% is captured and the picture is faithful.)*

--------------------------------------------------------------------------------

## Negative results & scientific integrity

> A low similarity or a low Rand index is a *result*, not a failure — but be
> careful about which conclusion it licenses.

-   **A low Rand index is not proof the embeddings are bad.** Before concluding
    "ESMC can't separate my classes," try `--layer -2`, an intermediate layer,
    or raise `--pca-components`; then run `esmc-embedding-layer-sweep` if labels
    exist. Report *what you tried*, not just the failing default.
-   **Do not upgrade a floor-level cosine to a relationship.** ~0.5 between two
    proteins is the baseline, not homology. Cosine over embeddings is **not** an
    alignment: it reports no residue correspondence and no significance. For
    homology or an alignment, use BLAST/MSA tools, not this skill.
-   **Do not read a mutation off the mean vector.** A near-1.0 cosine between a
    protein and its point mutant is expected and says nothing about whether the
    mutation matters. That is a different question — use
    `esmc-mutation-effect-scoring`.
-   **Always state the model and layer** behind every number you report. A
    cosine or a Rand index is meaningless without them.
-   **Report determinism honestly.** ESMC is a deterministic encoder;
    re-embedding the same sequence returns a bit-identical vector (self-cosine
    `1.000`). If two "identical" inputs disagree, something upstream differs — a
    different layer or model, or a whitespace/case difference in the sequence —
    find it before reporting.

--------------------------------------------------------------------------------

## Reference measurements

Embedding-specific anchors (this skill, `esmc-600m-2024-12`, mean-pooled output
embedding):

-   Self-cosine `1.000` (deterministic), point mutant `0.984`, unrelated real
    proteins `0.507`, scramble `-0.103`.
-   Lysozyme vs barnase variant families, KMeans k=2: adjusted Rand index
    `1.00`, silhouette `+0.863`.

Platform-wide anchors from the shared API reference §7 (see
[../references/esm-biohub-api.md](../references/esm-biohub-api.md)) — useful
cross-checks that the same fixtures behave consistently across ESM tools:

-   ESMC recovers 75/76 of ubiquitin's wild-type residues under masking, and
    ubiquitin's pseudo-perplexity is `1.05` versus `18.70` for its own scramble.
    That ~18x gap is the *mutation-scoring* echo of the embedding result above:
    the scramble is detectably non-protein-like by every ESM metric, which is
    why its cosine (`-0.103`) is the only pair to fall below the floor.

--------------------------------------------------------------------------------

## Pre-report reasoning checklist

Complete before writing any embedding report:

-   [ ] **State model + layer.** Named in the report, and identical across every
        vector compared.
-   [ ] **Read numbers from the JSON/CSV**, never from the `.npy`. Do not
        re-implement cosine, PCA, KMeans, or the Rand index by hand.
-   [ ] **Anchor every cosine to its matrix.** Quote `mean_cosine`, `min_cosine`
        and `max_cosine`; frame each pair as its rank within the run, not as an
        absolute score.
-   [ ] **Did you call ~0.5 "similar"?** If so, stop — that is the floor.
-   [ ] **Is the question really about a mutation?** If yes, this is the wrong
        skill; use `esmc-mutation-effect-scoring`.
-   [ ] **If clustering underperformed**, record which layers and
        `--pca-components` you tried before concluding the representation is
        empty.
-   [ ] **Confirm the figure exists** and is the one your numbers came from.
