---
name: esmc-protein-embeddings
description: >
  Turn protein sequences into ESMC vector embeddings and use them downstream.
  Use when the user wants to embed / vectorize / featurize proteins, get
  representations or latent vectors for sequences, compare how similar two
  proteins are, build a similarity matrix, cluster or group a set of sequences
  into families, project proteins with PCA, or produce features to train a
  downstream classifier or regressor on. Also use for questions about which
  ESMC hidden layer to embed from. Do not use when the user wants to score
  mutations or predict variant effects (use `esmc-mutation-effect-scoring`),
  benchmark every layer against a labelled task with a trained probe (use
  `esmc-embedding-layer-sweep`), explain what a model has detected in
  human-readable features (use `esmc-sae-feature-interpretation`), predict a
  3D structure (use `esmfold2-structure-prediction`), or design new sequences
  (use `esm3-protein-design`).
---

# ESMC Protein Embeddings

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If `.licenses/esmc_protein_embeddings_LICENSE.txt`
    does not already exist in the workspace root directory then (1) prominently
    notify the user to check the terms at
    https://biohub.org/acceptable-use-policy/ and https://biohub.ai/, then
    (2) create the file recording the notification text and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for and request
    this key if this skill looks relevant to the user's request.

## Overview

Embeds protein sequences with ESMC on the Biohub Platform and uses the vectors
downstream. One script, `scripts/embed.py`, with four subcommands:

-   `embed` — one sequence -> a mean-pooled `(D,)` vector or a per-residue
    `(L, D)` matrix.
-   `batch` — a FASTA -> an `(N, D)` matrix, ready to be a feature table.
-   `similarity` — a FASTA -> the pairwise cosine-similarity matrix.
-   `cluster` — a FASTA -> PCA + KMeans assignments, a scatter plot, and an
    adjusted Rand index when ground-truth labels are supplied.

Everything runs **remotely**. No model weights are downloaded and torch is never
imported.

**Do NOT use when:**

-   The user wants the *effect of a mutation* (LLR, deleteriousness, fitness).
    Embeddings are the wrong tool — use **`esmc-mutation-effect-scoring`**.
-   The user wants to *choose a layer empirically* by training a probe at every
    layer against a labelled task — use **`esmc-embedding-layer-sweep`**.
-   The user wants to know *what the model has detected* in interpretable terms
    — use **`esmc-sae-feature-interpretation`**.
-   The user wants a 3D structure (**`esmfold2-structure-prediction`**) or a new
    designed sequence (**`esm3-protein-design`**).
-   The user wants a sequence *alignment* or a homology search. Cosine over
    embeddings is not an alignment and reports no residue correspondence; use
    BLAST/MSA tools.

## Core Rules

-   **NEVER download model weights.** No `torch`, no `transformers`, no
    `huggingface_hub`, no `esm` PyPI package. All inference is remote.
-   **ALWAYS run scripts with `uv run --no-project`.** The `--no-project` flag is
    mandatory: without it `uv` walks up the directory tree, finds an unrelated
    `pyproject.toml`, and tries to build *that* project instead.
-   **ALWAYS embed everything you intend to compare with the same `--model` AND
    the same `--layer`.** Vectors from different models or different layers live
    in different spaces and their cosine is meaningless. If you re-run `batch`
    with a new layer, re-run `similarity` and `cluster` too.
-   **Do not compute cosine similarity, PCA, KMeans or the Rand index yourself;
    always use the script's output.** Do not eyeball the `.npy` and reason about
    the numbers — read the CSV/JSON the script writes.
-   **NEVER read cosine as a percentage identity.** Two *completely unrelated*
    proteins still score ~0.5 (measured below). The floor is high; only
    differences within one run mean anything.
-   Reachable ESMC models are `esmc-300m-2024-12` (D=960), `esmc-600m-2024-12`
    (D=1152, the default) and `esmc-6b-2024-12` (D=2560). Any other model name
    returns HTTP 403.
-   If this skill is used, ensure this is mentioned in the output.

## Choosing a layer

This is the single most consequential choice in this skill, and the default is
often not the best answer.

ESMC's hidden states have **`n_layers + 1` rows**: row `0` is the embedding
layer (raw token embeddings, essentially no context) and row `n_layers` is the
final block.

| Model | Blocks | Hidden rows | `--layer -2` resolves to | D |
|---|---|---|---|---|
| `esmc-300m-2024-12` | 30 | 0..30 | 29 | 960 |
| `esmc-600m-2024-12` | 36 | 0..36 | 35 | 1152 |
| `esmc-6b-2024-12` | 80 | 0..80 | 79 | 2560 |

-   Omitting `--layer` uses the **final output embedding**, which is tuned for
    the pretraining objective (masked-token prediction) — not necessarily for
    *your* task.
-   **Intermediate layers often beat the last layer.** In the ESM tutorial this
    skill is built from, layer 12 of the 30-block 300M model separates adenylate
    kinase structural classes far better in the top principal components than
    layer 30 does.
-   **`--layer -2` (the second-to-last hidden state) is the recommended blind
    default** when you have no labels to test against. The next best blind guess
    is roughly two thirds of the way through the network (e.g. `--layer 24` on
    the 600M model).
-   The **principled** way is to evaluate every layer on your task and pick the
    winner. That is exactly what **`esmc-embedding-layer-sweep`** does; it is
    cheap because one `mean_hidden_states()` request returns *all* layers
    mean-pooled in a **single** call. Prefer it whenever labels exist.
-   Negative `--layer` values are resolved locally before the request is sent, so
    `--layer -2` is safe on every model, including 6B.

## Utility Scripts

All examples are runnable from the repository root. `--output` names the primary
artifact; **sidecar files share its stem** (`x.npy` -> `x.json`;
`x.csv` -> `x.json` + `x.png`).

**1. `embed` — one sequence**

Mean-pooled `(D,)` by default; `--per-residue` gives `(L, D)` with BOS/EOS
already trimmed, so row `i` is residue `i+1`.

```bash
# Mean-pooled 1152-d vector from the default 600M model -> ubi.npy + ubi.json
uv run --no-project scripts/embed.py embed \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --output /path/to/out/ubi.npy

# Per-residue (76, 1152) matrix from the second-to-last hidden layer
uv run --no-project scripts/embed.py embed \
  --fasta /path/to/ubiquitin.fasta --layer -2 --per-residue \
  --output /path/to/out/ubi_per_residue.npy
```

**2. `batch` — a FASTA to an `(N, D)` feature matrix**

Row `i` of the `.npy` is the embedding of `ids[i]` in the `.json`. Requests run
in parallel (`--max-workers`, capped at 8).

```bash
uv run --no-project scripts/embed.py batch \
  --fasta /path/to/proteins.fasta --model esmc-600m-2024-12 --layer -2 \
  --max-workers 8 --output /path/to/out/features.npy
```

**3. `similarity` — pairwise cosine matrix**

Writes the full `N x N` matrix as CSV; `--top` controls only how many
best/worst pairs are summarised in the JSON and on stdout.

```bash
uv run --no-project scripts/embed.py similarity \
  --fasta /path/to/proteins.fasta --top 10 \
  --output /path/to/out/similarity.csv
```

**4. `cluster` — PCA + KMeans**

Projects with PCA (2 components by default, as in the ESM tutorial), clusters
the projection with KMeans, and writes `clusters.csv` (assignments + PC
coordinates), `clusters.json` (summary) and `clusters.png` (scatter). Pass
`--labels` (an `id,label` CSV) to score the clustering with the adjusted Rand
index.

```bash
uv run --no-project scripts/embed.py cluster \
  --fasta /path/to/proteins.fasta --clusters 3 --layer -2 \
  --labels /path/to/families.csv \
  --output /path/to/out/clusters.csv
```

## Interpreting the Output

**Cosine similarity** (measured on `esmc-600m-2024-12`, mean-pooled output
embedding — use these as your reference points):

| Pair | Cosine | Reading |
|---|---|---|
| ubiquitin vs itself | `1.000` | ESMC is deterministic: identical input -> identical vector |
| ubiquitin vs ubiquitin R42F | `0.984` | a single substitution barely moves the vector |
| ubiquitin vs lysozyme | `0.507` | **two unrelated real proteins** |
| ubiquitin vs scrambled ubiquitin | `-0.103` | composition preserved, biology destroyed |

The critical lesson: **~0.5 is the floor for unrelated proteins, not 0.** A
cosine of 0.6 does not mean "60% similar" and does not mean the proteins are
related. Judge a pair only against the *spread of the other pairs in the same
matrix* — the JSON reports `mean_cosine`, `min_cosine` and `max_cosine` for
exactly this purpose. Values below ~0 indicate a sequence that is not
protein-like at all (scrambles, random strings).

**Clustering** (`clusters.json`):

-   `adjusted_rand_index` (only with `--labels`): `> 0.8` the clustering
    recovers the known groups · `0.5–0.8` partial recovery · `< 0.2` no better
    than chance. It is 0 for random labelling and 1 for a perfect match.
-   `silhouette`: `> 0.5` well-separated clusters · `0.25–0.5` weak structure ·
    `< 0.25` the clusters are not really there. Reported even without labels, so
    it is your only quality signal in the unsupervised case.
-   `explained_variance_ratio`: if PC1+PC2 capture only a small fraction of the
    variance, the 2-D picture is hiding structure — raise `--pca-components`
    before concluding the sequences do not separate.
-   A low Rand index is **not** proof the embeddings are bad: try `--layer -2`,
    an intermediate layer, or run **`esmc-embedding-layer-sweep`** before giving
    up on the representation.

Always tell the user which model and which layer produced the numbers you
report. For the full reading protocol — the cosine floor, layer choice, and a
pre-report checklist — see [docs/interpretation-guide.md](docs/interpretation-guide.md),
and review the [worked example](docs/examples/clustering/report.md) before
writing a report.

## Common Mistakes

-   **Comparing vectors from different layers or models.** A cosine between a
    600M layer-35 vector and a 300M layer-12 vector is noise. The sidecar JSON
    records `model` and `layer` for every artifact — check them before combining
    files.
-   **Passing a negative layer straight to the API.** The Biohub service rejects
    a negative `ith_hidden_layer` with **HTTP 500 ("tuple index out of range")**;
    only `-1`, meaning "all layers", is special-cased, and even that is rejected
    for 6B. `scripts/embed.py` resolves negative `--layer` values locally, so
    always go through the script rather than hand-rolling a request.
-   **Treating the last layer as the best layer.** It is tuned for masked-token
    prediction, not for your task. See "Choosing a layer" above.
-   **Reading `mean` pooling as if it were per-residue.** The mean vector cannot
    tell you *where* two proteins differ; a single point mutation moves it by
    only ~0.016 cosine. To localise a difference, use `--per-residue`, or use
    **`esmc-mutation-effect-scoring`** for a proper position-wise score.
-   **Forgetting `--no-project`**, or writing `--output features` without the
    `.npy` extension (the script rejects it rather than let numpy silently
    rename the file and break the sidecar).

## Dependencies

-   **`uv`** — required to run every script.
-   **`credentials`** — the safe protocol for obtaining `BIOHUB_API_KEY`.
-   **`esmc-embedding-layer-sweep`** — pick the best layer with a labelled task
    instead of guessing. Use it before this skill whenever labels exist.
-   **`esmc-mutation-effect-scoring`** — variant effects and position-wise
    scores; use instead of embeddings for anything mutation-related.
-   **`esmc-sae-feature-interpretation`** — human-readable features when the
    question is "what did the model see?" rather than "how similar are these?".

## References

-   **Review the [worked example](docs/examples/clustering/report.md) first** —
    lysozyme vs barnase variant families clustered at a real adjusted Rand index
    of 1.00, plus the cosine-floor demonstration (unrelated proteins 0.51,
    near-identical mutant 0.98, scramble -0.10). Context and the free, replayable
    reproduce commands are in its [README](docs/examples/clustering/README.md).
-   **Read the [interpretation guide](docs/interpretation-guide.md) before
    interpreting** — the cosine floor (~0.5 is the floor for unrelated proteins,
    not 0) and how to choose a layer are the two traps it exists to prevent; it
    also carries a pre-report checklist.
-   [ESM/Biohub API reference](references/esm-biohub-api.md) — endpoints,
    reachable models, the daily credit model, and the verified gotchas shared by
    every ESM skill.
-   `references/citation.bib` — cite ESMC (Candido et al., 2026) whenever this
    skill contributes to a result.
