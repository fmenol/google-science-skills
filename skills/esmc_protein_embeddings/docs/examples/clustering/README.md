# Clustering & Similarity Example

A verified, end-to-end walk-through of the `cluster` and `similarity`
subcommands on the shared fixture proteins. Every number quoted in
[report.md](report.md) comes from the runs below; the figure is the real
`clusters.png` those runs wrote. This is the worked example to review before you
cluster or compare embeddings for a user.

All embeddings here are `esmc-600m-2024-12`, mean-pooled, **final output
embedding** (no `--layer`). Cosine and PCA only mean something within one such
matrix — see [Choosing a layer](#choosing-a-layer) below.

## The task

Two protein families, five sequences each: lysozyme and barnase, unrelated folds
of similar length. Each of the ten sequences is a **single** random substitution
away from its parent (e.g. `lysozyme_F38S`, `barnase_T6H`). One residue cannot
separate the families, so the *only* signal KMeans can cluster on is the protein
family itself. Recovering it is the test that the embedding carries real
biology.

## Examples

### 1. Recovering two families — `cluster`

-   **Report**: [report.md](report.md) (§1)
-   **Figure**: [clusters.png](clusters.png)
-   **Inputs**: `families.fasta` (10 variants),
    `families_labels.csv` (ground-truth family per id)
-   **Outputs**: `clusters.csv`, `clusters.json`
-   **Result**: adjusted Rand index **1.00** (perfect recovery), silhouette
    **+0.86**. PC1 alone carries **92.9%** of the variance and puts barnase on
    the left, lysozyme on the right — the two families never touch.
-   **Key point**: a perfect Rand index here is a property of *how far apart the
    two folds are*, not of the mutations. The variants inside each family pile
    up on top of each other; the between-family gap dwarfs the within-family
    scatter.

### 2. The cosine floor — `similarity`

-   **Report**: [report.md](report.md) (§2)
-   **Inputs**: `semantics.fasta` (ubiquitin, its `R42F` point
    mutant, lysozyme, and a scramble of ubiquitin)
-   **Outputs**: `similarity.csv`,
    `similarity.json`
-   **Result**: point mutant **0.984**, unrelated real protein (lysozyme)
    **0.507**, scramble **-0.103**.
-   **Key point**: **~0.5 is the floor for two unrelated real proteins, not 0.**
    A cosine of 0.6 is not "60% similar" and does not mean the proteins are
    related. Only a sequence that is not protein-like at all (the scramble)
    falls near or below zero.

## Reproduce (free, from cassettes)

Both commands replay from the recorded cassettes — no API key, no credits, and a
cache miss is a hard error rather than a live call:

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
EX=skills/esmc_protein_embeddings/docs/examples/clustering

uv run --no-project skills/esmc_protein_embeddings/scripts/embed.py cluster \
  --fasta "$EX/families.fasta" --model esmc-600m-2024-12 --clusters 2 \
  --labels "$EX/families_labels.csv" --output "$EX/clusters.csv" --max-workers 5

uv run --no-project skills/esmc_protein_embeddings/scripts/embed.py similarity \
  --fasta "$EX/semantics.fasta" --model esmc-600m-2024-12 --top 3 \
  --output "$EX/similarity.csv"
```

The sequences are the exact fixtures the eval records
(`evals/eval_esmc_protein_embeddings.py`), so replay hits every payload.

## Choosing a layer

This example uses the **final output embedding** because the two folds are so
far apart that any layer would separate them. That is the easy case. When the
classes are subtle — structural sub-classes of *one* enzyme, say — the last
layer is often *not* the best, and an intermediate layer wins. Reach for
`--layer -2` as a blind default, and read
[../../interpretation-guide.md](../../interpretation-guide.md) before you
conclude that "the embeddings don't separate my classes."

## Best-practices checklist

-   **Same model, same layer, every sequence.** Cosine between vectors from
    different models or layers is noise. The JSON sidecars record `model` and
    `layer` — check them before combining files.
-   **Judge a pair against the rest of its matrix**, never against an absolute
    scale. `similarity.json` reports `mean_cosine`, `min_cosine` and
    `max_cosine` for exactly this.
-   **Read the numbers from the JSON/CSV**, not from the `.npy`. Do not
    re-implement cosine, PCA, KMeans or the Rand index by hand.
-   **A low Rand index is not proof the embeddings are bad.** Try `--layer -2`,
    an intermediate layer, or the `esmc-embedding-layer-sweep` skill first.
-   **Always report which model and which layer produced the numbers.**
