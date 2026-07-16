# Embedding Report: Lysozyme vs Barnase Variant Families

*Produced with the `esmc-protein-embeddings` skill. Model `esmc-600m-2024-12`,
mean-pooled, final output embedding (no `--layer`).*

## Summary

ESMC embeddings **perfectly recover** the two protein families from sequence
alone. Ten sequences — five single-substitution variants of lysozyme and five of
barnase — cluster into their true families with an **adjusted Rand index of
1.00** and a silhouette of **+0.86**. A single principal component (PC1,
**92.9%** of the variance) is enough to separate the folds. The companion
`similarity` run establishes the reading scale for cosine: a one-residue mutant
sits at **0.984**, two *unrelated* real proteins at **0.507**, and only a
biologically destroyed scramble drops to **-0.103**. The headline lesson is that
**~0.5 is the cosine floor for unrelated proteins, not 0**, so cosine must be
read relative to the spread of a matrix, never as a percent identity.

## 1. Clustering: two families, perfectly separated

![PCA of ESMC mean embeddings; colour = KMeans cluster, marker = family label;
adjusted Rand index = 1.00](clusters.png)

**Interpretation:**

-   **Separation (PC1):** PC1 captures **92.9%** of the total variance and
    splits the two families cleanly — barnase variants at PC1 ~ -0.21, lysozyme
    variants at PC1 ~ +0.21. PC2 adds only **4.4%** (97.2% cumulative), so the
    2-D picture is faithful, not a projection artifact.
-   **Cluster to label agreement:** KMeans (k=2) assigns all five barnase
    variants to cluster 0 and all five lysozyme variants to cluster 1 — a
    one-to-one match with the ground-truth labels. Hence **adjusted Rand index =
    1.00** (0 = chance, 1 = perfect).
-   **Cluster quality:** silhouette **+0.863** (`> 0.5` = well-separated
    clusters). The within-family points are nearly coincident because each
    variant differs from its parent by a single residue; the between-family gap
    dwarfs that scatter.
-   **Why it is easy:** the only variable that co-varies with the labels is the
    protein family. A single point mutation cannot move a vector across the gap
    (see §2), so KMeans has nothing to latch onto *except* fold identity — and
    it finds it.

## 2. Similarity: the cosine floor

The `similarity` run over ubiquitin, its `R42F` point mutant, lysozyme, and a
scramble of ubiquitin fixes the scale on which every cosine in this skill must
be read (full matrix in `similarity.csv`):

| Pair | Cosine | Reading |
|---|---|---|
| ubiquitin vs ubiquitin `R42F` | **0.984** | one substitution barely moves the vector |
| ubiquitin `R42F` vs lysozyme | **0.556** | mutant is still just an unrelated protein to lysozyme |
| ubiquitin vs lysozyme | **0.507** | **two unrelated real proteins — the floor** |
| ubiquitin vs scrambled ubiquitin | **-0.103** | same composition, biology destroyed |

Matrix summary (`similarity.json`): `mean_cosine` 0.368, `min_cosine` -0.103,
`max_cosine` 0.984.

**Interpretation:**

-   **The floor is high.** Ubiquitin (76 aa) and lysozyme (129 aa) share no fold
    and no function, yet they score **0.507**. Read naively, "0.51" sounds like
    a real relationship; it is the *baseline* two arbitrary proteins produce.
    Never read cosine as percent identity, and never call 0.5-0.6 "similar."
-   **A point mutation is nearly invisible to the mean vector.** The `R42F`
    mutant keeps the cosine at **0.984** — mean-pooling averages one changed
    residue over 76, so the signal is ~0.016. To localise *where* two proteins
    differ, use `--per-residue`; for a proper position-wise mutation score use
    `esmc-mutation-effect-scoring`.
-   **Only non-protein input breaks the floor.** The scramble preserves
    amino-acid composition but destroys the sequence's biology, and it is the
    only pair that falls to ~0 or below (**-0.103**). Values near or under zero
    flag input that is not protein-like (scrambles, random strings) — a useful
    sanity check on a FASTA.
-   **Judge within the matrix.** The meaningful statement is not "0.507" in
    isolation but that 0.984 >> 0.507 >> -0.103 *within this one run*. That
    ordering is the result; the absolute values are only anchors.

## 3. Method & provenance

-   **Skill / script:** `esmc-protein-embeddings`, `scripts/embed.py` (`cluster`
    and `similarity`).
-   **Model / pooling / layer:** `esmc-600m-2024-12`, mean-pooled, final output
    embedding (`layer: null` in every sidecar). D = 1152.
-   **Clustering:** PCA to 2 components, then KMeans (k=2, `random_state=0`);
    adjusted Rand index against `families_labels.csv`.
-   **Inputs:** the exact eval fixtures — five point mutants each of lysozyme
    (seeds 100-104) and barnase (seeds 200-204) for clustering; ubiquitin, its
    `R42F` mutant, lysozyme, and `scramble(ubiquitin, seed=1)` for similarity.
-   **Cost:** zero. Every embedding replayed from the recorded cassettes.

## 4. Conclusion

For proteins from **distinct families**, ESMC embeddings are a strong,
label-free family detector: KMeans over the top principal components recovered
lysozyme vs barnase at a **perfect Rand index of 1.00**, driven almost entirely
by PC1. But the same run is a warning about *reading* cosine — unrelated real
proteins already score **~0.5**, so similarity is only interpretable **relative
to the rest of its matrix**, and a single point mutation is essentially
invisible to the mean-pooled vector. When your own classes are subtler than two
different folds and fail to separate at the default layer, do not conclude the
representation is empty: try `--layer -2` or an intermediate layer, or run
`esmc-embedding-layer-sweep` — see
[../../interpretation-guide.md](../../interpretation-guide.md).
