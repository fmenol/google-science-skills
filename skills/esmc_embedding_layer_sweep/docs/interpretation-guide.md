# Interpretation Guide

How to choose an ESMC layer and how to read a sweep. Read this before quoting a
winning layer or writing up a sweep. It collects the layer-selection heuristics,
the metric thresholds, and the honesty checks the skill depends on.

--------------------------------------------------------------------------------

## 1. Layer indexing (get this right first)

Row **0 is the embedding layer** — raw token embeddings, no attention,
essentially amino-acid composition. Row `k` is the output of transformer block
`k`. An `esmc-300m` model has 30 blocks, so **31 rows indexed 0..30**;
`esmc-600m` has 36 blocks (37 rows); `esmc-6b` has 80 blocks (81 rows). Always
report layers in this **0-based** convention. The source tutorial prints 1-based
numbers ("layer 25 of 37") — do not mix the two.

## 2. Reading MCC

Matthews correlation coefficient runs from -1 (total disagreement) through 0
(random) to +1 (perfect). It uses all four confusion-matrix cells, so unlike
accuracy it stays honest under class imbalance.

| MCC | Reading |
|---|---|
| > 0.8 | Strong — the property is linearly decodable from this layer. |
| 0.5–0.8 | Useful signal. |
| 0.2–0.5 | Weak. |
| ~0 | The probe learned nothing. |

**Read the shape of the curve, not just the peak.** The expected shape is a rise
through the early and middle layers, a broad plateau, and a **drop at the final
layer**. If the reported `best_layer` is far from the last layer, say so — that
is the finding.

## 3. Default layers when you cannot sweep

These are the source tutorial's priors, to use **only** when you have no labels
to measure against. They are starting points, not answers.

| Task | Layer to use |
|---|---|
| Physicochemical / functional properties (thermostability, solubility, expression) | **`[-2]`**, the second-to-last layer |
| High-level functional similarity (EC number, GO term, homology, family) | **~3/4 of the way through** the network |
| Residue–residue contacts | **late** layers |
| Anything | **not** the embedding layer (row 0) — it has seen no context |

The worked example lands squarely on the second row: separating two protein
families (a homology / family task) is decoded best at **layer 23 / 30 ≈ 77%
depth**, i.e. ~3/4 of the way through.

## 4. Why the last layer degrades

ESMC's final layer is specialised to its **masked-token training objective** —
predicting held-out residues — not to your downstream task. Its representation
is tuned to reconstruct amino acids, so task-relevant structure that peaked
earlier is partly discarded. This is why the last hidden state is a poor default
and why the sweep exists: in the source tutorial an intermediate layer reached
**MCC 0.935 while the final layer reached only 0.843**. In the worked example
the MCC is saturated, so the degradation shows up on log-loss instead:
**0.0055 at the final layer vs 0.0046 at layer 23**.

## 5. MCC saturation and the log-loss tie-break (the key caveat)

MCC is bounded above by **1.0**. On an easy or small dataset many layers hit
that ceiling at once, and MCC then has **no power to rank them**. A naive
`argmax` over the tied vector returns the *lowest* index — **layer 0, the
embedding layer** — the worst possible advice. The script guards against this:

-   It reports `mcc_saturated: true` and `n_layers_tied_at_best_mcc` whenever
    more than one layer ties at the top.
-   It breaks the tie with cross-validated **log-loss**, a strictly proper
    scoring rule that keeps improving after the classification is already
    perfect (it rewards a more confident, better-separated boundary), so it
    still discriminates when MCC cannot. `selection_metric` records which metric
    chose the winner.

**When you see `mcc_saturated: true`, relay it to the user.** The honest reading
is that the dataset **cannot** separate those layers; the log-loss winner is
*indicative, not definitive*. Recommend a harder or larger dataset (more
classes, closer families, more sequences) before committing to a layer. In the
worked example all 31 layers tie at MCC 1.000; layer 23 wins only on a 0.0046 vs
0.0056 log-loss gap.

## 6. Sanity checks

-   **Does the winner beat the embedding layer?** Row 0 is amino-acid
    composition with no context. Confirm `mcc[best] >= mcc[0]` and, when MCC
    ties, `log_loss[best] < log_loss[0]`. If row 0 genuinely wins, the task is
    decidable from composition alone and ESMC is adding nothing — reconsider the
    problem, not the layer.
-   **Is the winner crowned at row 0?** It should never be. A sweep that returns
    layer 0 is either mis-indexing the stack or blindly `argmax`-ing a saturated
    metric.
-   **Is the std band tiny?** With fewer than ~20 sequences the fold estimates
    are noisy and can pin to the ceiling. Check `mcc_std` and the `+/-1 std`
    band before quoting a winning layer.
-   **The winner does not transfer.** The best layer for EC classification is
    not the best layer for stability regression, and the best layer for ESMC
    600M is not row-for-row the best layer for ESMC 6B. Re-sweep when the task
    or the model changes.

## 7. Negative results & scientific integrity

A layer sweep can legitimately conclude **"this dataset cannot choose a
layer."** That is a real result, not a failure, and it must be reported as such:

-   **Saturated MCC** (Section 5) means every layer is tied; do not present a
    log-loss winner as if MCC had ranked it.
-   **A flat log-loss curve** with overlapping std bands means even the
    tie-break is inconclusive; say so.
-   **Row 0 winning** means the signal is composition, not context; the sweep
    has told you ESMC is unnecessary here.

Do not manufacture a decisive layer out of a saturated or tiny sweep. Report the
uncertainty, quote the real numbers, and recommend more or harder data. Cost
note (see the [shared API reference](../references/esm-biohub-api.md) §5): the
sweep is cheap — one billed request per sequence returns all layers at once — so
enlarging the dataset, not looping over layers, is the right way to buy
confidence.

## 8. Pre-report reasoning checklist

Before writing up a sweep, confirm you can answer each of these from the sweep's
own JSON output — never by eyeballing the curve yourself:

1.  **Best layer and depth.** What is `best_layer`, and what fraction of the
    network is it (`best_layer_depth_fraction`)? Is it far from the last layer?
2.  **Saturation.** Is `mcc_saturated` true? How many layers tied
    (`n_layers_tied_at_best_mcc`)? Which metric chose the winner
    (`selection_metric`)?
3.  **Final-layer degradation.** Does the best layer beat the last layer, on MCC
    (`best_minus_last_mcc`) or, when tied, on log-loss?
4.  **Embedding-layer sanity.** Does the best layer beat row 0? Is row 0 *not*
    the winner?
5.  **Confidence.** Is the std band tight? Are there ≥ ~20 sequences, or should
    the result be flagged as indicative?
6.  **Transfer.** Are you reusing this layer on a different task or model? If
    so, re-sweep.
