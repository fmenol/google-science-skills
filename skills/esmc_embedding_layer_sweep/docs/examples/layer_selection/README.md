# Example: Choosing an ESMC Layer (lysozyme vs barnase)

**Task**: which ESMC layer best separates two protein families?
**Model**: `esmc-300m-2024-12` (31 layer rows) **Verdict**: **layer 23 / 30
(~77% depth)** — *not* the last layer, *not* the embedding layer.

## Why this is a good example

This is the canonical worked example for the skill: a clean, reproducible sweep
that reproduces the source tutorial's headline finding — **the last layer is
frequently not the best layer** — on an independent 2-family dataset.

1.  **The best layer is in the deep-middle.** The sweep crowns **layer 23**
    (77% of the way through the network), matching the "~3/4 depth" rule of
    thumb for functional / homology-level tasks. The **final layer (30)
    degrades** and the **embedding layer (0) is never chosen**.
2.  **A metric-saturation trap, handled honestly.** These two families are so
    far apart that **all 31 layers hit MCC = 1.000**. A naive `argmax` over the
    tied MCC would return layer 0 (the embedding layer) — the worst advice a
    sweep can give. The script detects the tie (`mcc_saturated: true`) and
    breaks it with cross-validated **log-loss** instead.
3.  **The figure tells the story.** The flat MCC panel says "MCC can't rank
    these"; the log-loss panel shows the real shape — a rise, a broad
    deep-middle minimum at layer 23, and a climb back up at the final layer.

## Files

-   [report.md](report.md) — the full analysis with the real numbers.
-   [lysozyme_barnase_layer_sweep.png](lysozyme_barnase_layer_sweep.png) — the
    `plot` output: MCC and log-loss vs layer.

## Key takeaway

**Read the shape of the curve, not just the peak, and never trust a saturated
MCC.** When every layer ties, MCC has no power to rank them; fall back to
log-loss, sanity-check that the winner beats the embedding layer, and treat the
result as *indicative* until the dataset is hard enough to give MCC headroom.

## Reproduce (free, from cassettes — no API credits)

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
# Inputs are the eval's exact 2-family dataset (see
# evals/eval_esmc_embedding_layer_sweep.py): 8 lysozyme + 8 barnase point
# mutants, embedded with esmc-300m-2024-12.
uv run --no-project skills/esmc_embedding_layer_sweep/scripts/layer_sweep.py \
    embed-dataset --fasta dataset.fasta --labels labels.csv \
    --model esmc-300m-2024-12 --output out/embeddings.npy
uv run --no-project skills/esmc_embedding_layer_sweep/scripts/layer_sweep.py \
    sweep --embeddings out/embeddings.npy --n-splits 5 \
    --output-json out/sweep.json --output-csv out/sweep.csv
uv run --no-project skills/esmc_embedding_layer_sweep/scripts/layer_sweep.py \
    plot --sweep out/sweep.json --output lysozyme_barnase_layer_sweep.png
```
