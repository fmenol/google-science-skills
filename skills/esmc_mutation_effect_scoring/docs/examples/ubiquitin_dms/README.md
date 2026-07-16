# Example: Deep Mutational Scan (Ubiquitin) — POSITIVE

**Protein**: Human ubiquitin (P0CG48, 76 aa) **Run**: `scan` on
`esmc-600m-2024-12` **Verdict**: **A cleanly constrained protein the model
recognises perfectly**

## Why this is a good example

Ubiquitin is one of the most conserved eukaryotic proteins — a hard, well-
understood target, so the model's calls can be checked against known biology.
This scan is the **positive** worked example: it shows what a real, well-formed
protein looks like end to end, and it exercises every reading in the
interpretation guide.

1.  **Whole-sequence fitness is textbook.** Pseudo-perplexity **1.05** and
    **100% wild-type recovery** — ESM C predicts the correct residue at every one
    of the 76 positions under leave-one-out masking. Compare with the scramble
    control (18.70, 6.6%).
2.  **Constraint has real spread.** Mean entropy **0.23 bits** but ranging from
    **0.007** at the invariant initiator Met to **1.70** at a flexible loop
    proline — the model distinguishes the immovable core from the tolerant
    surface.
3.  **The LLR heatmap is saturated red** (scale ±15.10): nearly every
    substitution at nearly every position is disfavoured, the visual signature of
    a heavily constrained fold.
4.  **It contains a real trap.** G76, the functionally essential C-terminal
    conjugation residue, scores as one of the *most tolerant* positions — a clean
    illustration that evolutionary likelihood is not functional essentiality.

## What is in this folder

| File | What it is |
| :--- | :--- |
| `report.md` | The written analysis, quoting the real numbers. |
| `ubiquitin_llr_heatmap.png` | The (76, 20) LLR matrix, `heatmap` subcommand. |
| `ubiquitin_entropy.png` | Per-position constraint profile, `entropy` subcommand. |
| `ubiquitin_entropy.csv` | The entropy/fraction-deleterious table behind the plot. |

## How it was generated (replay — no API credits)

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
UBQ=MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
S=skills/esmc_mutation_effect_scoring/scripts/mutation_scoring.py

uv run --no-project $S scan --sequence $UBQ --model esmc-600m-2024-12 \
    --output-dir ./ubq_scan --top 10
uv run --no-project $S heatmap --scan-dir ./ubq_scan \
    --output ./ubiquitin_llr_heatmap.png
uv run --no-project $S entropy --scan-dir ./ubq_scan \
    --output ./ubiquitin_entropy.png
# and the five named variants used in the report:
uv run --no-project $S score-variants --sequence $UBQ --model esmc-600m-2024-12 \
    --variants I44A,K48R,L8A,G76A,M1P --output ./variants.json
```

These are the exact fixture inputs the eval records, so every call replays from
the cassettes.

## Key takeaway

A real protein produces **low pseudo-perplexity, high wild-type recovery, and
entropy with structure** — some positions pinned, others free. Use the
`summary.json` ranked lists to read it; then sanity-check surprising tolerant
calls (like G76) against the biology, because the model scores *likelihood*, not
*function*.
