# Example: Scramble Control (Shuffled Ubiquitin) — NEGATIVE

**Sequence**: ubiquitin's 76 residues, randomly shuffled (composition preserved,
fold destroyed) **Run**: `pseudo-perplexity` + `scan` on `esmc-600m-2024-12`
**Verdict**: **Not protein-like — the model correctly signals it has no idea**

## Why this is a critical example

This is the **negative control**, and it is the single most important sanity
check in the skill. It answers "what does a *not-a-protein* result look like?"
by feeding ESM C a sequence with the **exact amino-acid composition of ubiquitin**
but its residues shuffled, so no fold survives. If the model were merely counting
amino acids it would score this the same as real ubiquitin. It does the opposite.

| Metric | Ubiquitin | Its scramble | What it means |
| :--- | ---: | ---: | :--- |
| Pseudo-perplexity | **1.05** | **18.70** | 17.8× worse — approaching the ~20 random ceiling. |
| Wild-type recovery | 100% | **6.6%** | Collapses to near chance (1/20 = 5%). |
| Mean entropy (bits) | 0.23 | **4.12** | Hard against the 4.32-bit ceiling of total ignorance. |
| LLR heatmap scale | ±15.10 | **±3.63** | 4× flatter — no position strongly constrains anything. |
| Positions with all-substitutions-deleterious | 76 / 76 | **5 / 76** | Constraint is gone; the fold is what pinned every alternative. |

## The one real signal — and why it is a trap

The strongest single score anywhere in the scramble is **positive**: `R1M` at
**+3.63** — the model *prefers* a methionine at position 1. That is not a
discovered fitness gain; it is the **N-terminal initiator-Met positional prior**,
true of almost every protein and completely independent of the (destroyed) fold.
It is the textbook case of the interpretation guide's rule: **a positive LLR is a
flag to check against biology, not a discovery.**

## What is in this folder

| File | What it is |
| :--- | :--- |
| `report.md` | The written analysis, quoting the real numbers. |
| `scramble_llr_heatmap.png` | The (76, 20) LLR matrix — pale, low-contrast noise. |
| `scramble_entropy.png` | Per-position entropy, almost flat near the ceiling. |
| `scramble_entropy.csv` | The table behind the plot. |

## How it was generated (replay — no API credits)

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
UBQ=MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
SCR=RGITKEPQLHLQKGLLPNLKEEDYTVRIKEIMFSLKFITDRIGIQTIVGNSGLDLGVVQTQTPESAAQDDKTEKRL
S=skills/esmc_mutation_effect_scoring/scripts/mutation_scoring.py

# The pseudo-perplexity contrast (both sequences, one invocation):
printf '>ubiquitin\n%s\n>ubiquitin_scrambled\n%s\n' $UBQ $SCR > pppl.fasta
uv run --no-project $S pseudo-perplexity --fasta pppl.fasta \
    --model esmc-600m-2024-12 --output ./pppl.json
# The scramble's own scan + figures (masks already recorded by the pppl run):
uv run --no-project $S scan --sequence $SCR --model esmc-600m-2024-12 \
    --output-dir ./scr_scan --top 10
uv run --no-project $S heatmap --scan-dir ./scr_scan \
    --output ./scramble_llr_heatmap.png
uv run --no-project $S entropy --scan-dir ./scr_scan \
    --output ./scramble_entropy.png
```

The scramble is `scramble(UBIQUITIN, seed=0)` — the exact negative control the
eval records — so every call replays from the cassettes.

## Key takeaway

**Always check pseudo-perplexity first.** A value near the ~20 ceiling with mean
entropy near 4.3 bits means the sequence is not protein-like; report that and
stop. Do not mine per-position "constraint" out of what is really uniform noise,
and treat any lone positive LLR (like the initiator-Met prior here) as a flag to
verify, not a finding.
