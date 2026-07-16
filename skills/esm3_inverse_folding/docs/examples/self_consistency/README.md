# Worked Example: Self-Consistency on the Ubiquitin Backbone

A complete inverse-folding run on a real backbone, end to end: **design →
self-consistency → recovery**, plus the **decoy control** that proves the scorer
works. Read this before running the skill for the first time.

## The case

- **Target backbone:** human ubiquitin (P0CG48, 76 aa), folded with ESMFold2
  (`esmfold2-fast-2026-05`) to pLDDT 0.820 / pTM 0.781. Using an ESM-predicted
  native backbone lets us check recovery against a known answer.
- **Design:** 3 sequences with ESM3 (`esm3-open-2024-03`) at T=0.1, via
  `generate(track='sequence', coordinates=...)` — the `/inverse_fold` endpoint is
  unsupported for this API key (HTTP 422).
- **Validation:** re-fold each design and score scTM / scRMSD against the target.

## Why this is a good example

1. **It shows the whole loop, not just `design`.** The headline lesson of this
   skill is that an inverse-folding design is worthless without self-consistency
   scores. This example runs `design`, `self-consistency`, and `recovery` so you
   see what each adds.
2. **It carries its own negative control.** Because a perfect design saturates
   scTM at 1.000, that number alone cannot distinguish a working scorer from one
   that always returns 1.0. The example therefore scores a **shuffled-
   correspondence decoy** (scTM 0.078, scRMSD 15.3 Å) through the same function —
   the control that makes the 1.000 meaningful.
3. **It separates recovery from quality.** Recovery is 100% here only because the
   backbone is ubiquitin's own; the report is explicit that scTM, not recovery,
   is the quality score.

## Files

- [report.md](report.md) — the analysis, with every real number and the plot.
- [recovery.png](recovery.png) — per-position native-sequence recovery (the
  figure `recovery` renders).
- `self_consistency.json` /
  `self_consistency.csv` — ranked scTM / scRMSD / pLDDT.
- `design.json` — the unvalidated sequences + recovery + diversity.
- `recovery.json` — per-position agreement (76/76 conserved).

## Headline numbers

| Quantity | Value | Reading |
| :------- | :---- | :------ |
| Best-design scTM | **1.000** | Adopts the target fold (PASS, > 0.8) |
| Best-design scRMSD | **~0.00 Å** | CA traces superpose exactly |
| Decoy scTM (control) | **0.078** | Rejected — the scorer discriminates |
| Native recovery | **100.0%** (76/76) | vs 5.0% chance |
| Diversity at T=1.0 | 3/3 unique, 94.7% identity | Resample at moderate T, don't overheat |

## Reproduce (free — cassette replay, no key, no credits)

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
# 1. Fold ubiquitin and write ubiquitin_backbone.pdb (setup).
# 2. Then, on that backbone:
S=skills/esm3_inverse_folding/scripts/inverse_fold.py
uv run --no-project $S self-consistency --pdb ubiquitin_backbone.pdb \
    --num-samples 3 --temperature 0.1 --output self_consistency.json
uv run --no-project $S recovery --pdb ubiquitin_backbone.pdb \
    --num-samples 3 --temperature 0.1 --output recovery.json
```
