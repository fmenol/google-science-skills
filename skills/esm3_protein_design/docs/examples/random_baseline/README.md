# Example: Random-Sequence Baseline, the Foldability Floor (NEGATIVE)

**Objective**: establish what a *non-design* folds to, so a real design has
something to be measured against. **Input**: five uniformly-random amino-acid
strings (seed 0), length 60, folded with `esmfold2-fast-2026-05`. **Verdict**:
**random strings fold poorly (mean pTM 0.196) — this is the negative control**
that makes a de novo pTM interpretable.

## Why this is a good example

This is the honest negative control every de novo campaign needs, and it exists
to defeat one specific error:

1.  **ESMFold2 assigns non-zero confidence to any string.** The random 60-mers
    come back at pTM 0.148-0.261 (mean 0.196) and pLDDT 0.41-0.46 (mean 0.433) —
    *not* zero. So "my design has pTM 0.3, therefore it folds" is a false
    inference: you have to beat the floor, and the floor is not at zero.
2.  **The floor overlaps the de novo designs.** The best random control here
    (pTM **0.261**) actually exceeds the *mean* ESM3 design (0.259, see
    `denovo_design/`). At length 60 the design signal is weak enough that a
    single sample — design or control — proves nothing. This is why the skill
    compares *means over multiple samples*.
3.  **It is a real run, not an assumption.** These are actual ESMFold2 folds of
    actual random sequences (cassette replay, zero credits), reproduced from the
    skill's eval with the identical construction (`random_sequences(5, 60,
    seed=0)`).

## Key Takeaway

Never report a de novo pTM without the random floor beside it. A design is only
"folding" if it clears **~0.20 pTM / ~0.43 pLDDT** at length 60 — and clearing it
by a hair (as the unconditioned designs do) is a *modest* result to be reported
honestly, not a mini-protein. If you want daylight above the floor, condition the
design (the `motif_scaffold/` example reaches pTM 0.98).

## Provenance

`random_baseline.json` holds the five random sequences and their ESMFold2 pTM /
pLDDT. `random_floor.png` is rendered from those metrics and the de novo means.
Construction identical to `evals/eval_esm3_protein_design.py`. Fold model:
`esmfold2-fast-2026-05`.
