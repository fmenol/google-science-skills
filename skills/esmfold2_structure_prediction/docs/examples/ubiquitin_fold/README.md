# Example: Confident Single-Chain Fold (Ubiquitin)

**Input**: human ubiquitin (P0CG48, 76 aa), sequence only
**Model**: `esmfold2-fast-2026-05`, single sequence, num_loops 10,
num_sampling_steps 50, `--include-pae`
**Verdict**: **Confident, well-folded structure** — mean pLDDT **0.82** (0-1
scale, i.e. 82/100), pTM **0.78** (plausible fold)

## Why this is a good example

This is the reference **positive** case: a real, small, well-characterised
protein that ESMFold2 folds confidently from sequence alone.

1.  **The pLDDT scale, read correctly.** The mean pLDDT is `0.823`. On this API
    that is a **good** structure (82/100), not a poor one. The report states the
    scale explicitly and never treats 0.82 as if it were on the 0-100 scale.
    This is the single most common misread of ESMFold2 output.
2.  **Bands over a single mean.** 93.4 % of residues are in the confident-or-
    better bands; only a short 5-residue tail is low-confidence. The distribution,
    not just the average, is what makes the fold trustworthy.
3.  **A low-confidence region that is biologically correct.** The one dip
    (residues 72-76, mean pLDDT 0.56) is the mobile C-terminal `RLRGG` tail whose
    `LRGG` conjugation motif is genuinely flexible in solution. A low pLDDT there
    is the model *correctly* reporting a real flexible element — not an error.
4.  **PAE reads the topology.** The heatmap is a single dark block (one rigid
    domain) with a pale fringe for the flexible tail — exactly what pLDDT implies.

## Key Takeaway

Report the **scale**, the **band distribution**, and the **low-confidence regions
as ranges** — not just a single number. A confident fold still has flexible
parts, and naming them (mobile terminus, loop) is more useful than hiding them in
an average. Compare with the sibling `scramble_control` example to see what the
*absence* of a fold looks like on the same 76 residues.
