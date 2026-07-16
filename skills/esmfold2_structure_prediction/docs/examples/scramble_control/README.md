# Example: No Confident Fold (Scrambled Ubiquitin — Negative Control)

**Input**: ubiquitin's 76 residues, **shuffled** (same composition, fold
destroyed; seed 0)
**Model**: `esmfold2-fast-2026-05`, single sequence, num_loops 10,
num_sampling_steps 50
**Verdict**: **No confident fold** — mean pLDDT **0.48** (0-1 scale), pTM
**0.24** (unreliable). This is the correct, informative answer.

## Why this is a good example

This is the reference **negative** case, and its whole point is that **a low-
confidence result is a real result, not a failure.** The sequence has the
*identical amino-acid composition* to the real ubiquitin in the sibling
`ubiquitin_fold` example — only the order is scrambled — so the collapse in
confidence is attributable to the destroyed fold, nothing else.

1.  **Both numbers fall, together.** pLDDT drops 0.82 → 0.48 and pTM drops
    0.78 → 0.24 relative to the real sequence. A real protein must beat its own
    scramble on **both** metrics; this one does, by a wide margin.
2.  **The bands invert.** The real fold is 93.4 % confident-or-better; the
    scramble is **0 %** confident-or-better and 69.7 % disordered. The entire
    chain (1-76) is flagged as one low-confidence region.
3.  **This is how you know the pipeline reads the sequence.** If a shuffled
    sequence scored as well as the real one, the model would be responding to
    composition, not structure. The gap proves it is not.
4.  **What to do with it.** `analyze`'s overall verdict — *largely disordered or
    unfoldable, do not use for downstream structural analysis* — is the
    recommendation. Do not send this to docking, Foldseek, or MD, and do not
    invent secondary structure to describe.

## Key Takeaway

**Report the negative plainly.** When pTM is 0.24 and the whole chain is
disordered, the honest, useful output is "ESMFold2 has no confident structure for
this sequence" — not a hedged structural story. Use this same scramble comparison
whenever a *designed* or unfamiliar sequence matters: if your candidate scores
like this, the model is telling you it has no fold.
