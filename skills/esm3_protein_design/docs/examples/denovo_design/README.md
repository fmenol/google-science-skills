# Example: De Novo Design, Folded and QC'd (POSITIVE)

**Objective**: design a 60-residue de novo protein from scratch (no
conditioning). **Command**: `generate --length 60 --num-samples 5`.
**Verdict**: **designs beat the random floor, modestly** — the expected outcome
for a short, unconditioned prompt, and a demonstration of the QC contract.

## Why this is a good example

This is the canonical de novo campaign, and it shows how to report one *honestly*:

1.  **Every sample is folded.** The five raw ESM3 sequences are meaningless until
    ESMFold2 gives each a pTM and pLDDT. The report quotes those numbers; it never
    presents a bare sequence as a "design".
2.  **The result is read against the floor, not in absolute terms.** Mean design
    pTM 0.259 only means something next to the random-string floor of 0.196
    (n=5 each) — see the sibling `random_baseline/` example. The margin is real
    but small, exactly what a 60-residue unconditioned prompt should give.
3.  **The best design is still low-complexity.** `design_1`, top-ranked at pTM
    0.376, ends in a poly-Lys tail (`...EKEKKKKKQDKKK`). High pTM did not make it
    a clean sequence — the report flags it rather than laundering it.

## Key Takeaway

A modest pTM on a short, unconditioned design is **expected, not a bug**. Report
it relative to the random floor, quote the whole ranked set (not just the best),
inspect the top design for low-complexity, and — if the user needs a confident
fold — point them to conditioning (see the `motif_scaffold/` example, pTM 0.98).

## Provenance

`denovo.fasta` / `denovo.json` are the actual skill outputs (cassette replay,
zero credits). `denovo_vs_random.png` is rendered from those metrics and the
random-baseline metrics. Models: `esm3-open-2024-03` (design),
`esmfold2-fast-2026-05` (QC fold).
