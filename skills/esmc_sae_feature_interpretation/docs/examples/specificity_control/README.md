# Example: Specificity Control & Model Limitation

**Proteins**: Hen egg-white lysozyme (P00698, 129 aa), its **E35Q** point mutant,
and human **ubiquitin** (P0CG48, 76 aa) as the unrelated negative control.
**Question**: *Are SAE feature descriptions specific to a protein's biology, or
generic boilerplate? And what does a feature-set comparison actually measure?*
**Verdict**: The descriptions are **specific** (ubiquitin scores 0/10 on lysozyme
terms), and the Jaccard feature set is a fingerprint of **identity/fold/family** —
**not** of mutation effect.

## Why this is a critical example

The positive lysozyme case only *means* something if the same feature labels do
not fire on an unrelated protein. This example supplies that control and, in
doing so, exposes a limitation you must respect.

1.  **Specificity (the negative control).** Lysozyme's diagnostic terms
    (peptidoglycan / muramidase / lysozyme / glycoside / hydrolase / cell wall)
    hit **9/10** of lysozyme's top features but **0/10** of ubiquitin's.
    Ubiquitin instead activates an entirely ubiquitin-appropriate set, led by
    feature `3995` ("Ubiquitin-like domain detector"). Because the terms do *not*
    leak onto an unrelated protein, they are diagnostic rather than generic.
2.  **Jaccard separates related from unrelated.** A single active-site point
    mutation (E35Q, the classic activity-killer) gives
    `jaccard(lysozyme, E35Q) = 1.00` (50/50 shared), while
    `jaccard(lysozyme, ubiquitin) = 0.02` (2/98 shared). Unrelated proteins share
    almost nothing; the two features they *do* share are generic
    compartment/terminus signals, not biology.
3.  **The limitation: a feature set is not a mutation-effect score.** That
    `jaccard(lysozyme, E35Q) = 1.00` is the teaching point. The top-50 feature
    *set* is completely unchanged by a mutation known to abolish catalysis — even
    though the per-residue activations *do* move, and move most at residue 35.
    The feature set answers "is this the same kind of protein?", never "does this
    mutation matter?".

## Key Takeaways

**SAE feature descriptions are auto-generated hypotheses, not curated
annotations.** Read them carefully, not by keyword: ubiquitin's #2 feature `7865`
("ERAD–p97 UBL/UBX proteostasis") legitimately mentions peptide:N-glycanase
(NGLY1/PNGase), so a naive "glycan" keyword match would have wrongly flagged
ubiquitin as glycan-active. This is exactly why the eval rejected "glycan" and
"disulfide" as diagnostic terms — they leak. Trust `top_swissprot` and the full
`summary`, not the label alone.

**Match the tool to the question:**

-   To compare **identity / fold / family**, use this skill's `compare`
    (Jaccard). A high value means "same kind of protein".
-   To score a **mutation's functional effect**, use
    **`esmc-mutation-effect-scoring`**. A Jaccard of 1.00 for E35Q does **not**
    mean the mutation is harmless.
-   For a **ground-truth annotation** (GO/EC/InterPro), use
    **`esm3-function-prediction`** or UniProt — never an SAE label.

## Reproduce (replay, no credits)

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
S=skills/esmc_sae_feature_interpretation/scripts/sae_features.py
LYZ=KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL
MUT=KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFQSNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL
UBQ=MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG

uv run --no-project $S extract --sequence $UBQ --output results/ubiquitin.npz
uv run --no-project $S describe --features results/ubiquitin.npz \
    --limit 10 --output results/ubiquitin_described.json
uv run --no-project $S compare --sequence-a $LYZ --sequence-b $MUT \
    --top-n 50 --output results/cmp_lyz_mut.json   # -> jaccard 1.00
uv run --no-project $S compare --sequence-a $LYZ --sequence-b $UBQ \
    --top-n 50 --output results/cmp_lyz_ubq.json   # -> jaccard 0.02
uv run --no-project $S plot --features results/ubiquitin.npz \
    --feature-ids 3995,7865,3230 --output ubiquitin_feature_tracks.png
```

Files here: `report.md` (the analysis), `ubiquitin_feature_tracks.png` (Fig 1).
