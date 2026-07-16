# Example: Positive Case (Hen Egg-White Lysozyme)

**Protein**: Hen egg-white lysozyme (UniProt **P00698**, 129 aa)
**Question**: *What does ESMC see in this protein? Which motifs/domains does it
detect?*
**Verdict**: The model's top features are, hypothesis-for-hypothesis, cell-wall /
peptidoglycan hydrolase features — and the #1 feature's single strongest
database exemplar is lysozyme itself.

## Why this is a good example

This is the clean positive case: the SAE fires on features whose auto-generated
descriptions match the protein's real biology, and the output gives you the
tools to *check* that match rather than take it on faith.

1.  **The rankings agree with the biology.** 9 of the top 10 features by max
    activation carry descriptions mentioning peptidoglycan, muramidase,
    lysozyme, glycoside, hydrolase, or cell wall. The dominant domain (feature
    `2773`, "Cell wall glycan hydrolase cores") is high on *both* the
    max-activation and prevalence lists.
2.  **The trust check passes.** Feature `2773`'s strongest `top_swissprot`
    exemplar is **P00698 (18.985)** — the exact hen egg-white lysozyme entry the
    input corresponds to. When the top exemplar of the top feature is the protein
    itself, the hypothesis has earned real credibility.
3.  **Both rankings matter.** The prevalence list surfaces feature `4916`
    ("Secretory luminal ectodomains", on 128/129 residues but peaking at only
    0.397) — a broad compartment signal the max-activation list buries. It states
    what the protein *is* (a secreted, disulfide-rich ectodomain), not what it
    catalyses.

## Key Takeaway

**Even in a textbook-clean case, the descriptions are hypotheses.** The report
states every functional claim with a hedge ("the model associates…"), leans on
`top_swissprot` as the trust signal, and names the one top feature (`12315`,
"Mature secretory lumenal domains") that is a generic compartment label rather
than lysozyme-specific biology. The activation-track figure shows the full range
from a sharp local motif (feature `9784`, spiking at A31) to a broad domain
plateau (feature `4916`, on across the whole chain) — a concrete picture of why
you read both rankings.

## Reproduce (replay, no credits)

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
S=skills/esmc_sae_feature_interpretation/scripts/sae_features.py
LYZ=KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL

uv run --no-project $S extract --sequence $LYZ --output results/lysozyme.npz
uv run --no-project $S top-features --features results/lysozyme.npz \
    --limit 10 --output results/lysozyme_top.json
uv run --no-project $S describe --features results/lysozyme.npz \
    --limit 10 --output results/lysozyme_described.json \
    --report results/lysozyme_report.md
uv run --no-project $S plot --features results/lysozyme.npz \
    --feature-ids 9784,1924,2773,8581,4916 \
    --output lysozyme_feature_tracks.png
```

Files here: `report.md` (the analysis), `lysozyme_feature_tracks.png` (Fig 1).
