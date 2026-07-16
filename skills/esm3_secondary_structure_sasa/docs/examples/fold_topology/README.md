# Worked example — fold topology from sequence

This example runs `esm3_secondary_structure_sasa` on two proteins with
textbook-opposite folds and contrasts the predictions. It shows that ESM3
recovers the *topology class* of each protein from sequence alone — without ever
folding it — and that the SASA track reproduces the basic biophysics of a buried
hydrophobic core under a charged, solvent-exposed surface.

The two proteins:

-   **Ubiquitin** (human, P0CG48, 76 aa) — a **β-grasp** fold: an N-terminal
    β-hairpin, one prominent α-helix, and a mixed β-sheet packed against it.
    Mixed α/β.
-   **Myoglobin** (sperm-whale, P02185, 153 aa) — the classic **all-α globin**:
    eight helices (A–H), no β-strand.

They are the same fixtures the eval uses, so the calls are already recorded in
the cassette and replay for free.

## Reproduce (no key, no credits)

The `predict` and `plot` calls are served from the recorded eval cassette, so
this runs offline and spends nothing:

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
S=skills/esm3_secondary_structure_sasa/scripts/predict_ss_sasa.py
D=skills/esm3_secondary_structure_sasa/docs/examples/fold_topology

uv run --no-project $S predict \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --temperature 0.1 --output $D/ubiquitin_ss.json
uv run --no-project $S predict \
  --sequence VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG \
  --temperature 0.1 --output $D/myoglobin_ss.json

uv run --no-project $S plot --input $D/ubiquitin_ss.json --output $D/ubiquitin_ss.png
uv run --no-project $S plot --input $D/myoglobin_ss.json --output $D/myoglobin_ss.png
```

The committed figures were rendered from a one-record FASTA per protein, so their
titles read `ubiquitin` / `myoglobin`; the record id is the FASTA header, and
`--sequence` instead labels the record `query`. The predicted SS8 and SASA are
identical either way — only the label differs.

## Files

-   `report.md` — the contrastive analysis, with real numbers and both figures.
-   `ubiquitin_ss.json`, `myoglobin_ss.json` — the raw per-residue predictions.
-   `ubiquitin_ss.png`, `myoglobin_ss.png` — the SS8 ribbon over the SASA profile.
