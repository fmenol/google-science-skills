# Example: Positive Controls (Ubiquitin, Lysozyme, Carbonic Anhydrase II)

**Sequences**: human ubiquitin (76 aa), hen egg-white lysozyme (129 aa), human
carbonic anhydrase II (260 aa). **Task**: predict function from sequence alone
with the ESM3 function track, then verify every predicted accession against the
authoritative InterPro entry. **Verdict**: **all three correctly recognised; all
seven predicted accessions confirmed real.**

## Why this is a good example

These three proteins are deeply characterised, so their true domains are known
and we can grade ESM3's predictions against reality. They exercise the three
coverage regimes you will meet in practice:

1.  **Whole-chain single domain** (ubiquitin, coverage **0.99**): one domain,
    reported by 10 redundant annotations that collapse to a single span 1-75.
2.  **A domain embedded in a longer chain** (lysozyme, coverage **0.69**): the
    catalytic domain sits at residues 27-115 of 129, so coverage is partial *by
    design*, not by weakness.
3.  **A large whole-chain domain with rich keyword decoration** (CA2, coverage
    **0.98**): 17 annotations over a single 4-258 span, including the correct
    `zinc` / `lyase` catalytic keywords.

## What it teaches

-   **Read domains off `domain_architecture`, not `len(annotations)`.** Every one
    of these proteins has exactly **one** domain but 10-17 annotations.
-   **Verify accessions — and see them pass.** All seven predicted `IPRxxxxxx`
    accessions resolve to real InterPro entries whose names match the ESM3 labels
    exactly. `IPR001916` really is Glycoside hydrolase family 22 (the lysozyme
    family); `IPR001148` really is the alpha carbonic anhydrase domain.
-   **Keywords are noisier than entries.** Lysozyme's accession-bearing entries
    are perfect, yet the same run emits a spurious `aminoacyl trna` keyword. The
    report leads with the verified entries and flags the stray keyword.

## Key Takeaway

Even on textbook proteins, ESM3 emits redundant, nested labels and the occasional
keyword confabulation. The domain call is trustworthy **once the accession is
verified**; the raw label list is not. Count domains from `domain_architecture`,
lead with verified InterPro entries, and treat unverified keywords as hints only.

## Reproduce (replay, no credits)

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
S=skills/esm3_function_prediction/scripts/predict_function.py
uv run --no-project $S predict --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG --id ubiquitin --output ubiquitin_function.json
uv run --no-project $S plot --input ubiquitin_function.json --output ubiquitin_domains.png
# ...and likewise for lysozyme and ca2 (sequences in evals/fixtures.py).
```
