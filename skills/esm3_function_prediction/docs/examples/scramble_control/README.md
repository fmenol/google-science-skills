# Example: Negative Control (Scrambled Ubiquitin)

**Sequence**: human ubiquitin (76 aa) with its residues randomly shuffled
(`scramble(UBIQUITIN, seed=0)`) — same amino-acid composition, no fold.
**Task**: predict function from sequence alone. **Verdict**: **0 annotations,
coverage 0.00 — no recognisable function.**

## Why this is a good example

It shows what "no confident function" looks like, and proves the tool is not
simply pattern-matching amino-acid composition. Real ubiquitin and this scramble
contain the *identical* multiset of residues; the only thing removed is their
order — and with it, the fold. ESM3's answer changes completely:

| Sequence | Annotations | Coverage | InterPro |
| :--- | :--- | :--- | :--- |
| Real ubiquitin | 10 | 0.99 | IPR000626, IPR019956, IPR029071 |
| **Scrambled ubiquitin** | **0** | **0.00** | **none** |

## What it teaches

-   **An empty result is a finding, not a failure.** Zero annotations is exactly
    what a non-protein returns. Report it as "no recognisable function", never as
    a tool error, and never stretch to invent a domain.
-   **The empty result is informative because the confident result is real.** The
    contrast between 0.99 and 0.00 coverage on composition-identical sequences is
    the whole point: ESM3 reads the fold, not the amino-acid bag.
-   **Corroboration exists.** The same scramble raises ubiquitin's ESMC
    pseudo-perplexity from **1.05 to 18.70** — an independent signal that the
    shuffled sequence no longer reads as a protein at all.

## Key Takeaway

When ESM3 returns nothing, that *is* the answer: the sequence carries no
recognisable domain. It could be genuine dark matter (a real fold outside the
training data) or not a plausible protein (a scramble, a frameshift, an assembly
artefact). Distinguish the two with orthogonal evidence — but do not paper over
the empty result with a speculative annotation.

## Reproduce (replay, no credits)

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
S=skills/esm3_function_prediction/scripts/predict_function.py
uv run --no-project $S predict --sequence RGITKEPQLHLQKGLLPNLKEEDYTVRIKEIMFSLKFITDRIGIQTIVGNSGLDLGVVQTQTPESAAQDDKTEKRL --id scrambled_ubiquitin --output scrambled_function.json
uv run --no-project $S plot --input scrambled_function.json --output scrambled_domains.png
```
