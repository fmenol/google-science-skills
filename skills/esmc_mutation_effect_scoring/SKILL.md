---
name: esmc-mutation-effect-scoring
description: >
  Zero-shot prediction of how bad a protein mutation is, using the ESM C protein
  language model — no MSA, no structure, no experimental labels. Use when the
  user asks whether a variant is deleterious, pathogenic, destabilising or
  tolerated; wants a deep mutational scan (DMS) or LLR heatmap; wants to rank or
  triage point mutations (e.g. "is K48R bad?", "which residues can I safely
  mutate?", "score these variants"); wants to find mutation-tolerant positions
  for library design; or wants to know how protein-like / well-formed a sequence
  is (pseudo-perplexity). Do not use when the user wants to DESIGN or generate
  new sequences (use `esm3_protein_design`), wants a predicted 3D structure or a
  structure-based stability estimate (use `esmfold2_structure_prediction`), wants
  embeddings to feed a downstream model (use `esmc_protein_embeddings`), or wants
  a mechanistic account of what a region does (use
  `esmc_sae_feature_interpretation`).
---

# Zero-Shot Mutation Effect Scoring with ESM C

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If
    `.licenses/esmc_mutation_effect_scoring_LICENSE.txt` does not already exist
    in the workspace root directory then (1) prominently notify the user to
    check the terms at https://biohub.org/acceptable-use-policy/ and
    https://biohub.ai/, then (2) create the file recording the notification text
    and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for and request
    this key.

## Overview

Scores point mutations by asking ESM C a single question at every residue:
*"I have hidden this position — what belongs here?"*

For each position `i` the model is shown `sequence[:i] + '_' + sequence[i+1:]`
(**leave-one-out masking**) and returns a probability distribution over amino
acids at `i`, computed without ever seeing the wild-type residue there. From
that distribution we get:

-   **LLR** (log-likelihood ratio) = `log p(mutant) - log p(wild-type)` for all
    20 amino acids. Negative = the model prefers the wild type = predicted
    deleterious.
-   **Entropy** (bits) = how sure the model is about the position at all. Low =
    evolutionarily constrained.
-   **Pseudo-perplexity** = `exp(-mean log p(actual residue))` over the whole
    sequence. A whole-sequence fitness / "is this a real protein?" score.

This is the method of Meier et al. (2021), and it is a strong baseline on
ProteinGym without any task-specific training.

**This skill does NOT:**

-   Predict ΔΔG in kcal/mol. The scores are **unitless**. See *Interpreting the
    Output*.
-   Use an MSA, a structure, or any experimental data. It is sequence-only.
-   Score insertions, deletions, frameshifts or multi-residue variants. Single
    amino-acid substitutions only.
-   Design or generate sequences — that is `esm3_protein_design`.
-   Tell you whether a mutation is clinically pathogenic. It reports
    evolutionary plausibility, which correlates with, but is not, pathogenicity.

## Core Rules

-   **NEVER eyeball the LLR matrix, and NEVER compute mutation effects
    yourself.** Do not read `llr.npy` / `llr.csv` and reason about the numbers
    by hand, do not try to rank variants mentally, and do not estimate an effect
    from "chemical intuition" about the substitution. **ALWAYS** use the ranked
    output the script already produced: `most_deleterious_substitutions`,
    `best_tolerated_substitutions`, `most_constrained_positions` and
    `most_tolerant_positions` in `summary.json`, or the ranked table from
    `score-variants`. The script has already done the arithmetic correctly; any
    number you derive by hand is a number you can get wrong.
-   **Variants are 1-INDEXED**, the standard biology convention. `K48R` means
    the 48th residue, which must be a lysine. `score-variants` **validates the
    wild-type residue at every position and fails loudly if it does not match**.
    If it fails, do NOT "fix" it by shifting the index — stop and work out
    whether you have the wrong isoform, a construct with an uncleaved signal
    peptide or an affinity tag, or a 0-indexed list.
-   **NEVER report LLR as an energy.** It is unitless. Use it to **rank**, not
    to quantify. Never write "ΔΔG = -9.4 kcal/mol".
-   **ALWAYS `uv run --no-project`.** Never bare `python`/`python3`, never
    `pip install`. Without `--no-project`, `uv` walks up the directory tree,
    finds an unrelated `pyproject.toml` and tries to build that project instead.
-   **NEVER download model weights.** Everything runs remotely on the Biohub
    Platform. Do not install `torch`, `transformers` or the `esm` package.
-   **Budget the API calls.** A `scan` costs **one request per residue** (L+1
    total). A 76-residue protein is 77 requests (~30 s); a 500-residue protein
    is 501. For a handful of known variants use `score-variants`, which only
    masks the positions you asked about.
-   **Before interpreting results, read [`docs/interpretation-guide.md`](docs/interpretation-guide.md).**
    It covers the sign convention, the magnitude and entropy bands, the
    ranking-not-ΔΔG caveat, and the measured anchors you calibrate against.
-   **Review at least one worked example in [`docs/examples/`](docs/examples/)
    before writing a report** — including the negative
    [`scramble_control`](docs/examples/scramble_control/report.md) so you know
    what a "not protein-like" result looks like.
-   **Write the report using [`docs/report-templates.md`](docs/report-templates.md).**
    Every number in it must trace to `summary.json` or the `score-variants` JSON.
-   **If this skill is used, ensure this is mentioned in the output.**

## Utility Scripts

All commands are run from the skill directory (or with an absolute path to the
script). All output paths should be absolute, or relative to the user's project
— never relative to the skill directory.

### 1. `scan` — full deep mutational scan

Every residue mutated to every amino acid. Costs L+1 API requests.

```bash
uv run --no-project scripts/mutation_scoring.py scan \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --model esmc-600m-2024-12 \
  --output-dir ./ubiquitin_scan \
  --top 10
```

Writes into `--output-dir`:

| File | Contents |
|---|---|
| `llr.npy` | `(L, 20)` float32 LLR matrix. Columns are `ACDEFGHIKLMNPQRSTVWY`. |
| `llr.csv` | The same matrix in long form: `position, wt, mutant, variant, llr`. |
| `positions.csv` | Per position: entropy, fraction deleterious, mean LLR, best/worst substitution. |
| `summary.json` | **Read this one.** Pseudo-perplexity, WT recovery, entropy profile, and the `--top` N ranked positions and substitutions. |

`--fasta` may be used instead of `--sequence`.

### 2. `score-variants` — score specific variants

Only masks the positions you ask about, so it is cheap. Variants are
**1-indexed** and the stated wild-type residue **must** match the sequence.

```bash
uv run --no-project scripts/mutation_scoring.py score-variants \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --variants I44A,K48R,L8A,G76A \
  --output ./variants.json
```

Or from a file (one per line, `#` comments allowed):

```bash
uv run --no-project scripts/mutation_scoring.py score-variants \
  --fasta protein.fasta --variants-file candidates.txt --output ./variants.json
```

Prints a table ranked most- to least-deleterious and writes the same, plus each
variant's `rank_at_position` (where it sits among the 19 alternatives at that
residue) and the best substitution available at that position.

### 3. `pseudo-perplexity` — sequence fitness

"How protein-like is this?" Useful for triaging designs, checking a construct,
or comparing a variant sequence against wild type.

```bash
uv run --no-project scripts/mutation_scoring.py pseudo-perplexity \
  --fasta designs.fasta --output ./fitness.json
```

With `--fasta` it scores every record and ranks them best (lowest) to worst.

### 4. `heatmap` — the LLR matrix as a picture

Reuse a finished scan (**no API calls**), or pass `--sequence` to compute a
fresh one.

```bash
uv run --no-project scripts/mutation_scoring.py heatmap \
  --scan-dir ./ubiquitin_scan --output ./heatmap.png
```

Diverging colormap centred at 0: **red = deleterious, white = neutral, blue =
tolerated**; a dot marks the wild-type residue. Rows are ordered by descending
isoelectric point. Add `--vmax 6` to boost contrast on a heavily constrained
protein, where the default full-range scale can look uniformly red.

### 5. `entropy` — per-position constraint profile

```bash
uv run --no-project scripts/mutation_scoring.py entropy \
  --scan-dir ./ubiquitin_scan --output ./entropy.png
```

Writes the PNG plus a CSV alongside it. Dashed reference lines show the entropy
of a uniform choice over 2, 4, 8 and 16 amino acids, so you can see how much the
model has actually narrowed the field.

## Interpreting the Output

### LLR — is this mutation bad?

**LLR < 0 means predicted deleterious. More negative = worse. LLR of the
wild-type residue is exactly 0 by construction.**

As a starting heuristic:

| LLR | Reading |
|---|---|
| `> 0` | Model prefers the mutant over the wild type. Rare in a well-conserved protein; worth flagging. |
| `-1 to 0` | Effectively neutral — the model is indifferent. |
| `-3 to -1` | Mildly disfavoured. |
| `-6 to -3` | Clearly deleterious. |
| `< -6` | Strongly deleterious. |
| `< -10` | Severe; typically a buried-core or catalytic residue. |

**These bands shift with the protein.** They are a starting point, not a
verdict. LLR magnitudes are **not comparable across proteins**: a hyper-conserved
protein pushes every score down. In the ubiquitin scan above, the *mean* LLR
over all 1,444 substitutions is **-9.0** — a score of -6 is *better than
average* there, while in a tolerant protein -6 would be alarming. **Always
calibrate within the protein you are scoring**, using the ranked lists in
`summary.json` or `rank_at_position` from `score-variants`.

### Entropy — is this position constrained?

**Low entropy = evolutionarily constrained. High entropy = mutation-tolerant.**
Measured over the full 64-wide output distribution, so the ceiling is
`log2(20) = 4.32 bits` (total ignorance).

| Entropy (bits) | Reading |
|---|---|
| `< 0.1` | Essentially invariant. The model is certain. |
| `0.1 – 0.5` | Highly constrained. |
| `0.5 – 1.5` | Moderately tolerant. |
| `1.5 – 3.0` | Tolerant — a candidate site for library design. |
| `> 3.0` | Unconstrained; the model has no idea. On a *natural* protein this usually means a disordered or low-complexity region. |

Reference points measured on this API (`esmc-600m-2024-12`): ubiquitin has mean
entropy **0.23 bits** (min 0.007 at the initiator Met, max 1.70). A *scrambled*
ubiquitin sits at **4.12 bits** — near the 4.32 ceiling, i.e. the model correctly
signals it has no idea what belongs anywhere in a sequence that is not a protein.

### Pseudo-perplexity — is this a real protein?

`exp(-mean log p(residue | leave-one-out context))`. Ranges from 1 (perfect
prediction) to ~20 (uniform guessing over 20 amino acids).

| Pseudo-perplexity | Reading |
|---|---|
| `1 – 2` | Highly protein-like; a natural, well-conserved sequence. |
| `2 – 5` | Plausible protein. |
| `5 – 10` | Questionable — a poor design, a shuffled or chimeric sequence, or a fragment out of context. |
| `> 10` | Not protein-like. Approaching the random ceiling. |

Measured: **ubiquitin = 1.05** (and ESM C recovers **100%** of its residues under
leave-one-out masking). Its own **scramble = 18.70** — an **17.8×** separation on
identical amino-acid composition.

### Fraction deleterious

The fraction of the **19 non-wild-type substitutions** with LLR < 0, so 1.0 means
*every* substitution at that position is predicted harmful. (Note this differs
from the ESM tutorial, which divides by 20 and therefore caps at 0.95.)

**This metric saturates on conserved proteins and then tells you nothing.** In
the ubiquitin scan it is exactly **1.0 at all 76 positions** — a true result
(ubiquitin is one of the most conserved eukaryotic proteins) but a useless
discriminator. **When it saturates, fall back on entropy and on the LLR
ranking**, which still separate positions cleanly.

### What to actually report

1.  The ranked variants, most to least deleterious, with their LLRs.
2.  Whether each sits at a constrained (low-entropy) or tolerant position.
3.  The calibration: where the score falls relative to the rest of *this*
    protein.
4.  The caveat below, every time.

## Common Mistakes

-   **Treating LLR as ΔΔG.** It is a unitless log-probability ratio, not an
    energy, and it is not calibrated to any experimental scale. ESM zero-shot
    scores are **best used for RANKING** — "which of these ten variants is
    worst" — and are genuinely good at it. They are **not** absolute stability
    predictions. Never quote them in kcal/mol, and never claim a variant is
    "3 kcal/mol destabilising" from an LLR.

-   **Confusing evolutionary likelihood with functional essentiality.** These
    come apart, and the failure is not hypothetical. In the ubiquitin scan,
    **G76** — the C-terminal glycine that forms the isopeptide bond to substrate
    lysines, without which ubiquitin cannot be conjugated to anything and is
    functionally dead — is scored as the **second most tolerant position in the
    protein** (entropy 1.52 bits), and `G76C` comes out as the single best
    tolerated substitution anywhere in ubiquitin (LLR -0.14, essentially
    neutral). The model is not wrong about its own objective: ubiquitin occurs
    throughout the training data as polyubiquitin and as ubiquitin-fusion
    precursors, where position 76 is followed by *more residues*, so the
    likelihood there is genuinely diffuse. But the model is answering "what
    residue is likely here?", not "what residue does this protein need to work?"
    **Always sanity-check a surprising tolerant call against known biology, and
    treat the final and first few residues with particular suspicion — they have
    one-sided context.**

-   **Off-by-one on the position.** Variants are 1-indexed. The script validates
    the wild-type residue and refuses to score anything if it mismatches — it
    will even tell you when the residue you named appears at `position + 1`,
    which is the signature of a 0-indexing error. Do not work around this by
    shifting indices until it passes; find out why the numbering disagrees.

-   **Running `scan` on a large protein without thinking.** It is O(L) API
    requests. If the user only cares about five variants, use `score-variants`
    (six requests) instead of a 500-request scan.

-   **Re-running the scan to make a plot.** `heatmap` and `entropy` accept
    `--scan-dir` and then make **zero** API calls. Only pass `--sequence` to
    them if you have not already scanned.

## Dependencies

-   `uv` — environment management (required).
-   `credentials` — safe handling of `BIOHUB_API_KEY` (required).

## Related Skills

| If the user wants… | Use |
|---|---|
| To generate or design new sequences | `esm3_protein_design` |
| A predicted 3D structure, pLDDT, pTM | `esmfold2_structure_prediction` |
| To screen candidate binders against a target | `esmfold2_binder_screening` |
| Embeddings/representations for a downstream model | `esmc_protein_embeddings` |
| Which layer's embeddings to use | `esmc_embedding_layer_sweep` |
| A mechanistic read on what a region does | `esmc_sae_feature_interpretation` |
| A sequence back from a structure | `esm3_inverse_folding` |
| Functional annotation of a sequence | `esm3_function_prediction` |

## References

-   [Interpretation guide](docs/interpretation-guide.md) — how to read LLR,
    entropy, and pseudo-perplexity; the sign convention; magnitude calibration;
    the Negative Results & Scientific Integrity section; and the pre-report
    reasoning checklist.
-   [Report templates](docs/report-templates.md) — numbered scaffolds for a
    variant report and a deep-mutational-scan report.
-   [Worked example: ubiquitin DMS](docs/examples/ubiquitin_dms/report.md)
    (positive) — a cleanly constrained protein; pseudo-perplexity 1.05, 100%
    wild-type recovery, and the G76 likelihood-≠-essentiality trap.
-   [Worked example: scramble control](docs/examples/scramble_control/report.md)
    (negative) — what "not protein-like" looks like: pseudo-perplexity 18.70 vs
    ubiquitin's 1.05 on identical composition.
-   [ESM/Biohub API reference](references/esm-biohub-api.md) — endpoints,
    reachable models, request/response shapes, the credit model, and gotchas
    shared by every ESM skill.
-   `references/citation.bib` — cite Meier et al. (2021) for the zero-shot
    masked-marginal method, and Salazar et al. (2020) for pseudo-perplexity.
