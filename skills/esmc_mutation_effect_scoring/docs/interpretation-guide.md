# Interpretation Guide

This guide covers how to read the three numbers ESM C zero-shot scoring produces
— **LLR**, **entropy**, and **pseudo-perplexity** — how to calibrate them within
a protein, the mechanisms they can and cannot see, and the pre-report reasoning
checklist. Read this before writing any report. It is the deeper companion to the
*Interpreting the Output* section of `SKILL.md`; where the two overlap they use
the same measured anchors.

All reference numbers below were measured on `esmc-600m-2024-12`, the default
ESMC model, and are reproducible from the worked examples in `docs/examples/`.

--------------------------------------------------------------------------------

## The mental model (read this first)

Every score comes from one operation repeated per residue: hide position *i*,
show the model `sequence[:i] + '_' + sequence[i+1:]`, and read the probability
distribution it predicts over the 20 amino acids at *i* — **without the model
ever seeing the wild-type residue there** (leave-one-out masking). From that one
distribution:

-   **LLR** (per substitution) = `log p(mutant) - log p(wild-type)`. How much
    *less* likely the mutant is than the wild type at this position.
-   **Entropy** (per position) = how spread the distribution is. How unsure the
    model is about the position *at all*.
-   **Pseudo-perplexity** (per sequence) = `exp(-mean log p(actual residue))`.
    How predictable the whole sequence is — "is this a real protein?"

> [!IMPORTANT] These are **evolutionary likelihoods, not energies.** The model
> answers "what residue is likely here?", which correlates with — but is not —
> "what does this protein need to work?" Use the scores to **rank** variants, the
> one thing they are genuinely good at. Never convert an LLR to kcal/mol.

--------------------------------------------------------------------------------

## Signal Patterns Reference

Read a variant from the pair (its **LLR**, its **position entropy**). The two
together identify the mechanism; neither alone is enough.

| LLR | Position entropy | Reading | What to report |
| :-- | :-- | :-- | :-- |
| Strongly negative (`< -6`) | Low (`< 0.5`) | Substitution at a **constrained** residue — buried core, catalytic, or a tightly conserved motif. | Predicted deleterious, high confidence. Name the constraint if known. |
| Strongly negative (`< -6`) | High (`> 1.5`) | Rare and worth a second look. The position tolerates *variety* but not *this* change. | Deleterious, but check the ranked alternatives — a different substitution may be safe. |
| Mildly negative (`-3 to -1`) | Moderate (`0.5–1.5`) | The model is mildly disfavouring, but indifferent overall. | Likely tolerated; a subtle effect at most. |
| Near zero (`-1 to 0`) | High (`> 1.5`) | **Tolerant surface / loop position.** The model has no strong preference. | Safe candidate for mutation or library design. |
| **Positive** (`> 0`) | Any | The model prefers the *mutant*. On a conserved protein this is rare and suspicious — often a **positional prior** (e.g. N-terminal Met) rather than a fitness gain. | Flag it. Sanity-check against biology before claiming a beneficial mutation. |

> [!CAUTION] **Magnitude is only comparable within one protein.** A hyper-
> conserved protein pushes every LLR down; a tolerant one lifts them all. `-6` is
> *better than average* in ubiquitin (mean LLR **-9.0**) but alarming in a
> permissive protein. Always calibrate with the ranked lists in `summary.json`
> or `rank_at_position` from `score-variants`, never against a fixed cutoff.

--------------------------------------------------------------------------------

## LLR — is this mutation bad?

**LLR < 0 means predicted deleterious. More negative = worse. The LLR of the
wild-type residue is exactly 0 by construction** (it is `log p(wt) - log p(wt)`),
which is also what proves each matrix row is aligned to its residue.

Starting heuristic (shifts with the protein — see the caution above):

| LLR | Reading |
| :-- | :-- |
| `> 0` | Model prefers the mutant over the wild type. Rare in a well-conserved protein; worth flagging. |
| `-1 to 0` | Effectively neutral — the model is indifferent. |
| `-3 to -1` | Mildly disfavoured. |
| `-6 to -3` | Clearly deleterious. |
| `< -6` | Strongly deleterious. |
| `< -10` | Severe; typically a buried-core or catalytic residue. |

**Calibrate inside the protein.** In the ubiquitin scan the mean LLR over all
1,444 substitutions is **-9.0**, and the most damaging single change is
**I30W at -15.10**. Against that backdrop `K48R` at **-5.29** is the *least* bad
substitution available at position 48 (`rank_at_position` 1 of 19) yet is still
clearly deleterious — a conservative lysine→arginine swap the model dislikes
because K48 sits in a conserved patch. Report the **rank**, not the raw number in
isolation.

--------------------------------------------------------------------------------

## Entropy — is this position constrained?

**Low entropy = evolutionarily constrained. High entropy = mutation-tolerant.**
Measured over the full 64-wide output distribution, so the ceiling is
`log2(20) = 4.32 bits` (total ignorance over the 20 amino acids).

| Entropy (bits) | Reading |
| :-- | :-- |
| `< 0.1` | Essentially invariant. The model is certain. |
| `0.1 – 0.5` | Highly constrained. |
| `0.5 – 1.5` | Moderately tolerant. |
| `1.5 – 3.0` | Tolerant — a candidate site for library design. |
| `> 3.0` | Unconstrained; the model has no idea. On a *natural* protein this usually means a disordered or low-complexity region. |

**Anchors (`esmc-600m-2024-12`).** Ubiquitin has mean entropy **0.23 bits** — a
tightly constrained protein — ranging from **0.007 at the initiator Met (M1)** to
**1.70 at P19** in a flexible loop. A *scrambled* ubiquitin, by contrast, sits at
mean **4.12 bits**, hard against the 4.32 ceiling: the model correctly signals it
has no idea what belongs anywhere in a sequence that is not a protein. That gap —
0.23 vs 4.12 — is the entropy signature of "folded protein" vs "random string of
the same amino acids".

--------------------------------------------------------------------------------

## Pseudo-perplexity — is this a real protein?

`exp(-mean log p(residue | leave-one-out context))`, one number for the whole
sequence. Ranges from 1 (the model predicts every residue perfectly) to ~20
(uniform guessing over 20 amino acids).

| Pseudo-perplexity | Reading |
| :-- | :-- |
| `1 – 2` | Highly protein-like; a natural, well-conserved sequence. |
| `2 – 5` | Plausible protein. |
| `5 – 10` | Questionable — a poor design, a shuffled or chimeric sequence, or a fragment out of context. |
| `> 10` | Not protein-like. Approaching the random ceiling. |

**Anchors.** Ubiquitin scores **1.05** (and ESM C recovers **100%** of its
residues as the single most likely amino acid under masking). Its own **scramble
— identical amino-acid composition, fold destroyed — scores 18.70**, a **17.8×**
separation, with wild-type recovery collapsing to **6.6%**. Pseudo-perplexity is
the right tool for triaging a design, checking a construct for frame/assembly
errors, or asking whether a putative sequence is protein-like at all. It is a
whole-sequence fitness proxy, not a per-position score — do not use it to locate
*which* residue is the problem (use the entropy and LLR profiles for that).

--------------------------------------------------------------------------------

## Fraction deleterious — and its saturation trap

`positions.csv` and `summary.json` report, per position, the fraction of the
**19 non-wild-type substitutions** with LLR < 0. It reads as "how completely is
this position pinned down".

> [!WARNING] **This metric saturates on conserved proteins and then tells you
> nothing.** In the ubiquitin scan it is exactly **1.0 at all 76 positions** — a
> true result (every substitution at every position is disfavoured) but a useless
> discriminator: the tolerant loop P19 and the invariant M1 both read 1.0. When
> it saturates, **fall back on entropy and the LLR ranking**, which still
> separate positions cleanly (P19 entropy 1.70 vs M1 entropy 0.007). On the
> scramble, fraction-deleterious is *not* saturated (mean 0.60, only 5 of 76
> positions at 1.0) — a reminder that the metric only pins to 1.0 when a genuine
> fold constrains every alternative.

--------------------------------------------------------------------------------

## Calibration: how to compare scores fairly

Raw LLR magnitudes are not comparable across positions, let alone across
proteins. Two tools give you an honest comparison:

-   **`rank_at_position`** (from `score-variants`): where a substitution sits
    among the 19 alternatives at its own residue. Rank 1 = the least-bad change
    available there. This is the only calibration a single variant admits without
    a full scan.
-   **The ranked lists in `summary.json`** (from `scan`):
    `most_deleterious_substitutions`, `best_tolerated_substitutions`,
    `most_constrained_positions`, `most_tolerant_positions`. These have already
    done the arithmetic over the whole matrix. **Use them. Never re-rank the
    matrix by hand.**

--------------------------------------------------------------------------------

## Negative Results & Scientific Integrity

> [!CRITICAL] **"This sequence is not protein-like" and "the model has no
> preference here" are real, reportable findings.** Do not stretch a flat or
> noisy result into a mechanism.

-   **Value of the negative control.** A high pseudo-perplexity (say `> 10`) with
    entropy near the 4.32-bit ceiling is a clean, informative result: the model
    is telling you the sequence does not look like a protein. The scramble
    example (18.70, 4.12 bits) is what that looks like. Report it plainly; do not
    hunt for "constrained" positions in noise.
-   **A positive LLR is a flag, not a discovery.** The model preferring the
    mutant usually reflects a **positional prior**, not a fitness gain. The single
    strongest signal anywhere in the scrambled ubiquitin is `R1M` at **+3.63** —
    the model "wants" a methionine at position 1 because nearly every protein
    starts with the initiator Met, a fact independent of any fold. Never report a
    positive LLR as a beneficial mutation without checking biology first.
-   **Likelihood is not functional essentiality — and they come apart.** In the
    ubiquitin scan, **G76** — the C-terminal glycine that forms the isopeptide
    bond to substrate lysines, without which ubiquitin is functionally dead — is
    the **second most tolerant position in the protein** (entropy 1.52) and
    `G76C` is the single **best-tolerated substitution anywhere** (LLR -0.14). The
    model is not wrong about its own objective: ubiquitin appears throughout the
    training data as polyubiquitin and as ubiquitin-fusion precursors, where
    position 76 is followed by *more residues*, so the likelihood there is
    genuinely diffuse. But it is answering "what residue is likely here?", not
    "what does this protein need?" **Always sanity-check a surprising tolerant
    call against known biology, and treat the first and last few residues with
    particular suspicion — they have one-sided context.**
-   **Strict anti-speculation.** Do not invent a mechanism ("destabilises the
    core", "breaks a salt bridge") that the score alone cannot support. The score
    is a sequence likelihood; it does not localise a structural cause. If you want
    a structural claim, fold the variant (`esmfold2_structure_prediction`).
-   **Never quote an LLR as an energy.** It is a unitless log-probability ratio,
    not calibrated to any experimental scale. "3 kcal/mol destabilising" derived
    from an LLR is a fabrication.

--------------------------------------------------------------------------------

## Model Limitations

ESM C zero-shot scoring is **sequence-only, unitless, and single-substitution**.
It does NOT model:

-   **ΔΔG in physical units.** Scores are unitless log-likelihood ratios. Good
    for ranking, not for absolute stability. There is no kcal/mol.
-   **Insertions, deletions, frameshifts, or multi-residue variants.** Single
    amino-acid substitutions only. `score-variants` refuses anything else.
-   **Explicit structure, an MSA, or experimental data.** The model infers
    constraint from sequence patterns it learned in pre-training; it never sees a
    fold, an alignment, or a label.
-   **Binding partners, ligands, cofactors, membranes, or PTMs.** A residue that
    is critical only in complex with another molecule (an interface, an active
    site coordinating a metal) may look tolerant in isolation.
-   **Clinical pathogenicity.** The model reports evolutionary plausibility, which
    correlates with — but is not — pathogenicity. A benign-looking LLR does not
    clear a variant; a deleterious-looking one is a hypothesis, not a diagnosis.
-   **Function at the termini.** Positions with one-sided context (the first and
    last few residues) are systematically less reliable — see G76 above.

**Rules:**

-   If the position is a known **catalytic or binding residue** and the model
    calls it tolerant, state the discrepancy and defer to the biology.
-   If the variant is an **indel or affects more than one residue**, say the
    method does not apply — do not approximate it with a nearby substitution.
-   If the question is **stability in kcal/mol or a 3D consequence**, route to
    `esmfold2_structure_prediction`; this skill cannot answer it.

--------------------------------------------------------------------------------

## Pre-Report Reasoning Checklist

Complete this BEFORE writing the report. Ground every claim in a number the
script produced, not in chemical intuition about the substitution.

> [!CAUTION] **Never eyeball the LLR matrix or rank variants in your head.** Use
> `most_deleterious_substitutions` / `best_tolerated_substitutions` /
> `most_constrained_positions` / `most_tolerant_positions` from `summary.json`,
> or the ranked table and `rank_at_position` from `score-variants`. Any number
> you derive by hand is a number you can get wrong.

### 1. Sequence-level sanity

-   [ ] **Check pseudo-perplexity first.** If it is high (`> 10`) with mean
    entropy near 4.3 bits, the sequence is not protein-like — report that and
    stop, rather than interpreting per-position noise.
-   [ ] **Check wild-type recovery.** Near 100% confirms the model recognises the
    sequence and the BOS offset is correct; near chance (~5%) means either a
    non-protein input or a numbering/frame problem.

### 2. Calibrate within the protein

-   [ ] **Read the mean and range of LLR and entropy** from `entropy_summary` and
    the ranked lists. Establish what "bad" and "tolerant" mean *here* before
    judging any single variant.
-   [ ] **Locate each variant of interest in the ranking**, not on the absolute
    scale. Cite its `rank_at_position` and how it sits versus the protein's mean.

### 3. Per-variant reading

-   [ ] **Pair LLR with position entropy** for every variant (see the Signal
    Patterns table). A deleterious LLR at a low-entropy position is high-
    confidence; at a high-entropy position, check the alternatives.
-   [ ] **Flag every positive LLR** and check it against biology before calling
    anything beneficial.
-   [ ] **Sanity-check surprising tolerant calls** at functionally critical or
    terminal residues against known biology.

### 4. Integrity

-   [ ] **No energies.** No kcal/mol, no ΔΔG. Rank, do not quantify.
-   [ ] **No invented mechanism.** Do not assert a structural cause the score
    cannot support; route structural questions to the folding skill.
-   [ ] **State the caveat every time:** ESM C zero-shot scores are evolutionary
    likelihoods, best used for ranking, and are not clinical calls.

### 5. Report readiness

-   [ ] **Embed the figures** the run produced (LLR heatmap, entropy profile) and
    confirm each referenced file exists in the report folder.
-   [ ] **Quote real numbers** from `summary.json` / the variants JSON — never a
    value you reconstructed by hand.
-   [ ] **Answer the user's actual question** (is this variant bad? which
    positions are safe to mutate?), and note that the skill was used.

--------------------------------------------------------------------------------
