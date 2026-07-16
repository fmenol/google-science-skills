# Report Templates

Scaffolds for reporting ESM C zero-shot mutation-effect results. Fill every
bracketed field with a **real number from the run** — `summary.json` for a scan,
the variants JSON for `score-variants`. Never write a value you reconstructed by
hand. Read `interpretation-guide.md` before you start.

> [!IMPORTANT] **DO NOT** use the agent's default artifact directory for the
> report. Save it as `report.md` directly in the run's output directory (next to
> `summary.json` / the heatmap), and reference figures with **relative** paths
> (`heatmap.png`), so the folder is self-contained and portable.

There are two templates. Use **A** when the user named specific variants ("is
K48R bad?", "score these ten"). Use **B** when you ran a full deep mutational
scan ("which positions can I safely mutate?", "make me an LLR heatmap"). For a
negative or not-protein-like result, keep the same skeleton but let the finding
be the flatness — see `docs/examples/scramble_control/report.md`.

--------------------------------------------------------------------------------

## Template A — Specific-variant report (`score-variants`)

```markdown
# Mutation Scoring Report: {protein_name}

## 1. Summary of Findings
[One paragraph. Which variants are predicted deleterious, which tolerated, and
the single most important call. State that these are ESM C zero-shot
(evolutionary-likelihood) rankings, not energies or clinical calls. Lead with the
biology, not the raw numbers.]

## 2. Method & Context
- **Protein**: {name} ({accession}, {length} aa)
- **Model**: {esmc_model}  (leave-one-out masked-marginal scoring)
- **Variants scored**: {list}  (1-indexed; wild-type residue validated)
- **Calibration basis**: rank within the 19 alternatives at each position
  (`rank_at_position`); absolute LLR magnitudes are not comparable across
  positions or proteins.

## 3. Ranked Variants (most to least deleterious)
[Straight from the variants JSON, already ranked. Pair each LLR with the position
entropy so the reader sees whether the position is constrained.]

| Variant | LLR | Position entropy (bits) | rank@pos | Prediction |
| :------ | --: | ----------------------: | :------- | :--------- |
| {V}     | {x} | {e}                     | {r}/19   | {deleterious/tolerated} |

## 4. Interpretation
[For each variant of interest: read LLR + entropy together (see the Signal
Patterns table). Is it deleterious at a constrained position (high confidence) or
at a tolerant one (check alternatives)? Cite the best substitution available at
that position where relevant. Flag any positive LLR.]

## 5. Caveats
[The unavoidable ones: unitless (no ΔΔG / kcal/mol), sequence-only (no structure,
MSA, binding partner, ligand, or PTM), single substitutions only, evolutionary
likelihood != clinical pathogenicity. Note any residue whose function is not
captured by likelihood — e.g. a catalytic or interface residue.]

## 6. Conclusion
[Direct answer to the user's question, with the ranking and the confidence.
Note that the ESM C mutation-scoring skill was used.]
```

--------------------------------------------------------------------------------

## Template B — Deep mutational scan report (`scan`)

```markdown
# Deep Mutational Scan Report: {protein_name}

## 1. Summary of Findings
[One paragraph. How constrained is the protein overall (pseudo-perplexity, mean
entropy)? Where are the hard constraints and where is it permissive? State the
model and that these are zero-shot evolutionary likelihoods for ranking.]

## 2. Method & Context
- **Protein**: {name} ({accession}, {length} aa)
- **Model**: {esmc_model}  (leave-one-out masked-marginal, {length}+1 API calls)
- **Whole-sequence fitness**: pseudo-perplexity {pppl}; wild-type recovery
  {wt_recovery} — [a natural, well-recognised protein / a questionable sequence].
- **Constraint**: mean entropy {mean} bits (range {min}–{max}); mean LLR over the
  {19*length} substitutions {mean_llr}.

## 3. Most Constrained Positions (do not mutate)
[From `most_constrained_positions` in summary.json — lowest entropy. Annotate
with known biology where you can: buried core, catalytic, conserved motif.]

| Position | WT | Entropy (bits) | Note |
| :------- | :- | -------------: | :--- |

## 4. Most Tolerant Positions (candidates for mutation / library design)
[From `most_tolerant_positions` — highest entropy. Flag any that are functionally
critical despite high tolerance (likelihood != essentiality).]

| Position | WT | Entropy (bits) | Note |
| :------- | :- | -------------: | :--- |

## 5. Most / Least Damaging Substitutions
[From `most_deleterious_substitutions` and `best_tolerated_substitutions`. These
are the extremes of the whole (L, 20) matrix, already ranked.]

## 6. Figures
![LLR heatmap](heatmap.png)
*Fig 1: (L, 20) LLR matrix. Red = deleterious, blue = tolerated, white = neutral;
dot = wild type. Note the symmetric colour scale (±{vmax} LLR) — a saturated-red
matrix signals a heavily constrained protein.*

![Entropy profile](entropy.png)
*Fig 2: Per-position entropy. Low bars = evolutionarily constrained; dashed lines
mark uniform choice over 2/4/8/16 amino acids.*

## 7. Caveats
[Unitless (no kcal/mol), sequence-only, single substitutions, likelihood !=
pathogenicity, one-sided context at the termini. Note that fraction-deleterious
saturates at 1.0 on conserved proteins and is uninformative there — entropy and
the LLR ranking do the discriminating.]

## 8. Conclusion
[Which positions are safe to mutate, which are off-limits, and the confidence.
Answer the user's original question. Note that the ESM C mutation-scoring skill
was used.]
```

--------------------------------------------------------------------------------

## Filling the templates honestly

-   **Every number traces to a file.** Scan numbers come from `summary.json`;
    variant numbers from the `score-variants` output JSON. If a number is not in
    one of those, do not put it in the report.
-   **Rank, don't quantify.** Present variants in the order the script ranked
    them. Do not translate an LLR into a fold-change, a ΔΔG, or a percentage.
-   **Pair LLR with entropy.** A report that cites LLR without the position's
    entropy has thrown away half the signal.
-   **Let a negative result be negative.** If the sequence is not protein-like or
    the model has no preference, that *is* the finding. Do not manufacture a
    mechanism to fill section 4/5.
