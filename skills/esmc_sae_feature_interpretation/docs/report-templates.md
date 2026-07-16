# Report Templates

This file provides the report scaffold for an ESMC SAE feature-interpretation
run. Follow the structure so every report covers the same ground and carries the
hypotheses-not-facts disclaimer.

> [!IMPORTANT] **DO NOT** use the agent's default artifact directory for the
> report. Save it as a regular file named `report.md` in the run's output
> directory (e.g. `results/lysozyme/report.md`) in the workspace. Reference
> plots with **relative** paths (`lysozyme_feature_tracks.png`), and make sure
> every figure you embed actually exists in that folder.

> [!CRITICAL] Every report MUST state, near the top, that SAE feature
> descriptions are **auto-generated hypotheses, not curated annotations**, and
> must use hedged language ("the model associates...") throughout. Read
> [`interpretation-guide.md`](interpretation-guide.md) before filling this in.

--------------------------------------------------------------------------------

## 1. Single-protein feature-interpretation report

Use this for "what does the model see in this protein?" / "which motifs/domains
does ESM detect?" / "what features are active at residue N?".

```markdown
# SAE Feature Interpretation: {protein_name} ({accession}, {L} aa)

> **SAE feature descriptions below are auto-generated hypotheses, not curated
> annotations.** A firing feature means the model *associates* a concept with a
> region; it is not evidence the protein has that function. Verify against
> UniProt/InterPro before treating anything here as fact. Analysis produced with
> the `esmc-sae-feature-interpretation` skill.

## 1. Summary of Findings
[2-4 sentences. What does the model most strongly associate with this protein,
by max activation AND by prevalence? Name the dominant domain (the feature high
on both lists). State the single best trust signal (e.g. the top SwissProt
exemplar of the #1 feature is the protein itself). Hedge every functional claim.]

## 2. Run Context
- **Protein**: {protein_name} ({accession}), {L} residues
- **SAE model**: `{sae_model}`  (ESMC: `{esmc_model}`)
- **Normalization**: {normalize_features}
- **Distinct features fired**: {n} of 16,384 ({pct}% of codebook)

## 3. Top Features by Max Activation (motif-like / local)
[Table. These spike at a few residues: candidate catalytic residues, binding
motifs, or short conserved patterns. Cite the real max activation, prevalence,
1-based peak residue, and the hypothesised label.]

| # | Feature | Hypothesised label | Category | Max act. | Prevalence | Peak residue |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 1 | `{id}` | {label} | {category} | {max} | {prev}/{L} | {peak} |

## 4. Top Features by Prevalence (domain / family-like)
[Table, same columns. Call out any feature that appears HERE but NOT on the
max-activation list -- it is broadly on and easy to miss, and often states what
the protein *is*.]

## 5. Plots and Visual Analysis
[Embed the per-residue activation track(s). Give each panel a caption that says
whether the track is a sharp local spike (motif-like) or a broad plateau
(domain-like), and name the peak residue.]

![Feature tracks]({protein}_feature_tracks.png)
*Fig 1: Per-residue activation. [Which panels are local vs broad; peak residues.
Labels are auto-generated hypotheses.]*

## 6. Trust Assessment (per quoted feature)
[For each feature you lean on: does `top_swissprot` resemble this protein? Is the
observed peak above the feature's own `threshold`? Does the observed track match
the described `activation_pattern`? Flag any feature that fails these.]

## 7. What the model does NOT tell you here
[State the blind spots relevant to the user's actual question: no ground-truth
function (-> esm3-function-prediction), no mutation effect (->
esmc-mutation-effect-scoring), no structure (-> esmfold2-structure-prediction).
Name the one top feature that is generic/compartment-level rather than
protein-specific.]

## 8. Conclusion
[Narrative synthesis in hedged language. What the model associates with this
protein, how trustworthy that association is (via top_swissprot / the negative
control logic), and the one sentence of ground-truth caveat.]
```

--------------------------------------------------------------------------------

## 2. Two-protein comparison report (shared machinery)

Use this for "do these two proteins share internal machinery?" / "how similar
does ESM think these are?".

```markdown
# SAE Feature Comparison: {protein_A} vs {protein_B}

> **Feature labels are auto-generated hypotheses, not curated annotations.**
> Analysis produced with the `esmc-sae-feature-interpretation` skill.

## 1. Summary
[Jaccard value, the calibration-table reading, and the one-line verdict: same
kind of protein / shared domain-or-fold / largely distinct / unrelated.]

## 2. Run Context
- **Sequence A**: {protein_A} ({L_A} aa)
- **Sequence B**: {protein_B} ({L_B} aa)
- **SAE model**: `{sae_model}`;  **top-N**: {top_n}

## 3. Overlap
- **Jaccard(top-{top_n})** = {jaccard}  ({n_shared}/{n_union} features shared)
- **Reading**: [from the calibration table]

## 4. Shared Features (strongest first)
[Table: feature id, hypothesised label, max activation in A and in B. Are the
shared features biology-specific, or generic compartment/terminus signals?]

## 5. Distinctive Features
[What fires in A but not B, and vice versa -- the machinery that separates them.]

## 6. Interpretation & Caveats
[What the overlap means: identity / fold / family similarity. State explicitly
that the feature set is a fingerprint of *identity*, not of mutation effect: a
single substitution can leave Jaccard at 1.00. If the comparison is a point
mutant, direct functional questions to esmc-mutation-effect-scoring.]
```

--------------------------------------------------------------------------------

## Worked examples

Two complete, replay-generated reports show the scaffold in use:

-   **Positive** -- [`examples/lysozyme_features/report.md`](examples/lysozyme_features/report.md):
    lysozyme's top features are cell-wall / peptidoglycan hydrolase features; the
    #1 feature's strongest exemplar is lysozyme itself.
-   **Negative / limitation** --
    [`examples/specificity_control/report.md`](examples/specificity_control/report.md):
    ubiquitin scores 0/10 on lysozyme terms, and Jaccard separates a point mutant
    (1.00) from an unrelated protein (0.02) -- teaching feature specificity and
    the hypotheses-not-facts caveat.
