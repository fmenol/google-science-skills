# Report Templates

Templates for reporting an ESM3 design campaign. Follow the structure so every
report carries the QC that makes a design a design: the fold metrics, the
comparison to the floor, and (for scaffolds) the motif RMSD.

> [!IMPORTANT] **DO NOT** use the agent's default artifact directory for this
> report. Save it as `report.md` directly in the campaign's output directory
> (next to the `.fasta`, `.json`, and PDBs the skill wrote), and reference
> figures with relative paths (`filename.png`).

> [!CRITICAL] Read `docs/interpretation-guide.md` first. Never fill in a pTM,
> pLDDT, or motif RMSD from anywhere but the skill's metrics JSON. If a value is
> not in the JSON, it does not go in the report.

--------------------------------------------------------------------------------

## 1. De Novo / Conditioned Design Campaign

Use for `generate`, `design-with-ss`, `design-with-sasa`, `chain-of-thought`.

```markdown
# Design Report: {objective}

## 1. Summary of Findings
[What was designed, and the verdict in one paragraph. State the objective
(e.g. "a 60-residue de novo mini-protein"), the subcommand and conditioning
used, how many samples were drawn, the top design's pTM/pLDDT, and — for a de
novo run — how that compares to the random-string floor. Be honest about
whether this cleared the bar or is a rejected campaign.]

## 2. Design Setup
- **Objective**: {what the user asked for}
- **Subcommand / conditioning**: {generate | design-with-ss | ... ; the prompt}
- **Length**: {L}
- **Samples drawn**: {N}, temperature {T}
- **Models**: esm3-open-2024-03 (design) / esmfold2-fast-2026-05 (QC fold)

## 3. Ranked Designs (QC Metrics)
[The ranked table from the metrics JSON. Report the whole set, not just rank 1.]

| Rank | Design | pTM | mean pLDDT | min pLDDT | Notes |
| :--- | :----- | :-- | :--------- | :-------- | :---- |
| 1 | design_1 | {ptm} | {plddt} | {plddt_min} | [low-complexity? clean?] |
| ... | | | | | |

**Floor comparison** (de novo only): mean design pTM {x} vs random-control floor
{y} (n={N} each); mean design pLDDT {x} vs {y}. [Cleared / did not clear.]

## 4. Top Design
- **Sequence**: `{sequence}`
- **pTM / pLDDT**: {ptm} / {plddt}
- **Low-complexity / artifact check**: [inspected; clean, or flag the run/segment]
- **Figure**: ![QC](figure.png)
  *Caption: what the figure shows, with the real numbers.*

## 5. Interpretation
[Read the metrics against the calibration table and the floor. Is the pTM what a
prompt of this length and conditioning should give? Is the top design usable, or
low-complexity? If the campaign underperformed, say why (short, unconditioned)
and what to change (condition it / more samples / longer).]

## 6. Limitations
[These are computational designs, not validated proteins: no expression,
stability, solubility, or function is established. State length-sensitivity of
pTM and any artifact flags. Cross-reference the right sibling skill if the user's
real need is folding / inverse folding / function.]

## 7. Conclusion
[Verdict and concrete next step. Cite hayes2024simulating (ESM3).]
```

--------------------------------------------------------------------------------

## 2. Motif-Scaffolding Campaign

Use for `scaffold-motif`. The extra, non-negotiable content is the motif QC:
verbatim + CA RMSD.

```markdown
# Scaffold Report: {motif} into a {L}-residue protein

## 1. Summary of Findings
[What motif was grafted, from what source structure, into what scaffold, and
whether the graft succeeded. State the top design's pTM/pLDDT AND the motif CA
RMSD in the first paragraph — both certify the result.]

## 2. Design Setup
- **Motif**: {motif_sequence} ({n} residues), from {PDB} chain {c}, author
  residues {motif_range}
- **Scaffold**: {scaffold_length} residues, motif placed at {motif_start}-{end}
- **Samples drawn**: {N}, temperature {T}
- **Models**: esm3-open-2024-03 / esmfold2-fast-2026-05

## 3. Ranked Designs (QC Metrics)

| Rank | Design | pTM | mean pLDDT | Motif verbatim | Motif CA RMSD (Å) |
| :--- | :----- | :-- | :--------- | :------------- | :---------------- |
| 1 | design_1 | {ptm} | {plddt} | {true/false} | {rmsd} |
| ... | | | | | |

## 4. Motif Preservation (the two independent checks)
- [ ] **Sequence**: motif appears verbatim at the requested offset
  (`motif_verbatim` = {true/false}; grafted = `{substring}`).
- [ ] **Geometry**: motif CA RMSD = {rmsd} Å → [excellent < 1 / good 1-2 /
  FAILED > 3]. A verbatim motif with high RMSD is a failed graft.
- **Figure**: ![Per-residue pLDDT](scaffold_plddt.png)
  *Caption: per-residue confidence with the motif region shaded; note the fold
  quality across the whole design and within the motif.*

## 5. Interpretation
[Did the whole protein fold confidently (pTM/pLDDT) AND did the active-site
geometry survive (RMSD)? Both are required. Conditioning on a rigid motif
typically yields high pTM — say so, and contrast with the unconditioned floor.]

## 6. Limitations
[Geometry preserved != function validated. The design is computational; catalytic
activity, expression, and stability are untested. Function prediction is
`esm3-function-prediction`.]

## 7. Conclusion
[Verdict and next step. Cite hayes2024simulating (ESM3); cite kabsch1976 for the
RMSD superposition if the method is described.]
```

--------------------------------------------------------------------------------
