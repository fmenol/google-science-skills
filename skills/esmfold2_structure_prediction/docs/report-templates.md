# Report Templates

Scaffolds for reporting an ESMFold2 structure prediction. Fill the sections in
order; every number must come from the `analyze` report, and every figure must
exist in the report folder. Read `interpretation-guide.md` first.

> [!IMPORTANT] **DO NOT** use the agent's default artifact directory for this
> report. Save it as a regular file named `report.md` directly in the
> prediction's output directory (e.g. `ubiquitin_fold/report.md`) in the
> workspace, and reference plots with **relative** paths (`filename.png`).

> [!CRITICAL] **pLDDT from this API is on a 0-1 scale.** State the scale every
> time you quote a value ("0.82 on the 0-1 scale, i.e. 82/100"). The only
> already-rescaled value is the PDB B-factor column (0-100).

--------------------------------------------------------------------------------

## 1. Single-Chain Prediction Template

Use this for one protein folded with `fold`.

```markdown
# Structure Prediction Report: {name} ({length} aa)

## 1. Summary
[One paragraph. What was folded, with which model/settings, and the bottom-line
confidence verdict in plain words: is this a confident structure, a partial one,
or a negative result? Quote mean pLDDT (with scale) and pTM (with verdict).
Example: "ESMFold2 (esmfold2-fast-2026-05, single sequence) predicts a confident,
well-folded structure for X: mean pLDDT 0.82 on the 0-1 scale (82/100), pTM 0.78
(plausible fold)."]

## 2. Input & Run Settings
- **Sequence / source**: {name}, {length} aa [FASTA id or origin if known]
- **Model**: {esmfold2-fast-2026-05 | esmfold2-2026-05}  (MSA used: {yes/no})
- **Settings**: num_loops {N}, num_sampling_steps {N}, include_pae {yes/no}
- **Outputs**: PDB `{file}`, metrics `{file}`

## 3. Global Confidence
| Metric | Value | Verdict |
| :----- | :---- | :------ |
| Mean pLDDT (0-1) | {plddt_mean} ({x100}/100) | {overall verdict} |
| pTM | {ptm} | {confident / plausible / unreliable} |

- **Overall**: [quote the report's overall verdict verbatim and interpret it.]

## 4. Confidence Bands
[Fraction of residues in each band, from the report. Say what the distribution
means — e.g. "93% confident-or-better, one short low-confidence tail".]

| Band | pLDDT (0-1) | Residues | Fraction |
| :--- | :---------- | :------- | :------- |
| Very high | > 0.90 | {n} | {%} |
| Confident | 0.70-0.90 | {n} | {%} |
| Low | 0.50-0.70 | {n} | {%} |
| Disordered / very low | < 0.50 | {n} | {%} |

## 5. Low-Confidence Regions
[Every region from the report as a residue range. Identify the structural
element where you can — mobile terminus, loop, disordered linker. If none, say
so.]

- chain {A} {start}-{end} ({length} aa), mean pLDDT {x} — [interpretation]

## 6. Plots
![pLDDT vs residue]({plddt_plot}.png)
*Fig 1: Per-residue pLDDT with confidence bands. [What the trace shows: where it
is confident, where it dips, and why.]*

![PAE heatmap]({pae_plot}.png)
*Fig 2: Predicted Aligned Error. [One rigid domain or several? Any flexible
tails/hinges? Only include if the fold was run with --include-pae.]*

## 7. Domain / Topology Reading
[From the PAE heatmap: single rigid domain, or multiple domains with flexible
hinges? Give the rigid residue ranges. Reconcile with the pLDDT dips.]

## 8. Downstream Suitability & Limitations
[Which residue ranges are safe for docking / Foldseek / MD, and which
(low-confidence, disordered) must be excluded. State the model blind spots that
apply: confidence != validity/function; no affinity; single static conformation;
single-sequence unless an MSA was used. Note per-residue pLDDT is in the PDB
B-factor column (0-100).]

## 9. Conclusion
[Bottom line: confident structure, partial, or negative result — and what the
user can and cannot conclude from it. State that ESMFold2 (Biohub Platform) was
used.]
```

--------------------------------------------------------------------------------

## 2. Complex Prediction Template

Use this for a `fold-complex` run (proteins + DNA/RNA/ligands). It adds the
interface sections; keep everything from the single-chain template too.

```markdown
# Complex Prediction Report: {chains described}

## 1. Summary
[What complex was folded, and the bottom line: does the model predict a confident
interface? Quote iPTM with its verdict — that is the number that answers "do
these bind?" — and immediately state the affinity caveat: iPTM is confidence, not
a Kd. Example: "ESMFold2 predicts a confident interface between barnase and
barstar: iPTM 0.96 (confident interface). This reflects the model's confidence in
the geometry, not a measured affinity."]

## 2. Input & Run Settings
- **Entities**: chain A protein ({n} aa), chain B protein ({n} aa), [DNA/RNA/ligand ...]
- **Model / settings / outputs**: [as in the single-chain template]

## 3. Global & Interface Confidence
| Metric | Value | Verdict |
| :----- | :---- | :------ |
| Mean pLDDT (0-1) | {x} ({x100}/100) | {overall} |
| pTM | {x} | {verdict} |
| **iPTM** | {x} | {confident / possible / likely no binding} |

## 4. Per-Chain Confidence
[Mean pLDDT per chain from the report. Flag any chain the model folded poorly —
one confident chain next to a low-confidence partner must be treated separately.]

| Chain | Type | Length | Mean pLDDT (0-1) |
| :---- | :--- | :----- | :--------------- |
| A | protein | {n} | {x} |
| B | protein | {n} | {x} |

## 5. Interface Reading (PAE)
[The off-diagonal chain block IS the interface confidence. Quote the inter-chain
mean/min PAE from the report and reconcile with iPTM. Dark off-diagonal + high
iPTM = confident docking geometry.]

- PAE {A}-{B}: mean {x} Å, min {x} Å — [interpretation]

## 6. Plots
![pLDDT vs residue]({plddt_plot}.png)
*Fig 1: Per-residue pLDDT, chain boundaries marked.*

![PAE heatmap]({pae_plot}.png)
*Fig 2: PAE. On-diagonal blocks = each chain's internal rigidity; off-diagonal
blocks = inter-chain placement confidence.*

## 7. Interpretation & Limitations
[Does the model predict these associate? State plainly, with the caveats: iPTM is
not affinity and not evidence the complex forms in a cell (interpretation-guide
sections 6-7). Note which regions are reliable for downstream work.]

## 8. Conclusion
[Bottom line on the interface and what it does / does not license the user to
conclude. State that ESMFold2 (Biohub Platform) was used.]
```

--------------------------------------------------------------------------------

## 3. Worked Examples

Two filled-in reports live in `docs/examples/`; review at least one before
writing your own, including the negative one:

- `examples/ubiquitin_fold/report.md` — a confident single-chain fold (positive).
- `examples/scramble_control/report.md` — a scrambled sequence with no confident
  fold (negative): what "no structure" looks like, and why it is a real result.
