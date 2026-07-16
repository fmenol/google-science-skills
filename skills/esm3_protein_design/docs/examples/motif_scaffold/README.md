# Example: Motif Scaffolding, Geometry Preserved (POSITIVE)

**Objective**: graft the 1ITU renal-dipeptidase active-site motif (23 residues)
into a new 200-residue protein. **Command**: `scaffold-motif --motif-pdb 1ITU.pdb
--chain A --motif-range 124-146 --scaffold-length 200 --motif-start 73
--num-samples 2`. **Verdict**: **graft succeeded** — the whole protein folds
confidently and the motif geometry is preserved to sub-Ångström RMSD.

## Why this is a good example

This is the flagship workflow and the strongest illustration of "conditioning is
what makes ESM3 good", plus the two-part motif QC:

1.  **Conditioning transforms the result.** The *same* model that barely cleared
    the random floor on an unconditioned 60-mer (pTM 0.26, see `denovo_design/`)
    folds this 200-residue design to **pTM 0.981 / pLDDT 0.965** when given a
    rigid motif to build around.
2.  **Two independent motif checks, both required.** The motif must be (a)
    verbatim in the sequence and (b) geometrically preserved. Both designs kept
    the motif `VEGGHSIDSSLGVLRALYQLGMR` verbatim, and their CA traces superposed
    on the crystal motif at **0.26 Å and 0.28 Å** — the active-site geometry is
    essentially identical to 1ITU.
3.  **RMSD is the real certificate.** A verbatim motif with bad geometry would be
    a failed graft even at high pTM. Here both pass, so the graft is genuine.

## Key Takeaway

For scaffolding, report **both** the fold quality (pTM/pLDDT) **and** the motif CA
RMSD. Sequence identity in the motif is necessary but not sufficient; the
sub-Ångström RMSD is what says the functional site was held in place. Geometry
preserved is still not function validated — that is a separate, wet-lab or
`esm3-function-prediction` question.

## Provenance

`scaffold.fasta` / `scaffold.json` are the actual skill outputs (cassette replay,
zero credits). `scaffold_plddt.png` is rendered from the top design's PDB
(per-residue pLDDT) and JSON. The motif is author residues 124-146 of PDB 1ITU
chain A. Models: `esm3-open-2024-03` (design), `esmfold2-fast-2026-05` (QC fold).
