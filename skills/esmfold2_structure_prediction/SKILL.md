---
name: esmfold2-structure-prediction
description: >
  Predict an all-atom 3D structure from sequence with ESMFold2 — a single
  protein chain, or a molecular complex mixing proteins, DNA, RNA and
  small-molecule ligands (CCD code or SMILES). Use when the user has a
  SEQUENCE (wild type, mutant, designed, or de novo) and wants a structure,
  a PDB file, a confidence assessment (pLDDT / pTM / iPTM), or wants to know
  whether two chains are predicted to form a complex. Words that should
  trigger this skill: "fold this sequence", "predict the structure",
  "model this complex", "will these two proteins bind", "pLDDT", "iPTM".
  Do NOT use when: the user wants an EXISTING structure that somebody already
  solved or predicted — a UniProt ID goes to `alphafold-database-fetch-and-analyze`
  and a PDB ID goes to `pdb_database`; the user wants to find structural
  HOMOLOGS of a structure they already have (use `foldseek-structural-search`);
  or the user wants to rank MANY candidate binders against one target (use
  `esmfold2-binder-screening`, which is built for that batch job).
---

# ESMFold2: Structure Prediction

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If
    `.licenses/esmfold2_structure_prediction_LICENSE.txt` does not already exist
    in the workspace root directory then (1) prominently notify the user to
    check the terms at https://biohub.org/acceptable-use-policy/ and
    https://biohub.ai/, then (2) create the file recording the notification text
    and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for it. Never print
    the key.

## Overview

Runs ESMFold2 remotely on the Biohub Platform to predict an all-atom structure
**from sequence alone**, and writes:

*   a **PDB file** (per-residue pLDDT is in the B-factor column, rescaled 0-100),
*   a **metrics JSON** (per-residue pLDDT on the 0-1 scale, pTM, iPTM, per-chain
    breakdown, PAE saved beside it as `.npy`),
*   and, via `analyze`, a **confidence report** plus a pLDDT-vs-residue plot and
    a PAE heatmap.

It folds single chains **and** complexes: several proteins, DNA, RNA, and
small-molecule ligands given by PDB **CCD code** or **SMILES**.

**This skill does NOT:**

*   look up structures that already exist — for a UniProt ID use
    **`alphafold-database-fetch-and-analyze`**, for a PDB ID use
    **`pdb_database`**. Predicting a structure that has already been solved
    experimentally wastes compute and is less accurate than the real thing.
*   search for structural homologs — use **`foldseek-structural-search`** on a
    coordinate file (which this skill can produce for you first).
*   screen or rank a library of candidate binders — use
    **`esmfold2-binder-screening`**.
*   dock a ligand into a *known* pocket, or predict binding affinity. iPTM says
    how confident the model is in an interface, **not** how tightly it binds.
*   download model weights. Everything runs on the remote API.

## Core Rules

*   **pLDDT FROM THIS API IS ON A 0-1 SCALE, NOT 0-100.** A mean pLDDT of `0.82`
    is a *good* structure (equivalent to 82 in AlphaFold's convention), not a
    catastrophic one. NEVER report a raw value from the metrics JSON as though
    it were on the 0-100 scale, and never multiply by 100 without saying so. The
    metrics JSON carries `"plddt_scale": "0-1"` — read it. The PDB B-factor
    column is the one place the value is already rescaled to 0-100, by
    convention.
*   **Do not compute or eyeball confidence yourself; always use `analyze`.** It
    does the banding, the region-finding and the pTM/iPTM verdicts numerically.
    Do not average pLDDT by hand, do not read the PAE matrix by eye, and do not
    invent thresholds.
*   **ALWAYS `uv run --no-project`.** Without `--no-project`, `uv` walks up the
    directory tree, finds an unrelated `pyproject.toml`, and tries to build that
    project instead.
*   **NEVER pass an MSA to the default model.** `esmfold2-fast-2026-05` was not
    trained with MSAs and will **silently ignore** one — you get a
    single-sequence prediction that looks like an MSA-guided one. Only
    `esmfold2-2026-05` uses an MSA. `fold.py` enforces this for you: `--msa`
    auto-switches the model and refuses if you pin the fast model.
*   **A high pLDDT does not mean the sequence is a real or functional protein.**
    Confidence is not validity. Sanity-check designed sequences against a
    negative control if it matters.
*   **Before interpreting results, read [`docs/interpretation-guide.md`](docs/interpretation-guide.md).**
    It holds the pLDDT (0-1) / pTM / iPTM bands, how to read the PAE heatmap for
    domains and flexibility, when a fold is unreliable, and the negative-result
    rules. Do not judge confidence from memory.
*   **Review at least one worked example in [`docs/examples/`](docs/examples/)
    before writing a report**, including the negative `scramble_control` one, so
    you know what a confident fold *and* a no-confident-fold look like.
*   **Write the report using [`docs/report-templates.md`](docs/report-templates.md).**
*   If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All commands are run from the skill directory. Always pass **absolute** output
paths (or paths relative to the user's project), never paths relative to the
skill directory.

### 1. `fold` — a single protein chain

```bash
uv run --no-project scripts/fold.py fold \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --include-pae \
  --output-pdb /abs/path/ubiquitin.pdb \
  --output-metrics /abs/path/ubiquitin.json
```

`--fasta FILE` (with `--record ID` if the file holds more than one entry) can be
used instead of `--sequence`. Pass `--include-pae` whenever you intend to run
`analyze` — the PAE heatmap needs it.

Speed/quality knobs (defaults are the API's own, `--num-loops 20
--num-sampling-steps 100`):

```bash
# ~2 s for a 76 aa protein; good enough for triage and for scanning variants.
uv run --no-project scripts/fold.py fold --sequence MKT... \
  --num-loops 10 --num-sampling-steps 50 \
  --output-pdb /abs/path/x.pdb --output-metrics /abs/path/x.json
```

To use evolutionary information, pass an A3M alignment. This forces
`--model esmfold2-2026-05`, the only model that reads an MSA:

```bash
uv run --no-project scripts/fold.py fold --sequence MKT... --msa /abs/path/aln.a3m \
  --output-pdb /abs/path/x.pdb --output-metrics /abs/path/x.json
```

### 2. `fold-complex` — proteins + DNA + RNA + ligands

Convenience flags, one per chain. Each is `ID:VALUE`:

```bash
# Barnase + barstar, the textbook high-affinity pair.
uv run --no-project scripts/fold.py fold-complex \
  --protein A:AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIR \
  --protein B:KKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGCDITIILS \
  --include-pae \
  --output-pdb /abs/path/complex.pdb \
  --output-metrics /abs/path/complex.json
```

```bash
# Protein + a small molecule, by CCD code or by SMILES.
uv run --no-project scripts/fold.py fold-complex \
  --protein A:MKT... --ligand-ccd L:SAH \
  --output-pdb /abs/path/x.pdb --output-metrics /abs/path/x.json

uv run --no-project scripts/fold.py fold-complex \
  --protein A:MKT... --ligand-smiles 'L:CC(=O)Oc1ccccc1C(=O)O' \
  --output-pdb /abs/path/x.pdb --output-metrics /abs/path/x.json
```

**A homo-oligomer is ONE entity with several chain ids** — `A,B:SEQ`, not two
`--protein` flags:

```bash
uv run --no-project scripts/fold.py fold-complex --protein A,B:MNKIIIY... \
  --output-pdb /abs/path/dimer.pdb --output-metrics /abs/path/dimer.json
```

For anything richer (per-chain MSAs, non-canonical residues, covalent bonds),
use a JSON or YAML spec:

```yaml
# complex.yaml — RNase H1 homodimer bound to an RNA:DNA hybrid (PDB 4H8K)
sequences:
  - type: protein
    id: [A, B]              # one entity, two chains = homodimer
    sequence: MNKIIIYTDGGARGNPGPAGIGVVITDEKGN...
  - type: rna
    id: C
    sequence: CGACACCUGAUUCC
  - type: dna
    id: D
    sequence: GGAATCAGGTGTCG
```

```bash
uv run --no-project scripts/fold.py fold-complex --spec /abs/path/complex.yaml \
  --include-pae \
  --output-pdb /abs/path/rnaseh.pdb --output-metrics /abs/path/rnaseh.json
```

### 3. `analyze` — turn the numbers into a verdict

**This is the step that interprets confidence. Always run it; never eyeball the
raw arrays.**

```bash
uv run --no-project scripts/fold.py analyze \
  --metrics /abs/path/ubiquitin.json \
  --output-report /abs/path/ubiquitin_report.json \
  --output-plddt-plot /abs/path/ubiquitin_plddt.png \
  --output-pae-plot /abs/path/ubiquitin_pae.png
```

It prints a readable report and writes:

*   `--output-report`: mean/median pLDDT, the fraction of residues in each
    confidence band, low-confidence regions as **residue ranges**, per-chain
    statistics, inter-chain PAE, and the pTM/iPTM **verdict** strings.
*   `--output-plddt-plot`: pLDDT vs residue with the confidence bands shaded and
    chain boundaries marked.
*   `--output-pae-plot`: the PAE heatmap. For a complex, the off-diagonal blocks
    are the inter-chain confidence — dark means the model is sure how the chains
    sit relative to each other.

## Interpreting the Output

Use these thresholds, and take them from the report rather than recomputing:

| Quantity        | Reading                                                                             |
| --------------- | ----------------------------------------------------------------------------------- |
| **pLDDT (0-1)** | > 0.9 very high · 0.7-0.9 confident · 0.5-0.7 low · < 0.5 disordered                 |
| **pTM**         | > 0.8 confident fold · 0.5-0.8 plausible · < 0.5 unreliable                           |
| **iPTM**        | > 0.8 confident interface · 0.5-0.8 possible · < 0.5 likely no binding                |
| **PAE (Å)**     | low = the two residues are confidently placed *relative to each other*               |

Reference points measured on this API (`--num-loops 10 --num-sampling-steps 50`):

| System                        | pLDDT | pTM  | iPTM | Reading                        |
| ----------------------------- | ----- | ---- | ---- | ------------------------------ |
| Ubiquitin (76 aa, real)       | 0.82  | 0.78 | —    | a real, well-folded domain     |
| Ubiquitin, **sequence shuffled** | 0.48 | 0.24 | —    | same composition, no fold      |
| Barnase + barstar (real complex) | 0.95 | 0.97 | 0.96 | a genuine, tight interface     |

That second row is the useful calibration: a scrambled sequence collapses to
pLDDT 0.48 / pTM 0.24. **If a designed or unfamiliar sequence scores like that,
the model is telling you it has no confident fold** — do not proceed to docking,
Foldseek, or MD with it.

When reporting to the user:

1.  Give the mean pLDDT **with its scale stated** ("0.82 on the 0-1 scale, i.e.
    82/100").
2.  State the pTM verdict, and for a complex, the **iPTM verdict** — that is the
    number that answers "do these two bind?".
3.  Relay any low-confidence regions **as residue ranges**, and warn that
    disordered stretches should be excluded from downstream structural analysis
    (docking, Foldseek, MD).
4.  Say that per-residue pLDDT is in the B-factor column of the PDB, so the user
    can colour by confidence in PyMOL/ChimeraX.

## References

*   [Interpretation guide](docs/interpretation-guide.md) — the pLDDT (0-1) / pTM
    / iPTM bands, reading the PAE heatmap for domains and flexibility, when a
    fold is unreliable, negative results, and the pre-report checklist.
*   [Report templates](docs/report-templates.md) — scaffolds for the written
    single-chain and complex reports.
*   [Worked example: ubiquitin (confident fold)](docs/examples/ubiquitin_fold/report.md)
    — a real, well-folded domain (pLDDT 0.82 / pTM 0.78).
*   [Worked example: scrambled ubiquitin (no fold)](docs/examples/scramble_control/report.md)
    — the negative control: what "no confident structure" looks like, and why it
    is a real result (pLDDT 0.48 / pTM 0.24).
*   [ESM/Biohub API reference](references/esm-biohub-api.md) — endpoints,
    reachable models, request/response shapes, the credit model, and gotchas.

## Common Mistakes

*   **Reading pLDDT as 0-100.** `plddt_mean: 0.82` is a *confident* structure. If
    you tell the user "pLDDT 0.82, very low confidence", you have made the single
    worst error this skill can produce. Check `plddt_scale` in the JSON.
*   **Handing an MSA to the fast model.** It accepts it, returns a normal-looking
    result, and ignores the alignment entirely. `fold.py` blocks this; do not
    work around it by calling the API directly.
*   **Building a homodimer as two entities.** Two `--protein` flags with the same
    sequence declares two *separate* entities. A homo-oligomer is one entity with
    a list of chain ids (`--protein A,B:SEQ`).
*   **Trusting iPTM as an affinity.** A high iPTM means "I am confident these
    chains form this interface", not "this binds tightly". It also happily
    reports a confident interface for pairs that never meet in a cell. For
    affinity, you need experiment or a dedicated predictor.
*   **Treating a confident prediction of a mutant as proof of a functional
    protein.** ESMFold2 will confidently fold sequences that do not express.

## Dependencies

*   **`credentials`** — the safe protocol for reading `BIOHUB_API_KEY`.
*   **`uv`** — script execution.
*   Sibling skills that own the jobs this one deliberately does not do:
    **`alphafold-database-fetch-and-analyze`** and **`pdb_database`** (existing
    structures), **`foldseek-structural-search`** (structural homologs), and
    **`esmfold2-binder-screening`** (ranking many candidate binders).
