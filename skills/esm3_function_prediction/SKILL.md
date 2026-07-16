---
name: esm3-function-prediction
description: >
  Predict functional annotations — InterPro-style domains, families, sites and
  keywords — for a protein directly from its amino-acid sequence with the ESM3
  function track. No homology search, no BLAST, no database hit required. Use
  when a sequence has no detectable homologs, is a novel or de-novo design, or
  is metagenomic "dark matter", and InterPro/BLAST come up empty; or to get a
  first-pass domain architecture for an unannotated sequence. Do not use when
  the protein already has a UniProt accession or a known database entry — use
  the `interpro_database` skill for authoritative curated domains, or
  `quickgo_database` for curated GO terms. This skill PREDICTS annotations;
  those skills LOOK UP known ones.
---

# ESM3 Function Prediction

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If .licenses/esm3_function_prediction_LICENSE.txt does
    not already exist in the workspace root directory then (1) prominently notify
    the user to check the terms at https://biohub.org/acceptable-use-policy/ and
    https://biohub.ai/, then (2) create the file recording the notification text
    and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for it. Never print
    the key.

## Overview

ESM3's `function` track decodes functional annotations straight from the residue
sequence. The model was trained on InterPro annotations, so it emits InterPro
*keywords* ('ubiquitin domain', 'core', 'site') alongside full InterPro *entries*
that carry an accession ('Ubiquitin-like domain (IPR000626)').

Because nothing is looked up, this works precisely where homology-based
annotation fails: orphan sequences, de-novo designs, and metagenomic dark
matter. A scrambled sequence with no fold returns **nothing**, which makes an
empty result informative rather than a failure.

**This skill does NOT:**

-   Look anything up. It never queries InterPro, UniProt or GO. Every label is a
    model prediction. To retrieve *known* annotations for a protein that has a
    database entry, use `interpro_database` / `quickgo_database`.
-   Predict structure (use `esmfold2_structure_prediction`) or score mutations
    (use `esmc_mutation_effect_scoring`).
-   Return confidence scores. The function track gives labels and spans only —
    there is no per-annotation probability.

## Core Rules

-   **Read [`docs/interpretation-guide.md`](docs/interpretation-guide.md) before
    interpreting any result.** It covers 1-indexed inclusive ranges, IPR
    accession parsing, coverage bands, keyword noise, and the mandatory InterPro
    cross-check.
-   **Review a worked example in [`docs/examples/`](docs/examples/) before
    writing a report** — including the negative one
    ([`scramble_control/`](docs/examples/scramble_control/report.md)) so you know
    what "no confident function" looks like.
-   **Write the report with
    [`docs/report-templates.md`](docs/report-templates.md).**
-   **These are MODEL PREDICTIONS, not curated annotations.** ALWAYS label them
    as predictions in anything you show the user. ALWAYS report the caveat that
    ships in the `caveat` field of the output. NEVER present a predicted label as
    an established fact about the protein.
-   **ALWAYS cross-check a predicted InterPro accession against the real entry**
    using the `interpro_database` skill before relying on it. The script parses
    each accession into a structured `interpro_id` field precisely so you can
    hand it straight over. An accession that ESM3 invents looks exactly like one
    it gets right.
-   **Ranges are 1-indexed and INCLUSIVE.** `length == end - start + 1`. Residue
    `start` and residue `end` are both inside the annotation. Do not convert to
    0-indexed or half-open ranges yourself.
-   **Trust accession-bearing entries more than free-text keywords.** The
    keywords are noisier — real lysozyme reproducibly picks up a spurious
    `aminoacyl trna` label while its InterPro entries are exactly right.
-   **Use the Wrapper**: ALWAYS execute `scripts/predict_function.py`. Do not
    compute coverage, merge spans, or count domains yourself; always use the
    script's output.
-   **NEVER use `--temperature 0`.** It returns zero annotations for *every*
    sequence, including real ubiquitin (the script refuses it).
-   **Notification**: If this skill is used, ensure this is mentioned in the
    output, and state that the annotations are ESM3 predictions.

## Dependencies

This skill **predicts**; the following sibling skills **look up** what is already
known. They are complements, not alternatives — the intended workflow is
**predict, then verify**.

| Skill | Role | Use it to |
|---|---|---|
| `interpro_database` | Authoritative, curated | Resolve a predicted `interpro_id` (e.g. `IPR001916`) to the real InterPro entry: its true name, type, member signatures, GO mappings and species distribution. **This is how you confirm or refute a prediction.** |
| `quickgo_database` | Authoritative, curated | Resolve the GO terms behind a predicted keyword, and fetch curated GO annotations for the protein if it has a database entry. |

Recommended workflow:

1.  Run `predict` on the sequence (this skill).
2.  Take each `interpro_id` from the output and fetch the real entry with
    `interpro_database`. Does the entry's description match the protein you
    expect? Does its species distribution make sense?
3.  If the protein *does* have a UniProt accession, also fetch its curated
    annotations from `interpro_database` / `quickgo_database` and compare. Report
    agreements and disagreements explicitly.
4.  Report predictions and verified entries as two clearly separate categories.

If the sequence has no homologs at all (the case this skill exists for), step 2
is still mandatory: it tells you whether the *predicted family* is real biology,
even though the *protein* is unknown.

## Utility Scripts

All commands use `uv run --no-project`. This is mandatory — without
`--no-project`, `uv` walks up the directory tree, finds an unrelated
`pyproject.toml` and tries to build that project instead.

**1. `predict` — annotate one sequence**

```bash
uv run --no-project scripts/predict_function.py predict \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --id ubiquitin \
  --output ubiquitin_function.json
```

Or from a FASTA (first record):

```bash
uv run --no-project scripts/predict_function.py predict \
  --fasta orphan.fasta --output orphan_function.json --num-steps 16
```

Writes JSON with, for every annotation: `label`, `name` (accession stripped),
`start`, `end`, `length`, `interpro_id` (or `null`), and `kind`
(`interpro_entry` | `keyword`). Plus `interpro_ids` (all accessions),
`domain_architecture` (merged, non-overlapping spans sorted by position),
`covered_residues`, `coverage` (fraction of residues with any annotation), and
`caveat`.

**2. `batch` — annotate a whole FASTA**

```bash
uv run --no-project scripts/predict_function.py batch \
  --fasta metagenome_orfs.fasta \
  --output batch_function.json \
  --summary batch_summary.csv \
  --max-workers 8
```

Runs concurrently (capped at 8 workers). A record that fails is recorded with an
`error` field instead of sinking the batch. `--summary` writes one row per
record: id, length, n_annotations, coverage, InterPro accessions, domain
architecture, top labels.

**3. `plot` — domain-architecture diagram**

```bash
uv run --no-project scripts/predict_function.py plot \
  --input ubiquitin_function.json --output ubiquitin_domains.png
```

Draws every annotation as a labelled bar along the sequence axis, longest first.
InterPro entries are drawn in dark blue, free-text keywords in light blue.
Accepts a `batch` JSON too; select the record with `--record <id>`.

**Flags** (`predict` and `batch`): `--num-steps` (default 8; must be
`1..min(len, 100)`, clamped automatically; more steps → more annotations),
`--temperature` (default 1.0; **use 0.1–1.0**), `--model` (default
`esm3-open-2024-03`, the only ESM3 this key can reach).

## Interpreting the Output

**Coverage** — the fraction of residues carrying any annotation.

| Coverage | Reading |
|---|---|
| > 0.8 | ESM3 confidently recognises a domain spanning essentially the whole chain. Ubiquitin gives 0.99, carbonic anhydrase II 0.98. |
| 0.3 – 0.8 | A recognised domain inside a larger protein, or a partial/multi-domain hit. Lysozyme gives 0.70 (the mature chain's catalytic domain is 27–115 of 129). |
| < 0.3 | Weak. Few residues implicated; treat the labels as hints only. |
| 0 annotations | ESM3 does not recognise the sequence. This is a **real signal**: it is what a scrambled sequence returns. Either genuine dark matter, or not a plausible protein. Report it as "no recognisable function", never as a tool failure. |

**Annotation count is NOT domain count.** The labels are redundant and nested by
design — ubiquitin returns 10 annotations describing *one* domain
(`ubiquitin`, `ubiquitin like`, `ubiquitin domain`, three InterPro entries, and
two `site` spans). To count domains, read `domain_architecture`: ubiquitin
collapses to a single span, 1–75.

**Sites vs domains.** Short spans labelled `site` / `core` are predicted
functional sites. For ubiquitin these land on 30–34 and 55–61 — the
hydrophobic-patch region around Ile44 and the C-terminal approach. Long spans are
domains/families.

**Reference outputs** (verified on this API, reproducible across runs):

| Protein | Annotations | Coverage | Predicted InterPro |
|---|---|---|---|
| Ubiquitin (76 aa) | 10 | 0.99 | IPR000626, IPR019956, IPR029071 |
| Lysozyme (129 aa) | ~15 | 0.70 | IPR001916 (GH22), IPR023346 |
| Carbonic anhydrase II (260 aa) | ~17 | 0.98 | IPR001148, IPR023561 |
| Scrambled ubiquitin | **0** | 0.00 | — |

All six accessions are correct: GH22 is exactly the family of hen egg-white
C-type lysozyme, and IPR001148 is the alpha-class carbonic anhydrase domain. The
model gets the biology right — and still emits `aminoacyl trna` noise on
lysozyme. Both facts must inform how you report it.

## Common Mistakes

1.  **Reporting a prediction as a fact.** The single most damaging error. These
    labels never touched a database. Say "ESM3 predicts…", then verify the
    accession with `interpro_database` and say what that showed.
2.  **`--temperature 0` and `--temperature > 1`.** At 0 the function track
    collapses onto the null token and returns **zero annotations for every
    sequence, including real ubiquitin** — trivially misread as "this protein has
    no function". The script refuses it. Above 1.0 it confabulates: at
    temperature 1.5 a *scrambled* ubiquitin produced 36 annotations and an
    invented `NAD(P)-binding domain superfamily (IPR036291)`. Stay in 0.1–1.0.
3.  **Counting annotations as domains.** Ten annotations, one domain. Use
    `domain_architecture`, not `len(annotations)`.
4.  **Treating ranges as 0-indexed or half-open.** They are 1-indexed and
    inclusive; `length == end - start + 1`. Slicing `sequence[start:end]` in
    Python silently drops the first and last residue of the domain — you want
    `sequence[start - 1 : end]`.

## References

-   [Interpretation guide](docs/interpretation-guide.md) — how to read the
    output: ranges, accession parsing, coverage bands, keyword noise, the
    mandatory InterPro cross-check, negative results, and the pre-report
    checklist.
-   [Report template](docs/report-templates.md) — scaffold for the written
    report (standard and no-confident-function variants).
-   [Worked example: known proteins](docs/examples/known_proteins/report.md) —
    ubiquitin, lysozyme, CA2 predicted then verified against InterPro (positive).
-   [Worked example: scramble control](docs/examples/scramble_control/report.md)
    — what "no confident function" (0 annotations) looks like (negative).
-   [ESM/Biohub API reference](references/esm-biohub-api.md) — shared platform
    reference: endpoints, reachable models, credit model, gotchas.
-   [`references/citation.bib`](references/citation.bib) — cite Hayes et al.
    (2025) for the ESM3 function track; Blum et al. (2025) for InterPro when you
    resolve a predicted accession.
