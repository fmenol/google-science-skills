# ESM Skills — Documentation Parity Spec

Every ESM skill must reach the documentation standard set by
`alphagenome_single_variant_analysis` (the model-based analog). This spec defines
the artifacts and how to build them. Follow it exactly so the 12 skills read as
one consistent set.

## The standard (what AlphaGenome does)

- `docs/interpretation-guide.md` — how to read the outputs: thresholds, signal
  patterns, an explicit **"Negative Results & Scientific Integrity"** section,
  and a **pre-report reasoning checklist**.
- `docs/report-templates.md` — a numbered report scaffold the agent fills.
- `docs/examples/<case>/` — worked case studies, each with `README.md` (context),
  `report.md` (the analysis, quoting real numbers + embedding the figures), and
  the actual figure PNGs. The set MUST span a **positive** case and a
  **negative-result or model-limitation** case.
- `SKILL.md` mandates them: "review the worked example first", "read the
  interpretation guide before interpreting", "use the report template".

## Tiers

**Tier A** (interpretation-heavy; the agent produces a report and must judge
significance): `esmc_mutation_effect_scoring`, `esmfold2_structure_prediction`,
`esmfold2_binder_screening`, `esm3_function_prediction`,
`esmc_sae_feature_interpretation`, `esm3_protein_design`.
→ Full treatment: interpretation guide + report template + ≥2 worked examples
(≥1 positive, ≥1 negative/limitation).

**Tier B** (capability skills): `esmc_protein_embeddings`,
`esmc_embedding_layer_sweep`, `esm3_inverse_folding`, `esm3_guided_generation`,
`esm_protein_tracks`, `esm3_secondary_structure_sasa`.
→ Lighter: 1 worked example (`docs/examples/<case>/` with README + report + figure)
and a short "Interpreting the Output" already in SKILL.md. No separate report
template required, but an interpretation guide is welcome if the outputs need it.

## How to generate example data — FREE, from cassettes

The evals already recorded every skill's calls on canonical proteins. Regenerate
the artifacts by running your skill's CLI in **replay** mode (no key, no credits):

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay
uv run --no-project skills/<skill>/scripts/<script>.py <subcommand> \
    <the SAME inputs/params your eval uses> \
    --output skills/<skill>/docs/examples/<case>/result.json
```

Rules:
- **Use the exact fixture inputs + params your eval uses** (see
  `evals/eval_<skill>.py` and `evals/fixtures.py`: UBIQUITIN, LYSOZYME, BARNASE,
  BARSTAR, CA2, MYOGLOBIN, GFP). Those payloads are in the cassettes, so replay
  hits them for free. If you must issue a genuinely new call, keep it minimal —
  there is a 100-credit/day budget shared across everyone.
- Copy the figures your skill renders (PNGs) into the example folder.
- Quote **real numbers** from the output in `report.md`. Never invent values.
- The negative example must be a real run on a negative control (a scrambled
  sequence, a decoy binder, a disordered/low-confidence case) — the evals already
  compute these; reuse them.

## Required file layout per skill

```
skills/<skill>/
  references/
    citation.bib                     (exists)
    esm-biohub-api.md  -> symlink or 1-line pointer to the shared reference
  docs/
    interpretation-guide.md          (Tier A; optional Tier B)
    report-templates.md              (Tier A)
    examples/
      <positive_case>/
        README.md   report.md   *.png
      <negative_case>/             (Tier A: required)
        README.md   report.md   *.png
```

For the shared API reference, do NOT copy the file. Add a
`references/esm-biohub-api.md` that is a one-line pointer:
`See the shared ESM/Biohub API reference: ../../esm_common/references/esm-biohub-api.md`
and link that path from SKILL.md. (Skills are installed together in the bundle,
so the relative path resolves.)

## SKILL.md wire-up (all skills)

Add a **References** section listing the new docs with relative links, e.g.:

```markdown
## References
- [Interpretation guide](docs/interpretation-guide.md) — how to read the output.
- [Report template](docs/report-templates.md) — scaffold for the written report.
- [Worked example: <positive>](docs/examples/<case>/report.md)
- [Worked example: <negative>](docs/examples/<case>/report.md) — what "no effect" looks like.
- [ESM/Biohub API reference](references/esm-biohub-api.md)
```

Tier A skills additionally add these imperatives to **Core Rules**:
- "Before interpreting results, read `docs/interpretation-guide.md`."
- "Review at least one worked example in `docs/examples/` before writing a report,
  including a negative one so you know what 'no effect' looks like."
- "Write the report using `docs/report-templates.md`."

## Quality bar

- Prose in Google style, wrapped ~80 cols, matching the voice of the existing
  SKILL.md files.
- Every figure referenced in a `report.md` must exist in that folder.
- Every doc link in SKILL.md must resolve (the offline suite will check this).
- Interpretation guides must include the measured thresholds from the shared API
  reference §7 and the skill's own "Interpreting the Output" numbers.
