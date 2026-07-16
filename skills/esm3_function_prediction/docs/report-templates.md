# Report Templates

This file provides a scaffold for reporting ESM3 function-prediction results.
Follow the structure so every report separates **predictions** from **verified
entries** and never presents a decoded label as an established fact.

> [!IMPORTANT] **DO NOT** use the agent's default artifact directory for this
> report. Save it as a regular file named `report.md` directly in the sequence's
> output directory (e.g. `orphan_kinase/report.md`) in the workspace. Use
> relative paths for the embedded figure (`filename.png`), and generate that
> figure with `predict_function.py plot`.

> [!CRITICAL] Before filling this in, read `docs/interpretation-guide.md`, review
> the worked examples in `docs/examples/` (including the negative one so you know
> what "no confident function" looks like), and resolve **every** predicted
> `interpro_id` with the `interpro_database` skill.

--------------------------------------------------------------------------------

## 1. Standard Function-Annotation Report Template

Use this when ESM3 returns one or more annotations.

```markdown
# Function Prediction Report: {sequence_id} ({length} aa)

## 1. Summary of Findings
[One paragraph. State the predicted domain/family in plain language, where it
sits on the chain, and the residue coverage. Make the prediction status explicit:
"ESM3 predicts ...". State up front whether the InterPro cross-check CONFIRMED the
call. If 0 annotations, say so plainly and use the negative template below.]

## 2. Sequence Context
- **ID**: {sequence_id}
- **Length**: {length} aa
- **Provenance**: [orphan / de-novo design / metagenomic ORF / known protein]
- **Model / params**: esm3-open-2024-03, num_steps={n}, temperature={t}

## 3. Predicted Domain Architecture
[Read from `domain_architecture`, NOT len(annotations). One row per merged span.]

| Span (1-indexed, incl.) | Length | # annotations | InterPro accessions |
| :--- | :--- | :--- | :--- |
| {start}-{end} | {len} | {n} | {IPR ids} |

- **Coverage**: {coverage:.1%} of the chain -- [>0.8 whole-chain / 0.3-0.8
  partial or embedded / <0.3 weak].

## 4. Predicted Annotations
[List InterPro entries and free-text keywords SEPARATELY. The headline finding
rests on the accession-bearing entries; keywords are corroboration only.]

**InterPro entries (accession-bearing -- the reliable signal):**

| Label | Span | Accession |
| :--- | :--- | :--- |
| {label} | {start}-{end} | {IPRxxxxxx} |

**Free-text keywords (noisier -- corroboration only):**
{list}. [Flag any keyword with no supporting accession over the same span as
UNVERIFIED, e.g. a spurious `aminoacyl trna` on a lysozyme.]

## 5. InterPro Verification (mandatory)
[Resolve each accession with the `interpro_database` skill and record the verdict.
This is what turns a prediction into a defensible statement.]

| Predicted accession | Authoritative InterPro name | Type | Verdict |
| :--- | :--- | :--- | :--- |
| {IPRxxxxxx} | {name from interpro_database} | {domain/family/...} | CONFIRMED / MISMATCH / INVENTED |

- **Command used**: `uv run ./scripts/interpro_client.py fetch entry --source_db
  interpro --accession {IPRxxxxxx} --output {id}.jsonl`
- [If the protein has a UniProt accession, also fetch its curated annotations and
  compare; report agreements and disagreements explicitly.]

## 6. Domain-Architecture Figure
![Predicted domain architecture]({sequence_id}_domains.png)
*Fig 1: ESM3-PREDICTED domain architecture. Dark blue = InterPro entries
(verify each accession); light blue = free-text keywords. Ranges 1-indexed,
inclusive.*

## 7. Reliability Assessment
[Coverage band + its reading; entries-vs-keywords split; the outcome of the
InterPro cross-check. State that no per-annotation confidence score exists.]

## 8. Limitations & Caveats
[Reproduce the `caveat` field. State the function track does not look anything up,
returns no confidence score, and does not predict structure/mutations. Note any
stray keywords excluded from the headline finding.]

## 9. Conclusion
[Narrative synthesis: the predicted function, whether verification confirmed it,
and how confident the reader should be. If the sequence is a genuine orphan,
state that the cross-check validates the predicted FAMILY even though the protein
itself is uncharacterised.]
```

--------------------------------------------------------------------------------

## 2. Negative / No-Confident-Function Template

Use this when ESM3 returns **zero annotations** (coverage 0.00). An empty result
is a finding -- see `docs/examples/scramble_control/`.

```markdown
# Function Prediction Report: {sequence_id} ({length} aa) -- No Confident Function

## 1. Summary
ESM3 returned **0 predicted annotations** (coverage 0.00) for {sequence_id}. The
function track does not recognise a domain in this sequence. This is a real
signal, not a tool failure: it is the same result a scrambled sequence returns.

## 2. Sequence Context
- **ID / Length**: {sequence_id} / {length} aa
- **Provenance**: [orphan / de-novo / metagenomic / scrambled control]
- **Model / params**: esm3-open-2024-03, num_steps={n}, temperature={t}
  (temperature is in the valid 0.1-1.0 band, so the empty result is not an
  artefact of temperature=0).

## 3. Interpretation
[Two possibilities, both worth stating: (a) genuine dark matter -- a real protein
with no fold ESM3 recognises; (b) not a plausible protein -- a scramble, a
frameshift, an assembly artefact. Corroborate with any independent evidence,
e.g. structure prediction confidence or pseudo-perplexity vs a real protein.]

## 4. Figure
![Empty domain architecture]({sequence_id}_domains.png)
*Fig 1: The backbone carries no predicted annotation -- coverage 0.00.*

## 5. What This Does and Does Not Rule Out
[The empty result means ESM3 sees no recognised domain. It does NOT prove the
sequence is non-functional -- a novel fold outside training data would look the
same. Recommend orthogonal checks (structure prediction, HMM search) if function
matters.]
```

--------------------------------------------------------------------------------
