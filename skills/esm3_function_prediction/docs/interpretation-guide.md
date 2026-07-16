# Interpretation Guide

This guide covers how to read ESM3 function-track output: annotation ranges,
InterPro accession parsing, domain-architecture coverage, the difference between
accession-bearing entries and free-text keywords, the mandatory cross-check
against the authoritative InterPro entry, an explicit **"Negative Results &
Scientific Integrity"** section, and a **pre-report reasoning checklist**. Read
this before writing any function-annotation report.

> [!CRITICAL] **Every label here is a PREDICTION, not a curated fact.** The ESM3
> function track decodes InterPro-style annotations straight from the residue
> sequence. Nothing was looked up; no database was consulted. A predicted
> accession that the model invents looks *exactly* like one it gets right. You
> cannot rely on any accession until you have resolved it against the real
> InterPro entry with the `interpro_database` skill.

--------------------------------------------------------------------------------

## The output at a glance

`predict_function.py` writes one JSON record per sequence. The fields you will
interpret:

| Field | What it is |
| :--- | :--- |
| `annotations` | The list of predicted `[label, start, end]` triples, parsed. |
| `annotations[i].label` | The raw label ESM3 emitted, e.g. `Ubiquitin-like domain (IPR000626)` or the bare keyword `glycoside`. |
| `annotations[i].name` | The label with any `(IPRxxxxxx)` suffix stripped. |
| `annotations[i].start` / `.end` | **1-indexed, inclusive** residue range. |
| `annotations[i].length` | `end - start + 1`. |
| `annotations[i].interpro_id` | The parsed `IPRxxxxxx` accession, or `null`. |
| `annotations[i].kind` | `interpro_entry` (carries an accession) or `keyword`. |
| `interpro_ids` | Every distinct accession across all annotations. |
| `domain_architecture` | Merged, non-overlapping spans — the domain count. |
| `coverage` | Fraction of residues carrying any annotation, in `[0, 1]`. |
| `caveat` | The predictions-not-facts warning. Reproduce it to the user. |

--------------------------------------------------------------------------------

## Reading annotation ranges (1-indexed, inclusive)

Every `start`/`end` is **1-indexed and inclusive**. Residue `start` and residue
`end` are both *inside* the annotation, and `length == end - start + 1`. This is
the same convention InterPro and UniProt use, and the script never converts it
for you.

> [!CAUTION] **The slicing trap.** To pull a domain's residues in Python you
> want `sequence[start - 1 : end]`. Writing `sequence[start:end]` silently drops
> the domain's first *and* last residue. A "site" predicted at 30–34 is five
> residues (30, 31, 32, 33, 34), not four.

- **Do not** re-index to 0-based or half-open ranges yourself before reporting.
- When you quote a span to the user, quote it exactly as the model gave it
  (e.g. "residues 27–115").

--------------------------------------------------------------------------------

## Parsing IPR accessions

An InterPro accession is the literal string `IPR` followed by **exactly six
digits** (`IPR000626`). The function track embeds it *inside* the label, as
`Ubiquitin-like domain (IPR000626)`. The script pulls it out into the structured
`interpro_id` field precisely so you can hand it straight to the
`interpro_database` skill.

- An annotation with an accession has `kind == 'interpro_entry'` and a non-null
  `interpro_id`.
- An annotation without one has `kind == 'keyword'` and `interpro_id == null`.
- `interpro_ids` (top level) is the de-duplicated set of every accession in the
  record — this is your verification worklist.

--------------------------------------------------------------------------------

## Domain-architecture coverage

**Coverage** is the fraction of residues carrying any annotation. Read it off
`coverage`, and read the number of *domains* off `domain_architecture` (the
merged, non-overlapping spans) — never off `len(annotations)`.

| Coverage | Reading |
| :--- | :--- |
| **> 0.8** | ESM3 confidently recognises a domain spanning essentially the whole chain. Ubiquitin gives **0.99** (0.9868), carbonic anhydrase II **0.98** (0.9808). |
| **0.3 – 0.8** | A recognised domain inside a larger protein, or a partial / multi-domain hit. Lysozyme gives **0.69** (0.6899): the mature chain's catalytic domain is residues 27–115 of 129. |
| **< 0.3** | Weak. Few residues implicated; treat the labels as hints only. |
| **0 annotations** | ESM3 does not recognise the sequence. A **real signal**, not a tool failure — see "Negative Results" below. |

> [!NOTE] **Annotation count is NOT domain count.** The labels are redundant and
> nested by design. Real ubiquitin returns **10 annotations describing one
> domain**: three InterPro entries, three overlapping keyword spans
> (`ubiquitin like`, `ubiquitin domain`, `core`), and two short `site` spans.
> They all collapse to a single `domain_architecture` span, **1–75**. To count
> domains, read `domain_architecture`; ubiquitin has exactly one.

### Sites vs domains

Short spans labelled `site` or `core` are predicted functional sites; long spans
are domains or families. For ubiquitin the two `site` spans land on **30–34** and
**55–61**. Report them as predicted site spans, but weight them lightly: they are
keyword spans (see below), and short keyword spans are the noisiest thing the
track emits.

--------------------------------------------------------------------------------

## InterPro entries vs free-text keywords

The track emits two kinds of label, and they are **not** equally trustworthy.

- **InterPro entries** carry an accession (`kind == 'interpro_entry'`). These are
  the reliable signal. Verify the accession, then trust it.
- **Free-text keywords** carry no accession (`kind == 'keyword'`). Noisier. Use
  them only as corroboration of the accession-bearing entries.

> [!WARNING] **Keywords confabulate even when the domain call is correct.** Real
> hen egg-white lysozyme reproducibly picks up the spurious keyword **`aminoacyl
> trna`** (and `trna`, `donors`, `group of`) — lysozyme has nothing to do with
> tRNA — while its InterPro entries (`IPR001916`, `IPR023346`) are exactly right.
> The model gets the biology right at the accession level and still emits noise
> at the keyword level. Both facts must inform how you report it: lead with the
> verified entries, and either drop the stray keywords or flag them as unverified.

**Rule of thumb:** if a keyword has no supporting accession-bearing entry over
the same span, treat it as unverified and do not put it in the headline finding.

--------------------------------------------------------------------------------

## Predictions vs curated annotations — the mandatory cross-check

This is the step that separates a defensible report from a fabrication. The
function track never queried InterPro; it *reconstructed* these accessions from
sequence. Confirm each one against the real entry before you rely on it.

**How to verify** — hand each `interpro_id` to the `interpro_database` skill:

```bash
uv run ./scripts/interpro_client.py fetch entry \
    --source_db interpro --accession IPR000626 \
    --output ipr000626.jsonl
```

Then ask three questions of the authoritative entry:

1.  **Does it exist?** A 404 means ESM3 invented the accession. Discard it.
2.  **Does its name/type match** the protein you expect? A ubiquitin sequence
    predicting a kinase domain is a red flag even if the accession is real.
3.  **Does its species distribution make sense** for your sequence's origin?

The three worked proteins in `docs/examples/known_proteins/` were verified this
way. **Every predicted accession resolved to a real entry whose name matched the
ESM3 label exactly:**

| Predicted accession | InterPro name (authoritative) | InterPro type |
| :--- | :--- | :--- |
| `IPR000626` | Ubiquitin-like domain | Domain |
| `IPR019956` | Ubiquitin domain | Domain |
| `IPR029071` | Ubiquitin-like domain superfamily | Homologous superfamily |
| `IPR001916` | Glycoside hydrolase, family 22 | Family |
| `IPR023346` | Lysozyme-like domain superfamily | Homologous superfamily |
| `IPR001148` | Alpha carbonic anhydrase domain | Domain |
| `IPR023561` | Carbonic anhydrase, alpha-class | Family |

`IPR001916` (Glycoside hydrolase family 22) is *exactly* the family of hen
egg-white C-type lysozyme, and `IPR001148` (Alpha carbonic anhydrase domain) is
the correct domain for human CA2. That is a confirmed prediction — and you only
know it is confirmed because you looked it up. Report predictions and verified
entries as two clearly separate categories.

If the sequence has no homologs at all (the case this skill exists for), the
cross-check is *still* mandatory: it tells you whether the *predicted family* is
real biology, even though the *protein* itself is unknown.

For the GO terms behind a predicted keyword, and for curated GO annotations when
the protein does have a database entry, use the `quickgo_database` skill.

--------------------------------------------------------------------------------

## Temperature and decoding effects

The two degenerate ends of the temperature range produce misleading output:

- **`temperature = 0` returns ZERO annotations for every sequence, including real
  ubiquitin.** The function track collapses onto the null token. This reads
  trivially as "this protein has no function". The script refuses it outright.
- **`temperature > 1.0` confabulates.** At 1.5 a *scrambled* ubiquitin produced
  36 annotations and an invented `NAD(P)-binding domain superfamily (IPR036291)`.
  Stay in **0.1–1.0**; the default is 1.0.

`--num-steps` (default 8) controls how many decoding passes run; more steps
surface more annotations but do not change the domain call. It is clamped to
`1..min(len, 100)`.

--------------------------------------------------------------------------------

## Negative Results & Scientific Integrity

> [!CRITICAL] **An empty result is a finding, not a failure.** A sequence with no
> recognisable fold returns **zero annotations**. That is exactly what a
> scrambled sequence returns, and it is informative. Report it as "no
> recognisable function", never as a tool error, and never stretch to invent one.

- **Value of the empty result.** Reporting "ESM3 predicts no recognisable domain"
  is a valid scientific finding. It distinguishes genuine dark matter (or a
  non-protein) from a recognised fold. The negative control in
  `docs/examples/scramble_control/` is a scrambled ubiquitin: it drops from 10
  annotations / 0.99 coverage to **0 annotations / 0.00 coverage**. (For
  corroboration, the same scramble raises ubiquitin's pseudo-perplexity from
  **1.05 to 18.70** under ESMC — it no longer reads as a protein at all.)
- **Never present a prediction as a fact.** These labels never touched a
  database. Always write "ESM3 predicts…", then state what the `interpro_database`
  cross-check showed. Reproduce the `caveat` field verbatim to the user.
- **Do not over-read keywords.** A stray keyword (`aminoacyl trna` on lysozyme)
  is not evidence of anything. If it has no supporting accession, do not build a
  claim on it.
- **Confidence is not reported.** The function track returns labels and spans
  only — there is no per-annotation probability. Coverage and the
  entry-vs-keyword distinction are your only confidence proxies. Do not invent a
  numeric confidence.

### What ESM3 function prediction does NOT do

- **It does not look anything up.** No InterPro / UniProt / GO query happens.
- **It does not return a confidence score** per annotation.
- **It does not predict structure** (use `esmfold2_structure_prediction`) or
  **score mutations** (use `esmc_mutation_effect_scoring`).
- **It does not resolve isoforms, catalytic mechanism, or exact active-site
  residues** — the `site`/`core` spans are coarse predicted regions, not curated
  active sites.

--------------------------------------------------------------------------------

## Pre-Report Reasoning Checklist

Complete this BEFORE writing the report. Ground every claim in the actual JSON
output and the InterPro cross-check.

### 1. Read the output correctly

-   [ ] **Domains, not annotations**: report the domain count from
    `domain_architecture`, not `len(annotations)`.
-   [ ] **Ranges**: quote spans exactly as given (1-indexed, inclusive); confirm
    `length == end - start + 1`.
-   [ ] **Coverage**: state `coverage` and place it in the >0.8 / 0.3–0.8 / <0.3
    band, with the biological reading for that band.

### 2. Separate the reliable signal from the noise

-   [ ] **Entries vs keywords**: list `interpro_entry` labels and `keyword`
    labels separately. The headline finding rests on the entries.
-   [ ] **Flag stray keywords**: identify keywords with no supporting accession
    over the same span (e.g. `aminoacyl trna` on lysozyme) and mark them
    unverified.

### 3. Verify every accession (mandatory)

-   [ ] **Resolve each `interpro_id`** with the `interpro_database` skill.
-   [ ] **Record the verdict** per accession: confirmed (real entry, matching
    name/type), mismatched (real entry, wrong biology), or invented (404).
-   [ ] **Two categories**: present "ESM3 predictions" and "verified InterPro
    entries" as clearly distinct blocks.

### 4. Integrity check

-   [ ] **Predictions labelled as predictions** everywhere the user sees them.
-   [ ] **Caveat reproduced** from the `caveat` field.
-   [ ] **Empty result respected**: if 0 annotations, report "no recognisable
    function" and stop — do not manufacture a mechanism.

### 5. Report readiness

-   [ ] **Figure present**: the domain-architecture PNG exists in the report
    folder and is embedded.
-   [ ] **Narrative**: the report tells the sequence's story (what domain, where,
    verified against what), not a raw table dump.
-   [ ] **Query answered**: the report addresses the original question, including
    when the answer is "no confident function".

--------------------------------------------------------------------------------
