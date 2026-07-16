---
name: esmc-embedding-layer-sweep
description: >
  Chooses the optimal ESMC layer for a downstream supervised task by training a
  linear probe on every layer's mean-pooled embedding with stratified
  cross-validation. Use when the user asks "which layer should I use?", wants a
  layer sweep, asks whether the last layer is the right one, wants to know where
  in the network a property is encoded, or is building a classifier / linear
  probe on ESMC features and needs to pick the layer empirically. Requires a
  labelled dataset (a FASTA plus a class label per sequence). Do not use when the
  user simply wants embedding vectors for sequences and has no labels and no
  layer decision to make — use the `esmc-protein-embeddings` skill instead. Do
  not use for scoring point mutations (use `esmc-mutation-effect-scoring`) or for
  interpreting what individual features mean (use
  `esmc-sae-feature-interpretation`).
---

# ESMC Embedding Layer Sweep

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH.
2.  **User Notification**: If
    .licenses/esmc_embedding_layer_sweep_LICENSE.txt does not already exist in
    the workspace root directory then (1) prominently notify the user to check
    the terms at https://biohub.org/acceptable-use-policy/ and
    https://biohub.ai/, then (2) create the file recording the notification text
    and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for and request
    this key.

## Overview

ESMC is a transformer. Its layers form a hierarchy: early layers hold local
residue chemistry, middle layers hold motifs and evolutionary constraint, and
the **final layer is specialised to the masked-token training objective** rather
than to your task. **The last layer is frequently not the best layer.** In the
source tutorial (enzyme-class prediction over ~600 SwissProt enzymes with ESMC
600M) an intermediate layer reached **MCC 0.935 while the final layer reached
only 0.843**.

This skill treats layer choice as a hyperparameter and measures it. It embeds a
labelled dataset once, capturing **every layer in a single request per
sequence**, then trains an independent linear probe (StandardScaler + logistic
regression) on each layer's mean-pooled embedding under stratified
cross-validation, and reports which layer carries the most linearly-decodable
task signal.

**Do NOT use when:**

-   The user just wants embedding vectors and has no labels — use
    **`esmc-protein-embeddings`**.
-   The user wants variant effect scores — use **`esmc-mutation-effect-scoring`**.
-   The user wants to know what a feature *means* — use
    **`esmc-sae-feature-interpretation`**.
-   There are no labels at all. A sweep is supervised by definition; without
    labels there is nothing to optimise the layer against.

## Core Rules

-   **NEVER download model weights.** No `torch`, no `transformers`, no
    `huggingface_hub`, no `esm` PyPI package. Everything runs remotely on the
    Biohub API through the vendored `esm_biohub.py`.
-   **ALWAYS run scripts with `uv run --no-project`.** Without `--no-project`,
    `uv` walks up the directory tree, finds an unrelated `pyproject.toml`, and
    tries to build that project instead.
-   **NEVER loop over layers issuing one request per layer.** One
    `mean_hidden_states` call returns **all** `n_layers + 1` layers at once. A
    sweep costs N requests, not N x n_layers. `embed-dataset` already does this.
-   **Do not compute MCC, pick the best layer, or eyeball the curve yourself;
    always use the script's output.** The scripts do the cross-validation and
    the model selection.
-   **`--n-splits` must be <= the smallest class count.** Stratified CV needs at
    least one member of every class in every fold. The script enforces this and
    tells you the maximum.
-   **Never mix pooling conventions.** Features used at inference must come from
    the same call that produced the training features (`mean_hidden_states`).
    See Common Mistakes.
-   If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All commands run from the skill directory.

**1. Embed the dataset (one request per sequence, all layers)**

Takes a FASTA plus a labels CSV (`id,label`) and caches an
`(n_sequences, n_layers + 1, hidden_dim)` array to `.npy`, with a `.meta.json`
sidecar carrying the ids, labels and model. The sweep is then **re-runnable
offline** with no further API calls.

```bash
uv run --no-project scripts/layer_sweep.py embed-dataset \
  --fasta enzymes.fasta \
  --labels enzymes_labels.csv \
  --model esmc-600m-2024-12 \
  --output /path/to/out/embeddings.npy
```

`--model` is required (it fixes the layer count and dimension):
`esmc-300m-2024-12` (30 blocks, dim 960), `esmc-600m-2024-12` (36 blocks, dim
1152), `esmc-6b-2024-12` (80 blocks, dim 2560). Use `--id-column` /
`--label-column` if your CSV headers differ from `id` / `label`.

**2. Sweep every layer**

Trains one probe per layer under `StratifiedKFold` and scores with Matthews
correlation coefficient. Writes per-layer mean/std to both CSV and JSON.

```bash
uv run --no-project scripts/layer_sweep.py sweep \
  --embeddings /path/to/out/embeddings.npy \
  --n-splits 5 \
  --output-json /path/to/out/sweep.json \
  --output-csv /path/to/out/sweep.csv
```

**3. Plot the sweep**

MCC vs layer with a +/-1 std band, best and last layer marked. A second panel
plots cross-validated log-loss, which is what you read when the MCC curve is
flat.

```bash
uv run --no-project scripts/layer_sweep.py plot \
  --sweep /path/to/out/sweep.json \
  --output /path/to/out/sweep.png
```

**4. Train, score and save a probe at one layer**

`--layer` accepts a 0-based row index or a Python-style negative index, so the
`[-2]` rule of thumb can be written literally. Optionally classifies new
sequences.

```bash
uv run --no-project scripts/layer_sweep.py probe \
  --embeddings /path/to/out/embeddings.npy \
  --layer 25 \
  --n-splits 5 \
  --output-model /path/to/out/probe.joblib \
  --output-json /path/to/out/probe_metrics.json \
  --predict-fasta unknowns.fasta \
  --output-predictions /path/to/out/predictions.csv
```

## Interpreting the Output

For the full layer-selection heuristics, the saturation caveat, and a pre-report
checklist, read [`docs/interpretation-guide.md`](docs/interpretation-guide.md);
for a worked run on real data, see the
[layer-selection example](docs/examples/layer_selection/report.md).

**Layer indexing (get this right).** Row **0 is the embedding layer** — raw
token embeddings, no attention, essentially amino-acid composition. Row `k` is
the output of transformer block `k`. A 30-block model therefore has **31 rows**,
indexed 0..30. Report layers in this 0-based convention.

**MCC.** Ranges from -1 (total disagreement) through 0 (random) to +1 (perfect).
It uses all four confusion-matrix cells, so unlike accuracy it stays honest
under class imbalance. Read it as: > 0.8 strong, 0.5-0.8 useful, 0.2-0.5 weak,
~0 the probe learned nothing.

**Read the shape of the curve, not just the peak.** The expected shape is a rise
through the early and middle layers, a broad plateau, and a **drop at the final
layer** — the last layer is optimised for masked-token prediction, not for your
task, and routinely degrades. If the reported `best_layer` is far from the final
layer, say so explicitly: that is the whole point of the sweep.

**If you cannot afford a sweep, use these defaults** (from the source tutorial):

| Task | Layer to use |
|---|---|
| Physicochemical / functional properties (thermostability, solubility, expression) | **`[-2]`**, the second-to-last layer |
| High-level functional similarity (EC number, GO term, homology, family) | **~3/4 of the way through** the network |
| Residue-residue contacts | **late** layers |
| Anything | **not** the embedding layer (row 0) — it has seen no context |

These are priors, not answers. The optimal layer **does not transfer** reliably
across tasks or across model sizes; if you have labels, measure.

**`mcc_saturated: true` is a warning, not a result.** On easy or small datasets
every layer can hit MCC = 1.0 at once. MCC is bounded at 1.0, so it then has no
power to rank layers, and a naive `argmax` would hand back the lowest tied index
— which is **layer 0, the embedding layer**, the worst possible advice. When
this happens the script says so loudly and breaks the tie with cross-validated
**log-loss**, a strictly proper score that keeps improving after the
classification is already perfect (it rewards a more confident, better-separated
boundary). **Relay the saturation warning to the user.** The honest conclusion
is that the dataset cannot distinguish those layers; more or harder examples are
needed before committing.

**Sanity check the sweep against the embedding layer.** If row 0 genuinely beats
every transformer layer, the task is decidable from amino-acid composition alone
and ESMC is adding nothing — reconsider the problem, not the layer.

## Common Mistakes

1.  **Requesting one layer at a time.** `mean_hidden_states` returns every layer
    in a single request. Looping over layers multiplies the cost by ~30-80x for
    identical results. Never do it.
2.  **Off-by-one on the layer index.** Row 0 is the *embedding* layer, not
    transformer block 1. The tutorial prints 1-based layer numbers ("layer 25 of
    37"); this skill uses 0-based rows throughout. Do not mix them.
3.  **Mixing pooling conventions between training and inference.**
    `client.mean_hidden_states(seq, model)[L]` is pooled **server-side over all
    tokens, BOS and EOS included**. `client.embed(seq, model, layer=L)` pools
    **locally over residues only, BOS/EOS trimmed**. These are *different
    vectors for the same layer*. A probe trained on one and applied to the other
    silently degrades. `embed-dataset` and `probe --predict-fasta` both use
    `mean_hidden_states`; keep it that way.
4.  **Trusting a saturated or tiny sweep.** With fewer than ~20 sequences the
    fold estimates are noisy and can pin to the ceiling. Check `mcc_saturated`
    and the std band before quoting a winning layer.
5.  **Assuming the winning layer transfers.** The best layer for EC
    classification is not the best layer for stability regression, and the best
    layer for ESMC 600M is not row-for-row the best layer for ESMC 6B. Re-sweep
    when the task or the model changes.

## Dependencies

-   **`uv`** — environment management for every script invocation.
-   **`credentials`** — the safe protocol for obtaining `BIOHUB_API_KEY`.
-   **`esmc-protein-embeddings`** — sibling skill; use it when the user wants
    embeddings rather than a layer decision.
-   **`esmc-mutation-effect-scoring`**, **`esmc-sae-feature-interpretation`** —
    sibling ESMC skills for variant effects and feature interpretation.

## References

-   [Worked example: choosing a layer](docs/examples/layer_selection/report.md)
    — lysozyme vs barnase on ESMC 300M. The sweep crowns **layer 23 / 30**
    (~77% depth); the final layer degrades and MCC saturates, so log-loss breaks
    the tie. See the [example README](docs/examples/layer_selection/README.md)
    for context.
-   [Interpretation guide](docs/interpretation-guide.md) — the layer-selection
    heuristics, the MCC bands, the `[-2]` and ~3/4-depth defaults, the
    MCC-saturation caveat, and a pre-report checklist. Read it before quoting a
    winning layer.
-   [ESM/Biohub API reference](references/esm-biohub-api.md) — shared platform
    reference: endpoints, reachable models, the daily credit model, and the
    verified gotchas that apply to every ESM skill.
-   [`references/citation.bib`](references/citation.bib) — ESMC, ESM3, MCC, and
    the layer-interpretability literature.
