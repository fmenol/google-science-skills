# ESM / Biohub Platform API Reference

The shared reference for every ESM skill. It documents the endpoints, the
reachable models, the request/response shapes, the credit model, and the verified
gotchas. Each skill's own `SKILL.md` covers its task; this covers the platform
they all sit on.

All ESM skills talk to `https://biohub.ai` through the vendored `esm_biohub.py`
client. They do **not** import the `esm` PyPI package (which pulls torch + a fork
of transformers and downloads weights); every capability here is reachable over
plain JSON with no model download and no GPU.

---

## 1. Authentication & setup

- **Credential:** `BIOHUB_API_KEY`, sent as `Authorization: Bearer <key>`.
  Register at https://biohub.ai/developer-console/api-keys. Verify its presence
  with the `credentials` skill's safe protocol; never print it.
- **Transport:** `POST https://biohub.ai/api/v1/<endpoint>` with a JSON body.
  The client (`BiohubClient`) handles auth, rate limiting (cross-process file
  lock via `polite-http`), and retry/backoff.
- **Runner:** always `uv run --no-project`. Without `--no-project`, `uv` finds an
  unrelated parent `pyproject.toml` and fails.

## 2. Reachable models

| Model | Family | Notes |
|---|---|---|
| `esmc-300m-2024-12` | ESMC | 30 layers, dim 960 |
| `esmc-600m-2024-12` | ESMC | 36 layers, dim 1152 — default ESMC |
| `esmc-6b-2024-12` | ESMC | 80 layers, dim 2560; `ith_hidden_layer=-1` rejected |
| `esm3-open-2024-03` | ESM3 | the only reachable ESM3 (alias `esm3-sm-open-v1`) |
| `esmfold2-fast-2026-05` | ESMFold2 | default folding model (~2 s for 76 aa) |
| `esmfold2-2026-05` | ESMFold2 | slower; the only one that uses an MSA |

`esm3-medium-*`, `esm3-large-*`, and `esm3-*-multimer` return **HTTP 403** — the
key has no access. Any tutorial written against them is re-pointed at
`esm3-open-2024-03`.

## 3. Endpoints

Shapes below use L = number of residues, so tokenized tracks are length **L+2**
(a BOS token, the residues, an EOS token). **Residue *i* is at token index
*i+1*.**

### `encode` (ESMC) — free
Tokenize a sequence. `{"inputs": {"sequence": str}, "model": str}` →
`outputs.sequence` = list of L+2 ints. **Not billed** — returns 200 even when the
daily credit allowance is exhausted, so it is useless as a credit probe.

### `logits` (ESMC) — billed
Forward pass. Request carries `inputs.sequence` (the L+2 tokens from `encode`)
and a `logits_config`. Returns, per the flags set:

| Field | Shape | Enable with |
|---|---|---|
| `logits.sequence` | (L+2, 64) | `sequence: true` |
| `embeddings` | (1, L+2, D) | `return_embeddings: true` |
| `mean_embedding` | (1, 1, D) | `return_mean_embedding: true` |
| `hidden_states` | (n_layers+1, 1, L+2, D) | `return_hidden_states: true` |
| `mean_hidden_state` | (1, n_layers+1, D) | `return_mean_hidden_states: true` |
| `sae_outputs[name]` | sparse top-k (below) | `sae_config` set |

- The logit vocabulary is padded to **64**; index it by `eb.VOCAB[aa]`, never by
  position. `eb.AA20_IDX` gives the 20 canonical columns.
- **Layer 0 is the embedding layer**; layers `1..n_layers` are transformer
  blocks. `ith_hidden_layer=-1` means "all layers" and is rejected for ESMC 6B
  and every ESM3 model. Negative indices other than the client-resolved ones
  crash the server (HTTP 500) — use `eb.BiohubClient.resolve_layer`.
- Mean-pooled tensors are pooled **including** BOS/EOS; per-residue tensors are
  trimmed `[1:-1]` by the client.

### `fold` (ESMFold2) — billed
`{"sequence": str, "model": str, "num_loops": int, "num_sampling_steps": int,
"include_pae": bool, ...}`. Returns **at the top level** (not under `outputs`):
`coordinates` (L, 37, 3) atom37, `plddt` (L,), `ptm` (float), `pae` (L, L) when
requested. **pLDDT is on a 0–1 scale** — multiply by 100 for the conventional
scale / B-factors.

### `fold_all_atom` (ESMFold2) — billed
Fold a molecular complex. `all_atom_input.sequences` is a list of entities:

```json
{"type": "protein", "id": "A", "sequence": "MKT..."}
{"type": "protein", "id": ["A", "B"], "sequence": "..."}   // homodimer
{"type": "dna", "id": "C", "sequence": "GATC"}
{"type": "rna", "id": "D", "sequence": "GAUC"}
{"type": "ligand", "id": "L", "ccd": ["SAH"]}
{"type": "ligand", "id": "L", "smiles": "CC(=O)O"}
```

Returns `complex` (render with `eb.complex_to_pdb`), `plddt`, `ptm`, and
`interface_ptm` (iPTM, the primary interface-quality metric). Covalent bonds go
in `all_atom_input.covalent_bonds`.

### `generate` (ESM3) — billed
Iterative generation on ONE track. Top-level request:
`{"model", "inputs": {"sequence", "secondary_structure", "sasa", "coordinates"},
"track", "num_steps", "temperature", "schedule", "strategy", ...}`.

- `track` ∈ `sequence | structure | secondary_structure | sasa | function`.
- **Masking:** `'_'` masks a sequence position; `NaN` masks a coordinate
  (the client converts to JSON null); `None` masks a SASA entry.
- `num_steps` must be ≤ sequence length and is capped at 100 by the API.
- Returns `outputs.{sequence, coordinates, plddt, ptm, pae, secondary_structure,
  sasa, function}` per track. `function` is a list of `[label, start, end]` with
  **1-indexed inclusive** ranges; labels may embed an `IPRxxxxxx` accession.
- Predicting `secondary_structure` returns a genuine 8-state DSSP string
  (`GHITEBSC`); predicting `sasa` returns per-residue Å².

### `forward_and_sample` (ESM3) — billed
Single forward pass with sampling. Returns per-track `tokens`, `logprobs`,
`entropy`, `topk_logprobs`, `topk_tokens`. (No ESM skill uses this as its primary
path; it is redundant with the leave-one-out formulation in
`esmc-mutation-effect-scoring`.)

### SAE feature descriptions — no auth
`GET https://biohub.ai/esm/protein/api/v1alpha1/features/{index}` →
`label`, `summary`, `description`, `category`, `activation_pattern`,
`exemplar_protein_families`, `top_swissprot_activations`, `threshold`, and more.
These are auto-generated hypotheses, not curated annotations.

## 4. SAE output format

Over JSON, `sae_outputs[model_name]` is a **sparse top-k** object, not a dense
tensor:

```json
{"feature_indices": [[...], ...],  // (L+2, k) active feature ids per residue
 "values":          [[...], ...],  // (L+2, k) their activations
 "shape":           [L+2, 16384]}  // full (padded) codebook size
```

`eb.sae_features()` returns it trimmed of BOS/EOS as numpy. SAE model names are
`{esmc_model}-sae-layer{L}-k{k}-codebook{C}`; the released one is
`esmc-6b-2024-12-sae-layer60-k64-codebook16384`. `normalize_features=True` is
rejected for any 300M SAE.

## 5. The credit model (read this)

The account has a **hard allowance of 100 credits per day**, resetting at 00:00
UTC. It is a wall, not a rate limit.

- **HTTP 429 means two different things** — transient rate limiting *and*
  daily-credit exhaustion — distinguished only by the response body. The client
  raises `BiohubQuotaError` immediately on the latter (retrying cannot help).
- **Only `logits`, `fold`, `fold_all_atom`, `generate`, `forward_and_sample` are
  billed.** `encode` is free.
- **No usage endpoint** exists; you learn you are out only by hitting the wall.
- **Cost scales with the work:** embedding or folding one sequence ≈ 1 call;
  screening N candidates ≈ N; a deep mutational scan of an L-residue protein ≈ L
  (one masked pass per residue), so a 300-residue scan cannot finish in one day.

To keep evals repeatable, the client supports record/replay cassettes
(`BIOHUB_CASSETTE`, `BIOHUB_CASSETTE_MODE`); replay needs no key and no credits.

## 6. Verified gotchas

| Symptom | Cause / fix |
|---|---|
| `include_distogram=True` → HTTP 422 | Not implemented server-side. Use `interface_ptm`. |
| `/inverse_fold` → "does not support" | No reachable model supports it. Inverse-fold via `generate(track='sequence', coordinates=...)`. |
| pLDDT looks tiny (0.4–0.9) | It is 0–1, not 0–100. Scale by 100 for B-factors. |
| `pair_chains_iptm` is empty | Use the scalar `interface_ptm`. |
| off-by-one in per-residue scores | BOS offset: residue *i* is at token index *i+1*. |
| HTTP 403 on a model | The key lacks access (esm3-medium/large/multimer). |
| `uv` tries to build someone else's project | Missing `--no-project`. |

## 7. Reference measurements

Set the eval thresholds; useful as sanity anchors:

- Ubiquitin folds to pLDDT 0.82 / pTM 0.78.
- ESMC recovers 75/76 of ubiquitin's wild-type residues under masking.
- Inverse folding recovers ubiquitin's native sequence at 100 % identity.
- Barnase + barstar (a true high-affinity complex) gives iPTM 0.96.
- Ubiquitin pseudo-perplexity 1.05 vs 18.70 for its own scramble.
