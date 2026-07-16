# ESM Skills — Implementation Contract

**Every ESM skill MUST follow this document exactly.** It encodes hard-won facts
about the Biohub API that were verified against the live service. Deviating will
produce a skill that fails its eval.

---

## 1. Non-negotiables

1. **NEVER download model weights.** No `huggingface_hub`, no `transformers`, no
   `torch`, no `esm` PyPI package. Everything runs remotely on the Biohub
   Platform. A skill that imports torch is a failed skill.
2. **Use the vendored `esm_biohub.py`** in your skill's `scripts/` directory.
   Do not reimplement HTTP, auth, rate limiting, PDB writing, or numerics.
   It is a byte-identical copy of `skills/esm_common/esm_biohub.py`.
3. **`uv run --no-project` only.** Every script carries a PEP 723 header. Never
   `pip install`, never bare `python`/`python3`.
   **`--no-project` is mandatory, not optional.** Without it `uv` walks up the
   directory tree, finds whatever `pyproject.toml` the user happens to be nested
   under, and tries to build *their* project — which fails. Verified: running
   from this repo picks up an unrelated project two levels up and errors out.
   Every command in your SKILL.md and every example MUST read:

   ```bash
   uv run --no-project scripts/<name>.py <subcommand> ...
   ```
4. **File output, not stdout.** Stdout gets a short status line. Payloads go to
   `--output`. Large arrays go to `.npy`/`.csv`, never stdout.
5. **Required args, no silent defaults** for anything that caps results
   (`--limit`, `--num-samples`). Forces the agent to choose explicitly.
6. **Google Python style**: 2-space indents, Apache 2.0 header, `snake_case`.

---

## 2. Script skeleton

```python
# Copyright 2026 Google LLC
# ... (Apache 2.0 header, copy from esm_biohub.py) ...

"""One-line summary of what this script does."""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "polite-http",
#   "python-dotenv",
#   "numpy",
# ]
# ///

from __future__ import annotations

import argparse
import sys

from esm_biohub import BiohubClient, BiohubError   # vendored, same directory
import esm_biohub as eb


def cmd_something(args) -> None:
  client = BiohubClient()
  ...
  eb.write_json(result, args.output)


def main() -> None:
  parser = argparse.ArgumentParser(description='...')
  sub = parser.add_subparsers(dest='command', required=True)

  p = sub.add_parser('something', help='...')
  p.add_argument('--sequence', required=True, help='...')
  p.add_argument('--output', required=True, help='Output JSON file path')
  p.set_defaults(func=cmd_something)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
```

Add `matplotlib`, `scikit-learn`, `pandas` to the PEP 723 `dependencies` list
only in the scripts that actually need them.

---

## 3. Verified API facts (do not re-derive)

### Models this API key can reach
| Model | Notes |
|---|---|
| `esmc-300m-2024-12` | 30 layers, dim 960 |
| `esmc-600m-2024-12` | 36 layers, dim 1152 — **default ESMC** |
| `esmc-6b-2024-12` | 80 layers, dim 2560. `ith_hidden_layer=-1` is **rejected** |
| `esm3-open-2024-03` | the only ESM3 available (alias `esm3-sm-open-v1`) |
| `esmfold2-fast-2026-05` | **default folding model**, ~2 s for 76 aa |
| `esmfold2-2026-05` | slower, accepts MSAs |

`esm3-medium-*`, `esm3-large-*`, `esm3-*-multimer` return **HTTP 403** — the key
has no access. Tutorials that use them (`gfp_design`, `esm3_guided_generation`)
must be re-pointed at `esm3-open-2024-03`.

### Shapes and conventions
* Tokens are `L+2` (BOS … EOS). **Residue `i` lives at token index `i+1`.**
* `logits.sequence` → `(L+2, 64)`. The vocab is zero-padded to 64; always index
  by `eb.VOCAB[aa]`, never positionally. `eb.AA20_IDX` gives the 20 canonical
  columns in alphabetical order matching `eb.AA20`.
* `embeddings` → `(1, L+2, D)`; `mean_embedding` → `(1, 1, D)`.
* `hidden_states` → `(n_layers+1, 1, L+2, D)`; `mean_hidden_state` →
  `(1, n_layers+1, D)`. **Layer 0 is the embedding layer**, so a 36-layer model
  yields 37 rows.
* Mean-pooled tensors are pooled **including BOS/EOS** (server-side). Per-residue
  tensors must be trimmed `[1:-1]` yourself. `eb.embed()` already does this.
* **pLDDT is on a 0–1 scale**, not 0–100. Multiply by 100 for B-factors.
* `fold` → `coordinates (L, 37, 3)` atom37, `plddt (L,)`, `ptm` float,
  `pae (L, L)`.
* `fold_all_atom` → `complex` dict; use `eb.complex_to_pdb()`.
* The mask character is `'_'`. In `ESMProtein`-style prompts, `None` masks a SASA
  entry and `NaN` masks a coordinate.

### Traps that will silently corrupt results
* `include_distogram=True` → **HTTP 422, not implemented server-side.** The
  distogram-based losses from `binder_design.py` are unreachable. Use `iptm`.
* `pair_chains_iptm` comes back empty; use the scalar `interface_ptm`.
* The `inverse_fold` endpoint is **not supported by any reachable model**.
  Inverse folding must go through `generate(track='sequence', coordinates=...)`.
  (Verified: this recovers ubiquitin's native sequence at 100% identity.)
* `normalize_features=True` is rejected for ESMC **300M** SAE models.
* `num_steps` must be ≤ sequence length and is capped at 100 by the API.
* Errors are returned as HTTP status + a JSON `message`. `esm_biohub` turns these
  into `BiohubError` with the server's text plus a hint. Let them propagate.

### SAE
* Model name format: `{esmc_model}-sae-layer{L}-k{k}-codebook{C}`.
  The released one is `esmc-6b-2024-12-sae-layer60-k64-codebook16384`.
* Over JSON, `sae_outputs[name]` is **sparse top-k**:
  `{'feature_indices': (L+2, 64), 'values': (L+2, 64), 'shape': [L+2, 16384]}`.
  `eb.sae_features()` returns it trimmed, as numpy. No torch needed.
* Feature descriptions: `GET https://biohub.ai/esm/protein/api/v1alpha1/features/{idx}`
  (no auth required). Returns `label`, `summary`, `description`, `category`,
  `activation_pattern`, `exemplar_protein_families`, `top_swissprot_activations`,
  `threshold`, and more.

---

## 4. `esm_biohub` API you will use

```python
client = eb.BiohubClient()            # loads BIOHUB_API_KEY, rate-limits, retries

# ESMC
client.encode(seq, model)                        -> list[int]           (L+2)
client.embed(seq, model, layer=None, per_residue=False)  -> np.ndarray
client.mean_hidden_states(seq, model)            -> (n_layers+1, D)     1 request
client.sequence_logits(seq, model)               -> (L+2, 64)
client.sae_features(seq, sae_model)              -> (idx (L,k), val (L,k), C)
client.logits(tokens, model, **flags)            -> raw dict (escape hatch)

# ESMFold2
client.fold(seq, model, num_loops=, num_sampling_steps=, include_pae=)  -> dict
client.fold_all_atom([{'type':'protein','id':'A','sequence':...}, ...]) -> dict

# ESM3
client.generate(track, model, sequence=, coordinates=, secondary_structure=,
                sasa=, num_steps=, temperature=)                        -> dict
client.forward_and_sample(tokens, model, ...)                           -> dict

# helpers
eb.validate_sequence(s)          eb.read_fasta(path)      eb.write_json(obj, path)
eb.to_array(nested_json)         eb.log_softmax(x)        eb.softmax(x)
eb.shannon_entropy_bits(p)       eb.atom37_to_pdb(coords, seq, plddt)
eb.complex_to_pdb(complex_dict)  eb.parse_pdb_atom37(path, chain)
eb.kabsch_rmsd(a, b)             eb.tm_score(ca_a, ca_b)
eb.AA20  eb.AA20_IDX  eb.VOCAB  eb.MASK_CHAR  eb.ATOM37
eb.DEFAULT_ESMC  eb.DEFAULT_ESM3  eb.DEFAULT_ESMFOLD2  eb.DEFAULT_SAE_MODEL
```

`BiohubClient` is safe to share across threads for read-only calls; for batch
work use `concurrent.futures.ThreadPoolExecutor` (the rate limiter is
cross-process, so concurrency is already bounded). Keep `max_workers <= 8`.

---

## 5. SKILL.md structure (mandatory)

```markdown
---
name: {skill-name-with-hyphens}
description: >
  {What it does. When to use it — include the words a user would actually say.
  Then an explicit "Do not use when ..." clause naming the sibling skill that
  *should* be used instead.}
---

# {Title}

## Prerequisites
1.  **`uv`**: Read the `uv` skill and follow its Setup instructions.
2.  **User Notification**: If .licenses/{skill_dir}_LICENSE.txt does not exist in
    the workspace root then (1) prominently notify the user to check the terms at
    https://biohub.org/acceptable-use-policy/ and https://biohub.ai/, then
    (2) create the file recording the notification text and timestamp.
3.  **`BIOHUB_API_KEY`**: Required. Register at
    https://biohub.ai/developer-console/api-keys. You **MUST** use the safe
    credentials protocol in the `credentials` skill to check for it.

## Overview
{What it does, and explicitly what it does NOT do.}

## Core Rules
{Imperatives. NEVER/ALWAYS. Guard the known failure modes.}

## Utility Scripts
{Every subcommand, with a runnable example.}

## Interpreting the Output
{Teach the agent to read the numbers — thresholds, what is significant.}

## Common Mistakes
{2-4 real pitfalls.}
```

Rules to bake into every SKILL.md:
* "**Do not compute X yourself; always use the script's output.**" The agent must
  not do numeric reasoning it can get wrong.
* "If this skill is used, ensure this is mentioned in the output."
* Cross-reference sibling skills by name instead of reimplementing them.
* Include a **Dependencies** section naming sibling skills you rely on.

## 6. Interpretation thresholds (use these; they were measured)

| Quantity | Reading |
|---|---|
| pLDDT (0–1) | > 0.9 very high · 0.7–0.9 confident · 0.5–0.7 low · < 0.5 disordered |
| pTM | > 0.8 confident fold · 0.5–0.8 plausible · < 0.5 unreliable |
| iPTM | > 0.8 confident interface · 0.5–0.8 possible · < 0.5 likely no binding |
| LLR | < 0 predicted deleterious; more negative = worse |
| Entropy (bits) | low = evolutionarily constrained; high = tolerant |

Reference points measured on this API: ubiquitin folds to pLDDT 0.82 / pTM 0.78;
barnase+barstar (a true complex) gives iPTM 0.96.
