# ESM GPU Skills — Implementation Contract

**This is the sibling of [`../esm_common/SPEC.md`](../esm_common/SPEC.md), not a
replacement for it.** That document governs the twelve skills that call the
hosted Biohub API. This one governs the three skills that run **open ESM weights
on GPUs rented by the second through [Modal](https://modal.com)**.

Read `SPEC.md` first. Where this document is silent, `SPEC.md` still applies
(Google Python style, file output not stdout, `SKILL.md` structure).

---

## 1. Why a second class of skill exists

The API skills are bounded by what `https://biohub.ai` will do for you. Three of
those bounds are hard, and each blocks a cookbook tutorial:

| Bound | Consequence |
|---|---|
| **No gradients.** The API returns detached JSON. | `esmc_finetune.ipynb` (LoRA) and the optimisation half of `binder_design.py` are unreachable. |
| **100 credits/day.** A wall, not a rate limit. | `gfp_design.ipynb` warns some prompts "may require **thousands** of generations". At 1 credit each that is a 30-day campaign. |
| **Gated model tiers.** `esm3-medium/large/multimer` → HTTP 403. | Two tutorials are written against models the key cannot reach. |

Local weights on rented GPUs lift the first two completely. **They do not lift
the third** — see § 6, because getting this wrong is the most tempting mistake.

---

## 2. Why Modal, and the shape of a Modal skill

Modal runs Python functions on cloud GPUs you pay for by the second. There is no
cluster to stand up, no Kubernetes, no persistent infrastructure — a job spins
up a GPU, runs, returns its result, and the GPU is released. That makes these
skills portable: anyone with a Modal account can run them, with no site-specific
configuration.

Modal also draws the local/remote boundary for you, so there is **no
control-plane / payload split to maintain by hand**. One file per skill:

```
skills/<gpu_skill>/
  SKILL.md
  scripts/
    <name>.py       <- the Modal app: local entrypoint + GPU function/class
    esm_modal.py    <- vendored, byte-identical to esm_gpu_common/esm_modal.py
  references/citation.bib
```

Inside that one file:

* Code in an `@app.function(...)` or `@app.cls` runs **remotely on the GPU**,
  inside the container image. This is where torch, transformers and esm live.
* Code in `@app.local_entrypoint()` runs **locally** — it is the CLI. It touches
  only stdlib + the `modal` client, calls `.remote()` / `.map()` on the GPU
  function, and writes the returned results to a local `--out` directory.
* Heavy imports sit inside `with IMAGE.imports():` or inside function bodies, so
  importing the app **locally** (to launch it) never needs torch on your laptop.

Invoke with the Modal CLI, never `python`:

```bash
modal run scripts/<name>.py --arg value
```

---

## 3. The single most important rule

`evals/test_offline.py` enforces it: **no `torch` / `transformers` / `peft` /
`esm` / `numpy` / `pandas` / `matplotlib` import at module top level.** They may
only be imported inside a function body or a `with IMAGE.imports():` block, both
of which run remotely. The test parses each app's AST and fails on any top-level
heavy import, and it also imports every app under `uv run --with modal` to prove
it builds its image and constructs its functions with nothing but stdlib +
modal available.

Two consequences worth internalising:

* A Modal `@app.cls` with a `modal.parameter(...)` **cannot** use
  `from __future__ import annotations` — the future import turns the parameter's
  type annotation into a string, which Modal cannot decode. The binder-design
  app omits it deliberately; there is a comment saying so.
* Anything the GPU function returns must be **JSON/pickle-safe**. Structures come
  back as PDB strings, adapters and figures as base64 blobs, tensors as floats.
  The local entrypoint decodes them to files.

---

## 4. Verified facts about the models (do not re-derive)

### The transformers fork is mandatory

Checked against live PyPI, not assumed:

```
transformers 4.57.1  -> models.esmc: absent, models.esmfold2: absent
transformers 5.14.1  -> models.esmc: absent, models.esmfold2: absent
```

Stock transformers ships only legacy `esm` (ESM-2). The `biohub/*` Hub repos
carry no `modeling_*.py` and no `auto_map`, so `trust_remote_code` cannot
substitute. `esm_modal.py` installs `git+https://github.com/Biohub/transformers`
pinned to commit `ef32577f`, not `@main`. Pin commits — a designed protein
someone orders from a synthesis vendor must be reproducible.

### esm and the pinned fork conflict under one pip resolution

`esm`'s pyproject declares `transformers @ ...@main`. Asking pip for that AND a
commit-pinned transformers in one resolution is two direct URLs for one package,
and it fails with `ResolutionImpossible` (observed on a live GPU). So every image
installs `esm` first (dragging `transformers@main`), then force-reinstalls the
pinned commit with `--no-deps`. `esm_modal._pin_transformers()` does this.

### The weights are public and ungated

Probed against the Hub API. Every checkpoint used here returned `gated=false`:
`biohub/ESMC-300M/600M/6B`, `biohub/esm3-sm-open-v1`, `biohub/ESMFold2*`, and the
released SAEs. **An `HF_TOKEN` is not required.** A Modal Secret named
`huggingface` is attached when it exists (anonymous Hub pulls are IP-rate-limited),
referenced by name only — no token value appears in any file.

### VRAM, measured or quoted from the source

| Workload | VRAM | `--gpu` |
|---|---|---|
| ESMC-300M LoRA, batch 8, len 1024 | ~12 GB (measured) | `a100-40` |
| ESMC-600M LoRA | ~15 GB | `a100-40` |
| Binder design, minibinder vs ~115-aa target | **19.6 GB (measured)** | `a100-80` |
| Binder design, large antibody config | 27–51 GB (upstream) | `a100-80` / `h100` |
| ESM3-open design campaign | < 16 GB | `a100-40` (even `t4` works) |

There is no VRAM guard in the upstream folding code — it OOMs silently — so every
payload logs `torch.cuda.max_memory_allocated()`.

---

## 5. Skill contract

Every GPU skill:

1. **Is one Modal app file** with a GPU `@app.function`/`@app.cls` and an
   `@app.local_entrypoint()`. No heavy import at module top level.
2. **Defaults `--gpu` to the smallest flavor that fits** (§ 4). A bigger card is
   not faster for these jobs.
3. **Returns JSON/pickle-safe results**; the entrypoint writes them locally.
4. **Records provenance**: the GPU model, seed, and peak VRAM. ESMFold2 and ESM3
   sampling are not seedable across GPU architectures, so "which card" is part
   of the result.
5. **Sets `matplotlib.use('Agg')` before pyplot** and never calls `plt.show()`.
6. **Names its Modal Volumes and Secret by stable, generic names**
   (`esm-model-cache`, `esm-data-cache`, `huggingface`) — nothing tied to a
   person or an organisation.

---

## 6. What Modal does NOT unlock — read before promising anything

**`esm3-medium`, `esm3-large` and `esm3-*-multimer` remain unreachable.** The
only public ESM3 weights anywhere are the 1.4 B `esm3-sm-open-v1` — checked
directly against the Hub. **The 403 is a licensing decision, not a
serving-capacity one, and renting GPUs does not route around it.**

Two tutorials target those tiers:

* `gfp_design.ipynb` uses the 7 B variant. The campaign skill runs the identical
  protocol at 1.4 B, where pass rates through the two RMSD gates are materially
  lower. Volume is the compensation, not parity — the `SKILL.md` says so.
* `fold_invfold.py` hard-wires `esm3-medium-2024-08`; the capability is already
  covered by the API skill `esm3-inverse-folding` via
  `generate(track='sequence', coordinates=…)`.

Overpromising here is the single easiest way to make these skills dishonest.
