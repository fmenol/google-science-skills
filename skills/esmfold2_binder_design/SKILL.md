---
name: esmfold2-binder-design
description: >
  Design de novo protein binders (minibinders) against a target, by
  gradient-based optimisation through the ESMFold2 folding trunk on a rented GPU
  via Modal. Use when asked to "design a binder", "design a minibinder", "make
  something that binds X", "de novo binder design", or to reproduce the ESM
  binder-design protocol. Produces candidate sequences ranked by multi-critic
  interface confidence.
  Do not use to RANK binders you already have — that is
  `esmfold2-binder-screening`, which needs no GPU and runs in seconds. Do not
  use for general protein design without a target — that is `esm3-protein-design`.
---

# De novo binder design with ESMFold2 on a Modal GPU

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions.
2.  **User Notification**: If `.licenses/esmfold2_binder_design_LICENSE.txt`
    does not exist in the workspace root then (1) prominently notify the user to
    check the terms at https://biohub.org/acceptable-use-policy/ and
    https://biohub.ai/, then (2) create the file recording the notification text
    and timestamp. **This skill designs novel binding proteins — the acceptable
    use policy is not a formality here.**
3.  **Modal**: `pip install modal` and `modal token new`. This skill uses an
    80 GB-class GPU (`a100-80` by default); check your Modal plan covers it.
4.  **No `BIOHUB_API_KEY` needed.** Zero Biohub credits; the weights are public.

## Overview

This implements Algorithms 11–15 of the ESM world-model paper: a **150-step
gradient descent on a soft sequence-logit tensor**, backpropagated through the
ESMFold2 trunk, then a **multi-critic scoring stage** that folds each design with
several independent checkpoints and ranks by interface confidence.

The optimised objective is built entirely from the **distogram**:

| Loss | Weight | What it asks for |
|---|---|---|
| `inter_contact` | 0.5 | every target residue has ≥1 confident binder contact within 22 Å |
| `intra_contact` | 0.5 | the binder makes its own internal contacts (a real fold) |
| `glob` | 0.2 | radius of gyration near the compact-globular law `2.38·n^0.365` |
| `plm` (ESMC-6B pseudo-perplexity) | 0.15 | the sequence looks like a natural protein |

Cysteine is banned at initialisation and masked from all gradients.

**Why this needs a GPU and the API cannot do it.** The optimisation requires
gradients through the folding trunk. The hosted API returns detached JSON with no
autograd graph, and `include_distogram=True` returns HTTP 422 — the distogram
these losses are computed from is not even exposed.

**The app runs BioHub's own implementation**, fetched at a pinned commit
(`67838dc8`) rather than reimplemented, because the gradient handling is subtle
enough that a reimplementation would silently diverge from the published
protocol. BioHub's `binder_design.py` is itself a Modal script; its own
`modal.App` would collide with this one, so the GPU worker fetches the source,
strips the Modal layer, and drives the plain `ESMFold2Design` class directly.

The skill is one Modal app: a `Designer` class loads the models once per GPU
container (`@modal.enter`) and designs per seed (`@modal.method`); the
`local_entrypoint` fans seeds out with `.map` and pools the results.

## Core Rules

* **NEVER run a single seed and report the result as a designed binder.** This is
  a low-yield search. Upstream's own campaign filters thousands of candidates.
  `--num-seeds 1` is a smoke test.
* **NEVER rank by a single critic's iPTM.** Use `mean_iptm` across the hero
  critics. The multi-critic ensemble exists precisely to stop one over-confident
  checkpoint from crowning a bad design.
* **ALWAYS report a negative campaign as negative.** If the best `mean_iptm` is
  below 0.5 there is no credible binder, regardless of what rank 1 says. The
  entrypoint prints this explicitly — do not omit it.
* **NEVER present a design as validated.** iPTM is a model's confidence in its own
  prediction, not evidence of binding. Every output is a hypothesis requiring
  wet-lab testing.
* If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

### Design minibinders against a preset target

```bash
modal run scripts/design.py \
  --target-name pd-l1 --binder-name minibinder \
  --num-seeds 8 --out ./designs
```

Presets: `cd45`, `ctla4`, `egfr`, `pd-l1`, `pdgfr`. (Upstream's notebook says
`pdgfrb`; the working key is `pdgfr`.)

### Design against your own target

```bash
modal run scripts/design.py \
  --target-name my_target --target-sequence MKTAYIAKQRQISFVK... \
  --binder-name minibinder --num-seeds 16 --out ./designs
```

### Poll and collect

Modal streams logs live, and `.map` returns when all seeds finish. Results land
in `--out`: `ranked_designs.json` (pooled across seeds), `designs.fasta`, and the
predicted complex PDBs under `structures/`.

## Interpreting the Output

`mean_iptm` is the headline:

| mean iPTM | Reading |
|---|---|
| **> 0.8** | Confident interface. A genuine candidate to order. |
| **0.5–0.8** | Possible interface. Worth a larger panel, not on its own. |
| **< 0.5** | Likely no binding. Do not report as a hit. |

Also read:

* **`min_iptm` vs `max_iptm`** — wide spread means the critics disagree and the
  mean is hiding it. A design at mean 0.7 with min 0.3 is much weaker than one at
  mean 0.7 with min 0.65.
* **`trajectory`** — `inter_contact_loss` should fall over the 150 steps. A flat
  trajectory means the optimisation never engaged the target.
* **`peak_vram_gib`** — right-size the next run.

**Expected yield is low.** Most seeds produce nothing above 0.5. That is normal
for this method; "8 seeds, best mean iPTM 0.42, no credible binder" is a correct
and useful result.

### Measured on a live GPU

A minibinder designed against the `pd-l1` preset, one seed, A100-80GB:

| | |
|---|---|
| Model load | 275 s (ESMC-6B + 2 inversion + 4 hero critics) |
| 150 optimisation steps | 722 s, **19.6 GiB peak VRAM** |
| Loss | total **8.00 → 3.32**, inter-contact 5.16 → 2.28 |
| Result | **mean iPTM 0.898** (critics: 0.736 / 0.938 / 0.949 / 0.969) |
| Binder | 158 residues, **0 cysteines** (the protocol bans them) |

The 19.6 GiB is below the 27–51 GB upstream quotes for its large CD45 +
antibody-scFv config; a minibinder against a ~115-residue target is cheaper. The
critic spread (0.736 → 0.969) is exactly why the mean across an ensemble is the
ranking key rather than any single score.

## Limitations

* **Antibody / nanobody design is available but unverified.** The Modal image
  ships ANARCI + HMMER, so the CDR-annotation path works, and the CLI will accept
  `--is-antibody true` with an antibody scaffold — but only minibinder design has
  been checked end-to-end here. Treat antibody results with extra caution.
* **Results are not bitwise reproducible across GPU models.** LM dropout is
  forced on during inference and `torch.compile` is applied, so the same seed on
  different Modal GPU types can give different sequences. The card is recorded in
  the result.

## Common Mistakes

1.  **Treating rank 1 as a hit.** Rank is relative to the batch; an entire batch
    can be bad. Read the absolute iPTM.
2.  **Running one seed.** See Core Rules.
3.  **Passing `--target-sequence` with a preset name.** Upstream raises a
    `ValueError`; this CLI catches it before launch.

## Dependencies

* `esm_gpu_common` — shared Modal image/config and `SPEC_GPU.md`.
* `esmfold2-binder-screening` — the API-based scoring half. **Use it instead** if
  you already have candidates and only need them ranked.
* `esmfold2-structure-prediction` — to fold a design on its own.
* `pymol` — to visualise the complex PDBs this skill writes.
