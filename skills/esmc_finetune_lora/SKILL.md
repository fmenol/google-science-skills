---
name: esmc-finetune-lora
description: >
  Fine-tune ESMC on your own labelled protein data using LoRA (parameter-efficient
  fine-tuning), running on a rented GPU via Modal. Use when you have a labelled
  dataset — sequences plus a class or a measured value — and want a predictor that
  beats zero-shot: enzyme-class or localisation classification, stability or
  activity or expression regression, DMS fitness, any supervised protein property
  task. Also use when asked to "train on my data", "fine-tune ESM", "build a
  protein classifier", or "PEFT/LoRA an ESM model".
  Do not use when you have NO labels — for zero-shot variant effects use
  `esmc-mutation-effect-scoring`, and to pick a frozen layer for your own probe
  use `esmc-embedding-layer-sweep` (far cheaper, no GPU). Do not use to design new
  sequences — that is `esm3-protein-design` or `esmfold2-binder-design`.
---

# Fine-tuning ESMC with LoRA on a Modal GPU

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions.
2.  **User Notification**: If `.licenses/esmc_finetune_lora_LICENSE.txt` does not
    exist in the workspace root then (1) prominently notify the user to check the
    terms at https://biohub.org/acceptable-use-policy/ and https://biohub.ai/,
    then (2) create the file recording the notification text and timestamp.
3.  **Modal**: `pip install modal` and authenticate once with `modal token new`
    (or `modal setup`). A Modal account has a free monthly compute credit that
    covers small runs. Nothing else to install or configure.
4.  **No `BIOHUB_API_KEY` needed.** This skill never calls the hosted API and
    spends **zero Biohub credits**. The weights are public and ungated.

## Overview

This trains a **LoRA adapter plus a fresh prediction head** on a frozen ESMC
backbone. Roughly 0.1–1 % of parameters are trainable, so it fits on the
smallest GPU Modal offers and finishes in minutes.

It is the one ESM capability the hosted API **cannot** provide at all: the API
returns detached JSON with no autograd graph, so there is no way to compute a
gradient through the model remotely. This is why the tutorial it comes from
(`esmc_finetune.ipynb`) was previously out of scope.

The whole skill is one Modal app: `train_remote` runs on the GPU, the
`local_entrypoint` runs on your machine, launches it, and writes the adapter,
metrics and figures to a local directory.

**What it does NOT do.** It does not train the backbone (that is full
fine-tuning). It does not design sequences. It requires labels.

## Core Rules

* **NEVER interpret accuracy alone on an imbalanced dataset.** Read `macro_f1`
  and `classes_never_predicted` from `result.json`. The upstream tutorial's own
  conclusion after 1000 steps was that some classes "are being missed entirely"
  while accuracy looked respectable. If `classes_never_predicted` is non-empty,
  say so in your report.
* **NEVER compute metrics yourself from the predictions.** Use the values in
  `result.json`; they are computed once, in one place.
* **1000 steps is a smoke run, not a result.** Upstream recommends ~5000 steps
  for a meaningful fine-tune. Report the step count alongside any metric.
* Do not request a big GPU. ESMC-300M peaks at ~12 GiB; the default `a100-40` is
  ample and cheaper.
* If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

All commands run from the skill directory, with the Modal CLI.

### Reproduce the upstream tutorial (CARE enzyme classification)

`care:train` / `care:test` download the CARE benchmark (Zenodo record 14004425)
into a Modal volume the first time, then hit the cache forever after.

```bash
modal run scripts/finetune.py \
  --train-csv care:train --test-csv care:test \
  --steps 1000 --out ./results/care-1k
```

### Train on your own labelled data

A local CSV is uploaded to the worker automatically.

```bash
modal run scripts/finetune.py \
  --train-csv mydata_train.csv --test-csv mydata_test.csv \
  --seq-col sequence --label-col localisation \
  --steps 5000 --out ./results/mine
```

### Regression instead of classification

```bash
modal run scripts/finetune.py \
  --train-csv stability.csv --seq-col sequence --label-col ddG \
  --task regression --steps 5000 --out ./results/ddg
```

### Optional: faster Hub downloads

Anonymous HuggingFace pulls are rate-limited per source IP. To lift that, create
a Modal secret once (the token value never touches this repo):

```bash
modal secret create huggingface HF_TOKEN=hf_your_token_here
```

## Choosing a model

| `--model` | Params | When |
|---|---|---|
| `biohub/ESMC-300M` | 333 M | Default. Start here always. |
| `biohub/ESMC-600M` | 600 M | If 300M underfits a large dataset. |
| `biohub/ESMC-6B` | 6.35 B | Only with strong evidence scale is the bottleneck. Set a larger `--gpu` in the app; much slower. |

## Interpreting the Output

`result.json` is the source of truth. Key fields:

| Field | Reading |
|---|---|
| `test_metrics.accuracy` | Top-1. **Insufficient alone on imbalanced data.** |
| `test_metrics.macro_f1` | Treats every class equally. The honest headline for imbalanced sets. |
| `test_metrics.classes_never_predicted` | Non-empty ⇒ the model has collapsed onto majority classes. Report it. |
| `test_metrics.per_class_recall` | Which classes are actually learned. |
| `test_metrics.pearson_r` / `spearman_rho` | Regression. Spearman is the robust one for assay data with a skewed range. |
| `trainable_fraction` | Sanity check that LoRA engaged; should be well under 1 %. |
| `peak_vram_gib` | Right-size the next run. |
| `history[].val` | Whether validation was still improving at the end — if so, train longer. |

Rough guides for classification: macro-F1 **> 0.8** strong, **0.5–0.8** useful,
**< 0.5** weak — and if it is near `1/num_labels`, the model learned nothing.

**There is no upstream number to match.** The tutorial's notebook outputs are
stripped, so it reports no metric. Treat your first run as the baseline; do not
claim parity with a published result that does not exist.

## Common Mistakes

1.  **Reading accuracy on CARE and declaring success.** The label distribution
    is heavily skewed. Always pair it with macro-F1.
2.  **Training for 1000 steps and reporting it as a fine-tune.** That is ~1/6 of
    one epoch on CARE. It is a smoke test.
3.  **Forgetting the test set is a hard split.** The default CARE test file is
    the 30–50 % sequence-identity band — a deliberate generalisation test.
    Scores are lower than a random split and that is the point.
4.  **Assuming the adapter contains the model.** `adapter/` holds only the LoRA
    matrices and the head; reloading it needs the base checkpoint from the Hub.

## Notes from live verification

The dependency stack and the LoRA configuration were verified on a live GPU
(A100-40GB). Two things the upstream notebook gets wrong in a real environment,
handled here:

* **`target_modules=['out_proj']` fails** where ESMC is built on NVIDIA
  Transformer Engine — `attn.out_proj` is a TE `Linear`, not `torch.nn.Linear`,
  and PEFT rejects it. The app detects the build and routes it through
  `target_parameters` instead.
* **`lora_dropout=0.01` is incompatible with `target_parameters`** on peft
  0.17.1; it is coerced to 0 with a note.

Reference run (CARE, 30 steps, A100-40GB): 182 683 train / 1 846 val / 7 classes;
LoRA trainable **1.37 %**; **~12 GiB peak VRAM**; step-0 loss 1.9463 ≈ ln(7).

## Dependencies

* `esm_gpu_common` — the shared Modal image/config and the GPU-skill contract
  (`SPEC_GPU.md`). Read it before modifying this skill.
* `esmc-embedding-layer-sweep` — the cheap, no-GPU alternative when a frozen
  linear probe would do. **Try that first** if your dataset is small.
* `esmc-protein-embeddings` — for producing features without training.
