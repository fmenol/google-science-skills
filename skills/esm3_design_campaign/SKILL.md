---
name: esm3-design-campaign
description: >
  Run a large ESM3 motif-scaffolding campaign on rented GPUs via Modal — hundreds
  to thousands of generations with structural rejection sampling — to find designs
  that scaffold a functional motif into a genuinely novel fold. Use when asked to
  "design a novel GFP", "scaffold this active site", "reproduce the gfp_design
  protocol", or whenever a design task needs far more samples than a metered API
  allows. Reports campaign-wide yield and per-gate rejection statistics.
  Do not use for a handful of designs against the hosted API — that is
  `esm3-protein-design`, which needs no GPU. Do not use for binder design against
  a target protein — that is `esmfold2-binder-design`.
---

# High-throughput ESM3 design campaigns on Modal GPUs

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions.
2.  **User Notification**: If `.licenses/esm3_design_campaign_LICENSE.txt` does
    not exist in the workspace root then (1) prominently notify the user to check
    the terms at https://biohub.org/acceptable-use-policy/ and https://biohub.ai/,
    then (2) create the file recording the notification text and timestamp.
3.  **Modal**: `pip install modal` and `modal token new`. ESM3-open is 1.4 B and
    the default GPU is `a100-40`; even a `t4` runs it.
4.  **No `BIOHUB_API_KEY` needed.** Zero Biohub credits, unmetered inference.

## Overview

This generalises `gfp_design.ipynb` into a sharded campaign. The protocol per
attempt:

1. Build a prompt: masked everywhere except a few pinned motif residues, plus
   optional template structure tokens.
2. **Generate a structure.**
3. **Gate 1** — keep only if the constrained site is faithful (`site RMSD < 1.5
   Å`) **and** the backbone is genuinely different (`backbone RMSD > 1.5 Å`). The
   second half is the point: a design that merely reproduces the template is not
   a novel protein.
4. **Generate a sequence** for the accepted structure.
5. **Refold deterministically** (`num_steps=1, temperature=0`).
6. **Gate 2** — the constrained site must have survived; backbone novelty is now
   don't-care. The motif residues must also literally be present.

**Why this needs rented GPUs.** Each generation is a billed API call and the
Biohub account allows 100 per day. The notebook warns some prompts "may require
**thousands** of generations". A 1000-generation campaign is a 10-day exercise
against the API and minutes on Modal. The skill is one app: `campaign_shard` runs
on a GPU, and the `local_entrypoint` fans shards out with `.map` and pools them.

## A limitation you must state in any report

`gfp_design.ipynb` targets **`esm3-medium` (7B)**, and the paper's real campaign
used the 7B variant. **No ESM3 above the 1.4 B open model is published anywhere** —
checked directly against the HuggingFace Hub. The 403 on
`esm3-medium`/`large`/`multimer` is a licensing decision; renting GPUs does not
route around it.

So this skill runs the protocol at **1.4 B**. The interface is identical and the
protocol unmodified, but **pass rates through the two RMSD gates are materially
lower than the notebook's**. The honest compensation is volume. Do not present
results here as reproducing the published esmGFP campaign.

## Core Rules

* **NEVER report a pass rate without the denominator.** "12 passed" is
  meaningless; "12/2000 (0.6 %)" is a result.
* **ALWAYS read `gate_breakdown`.** It says which constraint is binding — the
  only actionable diagnostic when a campaign yields nothing:
  * mostly `structure_site_rmsd` → the motif is too hard to scaffold at 1.4 B.
  * mostly `structure_too_similar_to_template` → the prompt over-constrains the
    fold; pin fewer structure positions.
  * mostly `motif_not_preserved` → generation overwrites pinned positions; lower
    `--temperature`.
* **A zero-yield campaign is a legitimate result.** Report it, with the best
  achieved RMSD. Do not loosen thresholds after the fact and present the
  survivors as hits.
* **NEVER call a design validated.** These are hypotheses.
* If this skill is used, ensure this is mentioned in the output.

## Utility Scripts

### Reproduce the GFP protocol at campaign scale

```bash
modal run scripts/campaign.py \
  --preset gfp --num-generations 2000 --num-shards 8 --out ./gfp_campaign
```

The preset supplies the 1qy3 template, the chromophore motif
(`59:T,62:T,63:Y,64:G,93:R,219:E`), the pinned structure positions
(`55-69,93,219`) and length 229. **All positions are 0-indexed RESIDUE
positions**; the app applies the BOS `+1` shift internally. (The notebook's
`[56:71],[94],[220]` are already-shifted TOKEN indices — copying them literally
double-shifts and leaves motif residues 93 and 219 unpinned. A guard warns if a
scored motif residue is not pinned.)

### Scaffold your own motif

```bash
modal run scripts/campaign.py \
  --motif "45:H,49:D,132:S" --length 180 \
  --template-pdb 1abc --template-chain A --structure-positions "40-55" \
  --num-generations 1000 --num-shards 4 --out ./my_campaign
```

`--num-generations` is the whole-campaign total; it is split across
`--num-shards` GPU containers, each with its own seed.

## Interpreting the Output

| Field | Reading |
|---|---|
| `pass_rate` | Campaign yield. Single-digit percent is normal; sub-1 % common at 1.4 B. |
| `gate_breakdown` | Where attempts died. The primary diagnostic. |
| `rmsd_distributions.constrained_site_rmsd.best` | **Read this whenever yield is 0.** It distinguishes "just missed" from "nowhere near", which imply opposite next steps. |
| `refold_site_rmsd` | **< 1.5 Å** the motif survived design and refolding. The ranking key. |
| `backbone_rmsd` | **> 1.5 Å** the fold is novel relative to the template. |
| `identity_to_template` | Low identity plus a preserved site is the interesting outcome. |

### Measured on a live GPU — size your run from this

The GFP preset was run at three scales on A100-80GB, ~2 s per generation:

| Generations | Best constrained-site RMSD | Passed (gate: < 1.5 Å) |
|---|---|---|
| 20 | 3.08 Å | 0 |
| 150 | 2.58 Å | 0 |
| **600** | **1.54 Å** | **0** — missed by 0.04 Å |

Two things follow:

1. **Do not conclude the method fails from a small run.** At 600 generations the
   campaign is *at* the threshold and still improving; designs should begin
   passing in the low thousands — exactly what the notebook warns and what a
   100-credit/day API cap makes impossible.
2. **A zero-yield run is still informative** if you report the best RMSD. "0/600,
   best 1.54 Å" and "0/600, best 8 Å" are completely different findings.

### Why the site gate is hard (and not a bug)

Structure-token conditioning works as advertised — measured **17/17 pinned tokens
preserved** in the generated output. But ESM3 structure tokens describe **local**
backbone geometry, not global position, so a perfectly preserved token can still
sit ~3.8 Å from where the template puts it. The gate is a genuine test of whether
the global fold places the motif correctly, which is what gets harder at 1.4 B.
Do not "fix" a low pass rate by loosening the threshold.

## Common Mistakes

1.  **Running 50 generations and concluding the method does not work.** At a 1 %
    pass rate the expected yield is 0.5 designs. Scale first.
2.  **Loosening `--site-rmsd-max` until something passes.** The 1.5 Å gate is from
    the published protocol; changing it is a legitimate experiment but must be
    reported as a changed threshold, not a better result.
3.  **Claiming parity with esmGFP.** See the limitation above.
4.  **Ignoring `motif_not_preserved` failures.** They mean the design lacks the
    functional residues, whatever the RMSD says.

## Dependencies

* `esm_gpu_common` — shared Modal image/config and `SPEC_GPU.md`.
* `esm3-protein-design` — the API-based single-design sibling. **Use it first** to
  check a prompt is sane before spending GPU time on a campaign.
* `esm-protein-tracks` — to build and inspect prompt tracks from a PDB.
* `esmfold2-structure-prediction` — to independently fold a passing design.
* `pymol` — to visualise designs against the template.
