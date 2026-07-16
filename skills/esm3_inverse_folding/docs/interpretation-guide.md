# Interpreting Inverse-Folding Output

How to read what `scripts/inverse_fold.py` returns, and how to decide whether a
design is real. Read this before you report or rank any design. Do not recompute
the numbers — the script already computed them; interpret them.

The one rule that governs everything below: **a design is only as good as its
self-consistency scores.** ESM3 returns a fluent, confident-looking sequence for
*any* backbone, including one it has completely failed to solve. `design` alone
therefore tells you nothing about whether a sequence works. Always run
`self-consistency`.

--------------------------------------------------------------------------------

## scTM — the primary quality score

scTM (self-consistency TM-score) measures how well a design, re-folded by an
**independent** model (ESMFold2), reproduces the target backbone. Judge every
design on scTM **first**.

| scTM band | Verdict | What to say |
| :-------- | :------ | :---------- |
| **> 0.8** | Adopts the target fold | Report as a success. TM-score > 0.5 already implies the same fold; > 0.8 is a confident match. |
| **0.5 – 0.8** | Uncertain — fold only partially recovered | **Do not present as a success.** Say the fold is partially recovered and the design is unconfirmed. |
| **< 0.5** | Failed design | Say so plainly. The re-fold does not reproduce the target. |

Then break ties and sanity-check with the two secondary metrics:

- **scRMSD (CA, Å):** `< 2 Å` is the practical success bar in the design field;
  `> 5 Å` means the fold was not reproduced. Use it to rank designs that tie on
  scTM (lower is better) — the script already sorts this way.
- **Refold pLDDT (0–1):** confidence of the *re-fold*. `> 0.9` very high,
  `0.7–0.9` confident, `< 0.5` disordered. A design with high scTM but refold
  pLDDT `< 0.7` is a weak result — the fold matched, but the folder was not sure
  of it. (Note pLDDT is on a **0–1** scale here; multiply by 100 for the
  conventional 0–100 scale.)

**scTM saturates — so keep a decoy control in mind.** When a design re-folds onto
the target essentially perfectly (as on an ESM-predicted native backbone), scTM
legitimately hits 1.000. That alone cannot distinguish a working scorer from one
that always returns 1.0. The worked example scores a **shuffled-correspondence
decoy** through the same function and gets **scTM 0.078 / scRMSD 15.3 Å** — the
control that proves the 1.000 is a real match. See
[docs/examples/self_consistency/report.md](examples/self_consistency/report.md).

--------------------------------------------------------------------------------

## Recovery — a benchmark, not a quality score

Native-sequence **recovery** is the fraction of positions where a design matches a
known native sequence. It is only defined when the backbone came from a real
protein; for a de novo backbone there is no native sequence and recovery is
meaningless.

- **Chance is ~5%** (1 in 20 amino acids). The `recovery` plot draws this as a
  grey dotted baseline; an explicit random baseline measured on ubiquitin is
  **4.9%**.
- **Reference points:** ProteinMPNN / ESM-IF report ~50% recovery on native
  backbones. On an ESM-predicted backbone at T=0.1, ESM3 returns the native
  sequence essentially **exactly** — measured **100% on ubiquitin** (76/76
  positions), matching the shared API reference's anchor.
- **A design can be excellent at low recovery.** Recovery measures agreement with
  *one* native sequence, not fold quality. Many sequences fold to the same
  backbone, so a genuinely good de novo design can sit at 30% recovery. **Never
  present recovery as proof a design works — scTM is that score.**

Read the per-position plot for *where* recovery fails: green = matches native,
red = differs. Persistent red at the same positions across designs flags either a
hard-to-place region or a backbone defect there.

--------------------------------------------------------------------------------

## Diversity — resample at a moderate temperature

Each `generate` call is an independent stochastic decode. Diversity comes from
**drawing more samples** (`--num-samples`), not from a hotter decode.

- **mean pairwise identity** (across designs): high (`> 0.95`) means either the
  backbone strongly determines the sequence, or the temperature is too low for the
  diversity you asked for.
- At **very low T** the decodes can collapse to identical sequences: on the
  ubiquitin backbone, T=0.1 gives 1 unique sequence out of 3 (the backbone pins
  it); T=1.0 gives 3 unique at 94.7% identity.
- **Temperature ladder** (measured on ubiquitin's own predicted backbone):

  | T | Native recovery | Use |
  | :-- | :-------------- | :-- |
  | 0.1 | ~100% | Maximum fidelity. The default. |
  | 1.0 | ~98% | Adds real but modest diversity. |
  | 1.5 | ~90% | Aggressive. |
  | 2.0 | ~40% | **Broken** — barely above chance. Do not use. |

Past **T ≈ 1.5** you are destroying the design, not exploring it. To cover more
of the sequence space, raise `--num-samples` at a moderate temperature (≤ ~1.0).

--------------------------------------------------------------------------------

## Platform facts you will hit

- **`/inverse_fold` is unsupported for this API key.** Every reachable model
  (`esm3-open-2024-03`, `esmfold2-*`, `esmc-*`) answers **HTTP 422** *"Model '…'
  does not support the 'inverse_fold' endpoint."* The ESM SDK's `inverse_fold()`
  / `InverseFoldingConfig` path is unreachable — do not reach for it and do not
  "fix" it. Inverse folding goes through `generate(track='sequence',
  coordinates=...)` with the sequence track masked, which is what the script does.
- **A CA-only trace silently produces garbage.** ESM3 builds its structure
  representation from each residue's **N–CA–C frame**. Given only CA atoms the API
  does *not* error — it returns a fluent sequence at chance-level identity
  (measured **0.046** on ubiquitin). The script refuses such input; fix the input,
  do not work around the check.
- **Conditioning is sidechain-blind by default.** ESM3 reads only N/CA/C, so the
  script strips to the N/CA/C/O backbone. Verified against the live API: all-37-
  atom, N/CA/C/O, and N/CA/C conditioning return byte-identical designs at 100%
  recovery — so `--all-atom` is a no-op, and a high recovery number is real
  inverse folding, not the model reading sidechains.
- **pLDDT is 0–1**, not 0–100. Multiply by 100 for B-factors / the conventional
  scale.
- **Runner:** always `uv run --no-project`. Without `--no-project`, `uv` finds an
  unrelated parent `pyproject.toml` and fails.

--------------------------------------------------------------------------------

## Negative results & scientific integrity

- **An unvalidated design is not a result.** If you ran only `design`, you have
  candidate sequences and no evidence they fold. Do not rank, recommend, or hand
  them over. Run `self-consistency` first.
- **Report failures plainly.** scTM `< 0.5` is a failed design — say so. A partial
  recovery (0.5–0.8) is not a success; do not round it up in prose.
- **Do not sell recovery as quality.** For a de novo backbone, recovery is
  undefined; quoting a recovery number there is meaningless. For a native
  backbone, high recovery only means the design matches *that* sequence.
- **Keep the decoy control honest.** When scTM saturates at 1.000, state that the
  scorer was checked against a decoy (or that it should be) — a scorer that always
  returns 1.0 would look identical without it.

## Pre-report checklist

Before writing anything about a design:

- [ ] Did I run `self-consistency`, not just `design`? (If not, stop.)
- [ ] Is the best design's **scTM > 0.8**? If 0.5–0.8, I will call it uncertain;
      if < 0.5, I will call it failed.
- [ ] Is **scRMSD < 2 Å** and refold **pLDDT ≥ 0.7**? If not, I will flag the
      design as weak even at high scTM.
- [ ] If I quote recovery, is there a real native sequence, and did I state that
      recovery is a benchmark, not a quality score?
- [ ] Does the scTM number survive the **decoy control** (or did I note the
      control)?
- [ ] Did I get diversity from `--num-samples` at a moderate T, not from a decode
      hotter than ~1.5?
- [ ] Does the figure referenced in my report actually exist in the folder?
