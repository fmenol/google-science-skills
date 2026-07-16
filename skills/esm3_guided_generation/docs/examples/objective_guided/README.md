# Worked Example: Objective-Guided Generation

This is the canonical worked example for `esm3-guided-generation`. It runs the
two commands that matter — `generate` for a hard, cheap target and `compare`
for the honest guided-vs-unguided test — and reads their real output. The
analysis is in [report.md](report.md).

Both runs are the ones the eval records, on the smallest budgets that still
demonstrate the machinery (22 and 24 API calls). Read the report before you
launch anything larger: guided decoding is the most credit-hungry skill in the
ESM suite, and the cost arithmetic in the report is the thing to internalise.

## The two runs

### 1. `generate --objective no-cysteine` (L=48, 4 steps x 3 samples)

A cheap, sequence-only objective with a deterministically checkable target: a
converged run contains **zero** cysteines. Unguided ESM3 48-mers carry a mean of
1.67 cysteines (measured on this API; 83% have at least one), so hitting zero is
a real result, not a freebie. Figure: `nocys_trajectory.png`.

### 2. `compare --objective hydrophobicity` (L=48, 3 steps x 4 samples)

The honest test. `compare` runs three arms through identical machinery — guided
SVDD, the unguided ablation (the same loop with candidate selection switched
off), and one plain ESM3 generation — and reports Kyte-Doolittle GRAVY for each.
Because the unguided arm is an exact ablation of the *selection* step, any
difference is attributable to guidance alone. GRAVY is stochastic, so the report
also averages over a second seed. Figure: `gravy_compare_seed0.png`.

## Reproduce (free, no credits, no API key)

Every payload for these inputs is in the recorded cassettes. Replay them:

```bash
export BIOHUB_CASSETTE=$PWD/evals/cassettes
export BIOHUB_CASSETTE_MODE=replay          # a cache miss errors; never bills
SCRIPT=skills/esm3_guided_generation/scripts/guided_design.py

# Run 1 — no-cysteine. 22 API calls when run for real.
BIOHUB_CASSETTE_NS=nocys uv run --no-project "$SCRIPT" generate \
  --length 48 --objective no-cysteine \
  --num-decoding-steps 4 --num-samples-per-step 3 \
  --seed 0 --output nocys.json

# Run 2 — hydrophobicity, guided vs unguided. 24 API calls each seed.
BIOHUB_CASSETTE_NS=cmp_gravy_seed0 uv run --no-project "$SCRIPT" compare \
  --length 48 --objective hydrophobicity \
  --num-decoding-steps 3 --num-samples-per-step 4 \
  --seed 0 --output cmp_gravy_seed0.json
BIOHUB_CASSETTE_NS=cmp_gravy_seed1 uv run --no-project "$SCRIPT" compare \
  --length 48 --objective hydrophobicity \
  --num-decoding-steps 3 --num-samples-per-step 4 \
  --seed 1 --output cmp_gravy_seed1.json
```

Notes:

-   `BIOHUB_CASSETTE_NS` isolates each invocation's recording. The guided loop is
    stochastic and re-issues near-identical payloads across steps, so without a
    per-run namespace the replays would cross-talk. Use exactly the namespaces
    above (`nocys`, `cmp_gravy_seed0`, `cmp_gravy_seed1`).
-   `--no-project` is mandatory: without it `uv` walks up the tree, finds an
    unrelated `pyproject.toml`, and tries to build that project.
-   The figures in this folder are copies of the `.png` written next to each run's
    `--output` JSON (`nocys.png`, `cmp_gravy_seed0.png`).
