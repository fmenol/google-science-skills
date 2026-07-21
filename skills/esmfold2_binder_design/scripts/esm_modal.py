# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Shared Modal building blocks for the ESM GPU skills.

This module is VENDORED into every ESM **GPU** skill's `scripts/` directory,
exactly as `esm_biohub.py` is vendored into every ESM **API** skill.
`skills/esm_gpu_common/esm_modal.py` is canonical; the copies are kept
byte-identical (`evals/sync_common.py --check`).

Why the GPU skills exist
------------------------
The API skills are bounded by what `https://biohub.ai` will do for you: no
gradients, no fine-tuning, 100 credits/day. The GPU skills lift that bound by
running the **open-weight** ESM models on GPUs you rent by the second through
[Modal](https://modal.com). No cluster, no Kubernetes, no persistent
infrastructure -- a design job spins up an H100, runs, returns its result, and
the GPU is released.

How a Modal skill is shaped
---------------------------
Modal draws the local/remote boundary for you, so there is no manifest to build
and no control/payload split to maintain by hand:

* Code inside an `@app.function(...)` (or `@app.cls`) runs **remotely on the
  GPU**, inside the container image defined here. This is where torch,
  transformers and esm live.
* Code inside `@app.local_entrypoint()` runs **locally** -- it is the CLI. It
  only touches stdlib plus the `modal` client, calls `.remote()` on the GPU
  function, and writes the returned results to local disk.
* Heavy imports are deferred with `with IMAGE.imports():` so that merely
  importing the app locally (to launch it) never needs torch on your laptop.

`evals/test_offline.py` enforces exactly that: no torch/transformers/peft/esm
import at module top level in a skill script, only inside functions or an
`imports()` block.

Reproducibility
---------------
Both forks (`transformers` and `esm`) are pinned to a COMMIT, not `@main` as
upstream's own pyproject does. A designed protein someone orders from a
synthesis vendor must be reproducible, and `@main` moving under you silently
breaks that.
"""

from __future__ import annotations

import modal

# --------------------------------------------------------------------------
# Pinned dependency versions (verified on a live GPU, 2026-07)
# --------------------------------------------------------------------------

# The BioHub fork of transformers is MANDATORY, not a convenience:
# `transformers.models.esmc` and `.esmfold2` exist in NO PyPI release (checked
# against 4.57.1 and 5.14.1), and the biohub/* Hub repos ship no `modeling_*.py`
# and no `auto_map`, so `trust_remote_code` cannot substitute. Pinned to a
# commit; upstream pins `@main`.
TRANSFORMERS_COMMIT = 'ef32577f55da19a4989cd7b22e004dc43a4998cb'
ESM_COMMIT = '67838dc8ac76f4145613e6cb36c5f3d758542f7c'

TRANSFORMERS_FORK = (
    f'transformers @ git+https://github.com/Biohub/transformers.git'
    f'@{TRANSFORMERS_COMMIT}'
)
ESM_FORK = f'esm @ git+https://github.com/Biohub/esm.git@{ESM_COMMIT}'

# --------------------------------------------------------------------------
# GPU selection
# --------------------------------------------------------------------------
#
# Modal GPU strings, with the VRAM that decides whether a workload fits. These
# map to the flavors Modal exposes; pick the SMALLEST that fits -- a bigger card
# is not faster for these jobs, only dearer and sometimes slower to schedule.

GPUS: dict[str, str] = {
    # key -> Modal gpu= string
    't4': 'T4',            # 16 GB  -- ESM3-open inference only
    'l4': 'L4',            # 24 GB
    'a10g': 'A10G',        # 24 GB
    'l40s': 'L40S',        # 48 GB
    'a100-40': 'A100',     # 40 GB  -- ESMC LoRA
    'a100-80': 'A100-80GB',  # 80 GB -- binder design
    'h100': 'H100',        # 80 GB
    'h200': 'H200',        # 141 GB
}

# Reference VRAM footprints (measured on a live GPU where noted).
VRAM_GB: dict[str, int] = {
    't4': 16, 'l4': 24, 'a10g': 24, 'l40s': 48,
    'a100-40': 40, 'a100-80': 80, 'h100': 80, 'h200': 141,
}


def gpu_arg(key: str) -> str:
  """Translate a short GPU key into Modal's `gpu=` string.

  Args:
    key: One of `GPUS` (e.g. 'a100-80').

  Returns:
    The Modal GPU string (e.g. 'A100-80GB').

  Raises:
    ValueError: Unknown key, listing the valid ones.
  """
  if key not in GPUS:
    raise ValueError(f'unknown GPU {key!r}; known: {", ".join(sorted(GPUS))}')
  return GPUS[key]


# --------------------------------------------------------------------------
# Persistent storage: the HuggingFace weight cache
# --------------------------------------------------------------------------
#
# A Modal Volume outlives any single run, so the ESM weights (ESMC-6B is ~12.7
# GB, and binder design also pulls six ESMFold2 checkpoints) are downloaded once
# and reused. Mounted at /models, with HF_HOME pointed there.

MODELS_DIR = '/models'
DATA_DIR = '/data'

# create_if_missing keeps the first run self-bootstrapping: no manual
# `modal volume create` step.
models_volume = modal.Volume.from_name('esm-model-cache', create_if_missing=True)
data_volume = modal.Volume.from_name('esm-data-cache', create_if_missing=True)

VOLUMES = {MODELS_DIR: models_volume, DATA_DIR: data_volume}

# HuggingFace token. The ESM open weights are PUBLIC and ungated (checked 2026-07:
# biohub/ESMC-*, biohub/ESMFold2*, biohub/esm3-sm-open-v1 all report
# gated=false), so a token is NOT required. A Modal Secret named
# "huggingface" is attached when it exists -- anonymous Hub pulls are
# rate-limited per source IP, and a token lifts that. Create one with:
#   modal secret create huggingface HF_TOKEN=hf_...
# It is referenced by name only; no token value ever appears in this code.
HF_SECRET_NAME = 'huggingface'

_COMMON_ENV = {
    'HF_HOME': MODELS_DIR,
    'TORCH_HOME': f'{MODELS_DIR}/torch',
    'HF_XET_HIGH_PERFORMANCE': '1',
    'TOKENIZERS_PARALLELISM': 'false',
    # Keeps long optimisation loops alive on smaller cards by avoiding
    # fragmentation, rather than raw-capacity OOMs.
    'PYTORCH_CUDA_ALLOC_CONF': 'expandable_segments:True',
}


# --------------------------------------------------------------------------
# Container images
# --------------------------------------------------------------------------


def _pin_transformers(image: 'modal.Image') -> 'modal.Image':
  """Force-reinstall transformers at the pinned commit, with no deps.

  `esm`'s pyproject declares `transformers @ ...@main`. Asking pip for that AND
  a commit-pinned transformers in one resolution is two direct URLs for one
  package, which fails with `ResolutionImpossible` (observed on a live GPU).
  So every image installs `esm` first (dragging transformers@main), then this
  overwrites transformers with the pinned commit using --no-deps so pip does
  not re-resolve and re-conflict the rest of the graph.

  Args:
    image: The image after `esm` is installed.

  Returns:
    The image with transformers pinned.
  """
  return image.pip_install(
      TRANSFORMERS_FORK, extra_options='--no-deps --force-reinstall'
  )


def esm_image() -> 'modal.Image':
  """Image for the ESMC / ESM3 skills (LoRA fine-tuning, design campaigns).

  Deliberately lighter than the folding image: no flash-attn, no
  transformer-engine, no conda toolchain. ESMC loads and runs fine without them
  (the LoRA skill detects at runtime whether `out_proj` is a Transformer-Engine
  Linear or a stock `nn.Linear` and adapts), and skipping them keeps the image
  build to minutes rather than the ~half hour flash-attn compilation costs.
  """
  image = (
      modal.Image.debian_slim(python_version='3.12')
      .apt_install('git', 'build-essential')
      .pip_install(
          'torch==2.8.0',
          index_url='https://download.pytorch.org/whl/cu128',
      )
      .pip_install(
          ESM_FORK,
          'peft==0.17.1',
          'accelerate==1.10.1',
          'scikit-learn==1.7.2',
          'pandas==2.3.3',
          'matplotlib==3.10.7',
      )
      .env(_COMMON_ENV)
  )
  return _pin_transformers(image)


def esmfold2_image() -> 'modal.Image':
  """Image for the ESMFold2 binder-design skill.

  Mirrors upstream `binder_design.py`'s own Modal image: the CUDA toolchain,
  flash-attn, transformer-engine and xformers that the folding trunk needs, plus
  ANARCI + HMMER (conda-only) for antibody CDR annotation. This is the heavy one
  -- flash-attn and transformer-engine compile from source -- but Modal caches
  the built image, so the cost is paid once.
  """
  image = (
      modal.Image.micromamba(python_version='3.12')
      .apt_install('git', 'build-essential')
      .micromamba_install(
          'anarci>=2020.04.03',
          'hmmer=3.4',
          'cuda-version=12.8',
          'cuda-libraries-dev=12.8',
          'cuda-nvcc=12.8',
          'cmake',
          'ninja',
          channels=['conda-forge', 'bioconda'],
      )
      .pip_install(
          'torch==2.8.0',
          'triton==3.4.0',
          index_url='https://download.pytorch.org/whl/cu128',
      )
      .pip_install(
          'flash-attn==2.8.3',
          'transformer-engine[core-cu12,pytorch]==2.13.0',
          'xformers==0.0.32.post1',
          extra_options='--no-build-isolation',
      )
      .pip_install(
          'abnumber',
          ESM_FORK,
          'biotite==1.6.0',
      )
      .env({**_COMMON_ENV, 'XFORMERS_IGNORE_FLASH_VERSION_CHECK': '1'})
  )
  return _pin_transformers(image)


def hf_secrets() -> list['modal.Secret']:
  """The optional HuggingFace secret, as a list ready for `secrets=`.

  Returns an empty list if the named secret does not exist, so a run never
  fails merely because the (optional) token was not created.
  """
  try:
    return [modal.Secret.from_name(HF_SECRET_NAME)]
  except modal.exception.NotFoundError:
    return []


# --------------------------------------------------------------------------
# Small helpers shared by the payloads
# --------------------------------------------------------------------------


# The pip commands that install this environment on a LOCAL machine (for the
# `local` backend). Mirrors esm_image(): install esm first, then force-reinstall
# the pinned transformers with --no-deps. The binder-design skill additionally
# needs flash-attn / transformer-engine / ANARCI and is better run on Modal;
# see its SKILL.md.
LOCAL_INSTALL = (
    'pip install torch==2.8.0 '
    '--index-url https://download.pytorch.org/whl/cu128\n'
    f"pip install '{ESM_FORK}' peft==0.17.1 accelerate==1.10.1 "
    'scikit-learn==1.7.2 pandas==2.3.3 matplotlib==3.10.7\n'
    f"pip install --no-deps --force-reinstall '{TRANSFORMERS_FORK}'"
)


def gpu_available() -> bool:
  """True only if a local CUDA GPU is actually usable right now.

  "Verified available" means both that torch is importable in this interpreter
  AND that it reports a CUDA device -- not merely that a GPU exists on the box.
  Used to decide whether the `local` backend may run at all; a login node or a
  laptop returns False, and the skill routes to Modal instead.
  """
  try:
    import torch
  except ImportError:
    return False
  try:
    return bool(torch.cuda.is_available())
  except Exception:  # noqa: BLE001  (a broken driver must read as "no GPU")
    return False


def require_local_gpu() -> None:
  """Raise a helpful SystemExit unless a local CUDA GPU is verified usable.

  Called at the top of every skill's `local` backend, so the failure is a clear
  instruction (install these, or use Modal) rather than a NameError or an OOM
  ten minutes in.
  """
  try:
    import torch
  except ImportError:
    raise SystemExit(
        'The local backend needs torch and the ESM forks installed in THIS '
        'Python environment. Install them with:\n\n  ' + LOCAL_INSTALL +
        '\n\nOr run on a rented GPU with `modal run` instead (no local install '
        'needed).'
    )
  if not torch.cuda.is_available():
    raise SystemExit(
        'No local CUDA GPU is available (torch.cuda.is_available() is False). '
        'This machine cannot run the local backend. Use `modal run` to run on '
        'a rented GPU instead, or move to a machine with a GPU.'
    )


def choose_backend(requested: str) -> str:
  """Resolve a --backend value to 'local' or 'modal'.

  Args:
    requested: 'auto', 'local' or 'modal'.

  Returns:
    'local' if it will run in-process here, else 'modal'.

  Raises:
    SystemExit: 'local' was forced but no GPU is verified (via
      `require_local_gpu`), or an unknown value was passed.
  """
  if requested == 'modal':
    return 'modal'
  if requested == 'local':
    require_local_gpu()
    return 'local'
  if requested == 'auto':
    # Prefer a verified local GPU (free, no upload); otherwise Modal.
    return 'local' if gpu_available() else 'modal'
  raise SystemExit(f"unknown --backend {requested!r}; use auto|local|modal")


def default_gpu_note(gpu_key: str, need_gb: int) -> str | None:
  """Return a warning string if the chosen GPU is below a workload's floor.

  Args:
    gpu_key: The `GPUS` key the user selected.
    need_gb: The workload's minimum VRAM in GB.

  Returns:
    A human-readable warning, or None if the card is big enough.
  """
  have = VRAM_GB.get(gpu_key, 0)
  if have and have < need_gb:
    return (
        f'GPU {gpu_key!r} has ~{have} GB VRAM but this workload needs '
        f'~{need_gb} GB; expect an out-of-memory error. Pick a larger --gpu.'
    )
  return None
