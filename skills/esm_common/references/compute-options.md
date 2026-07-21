# Providing compute to the ESM skills

The ESM skills come in two families, and there are up to three ways to give them
compute. This page is the single reference; skills link here.

## 1. The API skills — a Biohub API key (the default)

The twelve retrieval/prediction skills (embeddings, folding, mutation scoring,
function prediction, structure prediction, …) run on Biohub's **hosted API**.
They need a **`BIOHUB_API_KEY`**, and nothing else — no GPU, no downloads.

**Get one** at https://biohub.ai/developer-console/api-keys, then **put it in
`~/.env`** using the safe protocol in the `credentials` skill (which never prints
the key). The skills read it from `~/.env` automatically.

```bash
# The credentials skill generates a command like this for you:
printf 'BIOHUB_API_KEY=%s\n' "<your-key>" >> ~/.env
```

That is all most users need. The rest of this page is the **alternative**, for
when you would rather not use the hosted API — for example to avoid the
100-credit/day cap, or to keep sequences on your own hardware.

### Running the API-skill *capabilities* without the hosted API

Everything the API skills compute (ESMC embeddings and logits, ESMFold2
structure, ESM3 generation) is also computable from the **open weights** on a
GPU. That is exactly what the three GPU skills do, and two of them cover the
heavy cases directly:

* **Folding / structure** → `esmfold2-structure-prediction` is API-based, but
  `esmfold2-binder-design` runs ESMFold2 on the open weights (Modal or a local
  GPU).
* **Design at scale** → `esm3-design-campaign` runs ESM3 generation on the open
  weights, unmetered.
* **Fine-tuning / gradients** → only possible on the open weights;
  `esmc-finetune-lora`.

For the lighter capabilities (a handful of embeddings, one fold) the hosted API
is simplest; reach for the GPU path when the API's limits are the problem.

## 2. The GPU skills — Modal, or a local GPU

The three GPU skills (`esmc-finetune-lora`, `esmfold2-binder-design`,
`esm3-design-campaign`) run the open weights on a GPU. They need **no Biohub key
and spend zero Biohub credits**. There are two backends.

### Modal (default; no GPU of your own required)

[Modal](https://modal.com) rents GPUs by the second and builds the container
image for you — nothing to install locally beyond the client.

```bash
pip install modal
modal token new                 # one-time browser auth
modal run scripts/<name>.py --...   # runs on a rented GPU
```

Optional, to lift the per-IP anonymous HuggingFace rate limit (the weights are
public, so a token is not required):

```bash
modal secret create huggingface HF_TOKEN=hf_your_token_here
```

The token value lives only in your Modal account; it never appears in any file.

### A local GPU — **only if one is verified available**

If you are on a machine that actually has a CUDA GPU (a workstation, or a GPU
container/pod), you can run the same code **in-process**, with no Modal and no
upload. Invoke the script with plain Python instead of `modal run`:

```bash
python scripts/<name>.py --...        # runs locally IF a GPU is verified
```

"Verified" is checked at runtime — `torch.cuda.is_available()` must be true and
the ESM dependencies must be importable. If they are not (e.g. a login node or a
laptop), the script does **not** silently fall back to CPU; it prints the exact
`modal run` command to use instead. Force the choice with `--backend
local|modal|auto` (default `auto` = local if verified, else Modal).

Installing the local dependencies (a one-time, heavy install — the pinned
`transformers` and `esm` forks):

```bash
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu128
pip install 'esm @ git+https://github.com/Biohub/esm.git@67838dc8ac76f4145613e6cb36c5f3d758542f7c' \
    peft==0.17.1 accelerate==1.10.1 scikit-learn==1.7.2 pandas==2.3.3 matplotlib==3.10.7
pip install --no-deps --force-reinstall \
    'transformers @ git+https://github.com/Biohub/transformers.git@ef32577f55da19a4989cd7b22e004dc43a4998cb'
```

`esmfold2-binder-design` additionally needs flash-attn, transformer-engine and
ANARCI/HMMER (conda); locally that is a substantial build, so **Modal is usually
the easier backend for binder design**. The lighter ESMC/ESM3 skills install
cleanly with the commands above.

## Which should I use?

| Situation | Use |
|---|---|
| A few embeddings / one fold, occasional | **Biohub API key** (§1) |
| Hitting the 100-credit/day cap | GPU skill on **Modal** (§2) |
| No GPU of your own | **Modal** |
| You have a CUDA GPU and the deps installed | **local** backend (§2) |
| Fine-tuning, or gradient-based design | GPU skill (Modal or local); the API cannot do gradients |
