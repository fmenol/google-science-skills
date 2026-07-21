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

"""LoRA fine-tuning of ESMC on a Modal GPU.

Run it with the Modal CLI (never `python` directly):

    modal run scripts/finetune.py --train-csv care:train --test-csv care:test \\
        --steps 1000

The `@app.function` runs remotely on the GPU; the `@app.local_entrypoint`
runs locally, launches it, and writes the adapter, metrics and figures to a
local `--out` directory. See `../SKILL.md`.

Ported from `esm/cookbook/tutorials/esmc_finetune.ipynb`, with the changes a
headless job needs and a notebook does not: it saves the adapter, stages the
dataset itself, writes figures to PNG (no `plt.show()`), and reports macro-F1
and per-class recall -- not just accuracy, which hides the class collapse the
tutorial itself ran into on the imbalanced CARE benchmark.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = ["modal"]
# ///

from __future__ import annotations

import base64
import json
import pathlib
import sys

import modal

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import esm_modal as em  # noqa: E402  (vendored, same directory)

app = modal.App('esm-finetune-lora')
IMAGE = em.esm_image()

# Heavy imports run only on the GPU side. This block executes inside the image,
# so importing this file locally to launch it never needs torch.
with IMAGE.imports():
  import os
  import time

  import matplotlib
  matplotlib.use('Agg')  # No display in a headless container.
  import matplotlib.pyplot as plt
  import numpy as np
  import pandas as pd
  import torch
  from peft import LoraConfig, get_peft_model
  from transformers import AutoTokenizer

# EC top-level class names, for readable reports. CARE task-1 labels are the
# integers 1-7; these say what they mean. (Module scope so both sides see it.)
EC1_NAMES = {
    1: 'Oxidoreductases', 2: 'Transferases', 3: 'Hydrolases', 4: 'Lyases',
    5: 'Isomerases', 6: 'Ligases', 7: 'Translocases',
}

# The CARE enzyme-classification benchmark (Yang et al., NeurIPS 2024 D&B), the
# dataset the upstream tutorial trains on. `care:train` / `care:test` resolve to
# these paths after a one-time download into the Modal data volume.
CARE_URL = 'https://zenodo.org/records/14004425/files/CARE_datasets.zip'
CARE_SPLITS = {
    'train': 'CARE_datasets/splits/task1/protein_train.csv',
    # The 30-50% sequence-identity band: deliberately the harder generalisation
    # setting, and the one the tutorial uses.
    'test': 'CARE_datasets/splits/task1/30-50_protein_test.csv',
}


# --------------------------------------------------------------------------
# GPU side
# --------------------------------------------------------------------------


@app.function(
    image=IMAGE,
    gpu=em.GPUS['a100-40'],
    volumes=em.VOLUMES,
    secrets=em.hf_secrets(),
    timeout=6 * 3600,
)
def train_remote(cfg: dict, train_bytes: bytes, test_bytes: bytes | None) -> dict:
  """Modal wrapper around `_train`. Runs on a rented GPU."""
  return _train(cfg, train_bytes, test_bytes)


def _train(cfg: dict, train_bytes: bytes, test_bytes: bytes | None) -> dict:
  """Fine-tune ESMC with LoRA and return the adapter, metrics and figures.

  This is the shared core: `train_remote` calls it on a Modal GPU, and the
  `local` backend calls it in-process on a verified local GPU. It uses the
  module-level heavy imports (torch, etc.), which are defined by the
  `IMAGE.imports()` block both remotely (in the container) and locally (when the
  deps are installed). Returns a JSON-safe dict; the adapter and PNGs are
  base64-encoded inside it so the caller can write them to disk with no shared
  filesystem.

  Args:
    cfg: All hyperparameters and column names.
    train_bytes: The training table as raw bytes, or b'' to use `cfg['train']`
      as a `care:` spec / volume path.
    test_bytes: The test table as raw bytes, or None.
  """
  import io
  import zipfile

  out = pathlib.Path('/tmp/esm-lora-out')
  out.mkdir(parents=True, exist_ok=True)
  torch.manual_seed(cfg['seed'])
  np.random.seed(cfg['seed'])
  is_reg = cfg['task'] == 'regression'

  if not torch.cuda.is_available():
    return {'status': 'error', 'error': 'no CUDA device on the Modal worker'}
  gpu_name = torch.cuda.get_device_name(0)
  print(f'device: {gpu_name}', flush=True)

  # --- data ---------------------------------------------------------------
  def _read_table(raw: bytes | None, spec: str) -> pd.DataFrame:
    if raw:
      return _load_frame(io.BytesIO(raw), spec)
    path = _resolve_care(spec)
    return _load_frame(path, path)

  train_df = _read_table(train_bytes, cfg['train'])
  raw_labels = train_df[cfg['label_col']].to_numpy()
  seqs = train_df[cfg['seq_col']].astype(str).tolist()

  label_names: list[str] = []
  label_to_idx: dict = {}
  if is_reg:
    num_labels = 1
    y = raw_labels.astype(np.float32)
  else:
    classes = sorted(set(raw_labels.tolist()))
    label_to_idx = {c: i for i, c in enumerate(classes)}
    num_labels = len(classes)
    label_names = [
        f'{c} {EC1_NAMES[c]}'
        if isinstance(c, (int, np.integer)) and c in EC1_NAMES else str(c)
        for c in classes
    ]
    y = np.array([label_to_idx[v] for v in raw_labels], dtype=np.int64)

  tr_idx, va_idx = _stratified_split(
      y if not is_reg else np.zeros(len(y), int),
      cfg['val_fraction'], cfg['seed'],
  )
  print(f'train {len(tr_idx)} | val {len(va_idx)} | labels {num_labels}',
        flush=True)

  # --- model --------------------------------------------------------------
  try:
    from transformers import ESMCForSequenceClassification
  except ImportError:
    return {
        'status': 'error',
        'error': 'transformers has no ESMCForSequenceClassification; the '
                 'BioHub fork did not install. See esm_gpu_common/SPEC_GPU.md.',
    }

  tokenizer = AutoTokenizer.from_pretrained(cfg['model'])
  model = ESMCForSequenceClassification.from_pretrained(
      cfg['model'], num_labels=num_labels,
      problem_type='regression' if is_reg else 'single_label_classification',
  ).to('cuda')

  # ESMC fuses qkv into `layernorm_qkv` and uses SwiGLU with fc1/fc2 as bare
  # nn.Parameter tensors, unreachable through target_modules -- hence
  # target_parameters (needs peft >= 0.17). `attn.out_proj` is environment-
  # dependent: with Transformer Engine present it is a TE Linear, not an
  # nn.Linear, and PEFT's target_modules rejects it. Detected, not assumed.
  te_out_proj = [
      n for n, mod in model.named_modules()
      if n.endswith('attn.out_proj') and not isinstance(mod, torch.nn.Linear)
  ]
  target_modules: list[str] = []
  target_parameters = [
      'layernorm_qkv.weight', 'ffn.fc1_weight', 'ffn.fc2_weight',
  ]
  if te_out_proj:
    target_parameters.append('attn.out_proj.weight')
  else:
    target_modules.append('out_proj')

  # peft's ParamWrapper (used for every target_parameters entry) rejects a
  # non-zero dropout. ESMC's fused weights can only be reached that way, so
  # dropout is coerced to 0 with a note rather than failing.
  lora_dropout = cfg['lora_dropout']
  if target_parameters and lora_dropout:
    print(f'NOTE: lora_dropout={lora_dropout} unsupported for fused params; '
          f'using 0.0', flush=True)
    lora_dropout = 0.0

  model = get_peft_model(model, LoraConfig(
      r=cfg['lora_rank'], lora_alpha=cfg['lora_alpha'], lora_dropout=lora_dropout,
      target_modules=target_modules, target_parameters=target_parameters,
      modules_to_save=['classifier'],
  ))
  model.print_trainable_parameters()
  n_trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
  n_total = sum(p.numel() for p in model.parameters())
  optimizer = torch.optim.AdamW(
      [p for p in model.parameters() if p.requires_grad], lr=cfg['lr'])

  # --- train --------------------------------------------------------------
  torch.cuda.reset_peak_memory_stats()
  rng = np.random.default_rng(cfg['seed'])
  perm = rng.permutation(len(tr_idx))
  history: list[dict] = []
  t0 = time.time()
  model.train()
  bs, max_len = cfg['batch_size'], cfg['max_seq_len']

  for step in range(cfg['steps']):
    start = (step * bs) % len(tr_idx)
    if start + bs > len(tr_idx):
      perm = rng.permutation(len(tr_idx))
      start = 0
    idx = tr_idx[perm[start:start + bs]]
    batch = _batch(tokenizer, seqs, y, idx, max_len, is_reg)
    with torch.autocast(device_type='cuda', dtype=torch.bfloat16):
      outputs = model(**batch)
      loss = outputs.loss
    loss.backward()
    optimizer.step()
    optimizer.zero_grad()

    rec = {'step': step, 'train_loss': float(loss.detach())}
    if step % 20 == 0 or step == cfg['steps'] - 1:
      print(f'step {step}/{cfg["steps"]} loss {rec["train_loss"]:.4f}',
            flush=True)
    if len(va_idx) and ((step + 1) % cfg['val_every'] == 0
                        or step == cfg['steps'] - 1):
      vt = y[va_idx]
      _, vpred = _evaluate(model, tokenizer,
                           [seqs[i] for i in va_idx], [y[i] for i in va_idx],
                           max_len, cfg['eval_batch_size'], is_reg)
      m = (_reg_metrics(vt, vpred) if is_reg
           else _clf_metrics(vt, vpred, num_labels))
      rec['val'] = {k: v for k, v in m.items() if k != 'confusion_matrix'}
    history.append(rec)

  train_seconds = time.time() - t0
  peak_gib = torch.cuda.max_memory_allocated() / (1024 ** 3)

  # --- test ---------------------------------------------------------------
  test_metrics = None
  if test_bytes or cfg.get('test'):
    tdf = _read_table(test_bytes, cfg.get('test', ''))
    tseqs = tdf[cfg['seq_col']].astype(str).tolist()
    if is_reg:
      ty = tdf[cfg['label_col']].to_numpy().astype(np.float32)
    else:
      keep = [i for i, v in enumerate(tdf[cfg['label_col']].tolist())
              if v in label_to_idx]
      tseqs = [tseqs[i] for i in keep]
      ty = np.array([label_to_idx[tdf[cfg['label_col']].iloc[i]] for i in keep],
                    dtype=np.int64)
    tloss, tpred = _evaluate(model, tokenizer, tseqs, ty, max_len,
                             cfg['eval_batch_size'], is_reg)
    test_metrics = (_reg_metrics(ty, tpred) if is_reg
                    else _clf_metrics(ty, tpred, num_labels))
    test_metrics['loss'] = tloss
    test_metrics['n'] = int(len(ty))

  # --- package outputs ----------------------------------------------------
  adapter_dir = out / 'adapter'
  model.save_pretrained(str(adapter_dir))
  tokenizer.save_pretrained(str(adapter_dir))
  zbuf = io.BytesIO()
  with zipfile.ZipFile(zbuf, 'w', zipfile.ZIP_DEFLATED) as zf:
    for f in adapter_dir.rglob('*'):
      if f.is_file():
        zf.write(f, f.relative_to(adapter_dir))

  figs = {'training_curves.png': _plot_history(history, is_reg)}
  if test_metrics and not is_reg:
    figs['confusion_matrix.png'] = _plot_confusion(
        np.array(test_metrics['confusion_matrix']), label_names)

  return {
      'status': 'ok',
      'task': cfg['task'], 'model': cfg['model'], 'gpu_name': gpu_name,
      'seed': cfg['seed'], 'steps': cfg['steps'], 'batch_size': bs, 'lr': cfg['lr'],
      'lora': {'rank': cfg['lora_rank'], 'alpha': cfg['lora_alpha'],
               'dropout': lora_dropout, 'target_modules': target_modules,
               'target_parameters': target_parameters},
      'trainable_params': n_trainable, 'total_params': n_total,
      'trainable_fraction': n_trainable / max(n_total, 1),
      'n_train': int(len(tr_idx)), 'n_val': int(len(va_idx)),
      'num_labels': num_labels, 'label_names': label_names,
      'train_seconds': train_seconds, 'peak_vram_gib': peak_gib,
      'test_metrics': test_metrics, 'history': history,
      'adapter_zip_b64': base64.b64encode(zbuf.getvalue()).decode(),
      'figures_b64': {k: base64.b64encode(v).decode() for k, v in figs.items()},
  }


# --------------------------------------------------------------------------
# GPU-side helpers (module scope, imported into the image with the app)
# --------------------------------------------------------------------------


def _care_root() -> pathlib.Path:
  """Where to cache CARE: the Modal data volume if mounted, else a local dir."""
  vol = pathlib.Path(em.DATA_DIR)
  if vol.is_dir():          # On Modal the volume is mounted at /data.
    return vol / 'CARE'
  return pathlib.Path.home() / '.cache' / 'esm-data' / 'CARE'


def _resolve_care(spec: str) -> str:
  """Resolve `care:<split>` to a path, downloading once into the cache."""
  if not spec.startswith('care:'):
    return spec
  split = spec.split(':', 1)[1]
  if split not in CARE_SPLITS:
    raise SystemExit(f'unknown CARE split {split!r}')
  root = _care_root()
  target = root / CARE_SPLITS[split]
  if target.is_file():
    return str(target)
  import urllib.request
  import zipfile
  root.mkdir(parents=True, exist_ok=True)
  zpath = root / 'CARE_datasets.zip'
  if not zpath.is_file():
    tmp = zpath.with_suffix('.partial')
    urllib.request.urlretrieve(CARE_URL, tmp)
    tmp.rename(zpath)
  with zipfile.ZipFile(zpath) as zf:
    zf.extractall(root)
  try:
    em.data_volume.commit()  # Persist for the next run (Modal only; no-op else).
  except Exception:  # noqa: BLE001
    pass
  if not target.is_file():
    raise SystemExit(f'CARE archive missing {CARE_SPLITS[split]}')
  return str(target)


def _load_frame(src, name: str) -> 'pd.DataFrame':
  suffix = str(name).lower()
  if suffix.endswith(('.parquet', '.pq')):
    return pd.read_parquet(src)
  sep = '\t' if suffix.endswith(('.tsv', '.tab')) else ','
  return pd.read_csv(src, sep=sep)


def _stratified_split(labels, fraction, seed):
  rng = np.random.default_rng(seed)
  train, val = [], []
  for cls in np.unique(labels):
    idx = np.flatnonzero(labels == cls)
    rng.shuffle(idx)
    n_val = int(round(len(idx) * fraction))
    val.extend(idx[:n_val])
    train.extend(idx[n_val:])
  return np.array(sorted(train)), np.array(sorted(val))


def _batch(tokenizer, seqs, labels, idx, max_len, is_reg):
  batch = tokenizer([seqs[i] for i in idx], return_tensors='pt',
                    padding=True, truncation=True, max_length=max_len)
  batch = {k: v.to('cuda') for k, v in batch.items()}
  dtype = torch.float32 if is_reg else torch.long
  batch['labels'] = torch.tensor([labels[i] for i in idx], dtype=dtype,
                                 device='cuda')
  return batch


def _evaluate(model, tokenizer, seqs, labels, max_len, bs, is_reg):
  model.eval()
  losses, preds = [], []
  with torch.no_grad():
    for start in range(0, len(seqs), bs):
      idx = list(range(start, min(start + bs, len(seqs))))
      batch = _batch(tokenizer, seqs, labels, idx, max_len, is_reg)
      with torch.autocast(device_type='cuda', dtype=torch.bfloat16):
        out = model(**batch)
      losses.append(float(out.loss))
      logits = out.logits.float()
      preds.append(logits.squeeze(-1).cpu().numpy() if is_reg
                   else logits.argmax(-1).cpu().numpy())
  model.train()
  return float(np.mean(losses)), np.concatenate(preds)


def _clf_metrics(y_true, y_pred, n_classes):
  cm = np.zeros((n_classes, n_classes), int)
  for t, p in zip(y_true, y_pred):
    cm[t, p] += 1
  tp = np.diag(cm).astype(float)
  support = cm.sum(1).astype(float)
  predicted = cm.sum(0).astype(float)
  recall = np.divide(tp, support, out=np.zeros_like(tp), where=support > 0)
  precision = np.divide(tp, predicted, out=np.zeros_like(tp), where=predicted > 0)
  denom = precision + recall
  f1 = np.divide(2 * precision * recall, denom, out=np.zeros_like(tp),
                 where=denom > 0)
  return {
      'accuracy': float(tp.sum() / max(cm.sum(), 1)),
      'macro_f1': float(f1.mean()),
      'per_class_recall': recall.tolist(),
      'per_class_precision': precision.tolist(),
      'per_class_support': support.astype(int).tolist(),
      'confusion_matrix': cm.tolist(),
      'classes_never_predicted':
          [int(i) for i in range(n_classes) if predicted[i] == 0],
  }


def _reg_metrics(y_true, y_pred):
  err = y_pred - y_true

  def _corr(a, b):
    return float(np.corrcoef(a, b)[0, 1]) if a.std() and b.std() else float('nan')

  rank = lambda v: np.argsort(np.argsort(v)).astype(float)  # noqa: E731
  return {
      'rmse': float(np.sqrt(np.mean(err ** 2))),
      'mae': float(np.mean(np.abs(err))),
      'pearson_r': _corr(y_true, y_pred),
      'spearman_rho': _corr(rank(y_true), rank(y_pred)),
  }


def _plot_history(history, is_reg):
  import io
  steps = [h['step'] for h in history]
  loss = [h['train_loss'] for h in history]
  win = min(20, max(1, len(loss) // 10))
  smooth = np.convolve(loss, np.ones(win) / win, mode='valid')
  fig, ax = plt.subplots(1, 2, figsize=(11, 4))
  ax[0].plot(steps, loss, alpha=0.25, lw=0.8, label='per step')
  ax[0].plot(steps[win - 1:], smooth, lw=2, label=f'{win}-step mean')
  ax[0].set(xlabel='step', ylabel='train loss', title='Training loss')
  ax[0].legend()
  vsteps = [h['step'] for h in history if 'val' in h]
  if vsteps:
    key = 'rmse' if is_reg else 'accuracy'
    ax[1].plot(vsteps, [h['val'][key] for h in history if 'val' in h], 'o-',
               label=key)
    if not is_reg:
      ax[1].plot(vsteps, [h['val']['macro_f1'] for h in history if 'val' in h],
                 's-', label='macro F1')
    ax[1].legend()
  ax[1].set(xlabel='step', title='Validation')
  fig.tight_layout()
  buf = io.BytesIO()
  fig.savefig(buf, format='png', dpi=150)
  plt.close(fig)
  return buf.getvalue()


def _plot_confusion(cm, names):
  import io
  with np.errstate(invalid='ignore', divide='ignore'):
    norm = cm / np.maximum(cm.sum(1, keepdims=True), 1)
  fig, ax = plt.subplots(figsize=(1.1 * len(names) + 3, len(names) + 2))
  im = ax.imshow(norm, cmap='Blues', vmin=0, vmax=1)
  ax.set_xticks(range(len(names)), names, rotation=45, ha='right')
  ax.set_yticks(range(len(names)), names)
  ax.set(xlabel='predicted', ylabel='true',
         title='Confusion matrix (row-normalised)')
  for i in range(len(names)):
    for j in range(len(names)):
      ax.text(j, i, f'{norm[i, j]:.2f}', ha='center', va='center', fontsize=8,
              color='white' if norm[i, j] > 0.5 else 'black')
  fig.colorbar(im, ax=ax, fraction=0.046)
  fig.tight_layout()
  buf = io.BytesIO()
  fig.savefig(buf, format='png', dpi=150)
  plt.close(fig)
  return buf.getvalue()


# --------------------------------------------------------------------------
# Local entrypoint (the CLI)
# --------------------------------------------------------------------------


def _build_cfg(a: dict) -> dict:
  """Assemble the config dict shared by both backends from a params mapping."""
  return dict(
      train=a['train_csv'], test=a['test_csv'], model=a['model'],
      seq_col=a['seq_col'], label_col=a['label_col'], task=a['task'],
      steps=a['steps'], batch_size=a['batch_size'], lr=a['lr'],
      max_seq_len=a['max_seq_len'], val_every=a['val_every'],
      val_fraction=a['val_fraction'], lora_rank=a['lora_rank'],
      lora_alpha=a['lora_alpha'], lora_dropout=a['lora_dropout'],
      seed=a['seed'], eval_batch_size=a['eval_batch_size'],
  )


def _finish(res: dict, out: str) -> None:
  """Write results locally and print the summary (shared by both backends)."""
  outdir = pathlib.Path(out)
  outdir.mkdir(parents=True, exist_ok=True)
  if res.get('status') != 'ok':
    print(f'FAILED: {res.get("error")}')
    raise SystemExit(1)
  _write_results(res, outdir)
  _summarise(res)


@app.local_entrypoint()
def main(
    train_csv: str,
    test_csv: str = '',
    model: str = 'biohub/ESMC-300M',
    seq_col: str = 'Sequence',
    label_col: str = 'EC1',
    task: str = 'classification',
    steps: int = 1000,
    batch_size: int = 8,
    lr: float = 1e-4,
    max_seq_len: int = 1024,
    val_every: int = 250,
    val_fraction: float = 0.01,
    lora_rank: int = 8,
    lora_alpha: int = 16,
    lora_dropout: float = 0.01,
    seed: int = 0,
    eval_batch_size: int = 16,
    out: str = './lora_results',
):
  """Fine-tune ESMC with LoRA on a Modal GPU and save the adapter locally.

  This is the **Modal** entrypoint (`modal run scripts/finetune.py ...`). To run
  on a verified local GPU instead, invoke the script directly:
  `python scripts/finetune.py ...` (see the `__main__` backend below).

  `--train-csv` / `--test-csv` accept a local file (uploaded to the worker), a
  path inside the mounted data volume, or `care:train` / `care:test` for the
  CARE enzyme benchmark.
  """
  cfg = _build_cfg(locals())
  train_bytes = _maybe_read(train_csv)
  test_bytes = _maybe_read(test_csv) if test_csv else None
  print(f'launching on Modal GPU ({em.GPUS["a100-40"]})...')
  _finish(train_remote.remote(cfg, train_bytes, test_bytes), out)


def _maybe_read(spec: str) -> bytes:
  """Read a local file to bytes; return b'' for a `care:` spec or volume path."""
  if not spec or spec.startswith('care:'):
    return b''
  p = pathlib.Path(spec)
  return p.read_bytes() if p.is_file() else b''


def _write_results(res: dict, outdir: pathlib.Path) -> None:
  import zipfile
  import io
  adapter_zip = base64.b64decode(res.pop('adapter_zip_b64'))
  with zipfile.ZipFile(io.BytesIO(adapter_zip)) as zf:
    zf.extractall(outdir / 'adapter')
  for name, b64 in res.pop('figures_b64', {}).items():
    (outdir / name).write_bytes(base64.b64decode(b64))
  (outdir / 'result.json').write_text(json.dumps(res, indent=2))
  print(f'\nwrote adapter, figures and result.json to {outdir}')


def _summarise(res: dict) -> None:
  print(f'task          {res["task"]}')
  print(f'model         {res["model"]}')
  print(f'gpu           {res["gpu_name"]}')
  print(f'trainable     {res["trainable_params"]:,} / {res["total_params"]:,} '
        f'({100 * res["trainable_fraction"]:.3f}%)')
  print(f'train time    {res["train_seconds"]:.0f}s | peak VRAM '
        f'{res["peak_vram_gib"]:.1f} GiB')
  tm = res.get('test_metrics')
  if tm and res['task'] == 'regression':
    print(f'TEST          rmse {tm["rmse"]:.4f}  pearson {tm["pearson_r"]:.3f}  '
          f'spearman {tm["spearman_rho"]:.3f}')
  elif tm:
    print(f'TEST          accuracy {tm["accuracy"]:.4f}  '
          f'macro-F1 {tm["macro_f1"]:.4f}  (n={tm["n"]})')
    never = tm.get('classes_never_predicted') or []
    if never:
      print(f'  WARNING: {len(never)} class(es) never predicted: {never}. '
            f'Accuracy is not a sufficient summary here -- read macro-F1 and '
            f'per-class recall in result.json.')


# --------------------------------------------------------------------------
# Local backend (the `local` path — verified GPU, no Modal)
# --------------------------------------------------------------------------


def _run_local(a: dict) -> None:
  """Fine-tune in-process on a verified local GPU."""
  em.require_local_gpu()
  cfg = _build_cfg(a)
  train_bytes = _maybe_read(a['train_csv'])
  test_bytes = _maybe_read(a['test_csv']) if a['test_csv'] else None
  print('running locally on the detected GPU...')
  _finish(_train(cfg, train_bytes, test_bytes), a['out'])


if __name__ == '__main__':
  import argparse

  # This path is the LOCAL backend: `python scripts/finetune.py ...` runs on a
  # verified local GPU. `modal run scripts/finetune.py ...` uses the Modal
  # entrypoint above instead. `--backend auto` (default) prefers a local GPU if
  # one is verified, else prints the exact `modal run` command to use.
  ap = argparse.ArgumentParser(description='ESMC LoRA fine-tuning (local GPU).')
  ap.add_argument('--train-csv', required=True)
  ap.add_argument('--test-csv', default='')
  ap.add_argument('--model', default='biohub/ESMC-300M')
  ap.add_argument('--seq-col', default='Sequence')
  ap.add_argument('--label-col', default='EC1')
  ap.add_argument('--task', default='classification',
                  choices=('classification', 'regression'))
  ap.add_argument('--steps', type=int, default=1000)
  ap.add_argument('--batch-size', type=int, default=8)
  ap.add_argument('--lr', type=float, default=1e-4)
  ap.add_argument('--max-seq-len', type=int, default=1024)
  ap.add_argument('--val-every', type=int, default=250)
  ap.add_argument('--val-fraction', type=float, default=0.01)
  ap.add_argument('--lora-rank', type=int, default=8)
  ap.add_argument('--lora-alpha', type=int, default=16)
  ap.add_argument('--lora-dropout', type=float, default=0.01)
  ap.add_argument('--seed', type=int, default=0)
  ap.add_argument('--eval-batch-size', type=int, default=16)
  ap.add_argument('--out', default='./lora_results')
  ap.add_argument('--backend', default='auto',
                  choices=('auto', 'local', 'modal'),
                  help='auto: local GPU if verified, else Modal.')
  args = vars(ap.parse_args())

  backend = em.choose_backend(args['backend'])
  if backend == 'local':
    _run_local(args)
  else:
    quoted = ' '.join(
        f'--{k.replace("_", "-")} {v!r}' for k, v in args.items()
        if k != 'backend' and v not in ('', None))
    raise SystemExit(
        'No local GPU verified, so this must run on Modal. Use:\n\n'
        f'  modal run scripts/finetune.py {quoted}\n')
