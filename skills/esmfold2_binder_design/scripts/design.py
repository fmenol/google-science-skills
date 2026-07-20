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

"""De novo binder design with ESMFold2 on a Modal GPU.

    modal run scripts/design.py --target-name pd-l1 --num-seeds 8

Why this loads upstream code instead of reimplementing it
--------------------------------------------------------
The protocol is Algorithms 11-15 of the ESM world-model paper: a 150-step
optimisation of a *soft* sequence-logit tensor through the ESMFold2 trunk, with
entropy-based contact losses off the distogram, a radius-of-gyration term, and
an ESMC-6B pseudo-perplexity regulariser whose gradient is separately
L2-normalised. Subtleties (`logits.grad` assigned directly, per-step model
sampling under forced 0.5 dropout, temperature-scaled LR) are load-bearing and
easy to get quietly wrong. So this runs **BioHub's own implementation**,
fetched at a pinned commit.

BioHub's `binder_design.py` is itself a Modal script, but its `modal.App` and
image are constructed at import time and would collide with this app's. So the
GPU worker fetches the source, truncates it at the `# ---- Modal ----` section
(everything above is Modal-free and contains the whole design protocol), and
execs it. The `ESMFold2Design` class it defines is a plain class whose
`load()` / `design()` this app drives directly.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = ["modal"]
# ///

# NOTE: no `from __future__ import annotations` here. `modal.parameter` reads
# the real annotation type off the class (`use_scaling_critics: bool`), and the
# future import would turn it into the string 'bool', which Modal cannot decode.
import json
import pathlib
import re
import sys

import modal

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import esm_modal as em  # noqa: E402  (vendored, same directory)

app = modal.App('esm-binder-design')
IMAGE = em.esmfold2_image()

ESM_COMMIT = em.ESM_COMMIT
UPSTREAM_URL = (
    f'https://raw.githubusercontent.com/Biohub/esm/{ESM_COMMIT}'
    f'/cookbook/tutorials/binder_design.py'
)

# Targets and scaffolds baked into the upstream protocol. Listed so the CLI can
# validate before spending GPU time. Upstream's key is 'pdgfr' (its notebook
# says 'pdgfrb').
PRESET_TARGETS = ('cd45', 'ctla4', 'egfr', 'pd-l1', 'pdgfr')
PRESET_BINDERS = ('minibinder', 'trastuzumab_framework_vhvl',
                  'atezolizumab_framework_vhvl', 'ocankitug_framework_vhvl')
ANTIBODY_BINDERS = tuple(b for b in PRESET_BINDERS if b != 'minibinder')


def strip_modal(src: str) -> str:
  """Remove BioHub's Modal layer so its protocol can be exec'd standalone.

  Truncates at the `# ---- Modal ----` header (or `def get_base_image(` /
  `app = modal.App(` as fallbacks) and drops any top-level `import modal`.
  Post-checks that no executable `modal.` reference survives; a surviving one
  means upstream was restructured and this needs updating rather than failing
  obscurely after model load.
  """
  markers = (r'^#\s*-+\s*Modal\s*-+\s*$', r'^def get_base_image\(',
             r'^app\s*=\s*modal\.App\(')
  cut = None
  for pat in markers:
    m = re.search(pat, src, re.MULTILINE)
    if m:
      cut = m.start()
      break
  if cut is None:
    raise RuntimeError('no Modal boundary in upstream binder_design.py; '
                       'upstream restructured -- update this skill')
  head = src[:cut]
  head = re.sub(r'^import modal\s*$', '', head, flags=re.MULTILINE)
  head = re.sub(r'^from modal import .*$', '', head, flags=re.MULTILINE)
  for i, line in enumerate(head.split('\n'), 1):
    if re.search(r'\bmodal\.', line.split('#', 1)[0]):
      raise RuntimeError(f'modal ref survived stripping at line {i}')
  return head


@app.cls(
    image=IMAGE,
    gpu=em.GPUS['a100-80'],
    volumes=em.VOLUMES,
    secrets=em.hf_secrets(),
    timeout=3600,
    scaledown_window=60,
)
class Designer:
  """Loads the ESMFold2 design models once per container, designs per call."""

  use_scaling_critics: bool = modal.parameter(default=False)

  @modal.enter()
  def load(self):
    """Fetch + strip upstream, then load ESMC-6B + inversion + hero critics."""
    import time
    import types
    import urllib.request

    with urllib.request.urlopen(UPSTREAM_URL, timeout=120) as resp:
      src = strip_modal(resp.read().decode('utf-8'))
    mod = types.ModuleType('binder_design')
    exec(compile(src, 'binder_design.py', 'exec'), mod.__dict__)  # noqa: S102
    self._mod = mod
    t0 = time.time()
    self.designer = mod.ESMFold2Design()
    self.designer.load(self.use_scaling_critics)
    print(f'models loaded in {time.time() - t0:.0f}s', flush=True)

  @modal.method()
  def design_one(self, cfg: dict) -> dict:
    """Run one design (one seed) and return a JSON-safe result."""
    import time
    import torch

    is_ab = {'true': True, 'false': False, 'auto': None}[cfg['is_antibody']]
    torch.cuda.reset_peak_memory_stats()
    t0 = time.time()
    best_sequences, trajectory, critic_results = self.designer.design(
        target_name=cfg['target_name'], binder_name=cfg['binder_name'],
        target_sequence=cfg['target_sequence'] or None,
        binder_sequence=cfg['binder_sequence'] or None,
        is_antibody=is_ab, seed=cfg['seed'], batch_size=cfg['batch_size'],
    )
    peak = torch.cuda.max_memory_allocated() / (1024 ** 3)

    rows = _serialise_critics(critic_results)
    ranked = _rank(rows)
    designs = []
    for i, s in enumerate(best_sequences):
      target, _, binder = s.partition('|')
      designs.append({'batch_idx': i, 'target': target, 'binder': binder})
    traj = {str(k): {kk: (float(vv) if isinstance(vv, (int, float)) else str(vv))
                     for kk, vv in v.items()}
            for k, v in trajectory.items()}
    return {
        'status': 'ok', 'seed': cfg['seed'],
        'target_name': cfg['target_name'], 'binder_name': cfg['binder_name'],
        'esm_commit': ESM_COMMIT, 'gpu_name': torch.cuda.get_device_name(0),
        'peak_vram_gib': peak, 'design_seconds': time.time() - t0,
        'designs': designs, 'ranked': ranked, 'trajectory': traj,
        'critic_scores': rows,
    }


# --------------------------------------------------------------------------
# GPU-side serialisation helpers
# --------------------------------------------------------------------------


def _serialise_critics(critic_results: list) -> list[dict]:
  """Flatten critic results, embedding predicted complexes as PDB strings."""
  rows = []
  for i, r in enumerate(critic_results):
    row = {k: v for k, v in r.items()
           if k not in ('complex', 'logits')
           and isinstance(v, (str, int, float, bool, type(None)))}
    for key in ('iptm', 'final_loss', 'distogram_iptm_proxy',
                'cdr_distogram_iptm_proxy'):
      if key in r and key not in row:
        try:
          row[key] = float(r[key])
        except (TypeError, ValueError):
          row[key] = None
    cx = r.get('complex')
    if cx is not None:
      try:
        row['pdb'] = cx.to_pdb_string()
      except Exception as exc:  # noqa: BLE001
        row['pdb'] = None
        print(f'WARNING: could not serialise complex {i}: {exc}', flush=True)
    rows.append(row)
  return rows


def _rank(rows: list[dict]) -> list[dict]:
  """Aggregate per-critic scores into one ranked row per design (mean iPTM)."""
  by: dict[int, dict] = {}
  for r in rows:
    b = r.get('batch_idx', 0)
    pair = r.get('designed_sequence') or ''
    target, sep, binder = pair.partition('|')
    d = by.setdefault(b, {'batch_idx': b, 'sequence': binder if sep else pair,
                          'target_sequence': target if sep else '',
                          'iptms': [], 'critics': []})
    if isinstance(r.get('iptm'), float):
      d['iptms'].append(r['iptm'])
    d['critics'].append(r.get('critic_name'))
  out = []
  for d in by.values():
    iptms = d.pop('iptms')
    d['mean_iptm'] = sum(iptms) / len(iptms) if iptms else None
    d['min_iptm'] = min(iptms) if iptms else None
    d['max_iptm'] = max(iptms) if iptms else None
    d['n_critics'] = len(iptms)
    out.append(d)
  out.sort(key=lambda r: (r['mean_iptm'] is None, -(r['mean_iptm'] or 0)))
  for rank, r in enumerate(out, 1):
    r['rank'] = rank
  return out


# --------------------------------------------------------------------------
# Local entrypoint (the CLI)
# --------------------------------------------------------------------------


@app.local_entrypoint()
def main(
    target_name: str,
    binder_name: str = 'minibinder',
    target_sequence: str = '',
    binder_sequence: str = '',
    is_antibody: str = 'false',
    seed: int = 0,
    num_seeds: int = 8,
    batch_size: int = 1,
    use_scaling_critics: bool = False,
    out: str = './designs',
):
  """Design de novo binders on Modal GPUs and collect them ranked by iPTM.

  One GPU container per seed (fanned out with Modal `.map`); every design from
  every seed is pooled and ranked together by mean interface iPTM.
  """
  _validate(target_name, target_sequence, binder_name, binder_sequence,
            is_antibody)

  seeds = list(range(seed, seed + num_seeds))
  cfgs = [dict(target_name=target_name, binder_name=binder_name,
               target_sequence=target_sequence, binder_sequence=binder_sequence,
               is_antibody=is_antibody, seed=s, batch_size=batch_size)
          for s in seeds]

  print(f'designing {num_seeds} seed(s) on {em.GPUS["a100-80"]} '
        f'(one container per seed)...')
  designer = Designer(use_scaling_critics=use_scaling_critics)
  results = list(designer.design_one.map(cfgs))

  outdir = pathlib.Path(out)
  outdir.mkdir(parents=True, exist_ok=True)
  (outdir / 'structures').mkdir(exist_ok=True)

  all_designs = []
  for res in results:
    if res.get('status') != 'ok':
      continue
    s = res['seed']
    (outdir / f'seed{s}_result.json').write_text(json.dumps(
        {k: v for k, v in res.items() if k != 'critic_scores'}, indent=2))
    for row in res.get('critic_scores', []):
      pdb = row.pop('pdb', None)
      if pdb:
        name = re.sub(r'[^A-Za-z0-9_.-]', '_',
                      f'seed{s}_design{row.get("batch_idx")}_'
                      f'{row.get("critic_name")}.pdb')
        (outdir / 'structures' / name).write_text(pdb)
    for d in res.get('ranked', []):
      d['seed'] = s
      all_designs.append(d)

  all_designs.sort(key=lambda r: (r.get('mean_iptm') is None,
                                  -(r.get('mean_iptm') or 0)))
  for i, d in enumerate(all_designs, 1):
    d['global_rank'] = i
  (outdir / 'ranked_designs.json').write_text(json.dumps(all_designs, indent=2))
  (outdir / 'designs.fasta').write_text(''.join(
      f'>rank{d["global_rank"]}_seed{d["seed"]} mean_iptm={d.get("mean_iptm")}\n'
      f'{d.get("sequence", "")}\n' for d in all_designs if d.get('sequence')))

  print(f'\ncollected {len(all_designs)} design(s) into {outdir}')
  print(f'{"rank":>4}  {"seed":>4}  {"mean iPTM":>9}')
  for d in all_designs[:10]:
    iptm = d.get('mean_iptm')
    print(f'{d["global_rank"]:>4}  {d.get("seed"):>4}  '
          f'{"n/a" if iptm is None else f"{iptm:>9.3f}"}')
  if all_designs and (all_designs[0].get('mean_iptm') or 0) < 0.5:
    print('\nNOTE: best mean iPTM is below 0.5 -- no credible binder. Report '
          'that, rather than presenting rank 1 as a hit.')


def _validate(tname, tseq, bname, bseq, is_ab):
  if tname in PRESET_TARGETS and tseq:
    raise SystemExit(f'{tname!r} is a preset target; omit --target-sequence.')
  if tname not in PRESET_TARGETS and not tseq:
    raise SystemExit(f'{tname!r} is not a preset ({", ".join(PRESET_TARGETS)}); '
                     f'--target-sequence is required.')
  for label, seq in (('--target-sequence', tseq), ('--binder-sequence', bseq)):
    if seq and '|' in seq:
      raise SystemExit(f'{label} contains "|": multi-chain is not supported.')
  if bname in ANTIBODY_BINDERS and is_ab == 'false':
    raise SystemExit(f'{bname!r} is an antibody scaffold but --is-antibody is '
                     f'false; that mismatch trips an upstream assertion.')
  if is_ab != 'false':
    # The Modal image ships ANARCI + HMMER, so the CDR-annotation path that was
    # impossible on the pip-only cluster image now runs. But antibody design
    # has NOT been verified end-to-end here -- only minibinder has. Allowed,
    # with a clear caveat, rather than blocked.
    print('WARNING: antibody design is available (ANARCI/HMMER are in the '
          'image) but UNVERIFIED in this skill; only minibinder design has '
          'been checked end-to-end. Treat results with extra caution.')
