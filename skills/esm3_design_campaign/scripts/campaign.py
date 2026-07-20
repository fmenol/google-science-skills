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

"""High-throughput ESM3 motif-scaffolding campaign on Modal GPUs.

    modal run scripts/campaign.py --preset gfp --num-generations 2000 \\
        --num-shards 8

Generalises `esm/cookbook/tutorials/gfp_design.ipynb` into a sharded rejection-
sampling campaign. The notebook runs its accept/reject loop once, by hand; each
generation is a billed API call, so at 100 credits/day a thousand-generation
campaign is a ten-day exercise. On rented GPUs it is minutes, and the volume the
protocol depends on becomes practical.

The protocol per attempt: prompt with a few fixed motif residues (+ optional
template structure), generate a structure, keep it only if the constrained site
is faithful AND the backbone is genuinely novel, generate a sequence, refold
deterministically, and check the site survived.

Honest caveat: `gfp_design.ipynb` targets `esm3-medium` (7B); the only public
ESM3 weights are the 1.4B open model, so pass rates are lower and volume is the
compensation, not parity. See SKILL.md.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = ["modal"]
# ///

from __future__ import annotations

import json
import pathlib
import sys

import modal

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import esm_modal as em  # noqa: E402  (vendored, same directory)

app = modal.App('esm-design-campaign')
IMAGE = em.esm_image()

with IMAGE.imports():
  import os
  import time

  import numpy as np

# ESM3 special structure tokens (esm/utils/constants/esm3.py).
STRUCTURE_MASK, STRUCTURE_EOS, STRUCTURE_BOS = 4096, 4097, 4098

# The GFP chromophore prompt from `gfp_design.ipynb`. A93R is deliberate.
# ALL POSITIONS ARE 0-INDEXED RESIDUES; the +1 BOS shift is applied internally.
# The notebook's structure window [56:71],[94],[220] is in already-shifted TOKEN
# indices -- transcribing it literally double-shifts and leaves motif residues
# 93 and 219 unpinned. These are the residue equivalents.
GFP_PRESET = {
    'motif': '59:T,62:T,63:Y,64:G,93:R,219:E',
    'template_pdb': '1qy3', 'template_chain': 'A',
    'structure_positions': '55-69,93,219', 'length': 229,
}


# --------------------------------------------------------------------------
# GPU side
# --------------------------------------------------------------------------


@app.function(
    image=IMAGE,
    gpu=em.GPUS['a100-40'],
    volumes=em.VOLUMES,
    secrets=em.hf_secrets(),
    timeout=3 * 3600,
)
def campaign_shard(cfg: dict) -> dict:
  """Run one shard of a design campaign and return its records + distributions.

  Args:
    cfg: motif, length, template, thresholds, num_generations, seed.

  Returns:
    A JSON-safe dict with per-attempt records, passing designs, the gate
    breakdown and the RMSD distributions over all attempts.
  """
  import torch
  from esm.models.esm3 import ESM3
  from esm.sdk.api import ESMProtein, GenerationConfig

  if not torch.cuda.is_available():
    return {'status': 'error', 'error': 'no CUDA device on the Modal worker'}
  torch.manual_seed(cfg['seed'])
  np.random.seed(cfg['seed'])
  motif = _parse_motif(cfg['motif'])
  motif_positions = sorted(motif)

  esm3 = ESM3.from_pretrained(cfg['model'], device=torch.device('cuda'))

  template_ca = None
  template_tokens = None
  template_seq = ''
  struct_positions: list[int] = []
  if cfg['template_pdb']:
    from esm.utils.structure.protein_chain import ProteinChain
    # ESMProtein.from_pdb takes a FILE PATH; RCSB accessions go through
    # ProteinChain.from_rcsb.
    if os.path.exists(cfg['template_pdb']):
      tmpl = ESMProtein.from_pdb(cfg['template_pdb'],
                                 chain_id=cfg['template_chain'])
    else:
      chain = ProteinChain.from_rcsb(cfg['template_pdb'], cfg['template_chain'])
      tmpl = ESMProtein.from_protein_chain(chain)
    template_seq = tmpl.sequence or ''
    template_ca = _ca(tmpl)
    template_tokens = esm3.encode(tmpl).structure
    struct_positions = _parse_positions(cfg['structure_positions'])
    unpinned = [p for p in motif_positions if p not in set(struct_positions)]
    if unpinned:
      # A scored-but-unpinned motif residue is free to move, so the site gate
      # can never pass. This is exactly how a position-index off-by-one shows.
      print(f'WARNING: motif residues {unpinned} are scored but not pinned; '
            f'check --structure-positions (0-indexed residues).', flush=True)

  prompt = _build_prompt(esm3, cfg['length'], motif, template_tokens,
                         struct_positions)

  torch.cuda.reset_peak_memory_stats()
  t0 = time.time()
  records = []
  for i in range(cfg['num_generations']):
    try:
      rec = _run_one(esm3, prompt, template_ca, motif, motif_positions, cfg, i,
                     ESMProtein, GenerationConfig)
    except Exception as exc:  # noqa: BLE001
      rec = {'index': i, 'passed': False, 'failed_gate': 'exception',
             'error': repr(exc)}
    records.append(rec)
    if rec.get('passed'):
      print(f'  [{i}] PASS site {rec.get("refold_site_rmsd", float("nan")):.2f}',
            flush=True)
    elif (i + 1) % 25 == 0:
      print(f'  {i + 1}/{cfg["num_generations"]}, '
            f'{sum(r.get("passed") for r in records)} passed', flush=True)

  elapsed = time.time() - t0
  passed = [r for r in records if r.get('passed')]
  gate_counts: dict[str, int] = {}
  for r in records:
    g = r.get('failed_gate', 'passed')
    gate_counts[g] = gate_counts.get(g, 0) + 1
  for r in passed:
    if template_seq and r.get('sequence'):
      r['identity_to_template'] = _identity(
          r['sequence'], template_seq[:len(r['sequence'])])

  def _dist(key):
    vals = sorted(r[key] for r in records
                  if isinstance(r.get(key), float) and r[key] == r[key])
    return None if not vals else {
        'n': len(vals), 'best': vals[0], 'median': vals[len(vals) // 2],
        'worst': vals[-1]}

  return {
      'status': 'ok', 'model': cfg['model'], 'seed': cfg['seed'],
      'gpu_name': torch.cuda.get_device_name(0),
      'num_generations': cfg['num_generations'], 'num_passed': len(passed),
      'pass_rate': len(passed) / max(cfg['num_generations'], 1),
      'gate_breakdown': gate_counts,
      'rmsd_distributions': {k: _dist(k) for k in (
          'constrained_site_rmsd', 'backbone_rmsd', 'refold_site_rmsd')},
      'template_pdb': cfg['template_pdb'],
      'motif': {str(k): v for k, v in motif.items()},
      'thresholds': {'site_rmsd_max': cfg['site_rmsd_max'],
                     'backbone_rmsd_min': cfg['backbone_rmsd_min']},
      'elapsed_seconds': elapsed,
      'peak_vram_gib': torch.cuda.max_memory_allocated() / (1024 ** 3),
      'passed_designs': passed,
  }


# --------------------------------------------------------------------------
# GPU-side helpers
# --------------------------------------------------------------------------


def _kabsch_rmsd(a, b):
  if a.shape != b.shape or len(a) == 0:
    return float('nan')
  finite = np.isfinite(a).all(1) & np.isfinite(b).all(1)
  a, b = a[finite], b[finite]
  if len(a) < 3:
    return float('nan')
  a = a - a.mean(0)
  b = b - b.mean(0)
  u, _, vt = np.linalg.svd(a.T @ b)
  d = np.sign(np.linalg.det(u @ vt))
  rot = u @ np.diag([1.0, 1.0, d]) @ vt
  return float(np.sqrt(((a @ rot - b) ** 2).sum() / len(a)))


def _ca(protein):
  c = protein.coordinates
  arr = c.detach().cpu().numpy() if hasattr(c, 'detach') else np.asarray(c)
  return arr[:, 1, :]


def _identity(a, b):
  if not a or not b or len(a) != len(b):
    return float('nan')
  return sum(x == y for x, y in zip(a, b)) / len(a)


def _parse_motif(spec):
  motif = {}
  for part in filter(None, (p.strip() for p in spec.split(','))):
    pos, aa = part.split(':')
    motif[int(pos)] = aa.strip().upper()
  return motif


def _parse_positions(spec):
  out = []
  for part in filter(None, (p.strip() for p in spec.split(','))):
    if '-' in part:
      lo, hi = part.split('-')
      out.extend(range(int(lo), int(hi) + 1))
    else:
      out.append(int(part))
  return sorted(set(out))


def _build_prompt(esm3, length, motif, template_tokens, structure_positions):
  from esm.sdk.api import ESMProtein
  import torch
  seq = ['_'] * length
  for pos, aa in motif.items():
    if pos >= length:
      raise SystemExit(f'motif position {pos} beyond length {length}')
    seq[pos] = aa
  prompt = esm3.encode(ESMProtein(sequence=''.join(seq)))
  if template_tokens is not None and structure_positions:
    struct = torch.full((len(prompt.sequence),), STRUCTURE_MASK, dtype=torch.long)
    struct[0] = STRUCTURE_BOS
    struct[-1] = STRUCTURE_EOS
    for pos in structure_positions:
      if 0 <= pos < length and pos + 1 < len(struct):
        struct[pos + 1] = template_tokens[pos + 1]
    prompt.structure = struct.to(prompt.sequence.device)
  return prompt


def _run_one(esm3, prompt, template_ca, motif, motif_positions, cfg, idx,
             ESMProtein, GenerationConfig):
  rec = {'index': idx, 'passed': False}
  n_masked = int((prompt.structure == STRUCTURE_MASK).sum().item())
  steps = max(1, min(n_masked, cfg['max_steps']))

  struct_gen = esm3.generate(prompt, GenerationConfig(
      track='structure', num_steps=steps, temperature=cfg['temperature']))
  decoded = esm3.decode(struct_gen)
  gen_ca = _ca(decoded)
  if template_ca is not None:
    n = min(len(gen_ca), len(template_ca))
    site = [p for p in motif_positions if p < n]
    rec['constrained_site_rmsd'] = (_kabsch_rmsd(gen_ca[site], template_ca[site])
                                    if site else float('nan'))
    rec['backbone_rmsd'] = _kabsch_rmsd(gen_ca[:n], template_ca[:n])
    if not (rec['constrained_site_rmsd'] < cfg['site_rmsd_max']):
      rec['failed_gate'] = 'structure_site_rmsd'
      return rec
    if not (rec['backbone_rmsd'] > cfg['backbone_rmsd_min']):
      rec['failed_gate'] = 'structure_too_similar_to_template'
      return rec

  seq_gen = esm3.generate(struct_gen, GenerationConfig(
      track='sequence', num_steps=steps, temperature=cfg['temperature']))
  designed = esm3.decode(seq_gen)
  rec['sequence'] = designed.sequence

  refolded = esm3.decode(esm3.generate(
      esm3.encode(ESMProtein(sequence=designed.sequence)),
      GenerationConfig(track='structure', num_steps=1, temperature=0.0)))
  if template_ca is not None:
    re_ca = _ca(refolded)
    n = min(len(re_ca), len(template_ca))
    site = [p for p in motif_positions if p < n]
    rec['refold_site_rmsd'] = (_kabsch_rmsd(re_ca[site], template_ca[site])
                               if site else float('nan'))
    if not (rec['refold_site_rmsd'] < cfg['site_rmsd_max']):
      rec['failed_gate'] = 'refold_site_rmsd'
      return rec

  rec['motif_preserved'] = all(rec['sequence'][p] == aa
                               for p, aa in motif.items()
                               if p < len(rec['sequence']))
  rec['passed'] = bool(rec['motif_preserved'])
  if not rec['motif_preserved']:
    rec['failed_gate'] = 'motif_not_preserved'
  return rec


# --------------------------------------------------------------------------
# Local entrypoint (the CLI)
# --------------------------------------------------------------------------


@app.local_entrypoint()
def main(
    num_generations: int,
    preset: str = '',
    model: str = 'esm3_sm_open_v1',
    motif: str = '',
    length: int = 0,
    template_pdb: str = '',
    template_chain: str = 'A',
    structure_positions: str = '',
    num_shards: int = 1,
    temperature: float = 1.0,
    max_steps: int = 20,
    site_rmsd_max: float = 1.5,
    backbone_rmsd_min: float = 1.5,
    seed: int = 0,
    out: str = './campaign',
):
  """Run a sharded motif-scaffolding campaign on Modal GPUs and pool results."""
  vals = dict(motif=motif, length=length, template_pdb=template_pdb,
              template_chain=template_chain,
              structure_positions=structure_positions)
  if preset == 'gfp':
    for k, v in GFP_PRESET.items():
      if not vals.get(k):
        vals[k] = v
    print('using the gfp_design.ipynb preset (1qy3 template, chromophore motif)')
  elif preset:
    raise SystemExit(f'unknown preset {preset!r}')
  if not vals['motif']:
    raise SystemExit('--motif is required (or --preset gfp)')
  if not vals['length']:
    raise SystemExit('--length is required (or --preset gfp)')

  per_shard = -(-num_generations // num_shards)
  if num_generations < 100:
    print('NOTE: upstream warns some prompts need thousands of generations; '
          'a run this small is a smoke test, not a campaign.')

  cfgs = [dict(model=model, num_generations=per_shard, seed=seed + s,
               temperature=temperature, max_steps=max_steps,
               site_rmsd_max=site_rmsd_max, backbone_rmsd_min=backbone_rmsd_min,
               **vals)
          for s in range(num_shards)]

  print(f'campaign: {num_generations} generations across {num_shards} shard(s) '
        f'= {per_shard} each, on {em.GPUS["a100-40"]}...')
  results = list(campaign_shard.map(cfgs))

  outdir = pathlib.Path(out)
  outdir.mkdir(parents=True, exist_ok=True)
  designs, total, gate_totals, best_site = [], 0, {}, None
  for res in results:
    if res.get('status') != 'ok':
      continue
    total += res.get('num_generations', 0)
    for g, n in (res.get('gate_breakdown') or {}).items():
      gate_totals[g] = gate_totals.get(g, 0) + n
    dist = (res.get('rmsd_distributions') or {}).get('constrained_site_rmsd')
    if dist and (best_site is None or dist['best'] < best_site):
      best_site = dist['best']
    (outdir / f'shard{res["seed"]}_result.json').write_text(
        json.dumps(res, indent=2))
    for d in res.get('passed_designs', []):
      d['seed'] = res['seed']
      designs.append(d)

  designs.sort(key=lambda d: d.get('refold_site_rmsd', float('inf')))
  (outdir / 'passed_designs.json').write_text(json.dumps(designs, indent=2))
  (outdir / 'all_designs.fasta').write_text(''.join(
      f'>seed{d.get("seed")}_design{d.get("index")} '
      f'site_rmsd={d.get("refold_site_rmsd")}\n{d["sequence"]}\n'
      for d in designs if d.get('sequence')))

  rate = len(designs) / total if total else 0.0
  print(f'\ncampaign yield: {len(designs)}/{total} passed ({100 * rate:.2f}%)')
  print(f'gate breakdown: {json.dumps(gate_totals)}')
  if not designs:
    print('No design passed both gates -- a legitimate result. The gate '
          'breakdown says which constraint binds.')
    if best_site is not None:
      print(f'Best constrained-site RMSD across the campaign: {best_site:.2f} A. '
            f'Report it -- it says whether more sampling would plausibly help.')
    return
  print(f'{"rank":>4}  {"seed":>5}  {"site RMSD":>9}  {"identity":>8}')
  for i, d in enumerate(designs[:10], 1):
    ident = d.get('identity_to_template')
    print(f'{i:>4}  {d.get("seed"):>5}  {d.get("refold_site_rmsd", 0):>9.2f}  '
          f'{"n/a" if ident is None else f"{ident:>8.1%}"}')
