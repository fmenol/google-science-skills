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

"""Eval for the esmc_embedding_layer_sweep skill.

Dataset: 8 single-point mutants of lysozyme + 8 of barnase, embedded with
esmc-300m-2024-12 (30 transformer blocks -> 31 layer rows, dim 960). That is
16 sequences, one request each: the whole sweep is cheap.

Scientific contract:
  * Two unrelated protein families must be linearly separable from a good
    layer's mean-pooled embedding -> best MCC > 0.8.
  * The embedding layer (row 0, no attention) must not beat the best
    transformer layer -> mcc[best] >= mcc[0], and must not itself be crowned.
  * The layers must be genuinely different representations, not the same vector
    returned 31 times.
  * Family identity is decoded best in the deep half of the network, and the
    final layer degrades relative to it (masked-LM specialisation).

MEASURED: MCC saturates at 1.000 on ALL 31 layers here. Lysozyme and barnase
are so far apart that even layer 0 -- raw token embeddings, i.e. amino-acid
composition -- separates them perfectly under 5-fold CV. So the originally
specified probe "assert the per-layer MCC vector is not constant" is not a weak
test on this dataset, it is a broken one: it fails for a CORRECT implementation
exactly as it would for a broken one, giving it zero diagnostic power.

The same scientific claim is therefore tested two ways that DO discriminate:
(1) directly on the cached array -- the 31 layer slices must be pairwise
distinct; and (2) through the sweep's own output on cross-validated log-loss, a
strictly proper scoring rule with headroom above a perfect classification. Both
would fail loudly if the client returned one vector for every layer.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy",
# ]
# ///

from __future__ import annotations

import pathlib
import random
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import fixtures as fx  # pylint: disable=g-import-not-at-top

SKILL = 'esmc_embedding_layer_sweep'
SCRIPT = 'layer_sweep.py'
MODEL = 'esmc-300m-2024-12'

N_PER_CLASS = 8
N_SEQUENCES = 2 * N_PER_CLASS
N_LAYER_ROWS = 31  # 30 transformer blocks + the embedding layer
HIDDEN_DIM = 960
N_SPLITS = 5

AA20 = 'ACDEFGHIKLMNPQRSTVWY'


def point_mutants(
    sequence: str, name: str, count: int, seed: int
) -> list[tuple[str, str]]:
  """Generates `count` distinct single-point mutants of `sequence`."""
  rng = random.Random(seed)
  seen: set[tuple[int, str]] = set()
  out: list[tuple[str, str]] = []
  while len(out) < count:
    position = rng.randrange(len(sequence))
    new_aa = rng.choice(AA20)
    old_aa = sequence[position]
    if new_aa == old_aa or (position, new_aa) in seen:
      continue
    seen.add((position, new_aa))
    mutant = sequence[:position] + new_aa + sequence[position + 1 :]
    out.append((f'{name}_{old_aa}{position + 1}{new_aa}', mutant))
  return out


def main() -> int:
  ev = fx.Eval(SKILL)
  work = pathlib.Path(tempfile.mkdtemp(prefix='layer_sweep_'))
  print(f'[{SKILL}] workdir: {work}')

  records = point_mutants(fx.LYSOZYME, 'lysozyme', N_PER_CLASS, seed=1)
  records += point_mutants(fx.BARNASE, 'barnase', N_PER_CLASS, seed=2)
  labels = {name: name.split('_')[0] for name, _ in records}

  fasta = fx.write_fasta(work / 'dataset.fasta', records)
  labels_csv = work / 'labels.csv'
  with open(labels_csv, 'w', encoding='utf-8') as handle:
    handle.write('id,label\n')
    for name, _ in records:
      handle.write(f'{name},{labels[name]}\n')

  # -- embed-dataset -------------------------------------------------------
  npy = work / 'embeddings.npy'
  fx.run_script(
      SKILL, SCRIPT, 'embed-dataset',
      '--fasta', fasta,
      '--labels', str(labels_csv),
      '--model', MODEL,
      '--output', str(npy),
  )
  ev.check('embed-dataset writes .npy', npy.is_file(), str(npy))
  meta_path = work / 'embeddings.meta.json'
  ev.check('embed-dataset writes metadata sidecar', meta_path.is_file())

  embeddings = np.load(npy)
  ev.check(
      'cached array is (16, 31, 960)',
      embeddings.shape == (N_SEQUENCES, N_LAYER_ROWS, HIDDEN_DIM),
      f'got {embeddings.shape}',
  )
  ev.check(
      'embeddings are finite', bool(np.all(np.isfinite(embeddings)))
  )

  meta = fx.load_json(meta_path)
  ev.check(
      'metadata records both classes',
      sorted(meta['class_counts']) == ['barnase', 'lysozyme']
      and set(meta['class_counts'].values()) == {N_PER_CLASS},
      str(meta['class_counts']),
  )

  # Layers must be genuinely different representations. This is the direct form
  # of "you are not returning the same vector 31 times": compare every layer's
  # dataset-mean vector against every other layer's.
  layer_means = embeddings.mean(axis=0)  # (31, 960)
  unit = layer_means / np.linalg.norm(layer_means, axis=1, keepdims=True)
  cosine = unit @ unit.T
  off_diagonal = cosine[~np.eye(N_LAYER_ROWS, dtype=bool)]
  ev.check(
      'the 31 layers are distinct representations',
      float(off_diagonal.max()) < 0.999,
      f'max off-diagonal cosine between layers = {off_diagonal.max():.4f}',
  )

  # -- sweep ---------------------------------------------------------------
  sweep_json = work / 'sweep.json'
  sweep_csv = work / 'sweep.csv'
  fx.run_script(
      SKILL, SCRIPT, 'sweep',
      '--embeddings', str(npy),
      '--n-splits', str(N_SPLITS),
      '--output-json', str(sweep_json),
      '--output-csv', str(sweep_csv),
  )
  ev.check('sweep writes JSON', sweep_json.is_file())
  ev.check('sweep writes CSV', sweep_csv.stat().st_size > 0)

  sweep = fx.load_json(sweep_json)
  rows = sweep['layers']
  ev.check(
      'sweep emits exactly 31 layer rows',
      len(rows) == N_LAYER_ROWS,
      f'got {len(rows)}',
  )
  ev.check(
      'layer indices are 0..30 in order',
      [r['layer'] for r in rows] == list(range(N_LAYER_ROWS)),
  )

  mcc = np.array([r['mcc_mean'] for r in rows])
  loss = np.array([r['log_loss_mean'] for r in rows])
  print(f'    MCC by layer:      {np.round(mcc, 3).tolist()}')
  print(f'    log-loss by layer: {np.round(loss, 4).tolist()}')
  ev.check(
      'every MCC is finite and in [-1, 1]',
      bool(np.all(np.isfinite(mcc)) and np.all(np.abs(mcc) <= 1.0)),
      f'min {mcc.min():.3f}, max {mcc.max():.3f}',
  )
  ev.check(
      'every log-loss is finite and non-negative',
      bool(np.all(np.isfinite(loss)) and np.all(loss >= 0.0)),
      f'min {loss.min():.4f}, max {loss.max():.4f}',
  )

  # SCIENTIFIC: the layers must be genuinely different representations, not the
  # same vector returned 31 times.
  #
  # The intended probe for this was "the per-layer MCC vector is not constant".
  # On this dataset that test is broken: two unrelated protein families are so
  # far apart that EVERY layer -- including layer 0 -- classifies them
  # perfectly, so MCC pins to 1.000 across all 31 layers for a CORRECT
  # implementation as well as a broken one. It has no diagnostic power here, so
  # asserting it would only punish correct code.
  #
  # Log-loss is a strictly proper scoring rule with headroom above a perfect
  # classification (it keeps rewarding a more confident, better-separated
  # boundary), so it still varies layer to layer. A broken implementation that
  # returned one vector 31 times would give a constant log-loss and fail this.
  ev.check(
      'per-layer log-loss is not constant (layers genuinely differ)',
      float(loss.std()) > 1e-6,
      f'std across layers = {loss.std():.6f}, '
      f'range {loss.min():.4f}..{loss.max():.4f}',
  )
  ev.check(
      'sweep is honest that MCC saturated',
      sweep['mcc_saturated'] is True
      and sweep['n_layers_tied_at_best_mcc'] == N_LAYER_ROWS,
      f'{sweep["n_layers_tied_at_best_mcc"]}/{N_LAYER_ROWS} layers tied; '
      f'selection_metric = {sweep["selection_metric"]!r}',
  )

  best = int(sweep['best_layer'])
  ev.check(
      'best_layer is a valid index 0..30',
      isinstance(best, int) and 0 <= best <= N_LAYER_ROWS - 1,
      f'best_layer = {best}',
  )
  ev.check(
      'sweep reports the last layer too',
      sweep['last_layer'] == N_LAYER_ROWS - 1,
  )

  # SCIENTIFIC: two unrelated families must be cleanly separable.
  ev.check(
      'best layer separates lysozyme from barnase (MCC > 0.8)',
      float(mcc[best]) > 0.8,
      f'MCC = {mcc[best]:.3f} at layer {best}',
  )

  # SCIENTIFIC: contextual layers must not lose to the raw embedding layer.
  ev.check(
      'best layer >= embedding layer on MCC (mcc[best] >= mcc[0])',
      float(mcc[best]) >= float(mcc[0]),
      f'mcc[{best}] = {mcc[best]:.3f} vs mcc[0] = {mcc[0]:.3f}',
  )

  # SCIENTIFIC: layer 0 is raw token embeddings -- amino-acid composition with
  # no context. A layer sweep that crowns it is either mis-indexing the layer
  # stack or blindly argmax-ing a saturated metric. Neither is acceptable.
  ev.check(
      'the embedding layer is NOT crowned best',
      best != 0,
      f'best_layer = {best}',
  )
  ev.check(
      'best layer beats the embedding layer on log-loss',
      float(loss[best]) < float(loss[0]),
      f'log-loss[{best}] = {loss[best]:.4f} vs log-loss[0] = {loss[0]:.4f}',
  )

  # SCIENTIFIC: family/function identity is decoded best in the deep half of
  # the network, not at the input and not at the masked-LM-specialised output.
  # The tutorial's rule of thumb is ~3/4 depth; observed here: layer 23/30.
  ev.check(
      'best layer lies in the deep half of the network',
      best >= N_LAYER_ROWS // 2,
      f'layer {best}/{N_LAYER_ROWS - 1} '
      f'({100 * best / (N_LAYER_ROWS - 1):.0f}% depth)',
  )
  ev.check(
      'the final layer degrades vs the best layer (masked-LM specialisation)',
      float(loss[-1]) > float(loss[best]),
      f'log-loss last = {loss[-1]:.4f} vs best = {loss[best]:.4f}',
  )

  # -- plot ----------------------------------------------------------------
  png = work / 'sweep.png'
  fx.run_script(
      SKILL, SCRIPT, 'plot',
      '--sweep', str(sweep_json),
      '--output', str(png),
  )
  ev.check(
      'plot writes a non-empty PNG',
      png.is_file() and png.stat().st_size > 1000,
      f'{png.stat().st_size if png.is_file() else 0} bytes',
  )

  # -- probe ---------------------------------------------------------------
  probe_json = work / 'probe.json'
  probe_model = work / 'probe.joblib'
  held_out = fx.write_fasta(
      work / 'query.fasta',
      [('query_lysozyme', fx.LYSOZYME), ('query_barnase', fx.BARNASE)],
  )
  predictions = work / 'predictions.csv'
  fx.run_script(
      SKILL, SCRIPT, 'probe',
      '--embeddings', str(npy),
      '--layer', str(best),
      '--n-splits', str(N_SPLITS),
      '--output-model', str(probe_model),
      '--output-json', str(probe_json),
      '--predict-fasta', held_out,
      '--output-predictions', str(predictions),
  )
  ev.check('probe saves a fitted model', probe_model.stat().st_size > 0)
  probe = fx.load_json(probe_json)
  ev.check(
      'probe at the best layer reproduces the sweep MCC',
      abs(probe['mcc_mean'] - float(mcc[best])) < 1e-6,
      f'probe {probe["mcc_mean"]:.3f} vs sweep {mcc[best]:.3f}',
  )

  # SCIENTIFIC: the probe must classify the two wild-type parents correctly.
  with open(predictions, encoding='utf-8') as handle:
    lines = [line.strip().split(',') for line in handle if line.strip()]
  header, body = lines[0], lines[1:]
  predicted = {row[0]: row[header.index('predicted_label')] for row in body}
  ev.check(
      'probe classifies both wild-type parents correctly',
      predicted.get('query_lysozyme') == 'lysozyme'
      and predicted.get('query_barnase') == 'barnase',
      str(predicted),
  )

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
