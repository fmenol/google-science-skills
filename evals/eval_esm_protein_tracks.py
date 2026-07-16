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

"""Eval for the esm_protein_tracks skill.

Ground truth is PDB 1UBQ — the classic 1.8 A ubiquitin structure, chain A.

Structural checks
  * `extract` writes a .json + .npz pair; coords are (76, 37, 3).
  * `inspect` writes a PNG plot and a summary JSON.
  * `build-prompt --keep 10-20` masks exactly 65 of 76 positions on every
    track, which pins down the 1-indexed-inclusive range conversion.

Scientific checks (exact ground truth — do not weaken)
  * The extracted sequence EQUALS the UBIQUITIN fixture, character for
    character. 1UBQ chain A is human ubiquitin; any mismatch means the PDB
    parser or the chain selection is broken.
  * The SS8 track contains BOTH H and E: ubiquitin is a beta-grasp fold, an
    alpha-helix (~23-34) packed against a mixed beta-sheet. Neither may be
    missing.
  * SASA is finite and non-negative everywhere, has a positive mean, and shows
    a real hydrophobic core (near-zero minimum) alongside fully exposed
    residues (a much larger maximum).

Only `extract --predict-function` would need BIOHUB_API_KEY, and this eval does
not use it: track extraction is pure local computation.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy",
# ]
# ///

from __future__ import annotations

import itertools
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from fixtures import Eval, UBIQUITIN, load_json, run_script  # noqa: E402

SKILL = 'esm_protein_tracks'

# Persistent so that re-runs reuse the cached 1UBQ download instead of hammering
# RCSB, and so the eval still works if RCSB is unreachable on a later run.
WORK = pathlib.Path('/tmp/esm_protein_tracks_eval')
CACHE = WORK / 'rcsb_cache'

PDB_ID = '1UBQ'
CHAIN = 'A'
LENGTH = 76

# `--keep 10-20` is 1-indexed and INCLUSIVE, so it keeps 11 residues (the 10th
# through the 20th) and masks the other 65. Getting this off by one is the most
# likely bug in the whole skill, so it is asserted on every track.
KEEP = '10-20'
KEPT = 11
MASKED = LENGTH - KEPT  # 65


def extract() -> pathlib.Path:
  """Runs `extract` on 1UBQ.

  `--cache-dir` points at a persistent directory, and the skill validates a
  cached file before trusting it. So the first run downloads 1UBQ from RCSB and
  every later run reads the cached copy — the eval keeps working with RCSB
  offline, without ever weakening the sequence-equality assertion below.
  """
  out = WORK / 'tracks.json'
  run_script(
      SKILL, 'tracks.py', 'extract',
      '--pdb-id', PDB_ID, '--chain', CHAIN,
      '--cache-dir', str(CACHE), '--output', str(out),
  )
  return out


def main() -> int:
  ev = Eval(SKILL)
  WORK.mkdir(parents=True, exist_ok=True)

  # ---------------------------------------------------------------- extract
  tracks_json = extract()
  tracks_npz = tracks_json.with_suffix('.npz')

  ev.check('extract writes tracks JSON', tracks_json.is_file(), str(tracks_json))
  ev.check('extract writes tracks NPZ', tracks_npz.is_file(), str(tracks_npz))

  tracks = load_json(tracks_json)
  with np.load(tracks_npz) as data:
    coords = data['coordinates']
    sasa = data['sasa']

  ev.check(
      'coordinates shape is (76, 37, 3)',
      coords.shape == (LENGTH, 37, 3),
      f'got {coords.shape}',
  )
  ev.check(
      'residue_ids has 76 entries',
      len(tracks['residue_ids']) == LENGTH,
      f'got {len(tracks["residue_ids"])}',
  )

  # ------------------------------------------------- SCIENTIFIC: sequence
  # 1UBQ chain A is human ubiquitin. This is a literal biological fact, so the
  # comparison is exact — no identity threshold, no alignment.
  sequence = tracks['sequence']
  ev.check(
      'extracted sequence EQUALS the UBIQUITIN fixture (76 aa, exact)',
      sequence == UBIQUITIN,
      f'len {len(sequence)}; {"exact match" if sequence == UBIQUITIN else sequence}',
  )

  # ------------------------------------------------------ SCIENTIFIC: SS8
  ss8 = tracks['secondary_structure']
  ev.check('SS8 track has length 76', len(ss8) == LENGTH, f'got {len(ss8)}')
  helix = ss8.count('H')
  strand = ss8.count('E')
  ev.check(
      "SS8 contains helix ('H') — ubiquitin's alpha-helix (~23-34)",
      helix > 0,
      f'{helix} H residues',
  )
  ev.check(
      "SS8 contains strand ('E') — ubiquitin's mixed beta-sheet",
      strand > 0,
      f'{strand} E residues',
  )
  ev.check(
      'SS8 uses only the SS8 alphabet (GHITEBSC) plus the mask char',
      set(ss8) <= set('GHITEBSC_'),
      f'alphabet used: {"".join(sorted(set(ss8)))}',
  )
  # The helix must be a contiguous element, not 'H' scattered through noise.
  longest_helix = max(
      (len(list(run)) for state, run in itertools.groupby(ss8) if state == 'H'),
      default=0,
  )
  ev.check(
      'the helix is a contiguous run of >= 8 residues (a real alpha-helix)',
      longest_helix >= 8,
      f'longest H run = {longest_helix}',
  )

  # ----------------------------------------------------- SCIENTIFIC: SASA
  ev.check(
      'SASA is finite everywhere (no NaN — 1UBQ is fully resolved)',
      bool(np.all(np.isfinite(sasa))),
      f'{int(np.isnan(sasa).sum())} NaN of {sasa.size}',
  )
  ev.check(
      'SASA is non-negative everywhere',
      bool(np.all(sasa >= 0.0)),
      f'min {float(np.min(sasa)):.2f}',
  )
  ev.check('SASA mean > 0', float(np.mean(sasa)) > 0.0,
           f'mean {float(np.mean(sasa)):.1f} A^2')
  # A folded protein has a hydrophobic core: some residues are fully occluded.
  ev.check(
      'buried residues exist (min SASA near 0) — the protein has a core',
      float(np.min(sasa)) < 5.0,
      f'min {float(np.min(sasa)):.2f} A^2',
  )
  # ...and a surface: the most exposed residue is far more accessible.
  ev.check(
      'exposed residues exist (max SASA >> min) — the protein has a surface',
      float(np.max(sasa)) > 100.0,
      f'max {float(np.max(sasa)):.1f} A^2',
  )
  ev.check(
      'the SASA summary counts both buried and exposed residues',
      tracks['summary']['sasa']['buried'] > 0
      and tracks['summary']['sasa']['exposed'] > 0,
      f'buried {tracks["summary"]["sasa"]["buried"]}, '
      f'exposed {tracks["summary"]["sasa"]["exposed"]}',
  )

  # ---------------------------------------------------------------- inspect
  summary_json = WORK / 'inspect.json'
  plot_png = WORK / 'inspect.png'
  run_script(
      SKILL, 'tracks.py', 'inspect',
      '--tracks', str(tracks_json),
      '--output', str(summary_json),
      '--plot', str(plot_png),
  )
  ev.check('inspect writes a summary JSON', summary_json.is_file())
  ev.check(
      'inspect writes a plot PNG',
      plot_png.is_file() and plot_png.stat().st_size > 5000,
      f'{plot_png.stat().st_size if plot_png.is_file() else 0} bytes',
  )
  report = load_json(summary_json)
  ss_summary = report['summary']['secondary_structure']
  ev.check(
      'inspect reports SS composition summing to ~100%',
      abs(ss_summary['helix_pct'] + ss_summary['strand_pct']
          + ss_summary['coil_pct'] - 100.0) < 0.5,
      f'{ss_summary["helix_pct"]}% H + {ss_summary["strand_pct"]}% E + '
      f'{ss_summary["coil_pct"]}% C',
  )
  ev.check(
      'inspect reports 0 residues with missing coordinates for 1UBQ',
      report['summary']['coordinates']['residues_with_no_atoms'] == [],
      str(report['summary']['coordinates']['residues_with_no_atoms']),
  )

  # ----------------------------------------------------------- build-prompt
  # This is the precise test of the 1-indexed-inclusive conversion.
  prompt_json = WORK / 'prompt.json'
  run_script(
      SKILL, 'tracks.py', 'build-prompt',
      '--tracks', str(tracks_json),
      '--keep', KEEP,
      '--output', str(prompt_json),
  )
  ev.check('build-prompt writes a prompt JSON', prompt_json.is_file())
  prompt = load_json(prompt_json)

  # -- sequence track: '_' is the mask character.
  pseq = prompt['sequence']
  ev.check('prompt sequence has length 76', len(pseq) == LENGTH, f'{len(pseq)}')
  ev.check(
      f'prompt sequence has exactly {KEPT} UNMASKED residues',
      sum(c != '_' for c in pseq) == KEPT,
      f'{sum(c != "_" for c in pseq)} unmasked',
  )
  ev.check(
      f"prompt sequence has exactly {MASKED} '_' characters",
      pseq.count('_') == MASKED,
      f'{pseq.count("_")} masked',
  )
  # And they are the RIGHT 11: --keep 10-20 (1-indexed inclusive) == seq[9:20].
  ev.check(
      'the kept residues are exactly UBIQUITIN[9:20] (1-indexed inclusive)',
      ''.join(c for c in pseq if c != '_') == UBIQUITIN[9:20],
      f'kept {"".join(c for c in pseq if c != "_")!r}, '
      f'want {UBIQUITIN[9:20]!r}',
  )

  # -- coordinate track: NaN is the mask value.
  with np.load(prompt_json.with_suffix('.npz')) as data:
    pcoords = data['coordinates']
  ev.check(
      'prompt coordinates shape is (76, 37, 3)',
      pcoords.shape == (LENGTH, 37, 3),
      f'got {pcoords.shape}',
  )
  finite_rows = np.any(np.isfinite(pcoords), axis=(1, 2))
  ev.check(
      f'coordinate track has exactly {KEPT} rows with any finite value',
      int(finite_rows.sum()) == KEPT,
      f'{int(finite_rows.sum())} finite rows',
  )
  ev.check(
      f'coordinate track has exactly {MASKED} all-NaN rows',
      int((~finite_rows).sum()) == MASKED,
      f'{int((~finite_rows).sum())} all-NaN rows',
  )
  ev.check(
      'the finite coordinate rows are positions 10-20 (0-indexed 9..19)',
      np.array_equal(np.where(finite_rows)[0], np.arange(9, 20)),
      f'{np.where(finite_rows)[0].tolist()}',
  )

  # -- SASA track: None (JSON null) is the mask value.
  psasa = prompt['sasa']
  ev.check('prompt SASA has length 76', len(psasa) == LENGTH, f'{len(psasa)}')
  ev.check(
      f'SASA track has exactly {KEPT} non-None entries',
      sum(v is not None for v in psasa) == KEPT,
      f'{sum(v is not None for v in psasa)} non-None',
  )

  # -- SS8 track: '_' is the mask character.
  pss8 = prompt['secondary_structure']
  ev.check(
      f"SS8 track has exactly {MASKED} '_' characters",
      len(pss8) == LENGTH and pss8.count('_') == MASKED,
      f'len {len(pss8)}, {pss8.count("_")} masked',
  )

  # -- and the script must report the counts it actually produced.
  ev.check(
      'build-prompt reports the masked count per track',
      prompt['masked_positions'] == {
          'sequence': MASKED, 'coordinates': MASKED,
          'secondary_structure': MASKED, 'sasa': MASKED,
      },
      str(prompt['masked_positions']),
  )

  # -- inverse-folding prompt: keep all coordinates, mask the whole sequence.
  if_json = WORK / 'prompt_inverse_folding.json'
  run_script(
      SKILL, 'tracks.py', 'build-prompt',
      '--tracks', str(tracks_json),
      '--keep', f'1-{LENGTH}',
      '--condition-on', 'coordinates',
      '--output', str(if_json),
  )
  inverse = load_json(if_json)
  ev.check(
      'an inverse-folding prompt masks the sequence entirely and keeps all '
      'coordinates',
      inverse['masked_positions']['sequence'] == LENGTH
      and inverse['masked_positions']['coordinates'] == 0,
      f'sequence masked {inverse["masked_positions"]["sequence"]}/{LENGTH}, '
      f'coordinates masked {inverse["masked_positions"]["coordinates"]}/{LENGTH}',
  )

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
