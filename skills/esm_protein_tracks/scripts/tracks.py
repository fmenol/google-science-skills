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

"""Extracts, inspects and masks ESM3's conditioning tracks from a structure.

ESM3 is promptable on five parallel, residue-aligned tracks: sequence, 3D
coordinates, secondary structure (SS8), solvent accessibility (SASA) and
function annotations. This script prepares those tracks from a real experimental
structure and builds masked prompts for `esm3-protein-design` and
`esm3-inverse-folding`. It never predicts a structure and never designs a
protein.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "biotite",
#   "matplotlib",
#   "numpy",
#   "polite-http",
#   "python-dotenv",
# ]
# ///

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any

import biotite.structure as bs
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402  (backend must be set first)
import numpy as np  # noqa: E402
from polite_http import http_client  # noqa: E402
from polite_http.http_client import HttpError  # noqa: E402

from esm_biohub import BiohubError  # noqa: E402  (vendored, same directory)
import esm_biohub as eb  # noqa: E402

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

RCSB_DOWNLOAD_URL = 'https://files.rcsb.org/download'
RCSB_QPS = 2.0  # RCSB is a public good; stay well under any plausible limit.

# biotite's `annotate_sse` implements P-SEA, which is an SS3 classifier. Its
# three states map onto the SS8 alphabet ('GHITEBSC') as follows. G/I/T/B/S are
# therefore NEVER produced here -- this is an approximation, not DSSP.
SS3_TO_SS8 = {'a': 'H', 'b': 'E', 'c': 'C'}

# SS8 letters grouped for reporting.
SS8_HELIX = 'GHI'
SS8_STRAND = 'EB'
SS8_COIL = 'TSC'

# Theoretical maximum solvent-accessible surface area per residue, in
# Angstrom^2. Tien et al. (2013) PLoS ONE 8(11):e80635, "Theoretical" column.
# Used ONLY to turn absolute SASA into a relative accessibility (RSA) for the
# buried/exposed call. The SASA *track* itself always stays in absolute A^2,
# which is what ESM3 was trained on.
MAX_ASA = {
    'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, 'C': 167.0,
    'E': 223.0, 'Q': 225.0, 'G': 104.0, 'H': 224.0, 'I': 197.0,
    'L': 201.0, 'K': 236.0, 'M': 224.0, 'F': 240.0, 'P': 159.0,
    'S': 155.0, 'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0,
}
RSA_BURIED = 0.25   # rsa <  0.25          -> buried
RSA_EXPOSED = 0.50  # rsa >  0.50          -> exposed; in between -> intermediate

# Canonical track names, plus the aliases we accept on the command line.
TRACKS = ('sequence', 'coordinates', 'secondary_structure', 'sasa')
TRACK_ALIASES = {
    'seq': 'sequence',
    'sequence': 'sequence',
    'coords': 'coordinates',
    'coordinate': 'coordinates',
    'coordinates': 'coordinates',
    'structure': 'coordinates',
    'ss': 'secondary_structure',
    'ss8': 'secondary_structure',
    'secondary_structure': 'secondary_structure',
    'sasa': 'sasa',
}

SS_METHOD = (
    'biotite.structure.annotate_sse (P-SEA); an SS3 approximation '
    '(helix/strand/coil) written into the SS8 alphabet as H/E/C. NOT DSSP: the '
    'states G, I, T, B and S are never emitted.'
)


# --------------------------------------------------------------------------
# Structure input
# --------------------------------------------------------------------------


def fetch_pdb(pdb_id: str, cache_dir: str) -> pathlib.Path:
  """Downloads a PDB entry from RCSB, caching it on disk.

  Args:
    pdb_id: A four-character PDB accession, e.g. '1UBQ'.
    cache_dir: Directory to cache the downloaded file in.

  Returns:
    Path to the cached .pdb file.

  Raises:
    BiohubError: If the identifier is malformed or the download fails and no
      cached copy exists.
  """
  code = pdb_id.strip().upper()
  if not code.isalnum() or len(code) != 4:
    raise BiohubError(
        f'{pdb_id!r} is not a 4-character PDB ID (e.g. 1UBQ). To read a local '
        'file use --pdb instead.'
    )

  cache = pathlib.Path(cache_dir).expanduser()
  cache.mkdir(parents=True, exist_ok=True)
  path = cache / f'{code}.pdb'
  # Trust the cache only if it still looks like a PDB. A previous run that was
  # served an error page or a truncated body must not poison every future run,
  # so re-download rather than parsing garbage. This also means that once a
  # structure is cached the skill keeps working with RCSB offline.
  if path.is_file() and _looks_like_pdb(path.read_text(encoding='utf-8')):
    print(f'Using cached structure: {path}')
    return path

  url = f'{RCSB_DOWNLOAD_URL}/{code}.pdb'
  client = http_client.HttpClient(RCSB_DOWNLOAD_URL, qps=RCSB_QPS, timeout=120.0)
  try:
    text = client.fetch_text(url)
  except HttpError as exc:
    status = getattr(exc, 'status_code', None) or getattr(exc, 'status', None)
    hint = ''
    if status == 404:
      hint = (
          f'\nHint: RCSB has no PDB-format file for {code}. Very large entries '
          'are mmCIF-only; download the .cif yourself and pass --pdb.'
      )
    raise BiohubError(f'Could not download {url} (HTTP {status}).{hint}') from exc

  if not _looks_like_pdb(text):
    raise BiohubError(
        f'{url} did not return a PDB file (no ATOM records). Refusing to cache '
        'it. Check the accession, or download the structure yourself and pass '
        '--pdb.'
    )
  path.write_text(text, encoding='utf-8')
  print(f'Downloaded {code} from RCSB -> {path}')
  return path


def _looks_like_pdb(text: str) -> bool:
  """True if the text carries at least one ATOM record."""
  return any(line.startswith('ATOM  ') for line in text.splitlines())


def structure_path(args) -> tuple[pathlib.Path, str | None]:
  """Resolves --pdb / --pdb-id into a local path plus the PDB ID if known."""
  if args.pdb and args.pdb_id:
    raise BiohubError('Pass either --pdb or --pdb-id, not both.')
  if args.pdb:
    path = pathlib.Path(args.pdb).expanduser()
    if not path.is_file():
      raise BiohubError(f'Structure file not found: {path}')
    return path, None
  if args.pdb_id:
    return fetch_pdb(args.pdb_id, args.cache_dir), args.pdb_id.strip().upper()
  raise BiohubError('One of --pdb or --pdb-id is required.')


# --------------------------------------------------------------------------
# Track computation
# --------------------------------------------------------------------------


def build_atom_array(sequence: str, coordinates: np.ndarray) -> bs.AtomArray:
  """Builds a biotite AtomArray from an atom37 array.

  Residues are renumbered 1..L so that `res_id - 1` indexes the track directly.
  This mirrors `ProteinChain.atom_array_no_insertions` in the ESM SDK, which is
  what its `sasa()` and secondary-structure helpers consume -- so the SS8 and
  SASA tracks stay exactly aligned with the coordinate track by construction,
  with no residue-key matching to get wrong.

  Args:
    sequence: Length-L one-letter sequence.
    coordinates: (L, 37, 3) atom37 coordinates; NaN marks an absent atom.

  Returns:
    An AtomArray holding only the atoms that are actually present.
  """
  present = np.all(np.isfinite(coordinates), axis=-1)  # (L, 37)
  atoms = []
  for res_idx, residue in enumerate(sequence):
    for atom_idx in np.where(present[res_idx])[0]:
      name = eb.ATOM37[atom_idx]
      atoms.append(
          bs.Atom(
              coord=coordinates[res_idx, atom_idx],
              chain_id='A',
              res_id=res_idx + 1,
              res_name=eb.AA1_TO_AA3.get(residue, 'UNK'),
              hetero=False,
              atom_name=name,
              element=name[0],
          )
      )
  if not atoms:
    raise BiohubError('The structure has no resolved backbone atoms.')
  return bs.array(atoms)


def annotate_ss8(atom_array: bs.AtomArray, length: int) -> str:
  """Returns the SS8 track: an SS3 approximation from biotite's P-SEA.

  Residues that P-SEA cannot classify (no CA atom) are left masked ('_').
  """
  track = [eb.MASK_CHAR] * length
  try:
    sse = bs.annotate_sse(atom_array)
  except Exception as exc:  # pylint: disable=broad-except
    print(
        f'[warning] Secondary-structure annotation failed ({exc}); the SS8 '
        'track will be fully masked.',
        file=sys.stderr,
    )
    return ''.join(track)

  # P-SEA classifies CA atoms, so the annotation is aligned to the residues that
  # carry one. Scatter it back into the full-length track rather than assuming a
  # 1:1 correspondence with every residue.
  carriers = atom_array.res_id[atom_array.atom_name == 'CA']
  if len(sse) != len(carriers):
    residues = bs.get_residues(atom_array)[0]
    if len(sse) != len(residues):
      raise BiohubError(
          f'annotate_sse returned {len(sse)} states for {len(carriers)} CA '
          f'atoms / {len(residues)} residues; cannot align the SS8 track.'
      )
    carriers = residues

  for res_id, state in zip(carriers, sse):
    track[int(res_id) - 1] = SS3_TO_SS8.get(str(state), eb.MASK_CHAR)
  return ''.join(track)


def residue_sasa(
    atom_array: bs.AtomArray, coordinates: np.ndarray, length: int
) -> np.ndarray:
  """Returns the SASA track: per-residue solvent-accessible area in Angstrom^2.

  Residues with no resolved atoms get NaN (which serialises to `null`, the mask
  value for the SASA track).
  """
  per_atom = bs.sasa(atom_array)
  # biotite yields NaN for any atom it excludes from the surface calculation
  # (e.g. an atom with no van der Waals radius). Those contribute nothing.
  per_atom = np.nan_to_num(np.asarray(per_atom, dtype=np.float64), nan=0.0)

  # res_id is 1-based, so bin 0 is always empty; `minlength` keeps the array
  # full length even when the trailing residues are unresolved.
  totals = np.bincount(
      atom_array.res_id, weights=per_atom, minlength=length + 1
  )[1 : length + 1]

  resolved = np.any(np.all(np.isfinite(coordinates), axis=-1), axis=1)  # (L,)
  totals[~resolved] = np.nan
  return totals


def predict_function(sequence: str, model: str) -> list[list[Any]]:
  """Asks ESM3 for the function track. Requires BIOHUB_API_KEY.

  Returns:
    A list of [label, start, end] annotations. `start`/`end` are 1-indexed and
    inclusive, matching the convention used everywhere else in this script.
  """
  client = eb.BiohubClient()
  out = client.generate('function', model=model, sequence=sequence, num_steps=1)
  annotations = (out.get('outputs') or {}).get('function') or []
  return [[str(a[0]), int(a[1]), int(a[2])] for a in annotations if len(a) >= 3]


# --------------------------------------------------------------------------
# Summaries
# --------------------------------------------------------------------------


def summarize(
    sequence: str,
    coordinates: np.ndarray,
    ss8: str,
    sasa: np.ndarray,
    function_annotations: list[list[Any]] | None,
) -> dict[str, Any]:
  """Builds the per-track summary block."""
  length = len(sequence)
  present = np.all(np.isfinite(coordinates), axis=-1)  # (L, 37)

  counts = {letter: ss8.count(letter) for letter in eb.SS8_VOCAB}
  counts[eb.MASK_CHAR] = ss8.count(eb.MASK_CHAR)
  helix = sum(ss8.count(c) for c in SS8_HELIX)
  strand = sum(ss8.count(c) for c in SS8_STRAND)
  coil = sum(ss8.count(c) for c in SS8_COIL)

  finite = sasa[np.isfinite(sasa)]
  rsa = relative_sasa(sequence, sasa)
  rsa_finite = rsa[np.isfinite(rsa)]

  # A residue whose backbone N, CA and C are all present.
  backbone = [eb.ATOM37.index(a) for a in ('N', 'CA', 'C')]
  backbone_complete = int(np.all(present[:, backbone], axis=1).sum())
  missing = [i + 1 for i in np.where(~present.any(axis=1))[0]]

  return {
      'length': length,
      'sequence': {
          'masked': sequence.count(eb.MASK_CHAR),
          'unmasked': length - sequence.count(eb.MASK_CHAR),
      },
      'coordinates': {
          'atoms_present': int(present.sum()),
          'residues_with_no_atoms': missing,
          'backbone_complete_residues': backbone_complete,
          'masked_residues': length - int(present.any(axis=1).sum()),
      },
      'secondary_structure': {
          'method': SS_METHOD,
          'helix_pct': round(100.0 * helix / length, 1) if length else 0.0,
          'strand_pct': round(100.0 * strand / length, 1) if length else 0.0,
          'coil_pct': round(100.0 * coil / length, 1) if length else 0.0,
          'masked': counts[eb.MASK_CHAR],
          'counts': {k: v for k, v in counts.items() if v},
      },
      'sasa': {
          'units': 'Angstrom^2',
          'mean': _round(np.mean(finite)) if finite.size else None,
          'median': _round(np.median(finite)) if finite.size else None,
          'min': _round(np.min(finite)) if finite.size else None,
          'max': _round(np.max(finite)) if finite.size else None,
          'buried': int((rsa_finite < RSA_BURIED).sum()),
          'intermediate': int(
              ((rsa_finite >= RSA_BURIED) & (rsa_finite <= RSA_EXPOSED)).sum()
          ),
          'exposed': int((rsa_finite > RSA_EXPOSED).sum()),
          'masked': int(np.isnan(sasa).sum()),
          'classifier': (
              f'relative accessibility (Tien 2013 maxima); buried rsa < '
              f'{RSA_BURIED}, exposed rsa > {RSA_EXPOSED}'
          ),
      },
      'function': {
          'annotations': len(function_annotations)
          if function_annotations is not None
          else None,
      },
  }


def relative_sasa(sequence: str, sasa: np.ndarray) -> np.ndarray:
  """Converts absolute SASA (A^2) to relative accessibility in [0, ~1]."""
  maxima = np.array(
      [MAX_ASA.get(residue, np.nan) for residue in sequence], dtype=np.float64
  )
  with np.errstate(divide='ignore', invalid='ignore'):
    return np.asarray(sasa, dtype=np.float64) / maxima


def _round(value: Any, digits: int = 2) -> float:
  return round(float(value), digits)


# --------------------------------------------------------------------------
# Serialisation
# --------------------------------------------------------------------------


def npz_path_for(json_path: str | pathlib.Path) -> pathlib.Path:
  return pathlib.Path(json_path).expanduser().with_suffix('.npz')


def write_tracks(payload: dict[str, Any], coordinates: np.ndarray,
                 sasa: np.ndarray, output: str) -> pathlib.Path:
  """Writes the .json + companion .npz pair and returns the .npz path."""
  npz = npz_path_for(output)
  payload['arrays'] = npz.name
  eb.write_json(payload, output)
  npz.parent.mkdir(parents=True, exist_ok=True)
  np.savez_compressed(npz, coordinates=coordinates, sasa=sasa)
  print(f'Success! Arrays written to: {npz}')
  return npz


def load_tracks(path: str) -> dict[str, Any]:
  """Loads a tracks/prompt .json plus its companion .npz.

  Returns:
    The JSON payload with `coordinates` and `sasa` restored as numpy arrays.
  """
  json_path = pathlib.Path(path).expanduser()
  if not json_path.is_file():
    raise BiohubError(f'Tracks file not found: {json_path}')
  with open(json_path, encoding='utf-8') as handle:
    payload = json.load(handle)

  # Prefer the sibling with the same stem, so that renaming the pair together
  # keeps working; fall back to the name recorded at write time.
  candidates = [npz_path_for(json_path)]
  recorded = payload.get('arrays')
  if recorded:
    candidates.append(json_path.parent / recorded)
  npz = next((p for p in candidates if p.is_file()), None)
  if npz is None:
    raise BiohubError(
        f'Companion array file not found (looked for '
        f'{", ".join(str(p) for p in candidates)}). It is written next to the '
        'JSON by `extract` / `build-prompt`; keep the two files together, and '
        'if you rename one, rename the other to match.'
    )
  with np.load(npz) as data:
    payload['coordinates'] = data['coordinates']
    payload['sasa'] = data['sasa']
  return payload


# --------------------------------------------------------------------------
# Range parsing (1-INDEXED, INCLUSIVE)
# --------------------------------------------------------------------------


def parse_ranges(spec: str, length: int, residue_ids: list[int],
                 index: str) -> np.ndarray:
  """Parses '10-20,45-50' into a boolean selection mask of length L.

  Ranges are ALWAYS 1-indexed and INCLUSIVE, the convention biologists use:
  '10-20' selects 11 residues, the 10th through the 20th. A bare '30' selects
  one residue.

  Args:
    spec: Comma-separated ranges.
    length: Track length L.
    residue_ids: Author residue numbers from the structure, in track order.
    index: 'position' -> ranges count 1..L along the extracted sequence.
      'residue-id' -> ranges are author residue numbers from the PDB file.

  Returns:
    A boolean array of shape (L,), True where selected.

  Raises:
    BiohubError: On malformed, inverted or out-of-range input.
  """
  selected = np.zeros(length, dtype=bool)
  if index == 'position':
    lookup = {i + 1: [i] for i in range(length)}
    label = 'track position'
    bounds = f'1..{length}'
  else:
    lookup = {}
    for i, res_id in enumerate(residue_ids):
      lookup.setdefault(int(res_id), []).append(i)
    label = 'author residue id'
    bounds = f'{min(residue_ids)}..{max(residue_ids)}' if residue_ids else 'n/a'

  for chunk in spec.split(','):
    chunk = chunk.strip()
    if not chunk:
      continue
    if '-' in chunk.lstrip('-'):
      # Split on the last '-' that is not a leading sign, so negatives error out
      # below rather than parsing silently.
      first, _, last = chunk.partition('-')
      try:
        start, end = int(first), int(last)
      except ValueError as exc:
        raise BiohubError(
            f'Malformed range {chunk!r}. Expected START-END, 1-indexed and '
            'inclusive, e.g. 10-20.'
        ) from exc
    else:
      try:
        start = end = int(chunk)
      except ValueError as exc:
        raise BiohubError(
            f'Malformed range {chunk!r}. Expected START-END or a single '
            'residue number.'
        ) from exc

    if start > end:
      raise BiohubError(
          f'Inverted range {chunk!r}: start {start} > end {end}.'
      )
    hits = 0
    for value in range(start, end + 1):
      for i in lookup.get(value, []):
        selected[i] = True
        hits += 1
    if hits == 0:
      # Only suggest --index residue-id when it would actually change anything,
      # i.e. when the file's author numbering is not simply 1..L.
      hint = ''
      if index == 'position' and residue_ids and residue_ids[0] != 1:
        hint = (
            f' Note this structure numbers its residues '
            f'{min(residue_ids)}..{max(residue_ids)}; pass --index residue-id '
            'to address that author numbering instead.'
        )
      raise BiohubError(
          f'Range {chunk!r} selects nothing. Ranges are 1-indexed inclusive '
          f'and are read as {label} (valid: {bounds}).{hint}'
      )
  if not selected.any():
    raise BiohubError(f'Selection {spec!r} is empty.')
  return selected


def parse_track_list(spec: str) -> list[str]:
  """Parses 'sequence,coords' into canonical track names."""
  names = []
  for chunk in spec.split(','):
    chunk = chunk.strip().lower()
    if not chunk:
      continue
    if chunk not in TRACK_ALIASES:
      raise BiohubError(
          f'Unknown track {chunk!r}. Choose from: {", ".join(TRACKS)}.'
      )
    canonical = TRACK_ALIASES[chunk]
    if canonical not in names:
      names.append(canonical)
  if not names:
    raise BiohubError('--tracks selected nothing.')
  return names


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------

SS_COLORS = {'H': '#c1443c', 'E': '#e8a33d', 'C': '#9aa0a6', '_': '#ffffff'}


def _ss_group(letter: str) -> str:
  if letter in SS8_HELIX:
    return 'H'
  if letter in SS8_STRAND:
    return 'E'
  if letter in SS8_COIL:
    return 'C'
  return '_'


def plot_tracks(sequence: str, ss8: str, sasa: np.ndarray, title: str,
                output: str) -> None:
  """Draws the SS8 ribbon over the SASA profile."""
  length = len(sequence)
  x = np.arange(1, length + 1)
  rsa = relative_sasa(sequence, sasa)

  fig, (top, bottom) = plt.subplots(
      2, 1, figsize=(max(9.0, length / 12.0), 4.6), sharex=True,
      gridspec_kw={'height_ratios': [1, 4], 'hspace': 0.12},
  )

  # -- SS8 ribbon ---------------------------------------------------------
  for i, letter in enumerate(ss8):
    group = _ss_group(letter)
    top.add_patch(
        plt.Rectangle(
            (i + 0.5, 0.0), 1.0, 1.0,
            facecolor=SS_COLORS[group],
            edgecolor='#c8ccd0' if group == '_' else 'none',
            hatch='///' if group == '_' else None,
            linewidth=0.3,
        )
    )
  top.set_xlim(0.5, length + 0.5)
  top.set_ylim(0.0, 1.0)
  top.set_yticks([])
  top.set_ylabel('SS8', rotation=0, ha='right', va='center', fontsize=9)
  for spine in top.spines.values():
    spine.set_visible(False)
  handles = [
      plt.Rectangle((0, 0), 1, 1, facecolor=SS_COLORS[k], label=v)
      for k, v in (('H', 'helix (H)'), ('E', 'strand (E)'), ('C', 'coil (C)'))
  ]
  if eb.MASK_CHAR in ss8:
    handles.append(
        plt.Rectangle((0, 0), 1, 1, facecolor='#ffffff', edgecolor='#c8ccd0',
                      hatch='///', label='masked/unknown')
    )
  top.legend(handles=handles, loc='lower left', bbox_to_anchor=(0, 1.05),
             ncol=4, frameon=False, fontsize=8)
  top.set_title(title, fontsize=10, loc='right', pad=18)

  # -- SASA profile -------------------------------------------------------
  values = np.asarray(sasa, dtype=np.float64)
  bottom.fill_between(x, 0.0, np.nan_to_num(values, nan=0.0), where=np.isfinite(values),
                      color='#4a7fb5', alpha=0.35, linewidth=0)
  bottom.plot(x, values, color='#2c5d8f', linewidth=1.0)

  buried = np.isfinite(rsa) & (rsa < RSA_BURIED)
  exposed = np.isfinite(rsa) & (rsa > RSA_EXPOSED)
  bottom.scatter(x[buried], values[buried], s=12, color='#2f5d3a', zorder=3,
                 label=f'buried (rsa < {RSA_BURIED})')
  bottom.scatter(x[exposed], values[exposed], s=12, color='#d06a2c', zorder=3,
                 label=f'exposed (rsa > {RSA_EXPOSED})')
  missing = ~np.isfinite(values)
  if missing.any():
    for i in np.where(missing)[0]:
      bottom.axvspan(i + 0.5, i + 1.5, color='#c8ccd0', alpha=0.5, linewidth=0)

  bottom.set_xlim(0.5, length + 0.5)
  peak = float(np.nanmax(values)) if np.any(np.isfinite(values)) else 1.0
  bottom.set_ylim(0.0, peak * 1.28)  # headroom so the legend clears the trace
  bottom.set_xlabel('Track position (1-indexed)')
  bottom.set_ylabel(r'SASA ($\AA^2$)')
  bottom.legend(loc='upper right', frameon=False, fontsize=8)
  bottom.spines['top'].set_visible(False)
  bottom.spines['right'].set_visible(False)

  path = pathlib.Path(output).expanduser()
  path.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(path, dpi=150, bbox_inches='tight')
  plt.close(fig)
  print(f'Success! Plot written to: {path}')


# --------------------------------------------------------------------------
# Subcommands
# --------------------------------------------------------------------------


def cmd_extract(args) -> None:
  """Extracts every track from a structure."""
  path, pdb_id = structure_path(args)
  sequence, coordinates, residue_ids = eb.parse_pdb_atom37(str(path), args.chain)
  length = len(sequence)

  atom_array = build_atom_array(sequence, coordinates)
  ss8 = annotate_ss8(atom_array, length)
  sasa = residue_sasa(atom_array, coordinates, length)

  function_annotations = None
  if args.predict_function:
    try:
      function_annotations = predict_function(sequence, args.model)
      print(f'Function track: {len(function_annotations)} annotation(s).')
    except BiohubError as exc:
      print(
          f'[warning] Function prediction failed, continuing without the '
          f'function track: {exc}',
          file=sys.stderr,
      )

  chain = args.chain or 'first chain in file'
  payload: dict[str, Any] = {
      'kind': 'esm_protein_tracks/tracks',
      'version': 1,
      'source': {
          'pdb_id': pdb_id,
          'path': str(path),
          'chain': args.chain,
      },
      'length': length,
      'sequence': sequence,
      'secondary_structure': ss8,
      'secondary_structure_method': SS_METHOD,
      'sasa': [None if np.isnan(v) else _round(v) for v in sasa],
      'residue_ids': [int(r) for r in residue_ids],
      'function_annotations': function_annotations,
      'masking': {
          'sequence': "'_' marks a masked residue",
          'coordinates': 'NaN (JSON null) marks a masked atom',
          'secondary_structure': "'_' marks a masked residue",
          'sasa': 'null marks a masked residue',
          'ranges': '1-indexed, inclusive',
      },
      'summary': summarize(sequence, coordinates, ss8, sasa,
                           function_annotations),
  }

  write_tracks(payload, coordinates, sasa, args.output)

  summary = payload['summary']
  print(
      f'Extracted {length} residues (chain {chain}) from {path.name}: '
      f'{summary["secondary_structure"]["helix_pct"]}% helix, '
      f'{summary["secondary_structure"]["strand_pct"]}% strand, '
      f'{summary["secondary_structure"]["coil_pct"]}% coil; '
      f'mean SASA {summary["sasa"]["mean"]} A^2.'
  )


def cmd_inspect(args) -> None:
  """Summarises a tracks file and plots the SS8 ribbon + SASA profile."""
  payload = load_tracks(args.tracks)
  sequence = payload['sequence']
  ss8 = payload['secondary_structure']
  coordinates = payload['coordinates']
  sasa = payload['sasa']

  summary = summarize(
      sequence, coordinates, ss8, sasa, payload.get('function_annotations')
  )
  source = payload.get('source') or {}
  title = source.get('pdb_id') or pathlib.Path(
      source.get('path', args.tracks)
  ).stem
  if source.get('chain'):
    title = f'{title} chain {source["chain"]}'

  plot = args.plot or str(
      pathlib.Path(args.output).expanduser().with_suffix('.png')
  )
  plot_tracks(sequence, ss8, sasa, title, plot)

  report = {
      'kind': 'esm_protein_tracks/inspection',
      'version': 1,
      'tracks_file': str(pathlib.Path(args.tracks).expanduser()),
      'plot': plot,
      'source': source,
      'summary': summary,
      'function_annotations': payload.get('function_annotations'),
  }
  eb.write_json(report, args.output)

  ss = summary['secondary_structure']
  sa = summary['sasa']
  print(
      f'{title}: {summary["length"]} residues | '
      f'SS8 {ss["helix_pct"]}% H / {ss["strand_pct"]}% E / {ss["coil_pct"]}% C'
      f' (masked {ss["masked"]}) | '
      f'SASA mean {sa["mean"]} A^2, buried {sa["buried"]} / '
      f'intermediate {sa["intermediate"]} / exposed {sa["exposed"]} '
      f'(masked {sa["masked"]}) | '
      f'residues with no coordinates: '
      f'{len(summary["coordinates"]["residues_with_no_atoms"])}'
  )


def cmd_build_prompt(args) -> None:
  """Masks the tracks and writes a prompt ready for ESM3."""
  if bool(args.keep) == bool(args.mask):
    raise BiohubError('Pass exactly one of --keep or --mask.')

  payload = load_tracks(args.tracks)
  sequence = payload['sequence']
  ss8 = payload['secondary_structure']
  coordinates = np.array(payload['coordinates'], dtype=np.float64, copy=True)
  sasa = np.array(payload['sasa'], dtype=np.float64, copy=True)
  residue_ids = [int(r) for r in payload['residue_ids']]
  length = len(sequence)

  spec = args.keep or args.mask
  selected = parse_ranges(spec, length, residue_ids, args.index)
  # `--keep` names what survives; `--mask` names what is destroyed.
  keep = selected if args.keep else ~selected

  conditioned = parse_track_list(args.condition_on)
  # A track the user did not name is emitted FULLY masked: same length, no
  # information. That is equivalent to omitting it, but keeps every prompt the
  # same shape and makes the masking counts self-documenting.
  masks = {
      name: (keep if name in conditioned else np.zeros(length, dtype=bool))
      for name in TRACKS
  }

  prompt_sequence = ''.join(
      residue if masks['sequence'][i] else eb.MASK_CHAR
      for i, residue in enumerate(sequence)
  )
  prompt_ss8 = ''.join(
      state if masks['secondary_structure'][i] else eb.MASK_CHAR
      for i, state in enumerate(ss8)
  )
  coordinates[~masks['coordinates']] = np.nan
  sasa[~masks['sasa']] = np.nan

  # Function annotations are ranges, not per-residue values: keep an annotation
  # only if it overlaps the kept region at all, and say so.
  function_annotations = payload.get('function_annotations')
  kept_function = None
  if function_annotations is not None:
    kept_ids = {i + 1 for i in np.where(keep)[0]}
    kept_function = [
        a for a in function_annotations
        if kept_ids & set(range(int(a[1]), int(a[2]) + 1))
    ]

  masked_counts = {
      'sequence': int(prompt_sequence.count(eb.MASK_CHAR)),
      'coordinates': int((~np.any(np.isfinite(coordinates), axis=(1, 2))).sum()),
      'secondary_structure': int(prompt_ss8.count(eb.MASK_CHAR)),
      'sasa': int(np.isnan(sasa).sum()),
  }
  kept_counts = {k: length - v for k, v in masked_counts.items()}

  prompt: dict[str, Any] = {
      'kind': 'esm_protein_tracks/prompt',
      'version': 1,
      'source': payload.get('source'),
      'tracks_file': str(pathlib.Path(args.tracks).expanduser()),
      'length': length,
      'selection': {
          'keep': args.keep,
          'mask': args.mask,
          'index': args.index,
          'convention': '1-indexed, inclusive',
          'conditioned_tracks': conditioned,
          'fully_masked_tracks': [t for t in TRACKS if t not in conditioned],
          'kept_positions': [i + 1 for i in np.where(keep)[0]],
          'kept_residue_ids': [residue_ids[i] for i in np.where(keep)[0]],
      },
      'masking': {
          'sequence': "'_' = masked",
          'coordinates': 'null (NaN in the .npz) = masked',
          'secondary_structure': "'_' = masked",
          'sasa': 'null = masked',
      },
      'masked_positions': masked_counts,
      'unmasked_positions': kept_counts,
      # These four keys map 1:1 onto BiohubClient.generate(...) keyword
      # arguments. `sequence`/`secondary_structure`/`sasa`/`coordinates`.
      'sequence': prompt_sequence,
      'secondary_structure': prompt_ss8,
      'sasa': [None if np.isnan(v) else _round(v) for v in sasa],
      'coordinates': _nan_to_none(coordinates),
      'function_annotations': kept_function,
      'secondary_structure_method': payload.get('secondary_structure_method'),
  }

  npz = npz_path_for(args.output)
  prompt['arrays'] = npz.name
  eb.write_json(prompt, args.output)
  npz.parent.mkdir(parents=True, exist_ok=True)
  np.savez_compressed(npz, coordinates=coordinates, sasa=sasa)
  print(f'Success! Arrays written to: {npz}')

  verb = 'keeping' if args.keep else 'masking'
  print(
      f'Prompt for {length} residues, {verb} {spec} '
      f'({args.index}, 1-indexed inclusive). Masked per track: '
      + ', '.join(f'{k} {v}/{length}' for k, v in masked_counts.items())
      + '.'
  )


def _nan_to_none(array: np.ndarray) -> Any:
  """Nested lists with NaN -> None, i.e. exactly what the ESM3 API expects."""
  out = np.asarray(array, dtype=np.float64).tolist()

  def walk(node):
    if isinstance(node, list):
      return [walk(x) for x in node]
    if isinstance(node, float) and np.isnan(node):
      return None
    return node

  return walk(out)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def add_structure_args(parser: argparse.ArgumentParser) -> None:
  parser.add_argument('--pdb', help='Path to a local .pdb / .cif structure.')
  parser.add_argument('--pdb-id', help='PDB accession to fetch from RCSB, e.g. 1UBQ.')
  parser.add_argument(
      '--chain',
      help='Author chain id to extract. Defaults to the first chain in the file.',
  )
  parser.add_argument(
      '--cache-dir',
      default='rcsb_cache',
      help='Where to cache structures downloaded from RCSB (default: ./rcsb_cache).',
  )


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Extract, inspect and mask ESM3's conditioning tracks (sequence, "
          'coordinates, SS8, SASA, function) from an experimental structure.'
      )
  )
  sub = parser.add_subparsers(dest='command', required=True)

  p = sub.add_parser(
      'extract', help='Pull every track out of a PDB/mmCIF structure.'
  )
  add_structure_args(p)
  p.add_argument(
      '--predict-function',
      action='store_true',
      help=(
          'Also fill the function track by calling ESM3 (needs BIOHUB_API_KEY). '
          'Every other track is computed locally and needs no key. For a proper '
          'function analysis use the esm3-function-prediction skill instead.'
      ),
  )
  p.add_argument('--model', default=eb.DEFAULT_ESM3,
                 help=f'ESM3 model for --predict-function (default: {eb.DEFAULT_ESM3}).')
  p.add_argument('--output', required=True, help='Output tracks JSON path.')
  p.set_defaults(func=cmd_extract)

  p = sub.add_parser(
      'inspect', help='Summarise a tracks file and plot SS8 + SASA.'
  )
  p.add_argument('--tracks', required=True, help='A tracks JSON from `extract`.')
  p.add_argument('--output', required=True, help='Output summary JSON path.')
  p.add_argument('--plot', help='Output PNG path (default: alongside --output).')
  p.set_defaults(func=cmd_inspect)

  p = sub.add_parser(
      'build-prompt',
      help='Mask the tracks into a prompt for ESM3. Ranges are 1-indexed '
           'and INCLUSIVE.',
  )
  p.add_argument('--tracks', required=True, help='A tracks JSON from `extract`.')
  p.add_argument(
      '--keep',
      help='Ranges to KEEP, everything else is masked. 1-indexed inclusive, '
           'e.g. "10-20,45-50".',
  )
  p.add_argument(
      '--mask',
      help='Ranges to MASK, everything else is kept. 1-indexed inclusive, '
           'e.g. "30-40".',
  )
  p.add_argument(
      '--index',
      choices=('position', 'residue-id'),
      default='position',
      help='Read ranges as 1..L track positions (default) or as author residue '
           'numbers from the structure file.',
  )
  p.add_argument(
      '--condition-on',
      default=','.join(TRACKS),
      help='Tracks the selection applies to (default: all four). Any track not '
           'listed is emitted FULLY masked, i.e. it conditions on nothing. '
           'Example for inverse folding: --condition-on coordinates.',
  )
  p.add_argument('--output', required=True, help='Output prompt JSON path.')
  p.set_defaults(func=cmd_build_prompt)

  args = parser.parse_args()
  try:
    args.func(args)
  except BiohubError as exc:
    print(f'Error: {exc}', file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
