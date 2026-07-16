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

"""Shared client for the Biohub Platform (ESM) API.

This module is VENDORED into every ESM skill's `scripts/` directory. It is the
single source of truth; `skills/esm_common/esm_biohub.py` is canonical and the
copies are kept byte-identical (see `evals/test_common_sync.py`).

Design notes
------------
* Talks raw JSON to https://biohub.ai/api/v1/*. Deliberately does NOT depend on
  the `esm` PyPI package, which requires torch + a fork of transformers and
  would download model weights. Everything here runs remotely.
* Rate limiting is cross-process (file-lock) via `polite-http`, so concurrent
  sub-agents on the same machine collectively respect the limit.
* All heavy numeric payloads come back as plain JSON nested lists; we convert to
  numpy. No torch anywhere.
"""

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import sys
import threading
from typing import Any, Iterable, Sequence

import dotenv
import numpy as np
from polite_http import http_client
from polite_http.http_client import HttpError

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

BASE_URL = 'https://biohub.ai'
API_KEY_NAME = 'BIOHUB_API_KEY'

# Default requests/second. The vendor's own SDK (ForgeBatchExecutor) opens with
# 32 concurrent requests and backs off on HTTP 429, so the service tolerates
# substantial concurrency. We stay deliberately more conservative and let
# polite-http handle 429 with exponential backoff. Override with BIOHUB_QPS.
DEFAULT_QPS = 5.0

# ESMC sequence vocabulary (esm/utils/constants/esm3.py). Padded to 64 logits.
SEQUENCE_VOCAB = [
    '<cls>', '<pad>', '<eos>', '<unk>',
    'L', 'A', 'G', 'V', 'S', 'E', 'R', 'T', 'I', 'D', 'P', 'K',
    'Q', 'N', 'F', 'Y', 'M', 'H', 'W', 'C', 'X', 'B', 'U', 'Z',
    'O', '.', '-', '|',
    '<mask>',
]
VOCAB = {tok: i for i, tok in enumerate(SEQUENCE_VOCAB)}
LOGIT_DIM = 64  # RegressionHead(d_model, 64); vocab is zero-padded to 64.

# The 20 canonical amino acids, in conventional alphabetical order.
AA20 = list('ACDEFGHIKLMNPQRSTVWY')
AA20_IDX = [VOCAB[a] for a in AA20]

MASK_CHAR = '_'  # MASK_STR_SHORT -> <mask> (id 32)
SEQUENCE_MASK_TOKEN = 32
SEQUENCE_BOS_TOKEN = 0
SEQUENCE_EOS_TOKEN = 2
CHAIN_BREAK_STR = '|'

# ESM3 structure-track special tokens (VQVAE_CODEBOOK_SIZE = 4096).
STRUCTURE_MASK_TOKEN = 4096
STRUCTURE_EOS_TOKEN = 4097
STRUCTURE_BOS_TOKEN = 4098

SS8_VOCAB = 'GHITEBSC'

# atom37 ordering (esm/utils/residue_constants.py).
ATOM37 = [
    'N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD',
    'CD1', 'CD2', 'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3',
    'NE', 'NE1', 'NE2', 'OE1', 'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2',
    'CZ3', 'NZ', 'OXT',
]
AA1_TO_AA3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN',
    'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
    'Y': 'TYR', 'V': 'VAL', 'X': 'UNK',
}

# --------------------------------------------------------------------------
# Model registry
# --------------------------------------------------------------------------

# n_layers = transformer blocks. Hidden-state tensors carry n_layers + 1 rows
# (index 0 is the embedding layer).
ESMC_MODELS = {
    'esmc-300m-2024-12': {'n_layers': 30, 'dim': 960},
    'esmc-600m-2024-12': {'n_layers': 36, 'dim': 1152},
    'esmc-6b-2024-12': {'n_layers': 80, 'dim': 2560},
}
ESM3_MODELS = {'esm3-open-2024-03', 'esm3-sm-open-v1', 'esm3-open'}
ESMFOLD2_MODELS = {'esmfold2-fast-2026-05', 'esmfold2-2026-05'}

DEFAULT_ESMC = 'esmc-600m-2024-12'
DEFAULT_ESM3 = 'esm3-open-2024-03'
DEFAULT_ESMFOLD2 = 'esmfold2-fast-2026-05'

# SAE codebooks are named '{esmc_model}-sae-layer{L}-k{k}-codebook{size}'.
DEFAULT_SAE_MODEL = 'esmc-6b-2024-12-sae-layer60-k64-codebook16384'
SAE_FEATURE_API = 'https://biohub.ai/esm/protein/api/v1alpha1/features'

# `-1` (all hidden layers at once) is rejected for ESMC 6B and every ESM3 model.
ALL_LAYERS_UNSUPPORTED = {'esmc-6b-2024-12'} | ESM3_MODELS


class BiohubError(RuntimeError):
  """Raised when the Biohub API returns an error or the request is invalid."""


class BiohubQuotaError(BiohubError):
  """Raised when the account's daily credit allowance is exhausted.

  This is a HARD stop, not a transient rate limit: retrying or lowering the
  request rate cannot fix it. The allowance resets at 00:00 UTC.
  """


# The service returns HTTP 429 for BOTH transient rate limiting and permanent
# daily-credit exhaustion. Only the response body distinguishes them, and the
# difference matters enormously: retrying a rate limit works, whereas retrying
# an exhausted quota just burns wall-clock time behind exponential backoff.
_QUOTA_MARKERS = ('credit limit', 'credits', 'quota', 'daily limit')

# Endpoints that cost credits. `encode` is free: it still returns 200 with the
# quota exhausted, so it is useless as a credit probe but safe to call freely.
BILLED_ENDPOINTS = frozenset(
    {'logits', 'fold', 'fold_all_atom', 'generate', 'forward_and_sample'}
)


# --------------------------------------------------------------------------
# Credentials
# --------------------------------------------------------------------------


def load_api_key(required: bool = True) -> str:
  """Loads BIOHUB_API_KEY from the environment or a .env file.

  Search order: process env, then ./.env walking up to the repo root, then
  ~/.env. The value is never printed.

  Args:
    required: If True, exit(1) with actionable instructions when missing.

  Returns:
    The API key, or '' when absent and `required` is False.
  """
  if os.environ.get(API_KEY_NAME):
    return os.environ[API_KEY_NAME]

  candidates = []
  here = pathlib.Path(__file__).resolve()
  for parent in [pathlib.Path.cwd(), *here.parents]:
    candidates.append(parent / '.env')
  candidates.append(pathlib.Path.home() / '.env')

  seen = set()
  for env_path in candidates:
    if env_path in seen or not env_path.is_file():
      continue
    seen.add(env_path)
    dotenv.load_dotenv(env_path, override=False)
    if os.environ.get(API_KEY_NAME):
      return os.environ[API_KEY_NAME]

  if required:
    print(
        f'Error: {API_KEY_NAME} is not set.\n\n'
        'Create an API token at https://biohub.ai/developer-console/api-keys\n'
        'then store it (your typing stays hidden):\n\n'
        f'  printf "Enter {API_KEY_NAME} (typing hidden): " && read -s val && '
        f'echo && echo "{API_KEY_NAME}=$val" >> ~/.env && echo "Saved."\n',
        file=sys.stderr,
    )
    sys.exit(1)
  return ''


# --------------------------------------------------------------------------
# Numeric helpers (numpy; no torch)
# --------------------------------------------------------------------------


def to_array(data: Any, dtype=np.float32) -> np.ndarray:
  """Converts a nested JSON list to a numpy array, mapping null -> NaN."""
  if data is None:
    raise BiohubError('Expected numeric payload but the API returned null.')
  return np.array(
      json.loads(json.dumps(data).replace('null', 'NaN')), dtype=dtype
  )


def log_softmax(x: np.ndarray, axis: int = -1) -> np.ndarray:
  """Numerically stable log-softmax."""
  x = np.asarray(x, dtype=np.float64)
  shifted = x - np.max(x, axis=axis, keepdims=True)
  return shifted - np.log(np.sum(np.exp(shifted), axis=axis, keepdims=True))


def softmax(x: np.ndarray, axis: int = -1) -> np.ndarray:
  """Numerically stable softmax."""
  return np.exp(log_softmax(x, axis=axis))


def shannon_entropy_bits(probs: np.ndarray, axis: int = -1) -> np.ndarray:
  """Shannon entropy in bits. Zero-probability terms contribute 0."""
  p = np.clip(np.asarray(probs, dtype=np.float64), 0.0, 1.0)
  with np.errstate(divide='ignore', invalid='ignore'):
    terms = np.where(p > 0, p * np.log2(p), 0.0)
  return -np.sum(terms, axis=axis)


def validate_sequence(sequence: str, allow_mask: bool = False,
                      allow_chainbreak: bool = False) -> str:
  """Validates and upper-cases a protein sequence.

  Raises:
    BiohubError: If the sequence is empty or has characters ESM cannot tokenize.
  """
  if not sequence or not sequence.strip():
    raise BiohubError('Empty protein sequence.')
  seq = ''.join(sequence.split()).upper()
  allowed = set(AA20) | set('XBUZO')
  if allow_mask:
    allowed.add(MASK_CHAR)
  if allow_chainbreak:
    allowed.add(CHAIN_BREAK_STR)
  bad = sorted(set(seq) - allowed)
  if bad:
    raise BiohubError(
        f'Sequence contains characters ESM cannot tokenize: {bad}. '
        f'Allowed: {"".join(sorted(allowed))}'
    )
  return seq


def read_fasta(path: str) -> list[tuple[str, str]]:
  """Reads a FASTA file into a list of (id, sequence) pairs."""
  records: list[tuple[str, str]] = []
  name, chunks = None, []
  with open(path, encoding='utf-8') as handle:
    for line in handle:
      line = line.strip()
      if not line:
        continue
      if line.startswith('>'):
        if name is not None:
          records.append((name, ''.join(chunks)))
        name, chunks = line[1:].split()[0] if len(line) > 1 else 'seq', []
      else:
        chunks.append(line)
  if name is not None:
    records.append((name, ''.join(chunks)))
  if not records:
    raise BiohubError(f'No FASTA records found in {path}')
  return records


def write_json(data: Any, output_file: str) -> None:
  """Writes JSON with indent=2 and prints a one-line success message."""
  path = pathlib.Path(output_file).expanduser()
  path.parent.mkdir(parents=True, exist_ok=True)
  with open(path, 'w', encoding='utf-8') as handle:
    json.dump(data, handle, indent=2, default=_json_default)
  print(f'Success! Data written to: {path}')


def _json_default(obj: Any) -> Any:
  if isinstance(obj, np.ndarray):
    return obj.tolist()
  if isinstance(obj, (np.floating, np.integer)):
    return obj.item()
  if isinstance(obj, (np.bool_,)):
    return bool(obj)
  raise TypeError(f'Not JSON serializable: {type(obj)}')


# --------------------------------------------------------------------------
# Structure writers (no biotite / no torch)
# --------------------------------------------------------------------------


def atom37_to_pdb(
    coordinates: np.ndarray,
    sequence: str,
    plddt: np.ndarray | None = None,
    chain_id: str = 'A',
) -> str:
  """Renders an atom37 coordinate array as a PDB string.

  Args:
    coordinates: (L, 37, 3). NaN marks an absent atom.
    sequence: Length-L one-letter sequence. '|' marks a chain break.
    plddt: Optional (L,) confidence in [0, 1]; written to the B-factor column
      scaled to 0-100.
    chain_id: Chain identifier for the first chain.

  Returns:
    A PDB-formatted string.
  """
  coordinates = np.asarray(coordinates, dtype=np.float64)
  if coordinates.ndim != 3 or coordinates.shape[1:] != (37, 3):
    raise BiohubError(
        f'Expected atom37 coordinates of shape (L, 37, 3), got '
        f'{coordinates.shape}.'
    )

  chain_order = [chain_id] + [c for c in _CHAIN_IDS if c != chain_id]

  # A '|' chain break may or may not occupy a coordinate row depending on the
  # endpoint; decide by comparing lengths rather than assuming.
  breaks = sequence.count(CHAIN_BREAK_STR)
  break_has_row = (
      breaks > 0 and coordinates.shape[0] == len(sequence)
  )

  lines: list[str] = []
  atom_serial = 1
  res_serial = 0
  chain_idx = 0
  coord_idx = 0

  for aa in sequence:
    if aa == CHAIN_BREAK_STR:
      lines.append('TER')
      chain_idx += 1
      res_serial = 0
      if break_has_row:
        coord_idx += 1
      continue
    if coord_idx >= coordinates.shape[0]:
      break
    res_serial += 1
    res_name = AA1_TO_AA3.get(aa, 'UNK')
    chain = chain_order[chain_idx % len(chain_order)]
    bfac = 0.0
    if plddt is not None and coord_idx < len(plddt):
      value = float(plddt[coord_idx])
      bfac = 0.0 if np.isnan(value) else value * 100.0

    for atom_idx, atom_name in enumerate(ATOM37):
      xyz = coordinates[coord_idx, atom_idx]
      if not np.all(np.isfinite(xyz)):
        continue
      lines.append(
          _pdb_line(
              record='ATOM  ',
              serial=atom_serial,
              atom_name=atom_name,
              res_name=res_name,
              chain=chain,
              res_seq=res_serial,
              xyz=xyz,
              bfactor=bfac,
              element=atom_name[0],
          )
      )
      atom_serial += 1
    coord_idx += 1

  lines.append('TER')
  lines.append('END')
  return '\n'.join(lines) + '\n'


def complex_to_pdb(complex_data: dict[str, Any]) -> str:
  """Renders a `fold_all_atom` `complex` payload as a PDB string.

  Handles proteins, nucleic acids and ligands (written as HETATM). B-factors
  carry per-token pLDDT scaled to 0-100.
  """
  required = ('sequence', 'atom_positions', 'atom_names', 'token_to_atoms',
              'chain_id')
  missing = [k for k in required if k not in complex_data]
  if missing:
    raise BiohubError(f'complex payload missing fields: {missing}')

  res_names = complex_data['sequence']
  positions = to_array(complex_data['atom_positions'], dtype=np.float64)
  atom_names = complex_data['atom_names']
  elements = complex_data.get('atom_elements') or [n[0] for n in atom_names]
  hetero = complex_data.get('atom_hetero') or [False] * len(atom_names)
  token_to_atoms = complex_data['token_to_atoms']
  chain_ids = complex_data['chain_id']
  plddt = complex_data.get('plddt')
  chain_lookup = (complex_data.get('metadata') or {}).get('chain_lookup', {})

  lines: list[str] = []
  atom_serial = 1
  res_serial: dict[str, int] = {}
  prev_chain = None

  for token_idx, (start, end) in enumerate(token_to_atoms):
    raw_chain = chain_ids[token_idx]
    chain = str(chain_lookup.get(str(raw_chain), str(raw_chain)))[:1] or 'A'
    if prev_chain is not None and chain != prev_chain:
      lines.append('TER')
    prev_chain = chain
    res_serial[chain] = res_serial.get(chain, 0) + 1

    res_name = str(res_names[token_idx])[:3].upper()
    bfac = 0.0
    if plddt is not None and token_idx < len(plddt):
      value = float(plddt[token_idx])
      bfac = 0.0 if np.isnan(value) else value * 100.0

    for atom_idx in range(int(start), int(end)):
      xyz = positions[atom_idx]
      if not np.all(np.isfinite(xyz)):
        continue
      lines.append(
          _pdb_line(
              record='HETATM' if hetero[atom_idx] else 'ATOM  ',
              serial=atom_serial,
              atom_name=str(atom_names[atom_idx]),
              res_name=res_name,
              chain=chain,
              res_seq=res_serial[chain],
              xyz=xyz,
              bfactor=bfac,
              element=str(elements[atom_idx]),
          )
      )
      atom_serial += 1

  lines.append('TER')
  lines.append('END')
  return '\n'.join(lines) + '\n'


_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'


def _pdb_line(
    *,
    record: str,
    serial: int,
    atom_name: str,
    res_name: str,
    chain: str,
    res_seq: int,
    xyz: Sequence[float],
    bfactor: float,
    element: str,
    occupancy: float = 1.00,
) -> str:
  """Formats one PDB ATOM/HETATM record with strict column alignment.

  Columns (1-indexed, per the PDB v3.3 spec):
    1-6 record, 7-11 serial, 13-16 atom name, 17 altLoc, 18-20 resName,
    22 chainID, 23-26 resSeq, 27 iCode, 31-38/39-46/47-54 x/y/z,
    55-60 occupancy, 61-66 tempFactor, 77-78 element.
  """
  name = atom_name.strip()
  # Names of 1-3 characters are indented one column; 4-character names are not.
  name4 = name[:4] if len(name) >= 4 else f' {name:<3s}'
  # Serial numbers wrap past 99999 (the field is only 5 wide).
  serial = serial % 100000
  res_seq = res_seq % 10000
  return (
      f'{record}{serial:>5d} {name4} {res_name:>3s} {chain[:1]}{res_seq:>4d}'
      f'    {xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}'
      f'{occupancy:>6.2f}{bfactor:>6.2f}          {element.strip()[:2]:>2s}'
  )


def parse_pdb_atom37(
    path: str, chain: str | None = None
) -> tuple[str, np.ndarray, list[int]]:
  """Parses a PDB/mmCIF file into (sequence, atom37 coordinates, residue ids).

  Only ATOM records for the 20 canonical residues are read; HETATM, waters and
  ligands are ignored. Missing atoms become NaN. Alternate locations other than
  'A' are skipped.

  Args:
    path: Path to a .pdb / .ent / .cif / .mmcif file.
    chain: Author chain id to extract. Defaults to the first chain encountered.

  Returns:
    (sequence, coordinates of shape (L, 37, 3), author residue numbers)

  Note:
    Residues are kept distinct by (chain, resSeq, insertion code), so insertion
    codes are handled correctly. The returned residue ids are bare ints, so a
    structure that uses insertion codes (e.g. antibody CDRs numbered 52, 52A,
    52B) yields DUPLICATE ids. Index by position when that matters.
  """
  suffix = pathlib.Path(path).suffix.lower()
  if suffix in ('.cif', '.mmcif'):
    records = _parse_mmcif_atoms(path)
  else:
    records = _parse_pdb_atoms(path)

  aa3_to_aa1 = {v: k for k, v in AA1_TO_AA3.items() if k != 'X'}
  by_residue: dict[tuple[str, int, str], dict[str, list[float]]] = {}
  order: list[tuple[str, int, str]] = []
  for chain_id, res_seq, icode, res_name, atom_name, xyz in records:
    if res_name not in aa3_to_aa1:
      continue
    if chain is not None and chain_id != chain:
      continue
    key = (chain_id, res_seq, icode)
    if key not in by_residue:
      by_residue[key] = {}
      order.append(key)
    by_residue[key].setdefault(atom_name, xyz)
    by_residue[key]['__resname__'] = res_name  # type: ignore[assignment]

  if not order:
    raise BiohubError(
        f'No standard protein residues found in {path}'
        + (f' for chain {chain!r}.' if chain else '.')
    )

  target_chain = chain if chain is not None else order[0][0]
  order = [k for k in order if k[0] == target_chain]

  sequence_chars: list[str] = []
  coords = np.full((len(order), 37, 3), np.nan, dtype=np.float64)
  res_ids: list[int] = []
  for i, key in enumerate(order):
    atoms = by_residue[key]
    res_name = atoms['__resname__']  # type: ignore[index]
    sequence_chars.append(aa3_to_aa1[res_name])  # type: ignore[index]
    res_ids.append(key[1])
    for atom_name, xyz in atoms.items():
      if atom_name == '__resname__' or atom_name not in ATOM37:
        continue
      coords[i, ATOM37.index(atom_name)] = xyz

  return ''.join(sequence_chars), coords, res_ids


def _parse_pdb_atoms(path: str):
  """Yields (chain, resseq, icode, resname, atom, xyz) from PDB ATOM records."""
  out = []
  with open(path, encoding='utf-8', errors='replace') as handle:
    for line in handle:
      if not line.startswith('ATOM'):
        continue
      altloc = line[16]
      if altloc not in (' ', 'A'):
        continue
      out.append((
          line[21],
          int(line[22:26]),
          line[26],
          line[17:20].strip().upper(),
          line[12:16].strip(),
          [float(line[30:38]), float(line[38:46]), float(line[46:54])],
      ))
  return out


def _parse_mmcif_atoms(path: str):
  """Yields (chain, resseq, icode, resname, atom, xyz) from an mmCIF atom_site."""
  out = []
  columns: dict[str, int] = {}
  in_loop = False
  with open(path, encoding='utf-8', errors='replace') as handle:
    for line in handle:
      stripped = line.strip()
      if stripped.startswith('_atom_site.'):
        columns[stripped.split('.', 1)[1].split()[0]] = len(columns)
        in_loop = True
        continue
      if not in_loop:
        continue
      if stripped.startswith('#') or stripped.startswith('loop_'):
        if columns:
          break
        continue
      if not (stripped.startswith('ATOM') or stripped.startswith('HETATM')):
        continue
      fields = stripped.split()
      if len(fields) < len(columns):
        continue

      def get(name, default=' '):
        idx = columns.get(name)
        return fields[idx] if idx is not None and idx < len(fields) else default

      if get('group_PDB') != 'ATOM':
        continue
      altloc = get('label_alt_id', '.')
      if altloc not in ('.', '?', 'A'):
        continue
      chain = get('auth_asym_id', get('label_asym_id', 'A'))
      seq_raw = get('auth_seq_id', get('label_seq_id', '0'))
      try:
        res_seq = int(seq_raw)
      except ValueError:
        continue
      icode = get('pdbx_PDB_ins_code', '?')
      icode = ' ' if icode in ('?', '.') else icode
      out.append((
          chain,
          res_seq,
          icode,
          get('label_comp_id').strip().upper(),
          get('label_atom_id').strip().strip('"'),
          [float(get('Cartn_x', 'nan')), float(get('Cartn_y', 'nan')),
           float(get('Cartn_z', 'nan'))],
      ))
  return out


def kabsch_rmsd(
    mobile: np.ndarray, target: np.ndarray
) -> tuple[float, np.ndarray]:
  """Superposes `mobile` onto `target` (Kabsch) and returns (RMSD, aligned).

  Args:
    mobile: (N, 3) coordinates to move.
    target: (N, 3) reference coordinates.

  Returns:
    (rmsd in Angstrom, the superposed mobile coordinates)
  """
  mobile = np.asarray(mobile, dtype=np.float64)
  target = np.asarray(target, dtype=np.float64)
  if mobile.shape != target.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
    raise BiohubError(
        f'kabsch_rmsd needs two (N, 3) arrays of equal shape; got '
        f'{mobile.shape} and {target.shape}.'
    )
  keep = np.all(np.isfinite(mobile), axis=1) & np.all(
      np.isfinite(target), axis=1
  )
  if keep.sum() < 3:
    raise BiohubError('Fewer than 3 shared finite atoms; cannot superpose.')
  mob, tar = mobile[keep], target[keep]

  mob_c = mob - mob.mean(axis=0)
  tar_c = tar - tar.mean(axis=0)
  u, _, vt = np.linalg.svd(mob_c.T @ tar_c)
  d = np.sign(np.linalg.det(vt.T @ u.T))
  rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
  aligned = (rot @ (mobile - mob.mean(axis=0)).T).T + tar.mean(axis=0)
  diff = aligned[keep] - tar
  rmsd = float(np.sqrt((diff**2).sum() / keep.sum()))
  return rmsd, aligned


def tm_score(mobile_ca: np.ndarray, target_ca: np.ndarray) -> float:
  """TM-score of two equal-length, already-corresponded CA traces."""
  mobile_ca = np.asarray(mobile_ca, dtype=np.float64)
  target_ca = np.asarray(target_ca, dtype=np.float64)
  _, aligned = kabsch_rmsd(mobile_ca, target_ca)
  length = len(target_ca)
  d0 = 1.24 * (max(length, 19) - 15) ** (1.0 / 3.0) - 1.8
  d0 = max(d0, 0.5)
  dist = np.linalg.norm(aligned - target_ca, axis=1)
  finite = np.isfinite(dist)
  return float(np.sum(1.0 / (1.0 + (dist[finite] / d0) ** 2)) / length)


# --------------------------------------------------------------------------
# API client
# --------------------------------------------------------------------------


class _Cassette:
  """Records and replays API responses so evals can re-run without credits.

  The Biohub account has a hard daily credit allowance, so a test suite that
  calls the live API on every run is not repeatable. Recording each response and
  replaying it thereafter is the standard fix (cf. VCR/betamax): the whole
  pipeline — payload construction, parsing, numerics, interpretation, file output
  — is still exercised end to end; only the network hop is served from disk.

  Layout: one JSON file per (endpoint, payload) key, holding a LIST of responses.
  Repeated identical calls (e.g. sampling N designs from one prompt) record N
  responses and replay them round-robin, so stochastic diversity survives.
  Writes are atomic, so parallel evals recording different keys cannot corrupt
  each other.
  """

  def __init__(self, path: str, mode: str, namespace: str = ''):
    base = pathlib.Path(path).expanduser()
    # An optional namespace isolates one logical interaction (e.g. a single CLI
    # invocation) into its own subdirectory. Without it, two invocations that
    # happen to issue an identical payload would share cassette entries; for a
    # stochastic multi-step routine (guided decoding) that cross-talk makes the
    # trajectory unreproducible. Skills whose assertions accept any valid
    # response leave it unset and benefit from cross-invocation cache reuse.
    self.dir = base / namespace if namespace else base
    self.mode = mode  # 'record' | 'replay' | 'auto'
    self.dir.mkdir(parents=True, exist_ok=True)
    self._lock = threading.Lock()
    self._cursor: dict[str, int] = {}
    self.hits = 0
    self.misses = 0

  @staticmethod
  def key(endpoint: str, payload: dict[str, Any]) -> str:
    blob = json.dumps(
        {'endpoint': endpoint, 'payload': payload},
        sort_keys=True,
        default=str,
    )
    return f'{endpoint}-{hashlib.sha256(blob.encode()).hexdigest()[:24]}'

  def load(self, key: str) -> dict[str, Any] | None:
    """Returns the next recorded response for `key`, or None if unrecorded."""
    path = self.dir / f'{key}.json'
    with self._lock:
      if not path.is_file():
        self.misses += 1
        return None
      try:
        entries = json.loads(path.read_text(encoding='utf-8'))
      except (OSError, ValueError):
        self.misses += 1
        return None
      if not entries:
        self.misses += 1
        return None
      i = self._cursor.get(key, 0)
      self._cursor[key] = i + 1
      self.hits += 1
      return entries[i % len(entries)]

  def save(self, key: str, response: dict[str, Any]) -> None:
    """Appends a response for `key`. Atomic, so concurrent recorders are safe."""
    path = self.dir / f'{key}.json'
    with self._lock:
      entries: list[Any] = []
      if path.is_file():
        try:
          entries = json.loads(path.read_text(encoding='utf-8'))
        except (OSError, ValueError):
          entries = []
      entries.append(response)
      tmp = path.with_suffix('.json.tmp')
      tmp.write_text(json.dumps(entries), encoding='utf-8')
      tmp.replace(path)


class BiohubClient:
  """Rate-limited JSON client for the Biohub Platform inference API.

  Environment:
    BIOHUB_API_KEY   Required (unless replaying a full cassette).
    BIOHUB_QPS       Requests/second (default 5).
    BIOHUB_CASSETTE  Directory of recorded responses. When set, enables
                     record/replay so evals do not burn the daily credit
                     allowance on every run.
    BIOHUB_CASSETTE_MODE
                     'auto' (default): replay a recorded response when one
                     exists, otherwise call the API and record it.
                     'replay': never call the API; a cache miss is an error.
                     'record': always call the API and (re-)record.
  """

  def __init__(
      self,
      token: str | None = None,
      base_url: str = BASE_URL,
      qps: float | None = None,
      timeout: float = 900.0,
  ):
    cassette_dir = os.environ.get('BIOHUB_CASSETTE')
    self.cassette = (
        _Cassette(
            cassette_dir,
            os.environ.get('BIOHUB_CASSETTE_MODE', 'auto').lower(),
            namespace=os.environ.get('BIOHUB_CASSETTE_NS', ''),
        )
        if cassette_dir
        else None
    )
    # A pure replay run needs no credentials at all.
    replay_only = self.cassette is not None and self.cassette.mode == 'replay'
    self.token = token or load_api_key(required=not replay_only)

    self.base_url = base_url.rstrip('/')
    if qps is None:
      qps = float(os.environ.get('BIOHUB_QPS', DEFAULT_QPS))
    self.qps = qps
    self._http = http_client.HttpClient(
        self.base_url,
        qps=qps,
        timeout=timeout,
        max_retries=6,
        # 429 is deliberately NOT retried by the transport: the service uses it
        # for daily-credit exhaustion as well as for rate limiting, and blindly
        # backing off through six attempts (up to ~8 minutes) on an exhausted
        # quota is pure waste. `post()` inspects the body and decides.
        retryable_status_codes=frozenset({500, 502, 503, 504}),
        default_headers={
            'Authorization': f'Bearer {self.token}',
            'Accept': 'application/json',
        },
    )
    self.tokens_used = 0
    self.billed_calls = 0

  # -- transport ----------------------------------------------------------

  def post(self, endpoint: str, payload: dict[str, Any]) -> dict[str, Any]:
    """POSTs JSON to /api/v1/{endpoint} and returns the parsed response.

    Raises:
      BiohubQuotaError: If the daily credit allowance is exhausted.
      BiohubError: On any other non-2xx response, carrying the server's message.
    """
    key = None
    if self.cassette is not None:
      key = _Cassette.key(endpoint, payload)
      if self.cassette.mode in ('auto', 'replay'):
        cached = self.cassette.load(key)
        if cached is not None:
          return cached
        if self.cassette.mode == 'replay':
          raise BiohubError(
              f'Cassette miss for {endpoint} (key {key}) and '
              'BIOHUB_CASSETTE_MODE=replay forbids calling the API.\n'
              'Re-record with BIOHUB_CASSETTE_MODE=auto once credits are '
              'available.'
          )

    url = f'{self.base_url}/api/v1/{endpoint}'
    try:
      data = self._http.fetch_json(url, method='POST', json_body=payload)
    except HttpError as exc:
      status = getattr(exc, 'status_code', None)
      body = _body_of(exc)
      if status == 429 and _is_quota_error(body):
        raise BiohubQuotaError(
            f'Biohub daily credit allowance exhausted (during {endpoint}).\n'
            f'  Server said: {_server_message(body)}\n'
            '  This is NOT a rate limit — retrying or lowering BIOHUB_QPS '
            'cannot help. The allowance resets at 00:00 UTC.\n'
            '  Re-run with BIOHUB_CASSETTE set to replay previously recorded '
            'responses without spending credits.'
        ) from exc
      raise BiohubError(_explain(exc, endpoint)) from exc

    if not isinstance(data, dict):
      raise BiohubError(f'{endpoint}: unexpected response type {type(data)}')
    for message in data.get('warning_messages') or []:
      print(f'[biohub warning] {message}', file=sys.stderr)
    self.tokens_used += int(data.get('tokens_used') or 0)
    if endpoint in BILLED_ENDPOINTS:
      self.billed_calls += 1
    if self.cassette is not None and key is not None:
      self.cassette.save(key, data)
    return data

  # -- ESMC ---------------------------------------------------------------

  def encode(self, sequence: str, model: str = DEFAULT_ESMC) -> list[int]:
    """Tokenizes a sequence server-side. Returns L+2 tokens (BOS ... EOS)."""
    data = self.post(
        'encode', {'inputs': {'sequence': sequence}, 'model': model}
    )
    return data['outputs']['sequence']

  def logits(
      self,
      tokens: Sequence[int],
      model: str = DEFAULT_ESMC,
      *,
      sequence: bool = False,
      return_embeddings: bool = False,
      return_mean_embedding: bool = False,
      return_hidden_states: bool = False,
      return_mean_hidden_states: bool = False,
      ith_hidden_layer: int = -1,
      sae_models: Iterable[str] | None = None,
      normalize_features: bool = True,
  ) -> dict[str, Any]:
    """Runs a forward pass and returns the requested tensors as JSON lists.

    Shapes (L = residues, so L+2 with BOS/EOS):
      logits.sequence      (L+2, 64)
      embeddings           (1, L+2, D)
      mean_embedding       (1, 1, D)
      hidden_states        (n_layers+1, 1, L+2, D)
      mean_hidden_state    (1, n_layers+1, D)
      sae_outputs[name]    {feature_indices, values: (L+2, k); shape: [L+2, C]}
    """
    if (
        ith_hidden_layer == -1
        and return_hidden_states
        and model in ALL_LAYERS_UNSUPPORTED
    ):
      raise BiohubError(
          f'{model} does not support ith_hidden_layer=-1 (all layers at once). '
          'Request a single layer, or use return_mean_hidden_states=True.'
      )
    sae_config = None
    if sae_models:
      models = list(sae_models)
      if normalize_features and any('300m' in m.lower() for m in models):
        raise BiohubError(
            'normalize_features=True is not supported for ESMC 300M SAE '
            'models. Pass normalize_features=False.'
        )
      sae_config = {'models': models, 'normalize_features': normalize_features}

    payload = {
        'model': model,
        'inputs': {'sequence': list(tokens)},
        'logits_config': {
            'sequence': sequence,
            'return_embeddings': return_embeddings,
            'return_mean_embedding': return_mean_embedding,
            'return_hidden_states': return_hidden_states,
            'return_mean_hidden_states': return_mean_hidden_states,
            'ith_hidden_layer': ith_hidden_layer,
            'sae_config': sae_config,
        },
    }
    return self.post('logits', payload)

  def resolve_layer(self, layer: int, model: str) -> int:
    """Maps a possibly-negative layer index onto the hidden-state stack.

    The stack has `n_layers + 1` rows: row 0 is the embedding layer and rows
    1..n_layers are the transformer blocks. Negative indices are interpreted
    Python-style against that stack, so `-1` is the final transformer layer and
    `-2` the second-to-last (a good blind default for downstream probes).

    This resolution MUST happen client-side. The API only special-cases
    `ith_hidden_layer=-1` (meaning "return every layer"); any other negative
    value crashes the server with HTTP 500 `tuple index out of range`.

    Raises:
      BiohubError: If the model is unknown or the index is out of range.
    """
    if model not in ESMC_MODELS:
      raise BiohubError(
          f'Cannot resolve a layer index for unknown model {model!r}. '
          f'Known ESMC models: {sorted(ESMC_MODELS)}'
      )
    n_rows = ESMC_MODELS[model]['n_layers'] + 1
    resolved = layer + n_rows if layer < 0 else layer
    if not 0 <= resolved < n_rows:
      raise BiohubError(
          f'Layer {layer} is out of range for {model}, which has {n_rows} '
          f'hidden-state rows (0 = embedding layer, '
          f'1..{n_rows - 1} = transformer blocks).'
      )
    return resolved

  def embed(
      self,
      sequence: str,
      model: str = DEFAULT_ESMC,
      *,
      layer: int | None = None,
      per_residue: bool = False,
  ) -> np.ndarray:
    """Embeds one sequence.

    Args:
      sequence: Protein sequence.
      model: An ESMC model name.
      layer: Hidden layer index (0 = embedding layer). Negative values index
        from the end, so -1 is the final transformer layer and -2 the
        second-to-last. None uses the model's final output embedding rather
        than a hidden state.
      per_residue: Return (L, D) with BOS/EOS trimmed instead of (D,).

    Returns:
      (D,) mean-pooled over residues, or (L, D) per-residue.

    Note:
      Mean pooling here EXCLUDES BOS/EOS, whereas `mean_hidden_states()` is
      pooled server-side INCLUDING them. The two therefore differ slightly
      (cosine ~0.9999). Pick one convention and use it for both training and
      inference — mixing them silently degrades a downstream probe.
    """
    sequence = validate_sequence(sequence)
    tokens = self.encode(sequence, model)

    if layer is None:
      out = self.logits(
          tokens,
          model,
          return_embeddings=per_residue,
          return_mean_embedding=not per_residue,
      )
      if per_residue:
        arr = to_array(out['embeddings'])[0]  # (L+2, D)
        return arr[1:-1]
      return to_array(out['mean_embedding']).reshape(-1)

    resolved = self.resolve_layer(layer, model)
    out = self.logits(
        tokens, model, return_hidden_states=True, ith_hidden_layer=resolved
    )
    hidden = to_array(out['hidden_states'])  # (n, 1, L+2, D)
    # The server may return the whole stack or just the requested row.
    row = resolved if hidden.shape[0] > 1 else 0
    arr = hidden[row, 0]  # (L+2, D)
    arr = arr[1:-1]
    return arr if per_residue else arr.mean(axis=0)

  def mean_hidden_states(
      self, sequence: str, model: str = DEFAULT_ESMC
  ) -> np.ndarray:
    """Returns (n_layers+1, D): every layer, mean-pooled. One request.

    Note:
      Pooled server-side INCLUDING BOS/EOS — see the note on `embed()`.
    """
    sequence = validate_sequence(sequence)
    tokens = self.encode(sequence, model)
    out = self.logits(tokens, model, return_mean_hidden_states=True)
    return to_array(out['mean_hidden_state'])[0]

  def sequence_logits(
      self, sequence: str, model: str = DEFAULT_ESMC
  ) -> np.ndarray:
    """Returns (L+2, 64) sequence logits. Index residue i at row i+1 (BOS)."""
    tokens = self.encode(validate_sequence(sequence, allow_mask=True), model)
    out = self.logits(tokens, model, sequence=True)
    return to_array(out['logits']['sequence'])

  def sae_features(
      self,
      sequence: str,
      sae_model: str = DEFAULT_SAE_MODEL,
      model: str | None = None,
      normalize_features: bool = True,
  ) -> tuple[np.ndarray, np.ndarray, int]:
    """Extracts sparse SAE features.

    Returns:
      (feature_indices (L, k), values (L, k), codebook_size) with BOS/EOS
      already trimmed.
    """
    if model is None:
      model = sae_model.split('-sae-')[0]
    sequence = validate_sequence(sequence)
    tokens = self.encode(sequence, model)
    out = self.logits(
        tokens,
        model,
        sae_models=[sae_model],
        normalize_features=normalize_features,
    )
    sae = (out.get('sae_outputs') or {}).get(sae_model)
    if sae is None:
      raise BiohubError(
          f'No SAE output returned for {sae_model!r}. Available keys: '
          f'{list((out.get("sae_outputs") or {}).keys())}'
      )
    indices = np.asarray(sae['feature_indices'], dtype=np.int64)[1:-1]
    values = np.asarray(sae['values'], dtype=np.float32)[1:-1]
    codebook = int(sae['shape'][1])
    return indices, values, codebook

  # -- ESMFold2 -----------------------------------------------------------

  def fold(
      self,
      sequence: str,
      model: str = DEFAULT_ESMFOLD2,
      *,
      num_loops: int = 20,
      num_sampling_steps: int = 100,
      include_pae: bool = False,
      msa: dict[str, Any] | None = None,
      lm_dropout: float = 0.3,
      lm_mask_pct: float | None = None,
      msa_max_depth: int | None = 1024,
      msa_column_mask_rate: float = 0.1,
      include_embeddings: bool = False,
  ) -> dict[str, Any]:
    """Folds a single chain. Returns coordinates (L, 37, 3), plddt, ptm, pae.

    pLDDT is on a 0-1 scale (multiply by 100 for the conventional scale).
    """
    if lm_mask_pct is None:
      lm_mask_pct = 0.1 if model == 'esmfold2-fast-2026-05' else 0.0
    payload = {
        'sequence': sequence,
        'msa': msa,
        'include_distogram': False,  # not implemented server-side (HTTP 422)
        'include_pae': include_pae,
        'include_pair_chains_iptm': False,
        'num_sampling_steps': num_sampling_steps,
        'num_loops': num_loops,
        'lm_dropout': lm_dropout,
        'lm_mask_pct': lm_mask_pct,
        'msa_max_depth': msa_max_depth,
        'msa_column_mask_rate': msa_column_mask_rate,
        'include_embeddings': include_embeddings,
        'model': model,
    }
    return self.post('fold', payload)

  def fold_all_atom(
      self,
      sequences: list[dict[str, Any]],
      model: str = DEFAULT_ESMFOLD2,
      *,
      num_loops: int = 20,
      num_sampling_steps: int = 100,
      include_pae: bool = False,
      covalent_bonds: list[dict[str, Any]] | None = None,
      lm_dropout: float = 0.3,
      lm_mask_pct: float | None = None,
      msa_max_depth: int | None = 1024,
      msa_column_mask_rate: float = 0.1,
      include_embeddings: bool = False,
  ) -> dict[str, Any]:
    """Folds a molecular complex (protein / dna / rna / ligand).

    Args:
      sequences: Entities, e.g.
        {'type': 'protein', 'id': 'A', 'sequence': 'MKT...'}
        {'type': 'dna', 'id': 'B', 'sequence': 'GATC'}
        {'type': 'ligand', 'id': 'L', 'ccd': ['SAH']}
        {'type': 'ligand', 'id': 'L', 'smiles': 'CC(=O)O'}
      model: An ESMFold2 model name.

    Returns:
      The raw response: `complex`, `plddt`, `ptm`, `interface_ptm`, `pae`.
    """
    if lm_mask_pct is None:
      lm_mask_pct = 0.1 if model == 'esmfold2-fast-2026-05' else 0.0
    all_atom_input: dict[str, Any] = {'sequences': sequences}
    if covalent_bonds:
      all_atom_input['covalent_bonds'] = covalent_bonds
    payload = {
        'all_atom_input': all_atom_input,
        'model': model,
        'include_distogram': False,  # not implemented server-side (HTTP 422)
        'include_pae': include_pae,
        'include_pair_chains_iptm': False,
        'num_sampling_steps': num_sampling_steps,
        'num_loops': num_loops,
        'lm_dropout': lm_dropout,
        'lm_mask_pct': lm_mask_pct,
        'msa_max_depth': msa_max_depth,
        'msa_column_mask_rate': msa_column_mask_rate,
        'include_embeddings': include_embeddings,
    }
    return self.post('fold_all_atom', payload)

  # -- ESM3 ---------------------------------------------------------------

  def generate(
      self,
      track: str,
      model: str = DEFAULT_ESM3,
      *,
      sequence: str | None = None,
      coordinates: Any = None,
      secondary_structure: str | None = None,
      sasa: list[float | None] | None = None,
      num_steps: int = 8,
      temperature: float = 1.0,
      top_p: float = 1.0,
      schedule: str = 'cosine',
      strategy: str = 'entropy',
      temperature_annealing: bool = False,
      condition_on_coordinates_only: bool = True,
      invalid_ids: Sequence[int] = (),
  ) -> dict[str, Any]:
    """Runs ESM3 iterative generation on one track.

    Args:
      track: 'sequence' | 'structure' | 'secondary_structure' | 'sasa' |
        'function'.
      sequence: Prompt with '_' marking positions to generate.
      coordinates: (L, 37, 3) with NaN marking unconditioned atoms.
      num_steps: Decoding steps. Must be <= L; the API caps this at 100.

    Returns:
      The raw response; `outputs` holds sequence / coordinates / plddt / ptm /
      pae / function / secondary_structure / sasa.
    """
    if track not in (
        'sequence', 'structure', 'secondary_structure', 'sasa', 'function'
    ):
      raise BiohubError(f'Unknown ESM3 track: {track!r}')

    coords = coordinates
    if coords is not None and isinstance(coords, np.ndarray):
      coords = _nan_to_none(coords)

    payload = {
        'model': model,
        'inputs': {
            'sequence': sequence,
            'secondary_structure': secondary_structure,
            'sasa': sasa,
            'coordinates': coords,
        },
        'track': track,
        'invalid_ids': list(invalid_ids),
        'schedule': schedule,
        'strategy': strategy,
        'num_steps': int(num_steps),
        'temperature': float(temperature),
        'temperature_annealing': temperature_annealing,
        'top_p': float(top_p),
        'condition_on_coordinates_only': condition_on_coordinates_only,
        'only_compute_backbone_rmsd': False,
    }
    return self.post('generate', payload)

  def forward_and_sample(
      self,
      tokens: Sequence[int],
      model: str = DEFAULT_ESM3,
      *,
      temperature: float = 1.0,
      top_p: float = 1.0,
      only_sample_masked_tokens: bool = True,
      topk_logprobs: int = 0,
      structure_tokens: Sequence[int] | None = None,
  ) -> dict[str, Any]:
    """Single ESM3 forward pass with sampling. Returns per-track statistics."""
    payload = {
        'model': model,
        'inputs': {
            'sequence': list(tokens),
            'structure': list(structure_tokens) if structure_tokens else None,
            'secondary_structure': None,
            'sasa': None,
            'function': None,
            'residue_annotation': None,
            'coordinates': None,
        },
        'sampling_config': {
            'sequence': {
                'temperature': float(temperature),
                'top_p': float(top_p),
                'only_sample_masked_tokens': only_sample_masked_tokens,
                'invalid_ids': [],
                'topk_logprobs': int(topk_logprobs),
            }
        },
    }
    return self.post('forward_and_sample', payload)


def _nan_to_none(array: np.ndarray) -> Any:
  """Converts a numpy array to nested lists with NaN -> None (JSON null)."""
  out = array.tolist()

  def walk(node):
    if isinstance(node, list):
      return [walk(x) for x in node]
    if isinstance(node, float) and np.isnan(node):
      return None
    return node

  return walk(out)


def _body_of(exc: HttpError) -> str:
  """Extracts the response body text from an HttpError."""
  raw = getattr(exc, 'body', None)
  if raw:
    return raw.decode('utf-8', 'replace') if isinstance(raw, bytes) else str(raw)
  return str(exc)


def _server_message(body: str) -> str:
  """Pulls the server's `message` field out of a JSON error body."""
  try:
    message = json.loads(body).get('message', body)
  except (ValueError, AttributeError):
    return str(body)[:600]
  if isinstance(message, (dict, list)):
    message = json.dumps(message)
  return str(message)[:600]


def _is_quota_error(body: str) -> bool:
  """True when a 429 body indicates credit exhaustion, not rate limiting."""
  lowered = _server_message(body).lower()
  return any(marker in lowered for marker in _QUOTA_MARKERS)


def _explain(exc: HttpError, endpoint: str) -> str:
  """Turns an HttpError into an actionable message including the server body."""
  status = getattr(exc, 'status_code', None)
  message_text = _server_message(_body_of(exc))

  hints = {
      401: 'Your BIOHUB_API_KEY is missing or invalid. Create one at '
           'https://biohub.ai/developer-console/api-keys',
      403: 'Your API key does not have access to this model. Keys are commonly '
           'limited to: esmc-300m/600m/6b-2024-12, esm3-open-2024-03, '
           'esmfold2-fast-2026-05, esmfold2-2026-05.',
      422: 'The request payload was rejected. Check the model name supports '
           'this endpoint and that the field names are correct.',
      429: 'Rate limited. polite-http already backs off; if this persists, '
           'lower BIOHUB_QPS.',
  }
  hint = hints.get(int(status)) if status else None
  message = f'Biohub API {endpoint} failed'
  if status:
    message += f' (HTTP {status})'
  message += f': {message_text}'
  if hint:
    message += f'\nHint: {hint}'
  return message
