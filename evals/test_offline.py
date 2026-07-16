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

"""Offline conformance + unit suite. Spends zero API credits.

Runs on every change. It catches the whole class of defects that does not need
the network — SPEC violations, import errors, CLI typos, broken numerics, bad
PDB geometry — so that scarce API credits are only ever spent on the scientific
assertions that genuinely require the model.

  uv run --no-project evals/test_offline.py
"""

# /// script
# requires-python = ">=3.11"
# dependencies = ["polite-http", "python-dotenv", "numpy"]
# ///

from __future__ import annotations

import hashlib
import pathlib
import re
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(
    0,
    str(pathlib.Path(__file__).resolve().parents[1] / 'skills' / 'esm_common'),
)

import numpy as np  # noqa: E402

import esm_biohub as eb  # noqa: E402
from fixtures import UBIQUITIN, Eval  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO / 'skills'

SKILLS = [
    'esmc_protein_embeddings',
    'esmc_embedding_layer_sweep',
    'esmc_mutation_effect_scoring',
    'esmc_sae_feature_interpretation',
    'esmfold2_structure_prediction',
    'esmfold2_binder_screening',
    'esm3_protein_design',
    'esm3_inverse_folding',
    'esm3_function_prediction',
    'esm3_guided_generation',
    'esm_protein_tracks',
    'esm3_secondary_structure_sasa',
]

# A skill that imports any of these is downloading model weights, which is the
# one thing the whole design forbids.
FORBIDDEN_IMPORTS = re.compile(
    r'^\s*(?:import|from)\s+(torch|transformers|huggingface_hub|esm\b|peft|modal)',
    re.MULTILINE,
)


def main() -> int:
  ev = Eval('offline')

  # ---------------------------------------------------------------- SPEC ---
  print(f'\n--- SPEC conformance (all {{len(SKILLS)}} skills) ---')
  canonical = (SKILLS_DIR / 'esm_common' / 'esm_biohub.py').read_bytes()
  want = hashlib.sha256(canonical).hexdigest()

  for skill in SKILLS:
    d = SKILLS_DIR / skill
    scripts = [p for p in (d / 'scripts').glob('*.py')
               if p.name != 'esm_biohub.py']
    skill_md = d / 'SKILL.md'

    ev.check(f'{skill}: SKILL.md exists', skill_md.is_file())
    ev.check(f'{skill}: has a CLI script', bool(scripts))
    ev.check(f'{skill}: citation.bib exists',
             (d / 'references' / 'citation.bib').is_file())

    vendored = d / 'scripts' / 'esm_biohub.py'
    ev.check(
        f'{skill}: vendored esm_biohub.py in sync',
        vendored.is_file()
        and hashlib.sha256(vendored.read_bytes()).hexdigest() == want,
    )

    for script in scripts:
      src = script.read_text(encoding='utf-8')
      bad = FORBIDDEN_IMPORTS.findall(src)
      ev.check(
          f'{skill}/{script.name}: no model-weight deps',
          not bad,
          f'found {bad}' if bad else 'no torch/transformers/hf/esm/peft/modal',
      )
      ev.check(
          f'{skill}/{script.name}: PEP 723 header',
          '# /// script' in src and '# ///' in src,
      )

    if skill_md.is_file():
      md = skill_md.read_text(encoding='utf-8')
      ev.check(f'{skill}: frontmatter name+description',
               md.startswith('---') and 'name:' in md and 'description:' in md)
      uv_cmds = re.findall(r'uv run(?! --no-project)[^\n]*\.py', md)
      ev.check(
          f'{skill}: every `uv run` uses --no-project',
          not uv_cmds,
          'ok' if not uv_cmds else f'{len(uv_cmds)} missing: {uv_cmds[:2]}',
      )

  # ----------------------------------------------------------------- CLI ---
  print('\n--- CLI contract: --help works for every subcommand ---')
  for skill in SKILLS:
    scripts = [p for p in (SKILLS_DIR / skill / 'scripts').glob('*.py')
               if p.name != 'esm_biohub.py']
    for script in scripts:
      proc = subprocess.run(
          ['uv', 'run', '--no-project', '--quiet', str(script), '--help'],
          capture_output=True, text=True, cwd=REPO, timeout=300,
      )
      ok = proc.returncode == 0
      subs: list[str] = []
      if ok:
        m = re.search(r'\{([a-z0-9,\-_]+)\}', proc.stdout)
        subs = m.group(1).split(',') if m else []
      ev.check(f'{skill}: {script.name} --help', ok,
               f'{len(subs)} subcommands: {",".join(subs)}' if ok
               else proc.stderr[-160:])

      for sub in subs:
        p2 = subprocess.run(
            ['uv', 'run', '--no-project', '--quiet', str(script), sub, '--help'],
            capture_output=True, text=True, cwd=REPO, timeout=300,
        )
        ev.check(f'{skill}: {script.name} {sub} --help', p2.returncode == 0,
                 p2.stderr[-120:] if p2.returncode else '')

  # ------------------------------------------------------------- numerics ---
  print('\n--- shared library: numerics ---')
  x = np.array([1.0, 2.0, 3.0, 4.0])
  ev.check('log_softmax normalised',
           abs(float(np.exp(eb.log_softmax(x)).sum()) - 1.0) < 1e-9)
  ev.check('softmax matches exp(log_softmax)',
           np.allclose(eb.softmax(x), np.exp(eb.log_softmax(x))))
  ev.check('log_softmax is shift-invariant',
           np.allclose(eb.log_softmax(x), eb.log_softmax(x + 1000.0)),
           'no overflow on large logits')
  ev.check('entropy(uniform over 4) == 2 bits',
           abs(float(eb.shannon_entropy_bits(np.full(4, 0.25))) - 2.0) < 1e-9)
  ev.check('entropy(delta) == 0 bits',
           abs(float(eb.shannon_entropy_bits(np.array([1.0, 0, 0, 0])))) < 1e-9)
  ev.check('AA20_IDX maps back through VOCAB',
           [eb.SEQUENCE_VOCAB[i] for i in eb.AA20_IDX] == eb.AA20)
  ev.check('vocab has 64-wide logit space', eb.LOGIT_DIM == 64)

  print('\n--- shared library: geometry ---')
  rng = np.random.RandomState(0)
  a = rng.rand(30, 3) * 10
  theta = 0.9
  rot = np.array([[np.cos(theta), -np.sin(theta), 0],
                  [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
  b = (rot @ a.T).T + np.array([4.0, -2.0, 7.0])
  rmsd, _ = eb.kabsch_rmsd(a, b)
  ev.check('kabsch undoes a rigid transform', rmsd < 1e-9, f'RMSD={rmsd:.2e}')
  ev.check('tm_score of identical traces == 1', abs(eb.tm_score(a, a) - 1.0) < 1e-9)
  far = a + rng.rand(30, 3) * 60
  ev.check('tm_score of unrelated traces is low', eb.tm_score(far, a) < 0.5,
           f'{eb.tm_score(far, a):.3f}')
  try:
    eb.kabsch_rmsd(a, rng.rand(5, 3))
    ev.check('kabsch rejects mismatched shapes', False)
  except eb.BiohubError:
    ev.check('kabsch rejects mismatched shapes', True)

  print('\n--- shared library: PDB geometry ---')
  ref = ('ATOM      1  N   MET A   1      27.340  24.430   2.614  '
         '1.00  9.67           N')
  line = eb._pdb_line(record='ATOM  ', serial=1, atom_name='N', res_name='MET',
                      chain='A', res_seq=1, xyz=[27.340, 24.430, 2.614],
                      bfactor=9.67, element='N')
  ev.check('PDB line matches deposited 1UBQ byte-for-byte', line == ref)

  coords = np.full((len(UBIQUITIN), 37, 3), np.nan)
  coords[:, :5, :] = rng.rand(len(UBIQUITIN), 5, 3) * 30
  plddt = np.linspace(0.4, 0.95, len(UBIQUITIN))
  pdb = eb.atom37_to_pdb(coords, UBIQUITIN, plddt)
  atoms = [ln for ln in pdb.splitlines() if ln.startswith('ATOM')]
  ev.check('atom37_to_pdb writes 5 atoms/residue',
           len(atoms) == 5 * len(UBIQUITIN))
  ev.check('all atoms on chain A', {ln[21] for ln in atoms} == {'A'})
  ev.check('residues numbered 1..L',
           sorted({int(ln[22:26]) for ln in atoms}) == list(range(1, 77)))
  bfacs = [float(ln[60:66]) for ln in atoms]
  ev.check('B-factor carries pLDDT scaled to 0-100',
           39.9 < min(bfacs) < 40.1 and 94.9 < max(bfacs) < 95.1,
           f'{min(bfacs):.1f}..{max(bfacs):.1f}')
  ev.check('NaN atoms are omitted, not written as NaN', 'nan' not in pdb.lower())

  print('\n--- shared library: PDB round-trip ---')
  work = REPO / 'evals' / '_artifacts'
  work.mkdir(parents=True, exist_ok=True)
  pdb_path = work / 'offline_roundtrip.pdb'
  pdb_path.write_text(pdb, encoding='utf-8')
  seq2, coords2, res_ids = eb.parse_pdb_atom37(str(pdb_path), chain='A')
  ev.check('round-trip recovers the sequence', seq2 == UBIQUITIN)
  ev.check('round-trip recovers residue ids', res_ids == list(range(1, 77)))
  ca = eb.ATOM37.index('CA')
  dev = float(np.nanmax(np.abs(coords2[:, ca] - coords[:, ca])))
  ev.check('round-trip preserves CA coords to 3dp', dev <= 0.0011,
           f'max dev {dev:.4f} A')

  print('\n--- shared library: validation + guards ---')
  ev.check('validate_sequence upper-cases and strips',
           eb.validate_sequence(' mqif vktl ') == 'MQIFVKTL')
  for bad in ('MQIF123', 'MQIF*', ''):
    try:
      eb.validate_sequence(bad)
      ev.check(f'validate_sequence rejects {bad!r}', False)
    except eb.BiohubError:
      ev.check(f'validate_sequence rejects {bad!r}', True)
  ev.check('validate_sequence allows mask when asked',
           eb.validate_sequence('MK_TL', allow_mask=True) == 'MK_TL')

  print('\n--- shared library: quota + cassette semantics ---')
  ev.check('BiohubQuotaError is a BiohubError',
           issubclass(eb.BiohubQuotaError, eb.BiohubError))
  ev.check('quota 429 body is recognised',
           eb._is_quota_error(
               '{"message":"You have exceeded your daily credit limit of 100 '
               'credits."}'))
  ev.check('plain rate-limit 429 is NOT treated as quota',
           not eb._is_quota_error('{"message":"Too many requests, slow down"}'))
  ev.check('encode is not a billed endpoint',
           'encode' not in eb.BILLED_ENDPOINTS
           and 'logits' in eb.BILLED_ENDPOINTS
           and 'fold' in eb.BILLED_ENDPOINTS)
  k1 = eb._Cassette.key('logits', {'a': 1, 'b': 2})
  k2 = eb._Cassette.key('logits', {'b': 2, 'a': 1})
  ev.check('cassette key is order-independent', k1 == k2)
  ev.check('cassette key separates endpoints',
           k1 != eb._Cassette.key('fold', {'a': 1, 'b': 2}))

  # ----------------------------------------------------------- documentation ---
  # Every doc link in a SKILL.md must resolve, every figure a report.md embeds
  # must exist, and the interpretation-heavy (Tier A) skills must carry an
  # interpretation guide and at least one worked example folder.
  print('\n--- documentation: links, figures, parity ---')
  tier_a = {
      'esmc_mutation_effect_scoring',
      'esmfold2_structure_prediction',
      'esmfold2_binder_screening',
      'esm3_function_prediction',
      'esmc_sae_feature_interpretation',
      'esm3_protein_design',
  }
  md_link = re.compile(r'\[[^\]]+\]\(([^)]+)\)')
  img_link = re.compile(r'!\[[^\]]*\]\(([^)]+)\)')

  def _strip_code_fences(text: str) -> str:
    out, fenced = [], False
    for line in text.splitlines():
      if line.lstrip().startswith('```'):
        fenced = not fenced
        continue
      if not fenced:
        out.append(line)
    return '\n'.join(out)

  shared_ref = SKILLS_DIR / 'esm_common' / 'references' / 'esm-biohub-api.md'
  ev.check('shared ESM/Biohub API reference exists', shared_ref.is_file())

  for skill in SKILLS:
    d = SKILLS_DIR / skill
    skill_md = d / 'SKILL.md'
    if not skill_md.is_file():
      continue

    # Local doc links in SKILL.md resolve (skip http(s), anchors, and any
    # example links that live inside fenced code blocks).
    dead = []
    for target in md_link.findall(
        _strip_code_fences(skill_md.read_text(encoding='utf-8'))):
      if target.startswith(('http://', 'https://', '#', 'mailto:')):
        continue
      rel = target.split('#', 1)[0]
      if not (d / rel).resolve().exists():
        dead.append(target)
    ev.check(f'{skill}: SKILL.md doc links resolve', not dead,
             'all resolve' if not dead else f'DEAD: {dead}')

    ev.check(f'{skill}: links the shared API reference',
             (d / 'references' / 'esm-biohub-api.md').is_file())

    # Every figure a WORKED EXAMPLE embeds must exist beside it. Only the
    # rendered reports under docs/examples/ are checked — report-templates.md is
    # a scaffold whose placeholder image names (`{plddt_plot}.png`) are
    # deliberate. Fenced code blocks are stripped so an example that quotes the
    # template inline doesn't trip the check either.
    examples_dir = d / 'docs' / 'examples'
    orphan_figs = []
    for md in examples_dir.rglob('*.md') if examples_dir.is_dir() else []:
      body = _strip_code_fences(md.read_text(encoding='utf-8'))
      for img in img_link.findall(body):
        if img.startswith(('http://', 'https://')):
          continue
        if not (md.parent / img.split('#', 1)[0]).exists():
          orphan_figs.append(f'{md.parent.name}/{md.name}->{img}')
    ev.check(f'{skill}: worked-example figures all exist', not orphan_figs,
             'all present' if not orphan_figs else f'MISSING: {orphan_figs}')

    # Every local markdown link inside docs/ must resolve too (guards against a
    # doc that points at a provenance file which was later removed).
    doc_dead = []
    for md in (d / 'docs').rglob('*.md') if (d / 'docs').is_dir() else []:
      body = _strip_code_fences(md.read_text(encoding='utf-8'))
      for target in md_link.findall(body):
        if target.startswith(('http://', 'https://', '#', 'mailto:')):
          continue
        if '{' in target or '}' in target:  # template placeholder
          continue
        if not (md.parent / target.split('#', 1)[0]).exists():
          doc_dead.append(f'{md.parent.name}/{md.name}->{target}')
    ev.check(f'{skill}: docs/ internal links resolve', not doc_dead,
             'all resolve' if not doc_dead else f'DEAD: {doc_dead}')

    doc_dir = d / 'docs'
    if skill in tier_a:
      ev.check(f'{skill}: has an interpretation guide (Tier A)',
               (doc_dir / 'interpretation-guide.md').is_file())
      ev.check(f'{skill}: has a report template (Tier A)',
               (doc_dir / 'report-templates.md').is_file())
      examples = list((doc_dir / 'examples').glob('*/report.md')) \
          if (doc_dir / 'examples').is_dir() else []
      ev.check(f'{skill}: has >=2 worked examples (Tier A)',
               len(examples) >= 2, f'{len(examples)} example(s)')
    elif skill != 'esm_common':
      examples = list((doc_dir / 'examples').glob('*/report.md')) \
          if (doc_dir / 'examples').is_dir() else []
      ev.check(f'{skill}: has a worked example (Tier B)',
               len(examples) >= 1, f'{len(examples)} example(s)')

  return ev.finish()


if __name__ == '__main__':
  sys.exit(main())
