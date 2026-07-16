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

"""Eval for the esm3_function_prediction skill.

Exercises the real CLI end to end against the live Biohub API.

Structural: the `predict` schema, 1-indexed inclusive ranges inside the
sequence, InterPro accessions parsed out of the labels, `batch` over a FASTA,
and `plot` writing a real PNG.

Scientific ground truth (measured against the live API; every value below was
reproduced across >= 4 independent runs before being asserted):

  UBIQUITIN -> 'ubiquitin' labels + IPR000626 / IPR019956 / IPR029071
  LYSOZYME  -> lysozyme / glycoside-hydrolase labels + IPR001916 / IPR023346
               (IPR001916 = Glycoside hydrolase family 22, which is exactly the
               family of hen egg-white C-type lysozyme)
  CA2       -> carbonic-anhydrase / zinc / lyase labels + IPR001148 / IPR023561
               (IPR001148 = Alpha carbonic anhydrase domain; human CA2 is an
               alpha-class carbonic anhydrase)
  scramble(UBIQUITIN, seed=0) -> ZERO annotations, reproducibly (6/6 runs at
               temperature 0.7 and 1.0). Both available negative-control
               assertions therefore hold, and we assert both: no
               'ubiquitin' label at all, AND coverage strictly below real
               ubiquitin's.

The keyword sets are checked for SPECIFICITY too: the lysozyme keywords must not
fire on CA2 and the CA2 keywords must not fire on lysozyme. A keyword generic
enough to match any protein would pass the positive check while proving nothing,
so the eval refuses to accept one.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

from __future__ import annotations

import pathlib
import re
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from fixtures import CA2  # noqa: E402
from fixtures import Eval  # noqa: E402
from fixtures import LYSOZYME  # noqa: E402
from fixtures import UBIQUITIN  # noqa: E402
from fixtures import load_json  # noqa: E402
from fixtures import run_script  # noqa: E402
from fixtures import scramble  # noqa: E402
from fixtures import write_fasta  # noqa: E402

SKILL = 'esm3_function_prediction'
SCRIPT = 'predict_function.py'

IPR_RE = re.compile(r'\bIPR\d{6}\b')

# Diagnostic keywords. Deliberately specific to each protein's biology: none of
# them is generic enough to match an arbitrary protein, and the specificity
# check below enforces that.
LYSOZYME_KEYWORDS = (
    'lysozyme', 'glycoside', 'glycosyl', 'hydrolase', 'muramidase',
    'cell wall', 'peptidoglycan',
)
CA2_KEYWORDS = (
    'carbonic anhydrase', 'lyase', 'zinc', 'carbonate dehydratase',
)


def labels_of(record: dict) -> list[str]:
  return [a['label'].lower() for a in record.get('annotations', [])]


def hits(record: dict, keywords) -> list[str]:
  """Labels of `record` matching any of `keywords` (case-insensitive)."""
  return [l for l in labels_of(record) if any(k in l for k in keywords)]


def interpro_ids(record: dict) -> set[str]:
  """InterPro accessions the skill parsed out, cross-checked against labels."""
  return set(record.get('interpro_ids') or [])


def main() -> int:
  ev = Eval(SKILL)
  tmp = pathlib.Path(tempfile.gettempdir()) / f'eval_{SKILL}'
  tmp.mkdir(parents=True, exist_ok=True)

  scrambled = scramble(UBIQUITIN, seed=0)

  ubi_json = tmp / 'ubiquitin.json'
  scram_json = tmp / 'scrambled.json'
  batch_json = tmp / 'batch.json'
  summary_csv = tmp / 'summary.csv'
  png = tmp / 'architecture.png'
  fasta = tmp / 'three.fasta'

  # ---------------------------------------------------------------- predict
  run_script(
      SKILL, SCRIPT, 'predict',
      '--sequence', UBIQUITIN, '--id', 'ubiquitin',
      '--output', str(ubi_json),
  )
  ubi = load_json(ubi_json)
  anns = ubi.get('annotations') or []

  ev.check(
      'predict: non-empty annotation list for ubiquitin',
      len(anns) > 0,
      f'{len(anns)} annotations',
  )
  ev.check(
      'predict: every annotation carries label/start/end',
      all(
          isinstance(a.get('label'), str) and a['label']
          and isinstance(a.get('start'), int)
          and isinstance(a.get('end'), int)
          for a in anns
      ),
      f'{len(anns)} checked',
  )

  n = len(UBIQUITIN)  # 76
  bad = [a for a in anns if not 1 <= a['start'] <= a['end'] <= n]
  ev.check(
      'predict: all ranges satisfy 1 <= start <= end <= 76 (1-indexed, incl.)',
      not bad,
      f'{len(anns)} in range' if not bad else f'out of range: {bad}',
  )
  ev.check(
      'predict: length == end - start + 1 (inclusive range arithmetic)',
      all(a['length'] == a['end'] - a['start'] + 1 for a in anns),
  )
  ev.check(
      'predict: sequence length reported correctly',
      ubi.get('length') == n,
      f'{ubi.get("length")}',
  )

  # InterPro accessions must be parsed OUT of the labels into a structured
  # field, and must agree with what the labels actually contain.
  parsed = interpro_ids(ubi)
  from_labels = {
      m for a in anns for m in IPR_RE.findall(a['label'])
  }
  ev.check(
      'predict: at least one InterPro accession parsed for ubiquitin',
      len(parsed) > 0,
      f'{sorted(parsed)}',
  )
  ev.check(
      'predict: parsed accessions match the accessions present in the labels',
      parsed == from_labels,
      f'parsed={sorted(parsed)} labels={sorted(from_labels)}',
  )
  ev.check(
      'predict: annotations with an accession are typed interpro_entry, '
      'others keyword',
      all(
          (a['kind'] == 'interpro_entry') == (a['interpro_id'] is not None)
          for a in anns
      ),
  )
  ev.check(
      'predict: domain architecture spans are non-overlapping and sorted',
      _spans_disjoint_sorted(ubi.get('domain_architecture') or []),
      f'{[(s["start"], s["end"]) for s in (ubi.get("domain_architecture") or [])]}',
  )
  ev.check(
      'predict: coverage is a fraction in [0, 1]',
      0.0 <= ubi.get('coverage', -1) <= 1.0,
      f'{ubi.get("coverage")}',
  )

  # ------------------------------------------------------------- scientific
  ubi_hits = hits(ubi, ('ubiquitin',))
  ev.check(
      'SCIENCE ubiquitin: a label contains "ubiquitin"',
      bool(ubi_hits),
      f'{sorted(set(ubi_hits))}',
  )
  ev.check(
      'SCIENCE ubiquitin: predicts IPR000626 (Ubiquitin-like domain)',
      'IPR000626' in parsed,
      f'{sorted(parsed)}',
  )

  # ------------------------------------------------------------------ batch
  write_fasta(
      fasta,
      [('ubiquitin', UBIQUITIN), ('lysozyme', LYSOZYME), ('ca2', CA2)],
  )
  run_script(
      SKILL, SCRIPT, 'batch',
      '--fasta', str(fasta),
      '--output', str(batch_json),
      '--summary', str(summary_csv),
      '--max-workers', '3',
  )
  batch = load_json(batch_json)
  records = {r['id']: r for r in batch.get('records', [])}

  ev.check(
      'batch: all 3 FASTA records predicted, none failed',
      len(records) == 3 and batch.get('n_failed') == 0
      and all('error' not in r for r in records.values()),
      f'{sorted(records)} n_failed={batch.get("n_failed")}',
  )
  csv_lines = (
      summary_csv.read_text(encoding='utf-8').strip().splitlines()
      if summary_csv.is_file() else []
  )
  ev.check(
      'batch: summary CSV has a header + one row per record',
      len(csv_lines) == 4 and csv_lines[0].startswith('id,'),
      f'{len(csv_lines)} lines',
  )

  lys = records.get('lysozyme', {})
  ca2 = records.get('ca2', {})

  lys_hits = hits(lys, LYSOZYME_KEYWORDS)
  ev.check(
      'SCIENCE lysozyme: a label is diagnostic of lysozyme biology',
      bool(lys_hits),
      f'{sorted(set(lys_hits))}',
  )
  ev.check(
      'SCIENCE lysozyme: predicts IPR001916 (Glycoside hydrolase family 22 — '
      'the family of hen egg-white C-type lysozyme)',
      'IPR001916' in interpro_ids(lys),
      f'{sorted(interpro_ids(lys))}',
  )

  ca2_hits = hits(ca2, CA2_KEYWORDS)
  ev.check(
      'SCIENCE CA2: a label is diagnostic of carbonic anhydrase biology',
      bool(ca2_hits),
      f'{sorted(set(ca2_hits))}',
  )
  ev.check(
      'SCIENCE CA2: predicts IPR001148 (Alpha carbonic anhydrase domain)',
      'IPR001148' in interpro_ids(ca2),
      f'{sorted(interpro_ids(ca2))}',
  )

  # Specificity: a keyword that matched everything would make the two checks
  # above worthless. Prove they discriminate.
  ev.check(
      'SCIENCE specificity: lysozyme keywords do NOT fire on CA2 or ubiquitin',
      not hits(ca2, LYSOZYME_KEYWORDS) and not hits(ubi, LYSOZYME_KEYWORDS),
      f'ca2={hits(ca2, LYSOZYME_KEYWORDS)} ubi={hits(ubi, LYSOZYME_KEYWORDS)}',
  )
  ev.check(
      'SCIENCE specificity: CA2 keywords do NOT fire on lysozyme or ubiquitin',
      not hits(lys, CA2_KEYWORDS) and not hits(ubi, CA2_KEYWORDS),
      f'lys={hits(lys, CA2_KEYWORDS)} ubi={hits(ubi, CA2_KEYWORDS)}',
  )

  # ------------------------------------------------------- negative control
  run_script(
      SKILL, SCRIPT, 'predict',
      '--sequence', scrambled, '--id', 'scrambled_ubiquitin',
      '--output', str(scram_json),
  )
  scram = load_json(scram_json)
  scram_ubi_hits = hits(scram, ('ubiquitin',))

  # Both available assertions hold (the scramble reproducibly returns zero
  # annotations), so we assert both rather than picking the weaker one.
  ev.check(
      'NEGATIVE CONTROL: scrambled ubiquitin yields NO "ubiquitin" label',
      not scram_ubi_hits,
      f'{scram.get("n_annotations")} annotations, ubiquitin labels='
      f'{scram_ubi_hits}',
  )
  ev.check(
      'NEGATIVE CONTROL: scrambled coverage strictly below real ubiquitin',
      scram.get('coverage', 1.0) < ubi.get('coverage', 0.0),
      f'scrambled {scram.get("coverage"):.3f} < real {ubi.get("coverage"):.3f}',
  )

  # ------------------------------------------------------------------- plot
  run_script(
      SKILL, SCRIPT, 'plot',
      '--input', str(ubi_json), '--output', str(png),
  )
  size = png.stat().st_size if png.is_file() else 0
  ev.check(
      'plot: writes a non-empty PNG',
      png.is_file() and size > 5000
      and png.read_bytes()[:8] == b'\x89PNG\r\n\x1a\n',
      f'{size} bytes',
  )

  print(f'\nArtifacts: {tmp}')
  return ev.finish()


def _spans_disjoint_sorted(spans: list[dict]) -> bool:
  """True if spans are sorted by start and pairwise non-overlapping."""
  previous_end = 0
  for span in spans:
    if span['start'] <= previous_end or span['start'] > span['end']:
      return False
    previous_end = span['end']
  return True


if __name__ == '__main__':
  sys.exit(main())
