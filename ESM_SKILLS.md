# ESM Skills — augmenting the Google DeepMind Science Skills

Twelve new skills that bring the **ESM protein world model** (ESMC, ESMFold2,
ESM3, and the ESMC sparse autoencoders) into the Science Skills bundle.

They are derived from the tutorials in [`esm/cookbook/tutorials/`](esm/cookbook/tutorials/)
and run entirely against the **Biohub Platform API** (`https://biohub.ai`).
**No model weights are ever downloaded.**

---

## 1. Why these skills, and how they fit

The existing bundle is overwhelmingly **retrieval**: 30-odd skills that look
things up in curated databases (UniProt, PDB, ChEMBL, ClinVar, AlphaFold DB…).
Its one *predictive* skill is `alphagenome-single-variant-analysis`, which runs a
model over the **non-coding genome**.

The ESM skills add the missing half: **prediction, design and interpretation over
proteins**. Concretely, they let an agent answer questions no database can:

| Question | Nothing in the bundle could answer it | Now |
|---|---|---|
| Is this *novel* mutation deleterious? | ClinVar only knows *observed* variants | `esmc-mutation-effect-scoring` |
| What does this structure look like? | AlphaFold DB only has *precomputed* UniProt entries | `esmfold2-structure-prediction` |
| Does my *designed* binder actually bind? | — | `esmfold2-binder-screening` |
| Design me a protein with this active site | — | `esm3-protein-design` |
| What does this dark-matter protein *do*? | InterPro needs a homology hit | `esm3-function-prediction` |

They also compose with the existing skills rather than duplicating them:
`esm3-function-prediction` hands a predicted `IPRxxxxxx` accession to
`interpro_database`; `esmfold2-structure-prediction` produces the coordinate file
that `foldseek-structural-search` and `pymol` consume.

---

## 2. The twelve skills

Each maps to a cookbook tutorial (or the reachable half of one).

| # | Skill | Source tutorial | What it does |
|---|---|---|---|
| 1 | `esmc-protein-embeddings` | `embed.ipynb` | Embed sequences with ESMC (300M/600M/6B); per-residue or mean-pooled, any layer. Similarity, clustering, PCA. |
| 2 | `esmc-embedding-layer-sweep` | `esmc_layer_sweep.ipynb` | Pick the *optimal layer* for a downstream property via a linear probe + stratified CV. Intermediate layers routinely beat the last one. |
| 3 | `esmc-mutation-effect-scoring` | `esmc_mutation_scoring.ipynb` | Zero-shot variant effects: full L×20 log-likelihood-ratio matrix, per-position entropy, DMS heatmap, pseudo-perplexity. No MSA, no structure, no labels. |
| 4 | `esmc-sae-feature-interpretation` | `esmc_sae_feature_interpretation.ipynb` | Mechanistic interpretability: decompose ESMC-6B's internals into ~16 384 sparse features and fetch each one's natural-language description. |
| 5 | `esmfold2-structure-prediction` | `esmfold2.ipynb` | All-atom structure prediction — single chains **and** complexes with DNA, RNA and small-molecule ligands (CCD or SMILES). |
| 6 | `esmfold2-binder-screening` | `binder_design.ipynb` (reachable half) | Fold target:candidate complexes and rank binders by interface confidence (iPTM, interface PAE, contacts). |
| 7 | `esm3-protein-design` | `esm3_generate.ipynb`, `gfp_design.ipynb` | Promptable multimodal design: de novo, motif scaffolding, secondary-structure- and SASA-conditioned generation. |
| 8 | `esm3-inverse-folding` | `snippets/fold_invfold.py` | Structure → sequence, validated by self-consistency (design → refold → scTM/scRMSD). |
| 9 | `esm3-function-prediction` | `snippets/esm3.py` (function track) | Predict InterPro-style domain annotations *from sequence alone* — works where BLAST/InterPro find nothing. |
| 10 | `esm3-guided-generation` | `esm3_guided_generation.ipynb` | Steer generation toward an objective (pTM, no-cysteine, hydrophobicity, radius of gyration) with soft value-based decoding. |
| 11 | `esm-protein-tracks` | `esmprotein.ipynb` | Prepare and mask the five ESM3 prompt tracks (sequence, coordinates, SS8, SASA, function) from a PDB. The input stage for #7/#8/#10. |
| 12 | `esm3-secondary-structure-sasa` | API capability | Predict per-residue secondary structure (true DSSP SS8) and solvent accessibility from sequence alone — no folding. |

Skill #12 is not from a tutorial; it was added after an audit of the full API
surface (see [§ API coverage](#api-coverage) below) found that ESM3's
`secondary_structure` and `sasa` generation tracks expose a distinct, reachable
capability that no other skill covered. It fills a gap the `esm-protein-tracks`
skill itself flagged: that skill can only *measure* secondary structure from a
known structure, and only as a 3-state (H/E/C) biotite approximation. This one
*predicts* genuine 8-state SS8 from sequence.

### Deliberately NOT built

| Tutorial | Why not |
|---|---|
| `esmc_finetune.ipynb` | LoRA fine-tuning needs gradients through the backbone: local weights + a CUDA GPU. **There is no remote-gradient endpoint.** Building it would violate the "no model downloads" constraint and could not be made to work. |
| `binder_design.py` (the *optimisation* half) | Backpropagates structural losses through the ESMFold2 trunk; needs an H100 with 27–51 GB VRAM and Modal. Its **scoring/selection half is reachable**, and is skill #6. |

Rather than ship two skills that could never pass an eval, these are documented
as out of scope. Skill #6's `SKILL.md` says so explicitly, and points users at
#7/#8 to *generate* candidates before ranking them.

---

## API coverage — is anything reachable left unbuilt?

The twelve skills were checked against the **complete** set of endpoints and
model capabilities the API key can reach, each verified with a live call:

| API capability | Endpoint | Covered by |
|---|---|---|
| ESMC embeddings / hidden states (all layers) | `logits` | #1, #2 |
| ESMC sequence logits (variant scoring) | `logits` | #3 |
| ESMC SAE features + feature descriptions | `logits`, `features/{id}` | #4 |
| ESMFold2 single-chain fold (coords, pLDDT, pTM, PAE) | `fold` | #5 |
| ESMFold2 complexes (protein/DNA/RNA/ligand, iPTM) | `fold_all_atom` | #5, #6 |
| ESM3 generation — sequence / structure tracks | `generate` | #7, #8 |
| ESM3 generation — function track | `generate` | #9 |
| ESM3 generation — **secondary_structure / sasa tracks** | `generate` | **#12** |
| ESM3 objective-guided decoding | `generate` + `fold` | #10 |
| Track extraction / prompt building | (local) | #11 |

Three further capabilities are technically reachable but **deliberately not
built**, because each is redundant with an existing skill and would add surface
without adding science:

- **`forward_and_sample` single-pass token statistics** (per-position entropy +
  top-k alternatives for ESM3). Verified working, but the useful signal — which
  residues are constrained, what the alternatives are — is exactly what
  `esmc-mutation-effect-scoring` already delivers, with a cleaner leave-one-out
  formulation.
- **ESMFold2 pair-pooled embeddings** (`include_embeddings`). Verified to return
  a `(L, 256)` structure-aware embedding, but `esmc-protein-embeddings` already
  covers the embedding use cases, and this niche variant has no distinct
  downstream task in the tutorials.
- **Function-annotation-conditioned design** (`generate` with
  `function_annotations` as input). The open ESM3 model's function conditioning
  is weak, so a dedicated skill would overpromise; `esm3-protein-design` already
  conditions on the sequence/structure/SS/SASA tracks that work well.

So: every reachable capability is either built or consciously declined with a
reason. If the API later exposes the gated ESM3 models (`medium`, `large`,
`multimer` — currently HTTP 403) or a gradient endpoint, that would open genuinely
new skills (multimer design, the binder-design optimisation loop); nothing more
is reachable on the current key.

## Documentation

Every skill is documented to the standard set by
`alphagenome_single_variant_analysis` — the bundle's other model-based skill. The
depth scales with how much interpretation the skill demands
([DOC_PARITY_PLAN.md](DOC_PARITY_PLAN.md) records the tiering):

- **Shared:** [`skills/esm_common/references/esm-biohub-api.md`](skills/esm_common/references/esm-biohub-api.md)
  — the endpoints, reachable models, response shapes, credit model and verified
  gotchas that every skill links to.
- **Tier A** (the six interpretation-heavy skills — mutation scoring, structure
  prediction, binder screening, function prediction, SAE interpretation, design):
  each ships a `docs/interpretation-guide.md` (thresholds, a "Negative Results &
  Scientific Integrity" section, a pre-report checklist), a
  `docs/report-templates.md`, and `docs/examples/` with **worked positive and
  negative cases** — real `report.md` analyses with real figures.
- **Tier B** (the six capability skills): each ships an interpretation guide and
  one worked example.

Every worked example was generated by **replaying the evals from cassettes —
zero credits** — so every figure and every number in the docs is real model
output on a real protein, including the negative controls (a scrambled sequence
folding worse, a decoy binder rejected, ubiquitin scoring 0/10 on lysozyme's SAE
features). The offline suite checks that every doc link resolves, every
worked-example figure exists, and every Tier-A skill carries the full set — so
the documentation cannot silently rot.

## 3. Architecture

### No weights, no torch

The `esm` PyPI package pulls `torch` **and a fork of `transformers`**
(`git+https://github.com/Biohub/transformers`), and its local paths call
`huggingface_hub.snapshot_download`. That is incompatible with the constraint.

Probing the live API showed it is unnecessary: `POST https://biohub.ai/api/v1/*`
with a bearer token returns **plain JSON** — embeddings, logits, coordinates and
even SAE features come back as nested lists. (Only the opt-in `return-bytes`
path uses zstd + pickled torch tensors; we simply never ask for it. SAE features
arrive as a sparse top-k dict: `{feature_indices, values, shape}`.)

So every skill is **stdlib + numpy**, installs in seconds, and never touches a
GPU.

### One shared client, vendored

[`skills/esm_common/esm_biohub.py`](skills/esm_common/esm_biohub.py) is the single
source of truth: auth, cross-process rate limiting, retry/backoff, the model
registry, numerics, PDB writing/parsing, Kabsch/TM-score.

Skills must be independently installable, so it is **copied** into each skill's
`scripts/` rather than imported across directories.
[`evals/sync_common.py`](evals/sync_common.py) performs the copy and
`--check` fails the eval suite if any copy has drifted.

Rate limiting uses [`polite-http`](https://pypi.org/project/polite-http/) — the
same library the existing GDM skills use — whose limiter is a **cross-process
file lock**, so concurrent sub-agents on one machine collectively respect the
limit.

### `uv run --no-project` is mandatory

Without `--no-project`, `uv` walks up the tree, finds whatever `pyproject.toml`
the user happens to be nested under, and tries to build *their* project. This
actually broke here (an unrelated project two directories up). Every documented
command uses `--no-project`.

---

## 4. What the API key can actually reach

Probed directly, not assumed:

| Model | Status |
|---|---|
| `esmc-300m-2024-12`, `esmc-600m-2024-12`, `esmc-6b-2024-12` | available |
| `esm3-open-2024-03` (= `esm3-sm-open-v1`) | available |
| `esmfold2-fast-2026-05`, `esmfold2-2026-05` | available |
| `esm3-medium-*`, `esm3-large-*`, `esm3-*-multimer` | **HTTP 403 — no access** |

Three tutorials (`gfp_design`, `esm3_guided_generation`, `fold_invfold`) are
written against `esm3-medium`/`esm3-large`. Every skill re-points at
`esm3-open-2024-03`.

Two further traps found by probing, both of which would silently waste a user's
time:

* **`include_distogram=True` → HTTP 422, not implemented server-side.** The
  distogram-based iPTM proxy from the binder-design paper is unreachable. The
  API returns a real `interface_ptm` instead, which is better.
* **The `/inverse_fold` endpoint is supported by *no* reachable model.** Inverse
  folding must go through `generate(track='sequence', coordinates=…)` — which
  works, and recovers ubiquitin's native sequence at 100 % identity.

---

## 5. The credit cap — read this before running anything

**The Biohub account has a hard allowance of 100 credits per day, resetting at
00:00 UTC.** This is not a rate limit you can back off around; it is a wall. It
shapes how these skills must be used, and it is the single most important
operational fact here.

Three findings, all verified against the live service:

* **HTTP 429 means two different things.** The service returns it both for
  transient rate limiting *and* for daily-credit exhaustion; only the response
  body distinguishes them. `polite-http` retries 429 by default, so an exhausted
  quota used to send every single call through six exponential backoffs — about
  **8 minutes of waiting per request, guaranteed to fail**. The client now
  inspects the body and raises `BiohubQuotaError` immediately (measured: 1.1 s).
* **`encode` is free.** It still returns 200 with the quota at zero. Only
  `logits`, `fold`, `fold_all_atom`, `generate` and `forward_and_sample` are
  billed. (So `encode` is useless as a credit probe — a 200 from it means
  nothing.)
* **There is no usage endpoint.** Remaining credits cannot be queried; you only
  discover exhaustion by hitting it.

### What this costs in practice

Billed calls are what matter. A few are cheap; one is not:

| Operation | Billed calls |
|---|---|
| Embed one sequence | 1 |
| Fold one protein / one complex | 1 |
| Predict function | 1 |
| Screen N binder candidates | N |
| **Deep mutational scan of an L-residue protein** | **L** (one masked pass per residue) |

A full DMS scan of a 300-residue protein is ~300 billed calls and **cannot
complete within a single day's allowance**. Scan a domain, not a proteome.

### Cassettes: pay once, re-run free

Because of the cap, an eval suite that calls the live API on every run is not
repeatable — which would make "every skill is verified" a claim you could only
check once. The client therefore supports record/replay (the VCR/betamax
pattern):

```bash
export BIOHUB_CASSETTE=evals/cassettes
export BIOHUB_CASSETTE_MODE=auto    # replay if recorded, else call + record
```

On replay the whole pipeline still runs end to end — payload construction,
parsing, numerics, interpretation, file output — and only the network hop is
served from disk. Replay needs **no API key and no credits**.

Recording is coordinated by [`evals/record.py`](evals/record.py), which:

* runs the evals **serially** (running them in parallel makes them race for one
  credit pool, so a burst exhausts it and skills fail for reasons unrelated to
  the skill under test — this is exactly what happened during development);
* orders them **cheapest-first**, so a truncated pass still banks the most green
  skills;
* is **resumable** — responses are written as they arrive, so if the quota dies
  mid-pass the credits already spent are kept, and the next run replays those
  and spends only on what is still missing. The pass converges even when the full
  suite costs more than one day's allowance.

```bash
uv run --no-project evals/record.py --wait     # block until reset, then record
uv run --no-project evals/record.py --replay   # verify from cassettes (free)
```

Cassettes also make one eval much cheaper than it looks: `scan`, `entropy`,
`score-variants` and `pseudo-perplexity` all re-issue the *same* leave-one-out
payloads, so only the first is billed and the rest are cache hits.

---

## 6. Evals: how "working" is defined

Every skill has `evals/eval_<skill>.py`, which drives its **real CLI** end-to-end
against the **live API**. Structural checks (files written, shapes, schemas) are
necessary but not sufficient — each skill must also pass a **scientific** check
with a real positive and, wherever possible, a **negative control**.

The controls are the point. A skill that merely returns *a* number always
"passes"; a skill that must rank a true binder above a decoy cannot fake it.

| Skill | Scientific assertion |
|---|---|
| embeddings | A point mutant of lysozyme embeds closer to lysozyme than ubiquitin does; KMeans recovers two protein families (ARI > 0.8) |
| layer sweep | Best layer reaches MCC > 0.8, and beats the raw embedding layer; per-layer scores are not constant |
| mutation scoring | LLR(wild-type) ≡ 0; **pseudo-perplexity(real) ≪ pseudo-perplexity(shuffled)**; X→Pro is scored far worse than conservative aliphatic swaps |
| SAE | Top-feature sets: `jaccard(lysozyme, lysozyme-mutant) > jaccard(lysozyme, ubiquitin)`; descriptions are non-empty and diagnostic |
| structure prediction | Ubiquitin folds to pLDDT > 0.7 / pTM > 0.6; **its own shuffled sequence folds strictly worse on both** |
| binder screening | **Barstar ranks #1 against barnase (iPTM > 0.8), above both a shuffled-barstar decoy and unrelated lysozyme** |
| protein design | Designed sequences fold with higher pTM than random amino-acid strings; grafted motifs are preserved (CA RMSD < 3 Å) |
| inverse folding | Native-sequence recovery ≫ chance; best design re-folds to the target backbone with scTM > 0.7 |
| function prediction | Ubiquitin → "ubiquitin"; lysozyme → glycoside-hydrolase terms; **shuffled ubiquitin does not** |
| guided generation | The `no-cysteine` objective yields **exactly zero cysteines** — a hard, deterministic outcome |
| protein tracks | Sequence extracted from PDB 1UBQ **equals the reference ubiquitin sequence exactly**; SS8 contains both helix and strand |

### Two tiers, because credits are scarce

**Tier 1 — offline (free, run on every change):**

```bash
uv run --no-project evals/test_offline.py     # 171 checks, 0 credits, ~10 s
```

Enforces the SPEC across all 12 skills (no weight-downloading imports; PEP 723
headers; every `uv run` carries `--no-project`; the vendored client is in sync;
every subcommand's `--help` parses) and unit-tests the shared library — numerics,
Kabsch/TM-score, PDB column geometry, the PDB round-trip, quota detection and
cassette semantics. This catches the entire class of defect that does not need
the network, so credits are only ever spent on assertions that genuinely require
the model.

**Tier 2 — live science (credit-gated, recorded once):**

```bash
uv run --no-project evals/record.py --wait     # serial, resumable, cheapest-first
uv run --no-project evals/record.py --replay   # re-verify from cassettes, free
```

### Verification status — all green

| Skill | Live eval | Headline scientific result |
|---|---|---|
| `esm_protein_tracks` | **34/34** | sequence from 1UBQ == reference ubiquitin exactly; SS8 has both H and E |
| `esm3_function_prediction` | **23/23** | ubiquitin→"ubiquitin", lysozyme→hydrolase; scramble does not |
| `esmc_protein_embeddings` | **14/14** | point mutant embeds at 0.984 vs 0.507 for an unrelated protein; clustering ARI 1.0 |
| `esmc_embedding_layer_sweep` | **26/26** | best layer 23/30 (77 % depth) beats the last layer; embedding layer never crowned |
| `esmc_sae_feature_interpretation` | **26/26** | 9/10 of lysozyme's top features name its real biology; ubiquitin 0/10 |
| `esm3_guided_generation` | **16/16** | no-cysteine objective → **0 cysteines**; guidance improves GRAVY (+0.38 vs −0.52) |
| `esm3_inverse_folding` | **18/18** | native recovery 100 %; best design re-folds at scTM 1.0; shuffled-correspondence decoy rejected |
| `esmfold2_structure_prediction` | **25/25** | ubiquitin pLDDT 0.82 / pTM 0.78; its scramble folds strictly worse |
| `esmfold2_binder_screening` | **21/21** | barstar ranks #1 vs barnase (iPTM 0.96), 52 interface contacts, above decoy + unrelated |
| `esm3_protein_design` | **19/19** | designed pTM 0.32 vs random 0.19; motif preserved at CA RMSD 0.23 Å |
| `esmc_mutation_effect_scoring` | **21/21** | pseudo-perplexity 1.05 (ubiquitin) vs 18.70 (its scramble); LLR(WT)≡0 |
| `esm3_secondary_structure_sasa` | **21/21** | ubiquitin helix placed at res 23-34; myoglobin 76% helix / 0% strand; hydrophobics buried (36 vs 101 Å²) |
| **Offline conformance + units** | **183/183** | zero credits — SPEC + library, runs on every change |

**12/12 live, 264 scientific checks, recorded across coherent passes that each fit
inside a single day's 100-credit allowance.** The full suite then replayed from
cassettes twice, deterministically, with no API key and no credits.

Two eval-robustness bugs were found and fixed during the replay verification, not
papered over:
* **Guided generation wasn't replay-deterministic.** Its N per-step candidates
  share one prompt (one cassette key), and a thread race plus cross-invocation
  key collisions made a different candidate win on replay, diverging the
  trajectory. Fixed by pairing the i.i.d. completions to slices in sorted order
  (deterministic winner) and namespacing each CLI invocation's cassette.
* **The designed-vs-random pTM test was flaky at n=2.** Server sampling is not
  seedable and random 60-mers occasionally fold as well as an unlucky pair of
  designs. A clean n=8 measurement (design mean 0.315 vs random 0.192 on *every*
  statistic — mean, median, best, pLDDT) set an honest n=5 threshold.

Reference values, measured against the live API: ubiquitin folds to pLDDT 0.82 /
pTM 0.78; ESMC recovers 75/76 of its wild-type residues; inverse folding recovers
the native sequence at 100 % identity; barnase + barstar gives iPTM 0.96.

---

## 6. Getting started

```bash
# 1. API token -> https://biohub.ai/developer-console/api-keys
printf "Enter BIOHUB_API_KEY (typing hidden): " && read -s val && echo && \
  echo "BIOHUB_API_KEY=$val" >> ~/.env && echo "Saved."

# 2. Fold a protein
uv run --no-project skills/esmfold2_structure_prediction/scripts/fold.py fold \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --output ubiquitin.json --pdb-output ubiquitin.pdb

# 3. Score every possible point mutation
uv run --no-project skills/esmc_mutation_effect_scoring/scripts/mutation_scoring.py scan \
  --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
  --output ubq_scan.json
```

Set `BIOHUB_QPS` to tune the shared rate limit (default 5).

---

## 7. Licensing

ESM models and code are MIT-licensed by Biohub. Use of the hosted Biohub
Platform is governed by the
[Acceptable Use Policy](https://biohub.org/acceptable-use-policy/). The platform
applies guardrails on sequences corresponding to controlled pathogens and toxins.
Each skill notifies the user of these terms on first use, following the same
`.licenses/` convention as the rest of the bundle.
