# Quickstart Tutorial - ESM skills

A clear, step-by-step tutorial for downloading, installing, and testing the **ESM protein world model skills** in this repository.

> **Key Architecture Highlights:**
> - **Zero Model Weights Downloaded:** No multi-gigabyte downloads.
> - **No GPU Required:** Core skills run 100% against the hosted Biohub Platform API (`https://biohub.ai`).
> - **Fast Setup:** Pure Python (standard library + `numpy`), executed cleanly via `uv`.

---

## 0. Prerequisites

Before starting, make sure you have:

1. **Python 3.10+**
2. **`uv` package runner** (installs in seconds):
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   # or on macOS via Homebrew:
   brew install uv
   ```
3. **Biohub API Token**:
   Generate a key from the [Biohub Developer Console](https://biohub.ai/developer-console/api-keys).

---

## 1. Download the Repository

Clone the repository to your machine:

```bash
git clone https://github.com/fmenol/google-science-skills.git
cd google-science-skills
```

---

## 2. Configure & Install the Skills

### Step 2A: Set Your API Key
Add your token to `~/.env` so all skills can access it automatically:

```bash
echo "BIOHUB_API_KEY=your_token_here" >> ~/.env
```

*(Alternatively, run `export BIOHUB_API_KEY=your_token_here` in your terminal session).*

### Step 2B: Install Skills into Your AI Assistant

* **For Google Antigravity:**
  ```bash
  mkdir -p ~/.gemini/config/skills/
  cp -r skills/esm* ~/.gemini/config/skills/
  ```

* **For Claude Code:**
  ```bash
  mkdir -p ~/.claude/skills/
  cp -r skills/esm* ~/.claude/skills/
  ```

> [!NOTE]
> **Standalone CLI Usage:**
> You do **not** need an AI assistant to run these skills! Every skill includes standalone PEP 723 CLI scripts that you can execute directly from your shell using `uv run --no-project`.

---

## 3. Test the Top 3 Skills by Usefulness

> [!IMPORTANT]
> **Zero-Credit Offline Testing:**
> Biohub accounts have a daily credit allowance (100 credits/day). The repository includes pre-recorded test cassettes (`evals/cassettes`) so you can test all skills **instantly with 0 credits and no network wait**:
> ```bash
> export BIOHUB_CASSETTE=evals/cassettes
> export BIOHUB_CASSETTE_MODE=replay
> ```

---

### #1. `esmfold2-structure-prediction` — 3D Structure Prediction

* **Why it's #1:** Generates high-accuracy 3D atomic coordinates (`.pdb`) for novel proteins and complexes from raw sequence in seconds, outputting per-residue pLDDT, overall pTM, and PAE.
* **Test Command** (Fold Ubiquitin, 76 amino acids):
  ```bash
  uv run --no-project skills/esmfold2_structure_prediction/scripts/fold.py fold \
    --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
    --output-metrics ubiquitin_metrics.json \
    --output-pdb ubiquitin.pdb
  ```
* **Output Files Generated:**
  - `ubiquitin.pdb`: 3D structure ready for PyMOL, ChimeraX, or Mol* (B-factor column holds pLDDT scaled 0–100).
  - `ubiquitin_metrics.json`: Confidence metrics (`pLDDT ~0.82`, `pTM ~0.78`).

---

### #2. `esmc-mutation-effect-scoring` — Zero-Shot Variant Effect Prediction

* **Why it's #2:** Evaluates novel mutations without requiring homologous MSAs or experimental labels. Perfect for triage, deep mutational scanning (DMS), and pathogenicity scoring.
* **Test Command** (Score the functional impact of the `K48R` point mutation in Ubiquitin):
  ```bash
  uv run --no-project skills/esmc_mutation_effect_scoring/scripts/mutation_scoring.py score-variants \
    --sequence MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG \
    --variants K48R \
    --output ubq_k48r.json
  ```
* **Output Summary:**
  - `variant: K48R | LLR: -5.293 | rank@pos: 1/19 | prediction: deleterious`
  - Negative Log-Likelihood Ratio (LLR) indicates that K48 is critical for ubiquitin function and substitution is deleterious.

---

### #3. `esm3-protein-design` — Multimodal De Novo Protein Design

* **Why it's #3:** ESM3's flagship generative capability. Generates brand-new protein sequences de novo, scaffolds functional motifs/active sites, or designs sequences conditioned on secondary structure (SS8) or SASA profiles. Automatically folds each design with ESMFold2 to verify foldability before reporting.
* **Test Command** (Generate 2 de novo 50-residue designs):
  ```bash
  uv run --no-project skills/esm3_protein_design/scripts/design.py generate \
    --length 50 \
    --num-samples 2 \
    --output-prefix design_ubq
  ```
* **Output Files Generated:**
  - `design_ubq.fasta`: Synthesizable amino acid sequences (canonical 20 AAs).
  - `design_ubq.json`: Designs ranked by pTM quality score.
  - `design_ubq_sample_0.pdb` & `design_ubq_sample_1.pdb`: Folded 3D atomic structures for each candidate.

---

## 4. Summary Table: Top 3 Skills

| Rank | Skill Name | CLI Script | Primary Function | Output Format |
| :--- | :--- | :--- | :--- | :--- |
| **1** | [`esmfold2-structure-prediction`](file:///Users/fmenol/Downloads/science-skills/skills/esmfold2_structure_prediction/SKILL.md) | `fold.py` | 3D atomic structure folding (monomers & complexes) | `.pdb`, `.json` |
| **2** | [`esmc-mutation-effect-scoring`](file:///Users/fmenol/Downloads/science-skills/skills/esmc_mutation_effect_scoring/SKILL.md) | `mutation_scoring.py` | Zero-shot variant effect & DMS scoring | `.json`, `.png` |
| **3** | [`esm3-protein-design`](file:///Users/fmenol/Downloads/science-skills/skills/esm3_protein_design/SKILL.md) | `design.py` | De novo sequence generation & motif scaffolding | `.fasta`, `.pdb`, `.json` |

---

## 5. How to Prompt Your AI Agent

Once installed into Antigravity or Claude Code, prompt naturally in chat:
- *"Fold the sequence `MQIFVK...` with ESMFold2 and save the PDB file."*
- *"Evaluate the fitness effect of mutation K48R on ubiquitin."*
- *"Design a 50-residue helical protein scaffold using ESM3."*
- *"Screen this list of candidate binders against target protein X."*
