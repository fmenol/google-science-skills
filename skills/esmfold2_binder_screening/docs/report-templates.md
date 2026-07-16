# Report Templates

Scaffold for reporting a binder screen. Fill every section with **real numbers
from the script output** — never invent an iPTM, a contact count, or a residue.
Read [interpretation-guide.md](interpretation-guide.md) first and work its
pre-report checklist before writing.

> [!IMPORTANT] **DO NOT** use the agent's default artifact directory. Save the
> report as `report.md` directly in the screen's output directory (e.g.
> `results/report.md`) in the workspace. Embed figures with relative paths
> (`filename.png`) and confirm each figure exists in that folder.

Two templates: a **screen report** (rank a panel of candidates) and a
**single-interface report** (deep-dive one target:binder pair). Most jobs use
the first.

--------------------------------------------------------------------------------

## 1. Screen Report Template

Use this to report `screen` (and any `rank`/filter) over a candidate panel.

```markdown
# Binder Screen: {n} candidates vs {target_name}

## 1. Summary & Shortlist
[One paragraph. State the target, how many candidates were screened, and the
verdict: which candidate(s) cleared the confident-interface band (iPTM > 0.8)
and are shortlisted to order, and which did not. If the shortlist is empty, say
so plainly — "ESMFold2 predicts no confident interface for any candidate" is a
complete result. State that iPTM is a structural-confidence prediction to test
experimentally, NOT an affinity. Note that this ESMFold2 binder-screening skill
was used.]

**Shortlist ({top_n}):** {candidate(s), or "none — no candidate cleared iPTM 0.8"}

## 2. Screen Setup
- **Target:** {target_name} ({L} aa) — chain A
- **Candidates:** {n} ({provenance: designs, homologs, mutant panel, ...})
- **Model / params:** {model}, num_loops {k}, num_sampling_steps {s}
- **Ranked by:** iptm   **Filters:** {--min-iptm, --max-isoelectric-point, ...}

## 3. Ranked Table
[Sorted by iPTM, best first. Quote real numbers from results.csv. Mark the
shortlist. Read confident_contacts, NOT the raw geometric count.]

| Rank | Candidate | iPTM | Verdict | Interface PAE (Å) | Binder pLDDT | Confident contacts | Selection score | Shortlisted |
| :--- | :-------- | :--- | :------ | :---------------- | :----------- | :----------------- | :-------------- | :---------- |
| 1 | {name} | {iptm} | {band} | {pae} | {plddt} | {n} | {score} | {✓/—} |
| … |

*Ranked on iPTM. Contacts are the PAE-gated `confident_interface_contacts`; the
raw geometric count is not a binding signal.*

## 4. Per-Candidate Verdicts
[For each candidate that matters — every shortlisted hit, plus any instructive
near-miss or decoy — one short block. Tie the verdict to iPTM AND interface PAE,
and say what to do.]

### {candidate} — {confident interface | possible | likely no binding}
- **iPTM {x}** ({band}); **interface PAE {y} Å**; **{n} confident contacts**.
- Binder folds: pLDDT {z} ({"yes" if >0.7}). [Folding ≠ binding — note if a
  confident fold does NOT bind.]
- [If a hit: name the interface residues from the per-residue output. If not:
  state the non-binder signature — high interface PAE, ~0 confident contacts.]
- **Action:** {order | deprioritise | drop}.

## 5. Interface Evidence (PAE heatmaps)
[Embed the heatmap for each candidate discussed. Describe its signature.]

![{target}:{candidate} PAE]({target}_{candidate}_pae.png)
*Fig N: {"Whole matrix dark — confident interface." | "Dark diagonals, pale
off-diagonal cross-chain blocks — folds but does not bind." | "Candidate's own
block also pale — neither folds nor binds."}*

## 6. Limitations
[Restate: iPTM is a confidence, not a Kd — no affinity, kinetics, or specificity.
No expression/solubility guarantee (note pI if the pI<6 minibinder filter is
relevant). ESMFold2 folds one static complex; dynamic/context-dependent binding
is out of scope. This skill ranks candidates; it does not design them.]

## 7. Conclusion & Next Steps
[Which to order and why, in one narrative. Every verdict is a hypothesis for the
wet lab. State the experimental test that would confirm the top hit.]
```

--------------------------------------------------------------------------------

## 2. Single-Interface Report Template

Use this for `analyze-interface` on one target:binder pair (a confirmed hit, or
a diagnosis of why a candidate fails).

```markdown
# Interface Analysis: {target} : {binder}

## 1. Verdict
[One line: iPTM {x} → {band}. Interface PAE {y} Å. {n} confident contacts.
State plainly whether ESMFold2 predicts a confident interface — a hypothesis to
test, not an affinity.]

## 2. Metrics
| Metric | Value | Reading |
| :----- | :---- | :------ |
| iPTM | {x} | {band} |
| interface PAE | {y} Å | {<5 confident / >20 none} |
| binder pLDDT | {z} | {folds?} |
| confident contacts | {n} | {gated} |
| geometric contacts | {m} | (not a binding signal alone) |
| selection_score | {s} | composite |
| isoelectric point | {pI} | {solubility proxy; pI<6 for minibinders} |

## 3. Interface Residues
[From the per-residue output. Target hotspots and binder hotspots by contact
count. Only claim residues the script lists.]
- **{target}:** {Res#(contacts), ...}
- **{binder}:** {Res#(contacts), ...}

## 4. PAE Heatmap
![PAE]({target}_{binder}_pae.png)
*Fig 1: {signature description}.*

## 5. Interpretation & Limitations
[What the interface means and does not. iPTM ≠ Kd. No specificity/affinity.
Next experimental step.]
```
