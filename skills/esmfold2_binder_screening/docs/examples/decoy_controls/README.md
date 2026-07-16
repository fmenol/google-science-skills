# Example: What "No Binding" Looks Like (decoy controls)

**Target:** barnase (110 aa) **Candidates:** barstar (true binder) + two decoys
**Verdict:** **Only barstar binds; both decoys land below iPTM 0.5 and are
filtered out**

This is the example to study, because telling a non-binder from a binder is the
entire job. The screen puts one true binder against two decoys chosen to break
the two metrics people wrongly trust:

- **lysozyme** — a well-folded, completely unrelated protein. Breaks the trap of
  reading a high **binder pLDDT** (it folds beautifully) as binding.
- **scrambled barstar** — barstar's amino-acid composition, shuffled so the
  fold is destroyed. Breaks the trap of reading **raw contact count** as binding
  (same residues, so it packs plenty of incidental contacts).

## How it was produced

Real output from this skill in replay mode (no credits), at the eval settings:

```bash
# rank the panel
uv run --no-project scripts/screen_binders.py screen \
  --target-fasta barnase.fasta --candidates candidates.fasta \
  --output screen/ --top-n 1 --num-loops 3 --num-sampling-steps 50

# per-decoy PAE heatmaps (same folds, deep-dived)
uv run --no-project scripts/screen_binders.py analyze-interface \
  --target-fasta barnase.fasta --binder-fasta lysozyme.fasta --output ... \
  --num-loops 3 --num-sampling-steps 50    # and again for the scrambled decoy

# the min-iPTM filter that drops both decoys
uv run --no-project scripts/screen_binders.py rank \
  --results screen/results.json --output shortlist.json \
  --min-iptm 0.8 --top-n 1
```

## Why this is a good example

The measured numbers separate cleanly, and each decoy teaches one trap:

1.  **iPTM ranks correctly.** barstar **0.964** » lysozyme **0.230** » scrambled
    **0.144**. The true binder is ~0.73 iPTM above the better decoy; both decoys
    sit firmly in the "likely no binding" band (< 0.5).
2.  **Folding is not binding (lysozyme).** Its binder pLDDT is **0.897** —
    better than many real designs — yet iPTM is 0.230 with **0** confident
    contacts. A confident fold with no confident interface is not a binder.
3.  **Raw contacts do not discriminate.** lysozyme shows **38** raw geometric
    contacts and the scrambled decoy **39** — a large fraction of barstar's 52 —
    but once PAE-gated they drop to **0** and **1**. Rank on raw contacts and
    you would have shortlisted a non-binder.
4.  **Interface PAE confirms the call.** barstar 2.65 Å vs both decoys > 22 Å.
5.  **The filter does the right thing.** `--min-iptm 0.8` keeps **1 of 3**
    (barstar) and drops both decoys — an honest, mostly-empty shortlist.

## Key Takeaway

A non-binder can fold well (lysozyme) and rack up raw contacts (both decoys).
Neither is binding. Trust **iPTM** and **interface PAE**, corroborate with
**PAE-gated** contacts, and never rank on binder pLDDT or geometric contacts.
When the honest answer is "none of these binds," an empty shortlist is the
correct result.
