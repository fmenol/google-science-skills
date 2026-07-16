# Example: A TRUE Binder (barnase : barstar)

**Target:** barnase (110 aa, chain A) **Binder:** barstar (89 aa, chain B)
**Verdict:** **Confident interface — iPTM 0.964, 52 confident contacts**

Barstar is barnase's natural intracellular inhibitor and one of the tightest
protein–protein complexes ever measured (Kd ~ 10⁻¹⁴ M). It is the calibration
point for the whole skill: the shared API reference says barnase + barstar
gives iPTM 0.96, and this worked example reproduces it.

## How it was produced

Real output from this skill in replay mode (no credits), at the eval settings:

```bash
uv run --no-project scripts/screen_binders.py analyze-interface \
  --target-fasta barnase.fasta --binder-fasta barstar.fasta \
  --output barnase_barstar/ --num-loops 3 --num-sampling-steps 50
```

## Why this is a good example

This is what an unambiguous, real, high-affinity interface looks like through
ESMFold2 — the positive control every screen should be read against:

1.  **iPTM 0.964** — top of the "confident interface" band (> 0.8), and stable
    at ~0.96 across every sampling setting tried.
2.  **Interface PAE 2.65 Å** — the confirming metric agrees: the model places
    the two chains confidently relative to each other.
3.  **52 confident contacts — and all 52 geometric contacts survive the PAE
    gate.** A true interface is not just close packing; it is close packing the
    model is confident about.
4.  **The interface is the *real* one.** The model recovers barstar's acidic
    recognition loop — Asp35 (the single biggest hotspot) and Asp39 — plugging
    barnase's active-site residues His102, Arg59 and Arg83. This is the textbook
    binding mode, not an artefact.
5.  **The PAE heatmap is dark everywhere**, including the off-diagonal cross-
    chain blocks. Compare it against the pale off-diagonal blocks of the decoys
    in `../decoy_controls/`.

## Key Takeaway

A real binder shows **iPTM > 0.8 AND low interface PAE AND confident (PAE-gated)
contacts**, all agreeing, with a biologically sensible interface. Read every
screen against this control. But note the discipline even here: iPTM 0.964 is a
**structural-confidence prediction**, not the Kd. We happen to know barnase :
barstar binds; for a novel design, a score like this is a strong hypothesis to
put in front of the wet lab — not a measured affinity.
