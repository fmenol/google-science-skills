# Interface Analysis: barnase : barstar

## 1. Verdict

ESMFold2 predicts a **confident interface**: **iPTM 0.964** (top of the > 0.8
band), **interface PAE 2.65 Å**, **52 confident contacts**. Barnase and barstar
are a real, ultra-high-affinity complex (Kd ~ 10⁻¹⁴ M), and this is the
calibration point for the skill — the shape of an unambiguous hit. Read as a
prediction it is a strong hypothesis to test; iPTM is a structural confidence,
**not** the Kd. This ESMFold2 binder-screening skill produced the analysis.

## 2. Metrics

| Metric              | Value    | Reading                                        |
| :------------------ | :------- | :--------------------------------------------- |
| iPTM                | 0.964    | Confident interface (> 0.8)                    |
| interface PAE       | 2.65 Å   | < 5 Å — chains confidently placed              |
| binder pLDDT        | 0.937    | Barstar folds well on its own                  |
| confident contacts  | 52       | PAE-gated — the count to trust                 |
| geometric contacts  | 52       | All 52 survive the gate (nothing spurious)     |
| pTM                 | 0.967    | Whole complex confident                        |
| complex pLDDT       | 0.952    | —                                              |
| selection_score     | 0.946    | Composite (0–1)                                |
| isoelectric point   | 4.37     | Acidic; passes a pI < 6 minibinder filter      |
| min cross-chain dist| 3.72 Å   | Tight side-chain packing across the interface  |

iPTM and interface PAE agree, and every geometric contact clears the PAE gate —
the three independent signals all point the same way. That agreement, not any
single number, is what makes this a confident call.

## 3. Interface Residues

From the script's per-residue output — 19 residues on each chain. The contact
counts recover the *real* binding chemistry:

- **barnase (chain A) hotspots:** Arg59 (6 contacts), His102 (6), Tyr103 (5),
  Ala37 (4), Ser38 (4), Gln104 (4), plus Lys27, Arg83, Ser85 — barnase's
  active-site and substrate-binding residues: His102 is the catalytic histidine;
  Arg59 and Arg83 are the substrate-phosphate-binding arginines.
- **barstar (chain B) hotspots:** **Asp35 (7 contacts — the single biggest
  hotspot)**, Gly31 (5), Pro27 (4), Tyr29 (4), Asn33 (4), Gly43 (4), Trp44 (4),
  with Asp39 also at the interface. Barstar's acidic recognition loop (Asp35 /
  Asp39) is exactly the region that plugs and neutralises barnase's active site.

The model did not merely place two chains near each other; it reconstructed the
biologically correct interface, acidic loop into active site.

## 4. PAE Heatmap

![barnase:barstar PAE](barnase_barstar_pae.png)

*Fig 1: The entire PAE matrix is dark, including the two off-diagonal
cross-chain blocks either side of the red chain boundary. Dark off-diagonal
blocks = the model is confident where barstar sits relative to barnase. This is
the "whole matrix dark" signature of a confident interface. Contrast the pale,
washed-out off-diagonal blocks of the non-binders in `../decoy_controls/`.*

## 5. Interpretation & Limitations

ESMFold2 predicts a confident, biologically sensible barnase:barstar interface —
consistent with the known complex. What it does **not** establish, even at iPTM
0.964:

- **Affinity.** iPTM is not a Kd, ΔG, or IC50. We know barnase:barstar is
  femtomolar from experiment, not from this number.
- **Specificity.** The score is for this pair; it does not rule out barstar
  binding elsewhere.
- **Expression / solubility / stability.** Predicted only indirectly (pI 4.37 is
  a favourable coarse proxy).

For a *novel* design scoring like this, the next step is experimental: express
the pair and measure binding (e.g. SPR/BLI for a Kd, or a pulldown for a
yes/no). A confident interface is where the wet lab starts, not where it ends.
