# Single-D classifier — models (slide, point 3)

Task: separate true D mesons from combinatorial background, candidate by candidate.
- D⁰ → K⁻π⁺: 2 tracks
- D⁺ → K⁻π⁺π⁺: 3 tracks

## Training sample (same for all four variants)

**Signal**
- Truth-matched D in forced-decay MC (Pythia8 A3, ND).
- D⁰ K↔π swapped hypothesis excluded.

**Background**
- Combinatorial candidates (truth 8–99) in minimum-bias MC (Pythia8 A3, ND).
- MC only. No data is used in this comparison.

**Selection** (same as the nominal analysis)
- 3 < pT < 30 GeV.
- Mass in the species range.
- L_xy ≥ 0.25 mm (preselection in the h5 production).
- D* veto for D⁺.

**Training and evaluation region**
- Signal window |m − m₀| < 55 MeV for D⁰ and < 45 MeV for D⁺.
- Background outside the window is used only for the mass-sculpting test.

**Weights**
- Background reweighted to the signal (pT, |η|) distribution, as in the nominal training.
- Classes balanced.

**Validation**
- 5-fold cross-validation by event, with the same hash fold as the nominal BDT.
- All quoted numbers are out-of-fold.
- 10 % of the training folds are used for early stopping.

## Inputs

**II. High-level: the 18 features of the nominal BDT "H"**
- Candidate: log pT, |η|, vertex χ², L_xy, impact parameter significance (signed and |·|), |z0| of the impact parameter.
- Pointing: −log₁₀(1 − cos θ_P).
- Decay angle: |cos θ*|.
- Daughter d0: min |d0| and max |d0|.
- Daughter z0 and η: max |z0| and max |η|.
- PID: dE/dx of the K, log p of the K, number of dE/dx hits of the K; max dE/dx of the π, min log p of the π.

**I. Low-level: the daughter tracks plus the candidate L_xy**

Per track there are 11 numbers:
- log pT, η
- Δφ and Δη with respect to the D direction
- log p
- dE/dx, plus a "missing dE/dx" flag
- number of dE/dx hits
- d0, z0
- role (K or π)

Notes:
- Charge enters through the role, because candidates are built with the right charge combination.
- The absolute azimuth φ is a detector symmetry, so it is replaced by Δφ.
- The vertex fit quantities (χ², cos θ*, pointing) are not given. The model has to infer them from the tracks.

## Models

| | BDT (XGBoost) | NN |
|---|---|---|
| high-level | 18 features; ≤ 1500 trees, depth 6, learning rate 0.05, subsample and colsample 0.8, min_child_weight 5, early stopping 50 (nominal settings) | **MLP** 18 → 128 → 128 → 64 → 1, SiLU, dropout 0.1; ≈ 27.5 k parameters (+1 input per feature with missing values) |
| low-level | flat vector: K, then π ordered by pT, plus L_xy; missing values handled natively | **Deep Sets**: shared track encoder φ: 11 → 64 → 64 → 64; pooling sum ⊕ max over tracks, ⊕ L_xy; head ρ: 129 → 128 → 64 → 1; ≈ 34 k parameters |

**NN training**
- Inputs: rank → Gauss transform, fitted on the training folds; missing values set to the median.
- Loss: weighted binary cross-entropy.
- Optimiser: AdamW, learning rate 2·10⁻³, weight decay 10⁻⁴, one-cycle schedule.
- Batch 4096, at most 30 epochs, early stopping with patience 5.

**Why Deep Sets for the low-level input**
- The tracks form an unordered set, and the two pions of D⁺ are interchangeable.
- A permutation-invariant encoder builds this symmetry in exactly. It does not have to learn it from data.
- It shares weights between tracks, with fewer parameters than an MLP on concatenated tracks.
- The same network handles 2 and 3 tracks.
- The K/π role flag keeps the asymmetric information.
- A transformer or graph network brings no benefit for sets of 2–3 elements.

## Metrics
- ROC and AUC (± spread over folds).
- ε_B at ε_S = 90 % (the nominal working point WP90).
- ε_S at the ε_B of the nominal BDT at its WP90, i.e. at fixed background.
- Score distributions for signal and background.
- Feature importance for the high-level inputs:
  - permutation ΔAUC for both models;
  - mean |SHAP| for the BDT.
- Mass sculpting: background ε at WP90 inside the window / outside the window, and the background mass shape before and after the cut.

**Reference:** the production BDT "H" is applied to the same candidates, using the fold-matched model. It was trained against data sidebands, not MC background, so it is not optimised for this test sample.
