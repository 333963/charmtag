# PREZ-2 — single-D classifier: BDT vs NN × high-level vs low-level inputs (MC, OOF, signal window)

Signal: forced-decay MC, truth-matched D; background: minimum-bias MC combinatorics (truth 8–99); |m − m₀| < window; background reweighted to signal (pT, |η|), classes balanced. WP90 = ε_S 90 %. „ε_S @ ε_B(nom)” = signal efficiency at the background efficiency of the nominal BDT H at its WP90. Sculpting = ε_B(in window)/ε_B(outside window) at WP90 (1 = mass-independent).

## D0  (signal 200101, background 1329197 candidates in window)

| model | size | AUC | ε_B @ ε_S 90 % | rejection gain vs nominal | ε_S @ ε_B(nom) | sculpting |
|---|---|---|---|---|---|---|
| BDT (XGBoost), high-level | 820 trees | 0.9438 ± 0.0005 | 0.1778 | 1.02 | 0.902 | 1.01 |
| NN (MLP), high-level | 27521 params | 0.9459 ± 0.0005 | 0.1715 | 1.06 | 0.906 | 1.01 |
| BDT (XGBoost), low-level | 1367 trees | 0.9476 ± 0.0005 | 0.1623 | 1.12 | 0.911 | 1.00 |
| NN (Deep Sets), low-level | 34049 params | 0.9490 ± 0.0003 | 0.1571 | 1.16 | 0.914 | 1.00 |
| BDT H (nominal, trained vs data SB) | – | 0.9427 ± 0.0005 | 0.1815 | 1.00 | 0.900 | 1.02 |

## Dplus  (signal 135929, background 879173 candidates in window)

| model | size | AUC | ε_B @ ε_S 90 % | rejection gain vs nominal | ε_S @ ε_B(nom) | sculpting |
|---|---|---|---|---|---|---|
| BDT (XGBoost), high-level | 1108 trees | 0.9713 ± 0.0002 | 0.0719 | 1.02 | 0.902 | 1.18 |
| NN (MLP), high-level | 27521 params | 0.9729 ± 0.0002 | 0.0673 | 1.09 | 0.907 | 1.18 |
| BDT (XGBoost), low-level | 1498 trees | 0.9788 ± 0.0002 | 0.0468 | 1.57 | 0.932 | 1.16 |
| NN (Deep Sets), low-level | 34049 params | 0.9811 ± 0.0002 | 0.0407 | 1.80 | 0.939 | 1.14 |
| BDT H (nominal, trained vs data SB) | – | 0.9710 ± 0.0002 | 0.0734 | 1.00 | 0.900 | 1.20 |

