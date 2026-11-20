# PREZ 19.11 — slajd 6 „Train on pairs of D mesons” + szkice topologii (projekt A)
Generowane: `make_p6.py` (lokalnie, .venv_h5). Wejścia: generator Pythia 8 A3 `AFP_ddpred/gen/A3/gen_*.tsv.gz` (500 plików, prawda),
`pn_scores.npz` = PAIRNET v1 (EOS charm_q/pairnet/res, 17.11, tylko MC, OOF; kopia lokalna), `PAIRNET/runs/v3/sbtrain/sbtrain.json` (σ_μ/σ_nom).
Liczby: `p6_numbers.json`.
- p6_schematic_topologies — schematy: flavour creation (back-to-back), gluon splitting (kolinearne), DPS.
- p6_gen_dphi / p6_gen_deta / p6_gen_dphi_deta_os / p6_gen_pt1pt2 — topologie par D0/D⁺ (kanał analizy, pT 3–30, |η|<2,5) na poziomie generatora; OS/LS jako proxy (generator nie zapisuje procesu partonowego), linia przerywana = pary niezależne (granica DPS).
- p6_pnet_arch — architektura sieci par (Deep Sets; wejścia nóg jak v1: 18 cech H + gatunek + logit H; v2/v3: noga D⁺ z cechami Dalitza).
- p6_pnet_score, p6_pnet_roc — score i ROC na parach MC (v1, tło z MC), gwiazdka = kategorie nominalne.
- p6_pnet_sigma — oczekiwane σ_μ/σ_nom z SBTRAIN v3 (tło z pasm danych 428770, sygnał MC×SF): pnetd 0,79.
