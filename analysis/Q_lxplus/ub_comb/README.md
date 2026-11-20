# ub_comb — wspólne μ par D–D z próbek bez triggera dżetowego (projekt A, 19.11.2026)

`ubcomb_run.py`: 427929 (MB, pary UB427, N_D* NDS) + podzbiór 428770 (maska z bitów trigPassed; normalizacja jak UBTRIG: n_pair/n_real × r_S,
tło tylko z S, n_pair × N_NDS/N_biegu) → wspólne μ = min Σ profili NLL (próbka × kanał). Warianty: strict (łańcuchy bez obciążenia),
nonj12 (= klasy UBTRIG MBTS_noj20 ⊕ other; NIE czyste: sptrk_pt4/6/8, AFP, łańcuchy spoza listy), nonj12_noafp.
Kontrole: 427929 odtwarza UB427 (0,527 [0,175; 0,920]); 428770 całość w skali NDS = 1,50 (STOP inaczej).
Systematyki raportowane (budżet przejściowy): N_D* i model wag per run, mc_stat i tune wspólne, statystyka × f_jk (z UB427JK gdy jest).
lxplus: `./submit_uc.sh check|submit|status` (batch charm_ubcomb, EOS charm_q/ub_comb); etap ORCH `stages/ubcomb.sh`.
Kopie: `btrig_common.py` (add_paths bez charmbkg), `dstar_refit_aa.py` (AA_merge; nie charmtag z PF/src).
