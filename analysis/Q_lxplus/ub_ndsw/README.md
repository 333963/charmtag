# ub_ndsw — szerokość piku D* w klasach triggera 428770 vs 427929 (19.11.2026, projekt A)

**Pytanie.** W 427929 fit Δm daje skalę σ piku D* 1,40 × MC (UB427_NDS, `PF2_lxplus/runs/lxplus/unblind/r427929/nds/NDS_427929.md`), w 428770 nominał 1,02.
Ten sam detektor. Czy szeroki pik jest prawdziwą własnością zdarzeń minimum bias, czy tło wciągnięte w sygnał w małej próbce?
Od tego zależy N_D*(427929): skala wolna 2046 → μ = 0,53; skala = MC 1417 → μ ≈ 0,77.

**Pomiar.** Ten sam fit (kopie `dstar_refit_nds.py`, `nds_var.py` z ub427_nds; zmiana: parametr `s0` i zwracane `nll`, granice skali 0,5–2,2)
na D* 428770 (h5 v1.1 `*_none.h5`, próbka nominału 05.11) w klasach triggera z tabeli `ub_trig/trg_*.npz` (lookup po (run, event)):
MBTS = MBTS_j20 ∪ MBTS_noj20, J12only, other, osobno MBTS_j20 / MBTS_noj20, całość; oraz 427929 (h5 po filtrze UB427) w tym samym kodzie.
Skala σ z profilu NLL (siatka 0,60–2,00, ΔNLL = 0,5), biny pT: całość 3–30, 3–5, 5–8, 8–30, 3–8.
Kontrola: całość 428770 musi odtworzyć nominał (wpisy w oknie ≈ 482979, N ≈ 4436, skala ≈ 1,02); 427929 musi odtworzyć 1,40.

**Kryterium (zapisane 19.11 przed wynikiem).** σ_s = połowa przedziału 68 % skali MBTS(428770), całość 3–30 GeV.
- s(MBTS) w 2σ_s od 1,02 i NIE w 2σ_s od 1,40 → **szerokość w 427929 nieprawdziwa** (tło wciągnięte w sygnał). Konsekwencja dla `budzet_427929.py`:
  nominał N_D* ze stałą skalą σ (1417 ± …, wariant `sig_scale_fixed`), μ(427929) ≈ 0,53 × 2046/1417 ≈ 0,77, a syst NDS liczony bez wariantu
  „skala wolna” (kategoria `sig` zdominowana dziś przez 628 z różnicy wolna/stała). **Zmiana nominału = osobna decyzja.**
- s(MBTS) w 2σ_s od 1,40 i NIE w 2σ_s od 1,02 → **szerokość prawdziwa** (zdarzenia MB mają szerszy pik); μ(427929) = 0,53 zostaje,
  różnica wolna/stała przestaje być syst (stała skala z 428770 nie dotyczy MB) → syst N_D* maleje.
- inaczej (zgodne z oboma albo z żadnym) → **nierozstrzygnięte**; syst N_D* 427929 zostaje (±0,22 na μ).
Bin 3–8 GeV (tam 427929 ma 1,5) — pomocniczo, nie zmienia werdyktu. J12only/całość 428770 — kontrola (oczekiwane ≈ 1,0).

**Uruchomienie.** `submit_nw.sh check|submit|status` (1 job, batch `charm_ubndsw`, h5 czytane wprost z EOS), wynik EOS `charm_q/ub_ndsw`.
Etap ORCH `ORCH/stages/ubndsw.sh` → `PF2_lxplus/runs/lxplus/unblind/r427929/ndsw`.
Test lokalny: `ndsw.py --selftest --ref ndstar_full_428770.json --out DIR` (syntetyczne: MBTS s = 1,4, J12 s = 1,0) — 19.11 OK: MBTS 1,49 [1,39; 1,60], J12only 0,96 [0,92; 1,01].
