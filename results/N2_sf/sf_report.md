# N2 — SF dane/MC z tagu D* (D0) i fitów masy (D⁺), ½ lumi ∩ podpróbki S (2026-11-04)

Wejście: score `L_charmtag/runs/v11_full/scores/{D0_H,Dplus_H}` (½ lumi) × podpróbki `S_bkggen/sub` (D0: 13 % wierszy; D⁺: pełne score 33 M). Tag D*↔D0 po bitach (d0, z0) obu córek (100 % dopasowań; D0 w D* nie jest refitowany). Fit Δm: odporny `dstar_refit.fit_passfail` (V): 2 Gaussy z MC truth, tło u^a e^{−bu}(1+c1u+c2u²), jednoczesny pass/fail, bez osobnego fitu „bez cięcia”.

## D0 (tag D*): ε_data (fit) / ε_MC (prawda fd) → SF

> **KOREKTA 13.11.2026.** SF i N_D* w tabeli przeliczone nowym `fit_passfail` (12.11: tło próg × pol2, wszystkie starty,
> przedział z profilu, flaga `determined`) na tych samych histogramach `out/dstar_hist_*.npz`; stary fit był źle uwarunkowany (10/20 fitów valid=False).
> Błędy: **± stat z HESSE** (⊕ błąd ε_MC prawdy), jak dotąd; przedziały profilu 68 %: `NUM_audit_2026-11-12/passfail/compare_tables.md` (sekcja N2),
> surowe `…/passfail/out/real_N2_V_T4.jsonl`. Stare wartości: `sf_report.md.przed_2026-11-13`; `out/sf_D0.json` NIE zmieniony (nadal stary fit).
> **Wniosek zmieniony:** 3–5 i 5–8 GeV w całości nieokreślone (low_signal: N_D* < 3σ) — liczby w nawiasach nie są pomiarem; 8–30 i zintegr. określone.
> Syst. i closure MC poniżej przeliczone 13.11 (cz. 2) — zob. noty przy nich.

| pT | WP80 | WP70 | WP50 | WP65 | WP90 | N_D* (fit) |
|---|---|---|---|---|---|---|
| 3–5 | nieokr. (low_signal) (0,72±0,34) | nieokr. (low_signal) (0,44±0,22) | nieokr. (low_signal) (0,27±0,16) | nieokr. (low_signal) (0,51±0,26) | nieokr. (low_signal) (0,62±0,28) | (197–204±96) |
| 5–8 | nieokr. (low_signal) (0,76±0,33) | nieokr. (low_signal) (0,75±0,39) | nieokr. (low_signal) (1,14±0,67) | nieokr. (low_signal) (0,74±0,39) | nieokr. (low_signal) (0,77±0,34) | (64–141±36–63) |
| **8–30** | **1,05±0,19** | **1,08±0,21** | **1,18±0,24** | **1,20±0,25** | **1,05±0,17** | 118–141±27–35 |
| zintegr. | 0,80±0,22 | 0,67±0,19 | 0,69±0,21 | 0,74±0,21 | 0,74±0,19 | 363–374±104 |
| wrzesień (zintegr., v2) | 0,931 | 0,915 | 0,882 | ~0,91 | 0,938 | — |
Syst (wariant fitu stały μ/σ, zliczenie z odjęciem tła, closure MC; ⊕ jak `n2_fit.py`): 8–30 ≤0,21 (0,02/0,03/0,07/0,21/0,06 @80/70/50/65/90; WP65 z sb_count −0,19, reszta ≤0,07);
5–8 0,09–0,33, 3–5 0,58–1,56 — nie do użycia (tam także nominał i wariant stały μ/σ nieokreślone, low_signal).
> **KOREKTA 13.11.2026 (cz. 2).** Linia syst. przeliczona nowym `fit_passfail` (12.11) i nowym `sb_count` (NLL nasycone; w części binów inne minimum fitu boczne)
> na `out/dstar_hist_*.npz`; było: 8–30 ≤0,15 (0,07/0,15/0,13/0,14/0,07), 3–5 ≥0,7. Tabela per bin: `NUM_audit_2026-11-12/passfail_rest/compare_tables_rest.md`
> (sekcja N2 syst.). Wniosek: 8–30 syst. nadal mniejsza od stat. (±0,17–0,25) — SF 8–30 zgodne z 1 bez zmian; 3–5 bez zmian (nie do użycia).
- Artefakt V (stary fit: 0,41/0,80 @3–5): fit z 04.11 dawał 0,57±0,27 @80; **po korekcie 13.11 bin 3–5 (i 5–8) jest nieokreślony** (low_signal) —
  ani deficytu, ani jego braku nie da się tu stwierdzić (statystyka 13 % × niskie S/B: 3–5 ~200±96 D*, 5–8 64–141). 8–30: SF 1,05–1,20 ± 0,17–0,25,
  zgodne z MC (SF = 1) w ≤ 0,8σ i z wrześniem (zintegr.) w ≤ 1,2σ (dawniej 0,97–1,14, „w granicach 1σ”).
- Closure MC (ten sam fit na fd z tłem; ε_truth = zliczenie prawdy w 143–148,5 MeV): ε_fit/ε_truth (± σ_HESSE/ε_truth) 3–5 **1,004–1,016** (±0,008–0,014, pull +0,3…+2,1),
  5–8 **0,979–1,003** (±0,006–0,014, pull −1,5…+0,4), 8–30 **0,957–0,989** (±0,008–0,019, pull −2,3…−1,3; wszystkie 5 WP poniżej 1, WP w binie skorelowane).
  Średnia ważona 0,998 ± 0,003. **Wniosek zmieniony:** „metoda nieobciążona” nie jest uzasadniona dla 8–30 — obciążenie −1…−4 % (1,3–2,3σ na WP);
  w 3–5 +0,4…+1,6 %. Wielkość ≪ stat. SF danych (±0,17–0,25 w 8–30) i jest już w syst. (ΔSF_clos). Część efektu może być definicyjna
  (fit liczy cały sygnał, prawda tylko okno 143–148,5 MeV) — nie rozdzielone.
> **KOREKTA 13.11.2026 (cz. 2).** Closure przeliczone nowym `fit_passfail` (12.11); dawny zapis „0,99–1,02 (8–30), 0,96–1,07 (3–5)” nie zgadzał się nawet z
> `out/sf_D0.json` (stary fit: 8–30 0,928–0,986, 3–5 0,998–1,017). Wszystkie 15 fitów MC valid i określone. Tabela: `NUM_audit_2026-11-12/passfail_rest/compare_tables_rest.md`.

## D⁺ (inkluzywne fity masy; 2 Gaussy z MC + exp(pol2))
Fit „bez cięcia” w 3–5 i 5–8 niestabilny (S/B ≪ 1 %) → **tylko metoda względna do WP90** (jak wrzesień):
| pT | rel80 | rel70 | rel50 | rel65 | błąd |
|---|---|---|---|---|---|
| 3–5 | 1,24 | 1,36 | 1,45 | 1,37 | ±0,3 (N90 2248±593) |
| 5–8 | 1,19 | 1,24 | 1,13 | 1,21 | ±0,14 |
| **8–30** | **1,01** | **0,98** | **0,91** | **0,97** | ±0,05 |
| wrzesień v2 (cum.) | 1,00 | 1,00 | 0,93 | ~0,99 | — |
SF_abs(WP) = rel(WP) × SF(WP90); SF(WP90) z D0 (transfer, +5 % syst). 8–30 zgodne z v2 (0,91 vs 0,93 @50).

## Werdykt i v3
- Lokalna statystyka nie poprawia wrześniowych SF; **v3 = wartości v2 + jawne niepewności** (`unc`: 0,07–0,10) + dodane klucze 0,65/0,3 (kategorie par). Plik `sf_for_pairs_v3.json`, override `Q_lxplus/pf_pairs/cfg/sf_v3.yaml` (NIE wysłane).
- Wpływ na Z par: Y ∝ SF_D0·SF_D⁺ per kategoria → ±0,08–0,13 na Y → **ΔZ ≈ ±8–13 %** (systematyka, nie przesunięcie centralne).
- Docelowo: `v` (pipeline) → `dstar_nom` histogramy pełnej lumi → `n2_fit.py` czyta ten sam format (h_<bin>_<wp>_pass/fail) → SF per bin z błędem ~±0,05 @8–30, ~±0,1 @3–5.

## v4 (2026-11-05, po B2b — EV pełna lumi) — `n2_v4.py` → `sf_for_pairs_v4_proposed.json`
> Adnotacja 13.11.2026: wejścia v4 (SF B2b z fitu pass/fail D*, SF95 5–8 „0,82–1,00” z V) pochodzą ze starego `fit_passfail` i nie zostały przeliczone nowym (12.11).
D0 (species 10) = średnia ważona binów (wagi N_D* B2b: 3–5 0,17 / 5–8 0,33 / 8–30 0,50; wynik niewrażliwy na wagi: ±0,01 przy 1:1:1 i 2:1:1):
8–30 SF abs B2b 1,00/0,98/0,96/0,96/0,92 ±0,03–0,04; 5–8 SF95 0,91±0,10 (syst pokrywający 0,82–1,00) × SF_rel sW; 3–5 0,91±0,15.
**v4 D0: 0,949/0,935/0,939/0,937/0,931 ±0,047–0,049 @90/80/70/65/50** (v3: 0,938/0,931/0,915/0,911/0,882 ±0,07–0,10) — zgodne w 1 σ, błąd ~2× mniejszy; SF95(D0) = 0,956±0,047 (v3 zakłada 1 — różnica to deficyt 5–8/WP95).
D⁺ (4) = v3. Nominal unblind nadal v3 (przyjęte); v4 kandydat do `pairs.sf_json`.
