# T-FIN — pary D–D po poprawce weta wspólnych śladów (PF) — {{mode}}

Wygenerowano {{generated}} z `{{pf}}`{{unblind_src}}. Odniesienia: ½ lumi (Mac) 1,65σ · wrzesień v1.0 2,63σ · pełna lumi tilt+rho (V, stare weto) 2,27σ.

## 1. Weto wspólnych śladów (pairs_counts.json)
track_id `{{track_id}}`: par przed wetem {{pairs_before_veto}}, usuniętych {{pairs_shared_track}} ({{shared_pct}} %), w histogramach {{pairs_kept}} (region {{region}}).

## 2. Tło i czułość — ślepo (compare_pf)
| kanał | SB×SB nominal | SB×SB PF | Δ | Z_A nominal | Z_A PF | σ_μ nominal | σ_μ PF | dev/ndf nom | dev/ndf PF | ρ PF | tilt PF |
|---|---|---|---|---|---|---|---|---|---|---|---|
{{blind_rows}}
Z_A łączne (kwadratura): nominal {{Zq_nominal}} → PF **{{Zq_pf}}**.

## 3. Toys (PF, ślepo, μ=1)
| kanał | σ(pull) | bias/σ_μ | pokrycie 68 % | n_toys | G1 | G2 | G3 |
|---|---|---|---|---|---|---|---|
{{toys_rows}}

## 4. Systematyki po PF (U + Y + S8 + N)
| kanał | σ_stat | kształt+tło+κ | SF | normalizacja | exp/pol2 (Y) | SALAD (S8) | suma syst | syst/stat | G7 |
|---|---|---|---|---|---|---|---|---|---|
{{syst_rows}}
Źródła: `pairs/syst/*.json` ({{n_syst_json}} plików), `syst_bkg_model_Dplus.json` (Y), `bkg_shape_salad_*.json` (S8). Brakujące → „—” i G7 = ?.

## 5. Mixing (S5 nom_mix)
{{mix_line}}

## 6. Checklist GO/NO-GO (automatyczna część; reszta ręcznie w unblinding_checklist.md)
| # | kryterium | D0D0 | D⁺D0 | D⁺D⁺ |
|---|---|---|---|---|
{{gate_rows}}
{{g6_line}}

Werdykt automatyczny: **{{verdict}}** (G1–G4, G6, G7, G11 wszystkie kanały; G5, G8–G10, G12 ręcznie).

## 7. Wynik odsłonięty
{{unblind_section}}

## 8. Porównanie z wcześniejszymi liczbami
| | Z | uwagi |
|---|---|---|
| ½ lumi (Mac, 23/47) | 1,65 | stare weto |
| wrzesień v1.0 (35 plików) | 2,63 (toys 2,8) | stare weto, inny klasyfikator |
| pełna lumi tilt+rho (V) ślepo | 2,27 | stare weto |
| **PF ślepo (Asimov)** | **{{Zq_pf}}** | poprawne weto |
| **PF odsłonięte** | **{{Z_obs}}** | {{Z_obs_note}} |
