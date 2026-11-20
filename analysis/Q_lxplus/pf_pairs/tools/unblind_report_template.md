# Pary D–D (D⁰D⁰, D⁺D⁰, D⁺D⁺) — wynik po odsłonięciu, pełna lumi 2022 (run 428770), weto śladów d0z0

{{header_note}}
Wygenerowano {{generated}} · wyniki `{{ub}}` · ślepy bieg `{{pf}}` · normalizacja `{{norm_mode}}` · N_D* (dane) = {{ndstar}} ({{ndstar_method}}; stary `fit_dm`: {{ndstar_old}}).

## 1. Wynik
| kanał | μ̂ (dane/Pythia) | 68 % (profil) | σ_stat (Hesse) | σ_syst (U) | Z obs | Z oczek. (ślepy Asimov) | p lok. | dev/ndf | S oczek. (pełna płaszczyzna) |
|---|---|---|---|---|---|---|---|---|---|
{{main_rows}}

**Wspólne μ = {{mu_common}} {{mu_common_68}}, Z = {{Z_comb}}σ (p = {{p_comb}})**; oczekiwane (ślepo, Asimov) {{Z_exp_comb}}σ.
Z systematykami kształtu profilowanymi: μ = {{mu_comb_syst}}, Z = {{Z_comb_syst}}σ. Wspólna systematyka μ (U, total_common): ±{{syst_common}}.
{{kappa_lines}}

## 2. Kontrole stabilności
| kanał | bias toys @fit | μ po korekcie biasu | σ(pull) | pokrycie 68 % | μ bez kategorii z tagiem D* | κ=channel: μ (Z) |
|---|---|---|---|---|---|---|
{{toy_rows}}

## 3. Rzuty m₁, m₂ (suma po kategoriach score) — dane vs fit
{{figures}}

## 4. Budżet systematyk μ (U, po PF)
| kanał | σ_stat | kształty | SF | normalizacja | κ, kategorie | suma syst | μ (syst prof.) | Z (syst prof.) |
|---|---|---|---|---|---|---|---|---|
{{syst_rows}}
{{extra_syst}}

## 5. Porównanie z wrześniem (v1.0, 27.09: 3,6σ, stare weto śladów — nieskuteczne dla D⁺)
| kanał | pary oczek. (Pythia, pełna płaszczyzna) | wrzesień μ / Z | teraz μ / Z |
|---|---|---|---|
{{cmp_rows}}
Wrzesień łącznie: μ = {{mu_sept}}, Z = {{Z_sept}}σ → teraz μ = {{mu_common}}, Z = {{Z_comb}}σ.

## 6. Uwagi
- Ślepy bieg `pf_pairs` liczył N_D* starym `fit_dm`; węzeł `unblind` (src 2.0.2) liczy `dstar_refit` — μ skaluje się liniowo z 1/N_D*; Z nie zależy od N_D*.
- Toys przy μ̂ (100/kanał), kontrola bez tagu D* i κ=channel to część węzła `fit2d --unblind` (bez zmian konfiguracji po decyzji; hash cfg w UNBLIND_OK).
