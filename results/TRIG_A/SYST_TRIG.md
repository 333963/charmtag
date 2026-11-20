# SYST_TRIG — systematyka „trigger / selekcja zdarzeń” w budżecie μ (PRERejestracja 19.11.2026, przed wynikiem UBTRIG)

Projekt A (pary D–D). Spisane **zanim** poprawiony `ubtrig` (μ par per klasa triggera, run 428770) dał wynik; jedyne znane już liczby to UB427
(18.11: μ₄₂₇₉₂₉ = 0,53 [0,17; 0,92] vs μ₄₂₈₇₇₀ = 1,50 [1,10; 1,90], Δ = −0,97, −1,6σ / −1,5σ z N_D*), emulacja J12 (TRIGCAL) i udziały N_D*
klas (lokalnie 19.11: J12only 0,769 ± 0,019, MBTS_j20 0,135 ± 0,015, MBTS_noj20 0,010 ± 0,004, other 0,084 ± 0,012). Reguły poniżej
są zamrożone; licznik: `TRIG_A/syst_trig.py` (bez parametrów do strojenia po wyniku). Zmiana reguł po wyniku = tylko osobną decyzją, z adnotacją.

## 1. Co mierzymy
Obciążenie triggerem = względna nadwyżka μ nominału (428770, 77 % N_D* w zdarzeniach tylko z L1_J12) nad μ w zdarzeniach nieobciążonych:

  d ≡ (μ_nom − μ_czyste) / μ_nom      (d > 0 ⇔ trigger zawyża μ; korekta μ_corr = μ·(1 − d))

Model: klasa J12 ma [pary/D*] × κ względem MB przy tych samych (n_trk, lead_pt) (κ zdefiniowane PO wagach w2 — TRIGMC/TRIGCAL),
reszta κ = 1 ⇒ μ_nom = (r_J12·κ + 1 − r_J12)·μ_prawd., d = r_J12(κ − 1)/(1 + r_J12(κ − 1)).

Dwa niezależne pomiary d (w obu skala N_D* wypada albo jest w błędzie):
- **d₁ (UB427, między runami):** d₁ = (μ₄₂₈₇₇₀ − μ₄₂₇₉₂₉)/μ₄₂₈₇₇₀, σ² = (f_jk·σ₄₂₇₉₂₉)² + σ₄₂₈₇₇₀² + (μ₄₂₇₉₂₉·ε_N7)² + (μ₄₂₈₇₇₀·ε_N8)², podzielone przez μ₄₂₈₇₇₀.
  σ₄₂₇₉₂₉ = profil 68 % (bez korelacji zdarzeń) × **f_jk** = σ_jk/σ_profil z UB427JK (`Q_lxplus/ub427_jk`); dopóki brak — f_jk = 1,18 (428770, jk_ev200),
  oznaczone „założone”. σ₄₂₈₇₇₀ = 0,46 (już z ×1,18). ε_N = względny błąd stat. N_D* (0,237 / 0,151); kształt N_D* wspólną metodą — pomijany w Δ.
- **d₂ (UBTRIG, w runie 428770, klasy rozłączne):** d₂ = r_J12·(μ_J12 − μ_reszta)/μ_nom; μ_reszta = średnia ważona klas ≠ J12only z dopasowanym μ;
  σ klasy = (σ₆₈ × 1,18) ⊕ μ_c·δr_c/r_c (1,18 = czynnik zdarzeń z 428770; klasy to podzbiory tego runu). Skala μ (μ_nom = 1,19 vs 1,50) wypada w stosunku.
  Kontrola (informacyjnie): to samo z resztą „czystą” = MBTS_noj20 ⊕ other (bez HLT_j20).
- **d (wynik):** średnia ważona d₁, d₂ (traktowane jako niezależne: d₁ — inny run, d₂ — rozłączne klasy; wspólny tylko mianownik μ_nom).
  Zgodność |d₁ − d₂|/σ > 2 → flaga „niezgodne” (raport, wynik = d z większym σ_d = ostrożnie: σ_d → |d₁ − d₂|/2 ⊕ σ_d).
  Z_trig = d/σ_d.

Emulacja (TRIGCAL, kalibracje z χ²/ndf < 2: em 27,5/4 i all 40/9): κ k=0..2 = 1,146–1,388; **κ_emu = 1,24 ± 0,10** (środek ± połowa zakresu)
⇒ d_emu = r_J12(κ_emu − 1)/(1 + r_J12(κ_emu − 1)) = 0,156 (zakres 0,10–0,23) przy r_J12 = 0,769. Klasa MBTS_j20 (HLT_j20, też dżetowa) — bez kalibracji, κ = 1 (uwaga w raporcie).

## 2. Bramki ważności (STOP → wiersz „brak”, nic nie wchodzi do budżetu)
- UBTRIG: wszystkie r_c określone i skończone, Σ r_c ∈ [0,95; 1,05]; `nominal_ok`; Σ klas = pairs_full we wszystkich kanałach; nogi bez klasy = 0;
  fit wspólny μ dla J12only i ≥ 1 innej klasy. (Bieg 18.11 z r_c = NaN → BRAK.)
- UB427: μ₄₂₈₇₇₀ odtworzone (|μ − 1,50| < 0,01), σ skończone.
- Bez ważnego UBTRIG: tylko d₁ informacyjnie, decyzja odłożona (status „niekompletne”).

## 3. Reguły decyzji (μ = 1,50, zapis bieżący)
| warunek | status | wiersz syst. „trigger” | inne wiersze |
|---|---|---|---|
| \|Z_trig\| < 2 i d zgodne z d_emu (\|d − d_emu\|/σ_d ≤ 2) | brak istotnego efektu | **−μ·d_emu** (jednostronnie w dół; ≈ −0,23 przy d_emu = 0,156) | „model wag ±0,24” i κ_cls bez zmian |
| \|Z_trig\| < 2, d niezgodne z d_emu | emulacja odrzucona | −μ·max(\|d\|, σ_d) jednostronnie w dół jeśli d ≥ 0, inaczej ±μ·σ_d | bez zmian |
| 2 ≤ Z_trig < 3 | wskazówka | wiersz łączny „selekcja zdarzeń (trigger + wagi)” = max(μ·d, √(0,24² + (μ·d_emu)²)), jednostronnie w dół; **zastępuje** „model wag ±0,24” | κ_cls bez zmian |
| Z_trig ≥ 3 | efekt | jak wyżej + wariant μ_corr = μ·(1 − d) ± μ·σ_d w raporcie | zmiana nominału **wyłącznie osobną decyzją** |
| Z_trig ≤ −2 | nadwyżka w czystych | ±μ·\|d\| (dwustronnie), raport | bez zmian |

## 4. Podwójne liczenie
- **Model wag ±0,24** (w2 ↔ w1, mapa (n_trk, lead_pt)) opisuje rozkład zdarzeń; κ jest zdefiniowane przy ustalonych (n_trk, lead_pt) (po w2), więc w
  przypadku |Z| < 2 wiersze są rozłączne i dodawane w kwadraturze. Pomiar d₁ (między runami) zawiera też różnicę modelowania c_w obu runów — dlatego
  przy Z ≥ 2 nie dodaje się d do wiersza wag, tylko oba zastępuje jeden wiersz „selekcja zdarzeń” (max z pomiaru i z kwadratury).
- **κ_cls ±0,06** (korelacja ε klasyfikatora między nogami, MC fd_nd) — niezależne od triggera (ten sam klasyfikator we wszystkich klasach) → bez zmian.
- **N_D*** — d₂ niezależne od N_D* (stosunki w runie); d₁ zawiera stat. N_D* obu runów; wiersz N_D* ± 0,33 bez zmian.
- **Tune ±0,17, mc_stat ±0,46** — bez zmian (MC normalizacja wspólna dla wszystkich klas/runów, wypada z d).

## 5. Wyjście
`syst_trig.py --ubtrig ubtrig.json --ub427 ub427.json [--ub427jk ubjk.json]` → `SYST_TRIG.json` + tabela md: d₁, d₂, d, Z_trig, d_emu, status, wiersz.
Wynik do przeglądu przed wpisem do `UNB_WYNIK.md` / `VER_budget`.
