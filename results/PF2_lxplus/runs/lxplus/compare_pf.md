# PF2 — powtórka par z poprawnym wetem wspólnych śladów (ŚLEPO)

- track_id: `d0z0`; pary przed wetem 296782, usunięte (wspólny ślad) 226944 = 76.5 %; w histogramach 56657 (region blind)

| kanał | SB×SB nominal | SB×SB PF | zmiana | Z_A nominal | Z_A PF | σ_μ nominal | σ_μ PF | dev/ndf nominal | dev/ndf PF |
|---|---|---|---|---|---|---|---|---|---|
| D0D0 | — | — | — | 1.38 | 1.90 | 0.79 | 0.57 | 0.79 | 0.65 |
| DpD0 | — | — | — | 1.60 | 2.71 | 0.69 | 0.41 | 1.79 | 0.94 |
| DpDp | — | — | — | 0.47 | 1.10 | 2.72 | 0.90 | 3.05 | 1.07 |

Z łączne (kwadratura): nominal 2.1668303597163074 → PF 3.4866050989078103; odniesienia: ½ lumi 1,65σ, wrzesień 2,63σ, pełna lumi tilt+rho (V) 2,27σ — wszystkie ze starym wetem.

- ρ/tilt (PF): {'D0D0': None, 'DpD0': None, 'DpDp': None} — oczekiwane ρ≈0 po wecie (MC: 0,24→0,01 w D⁺D⁺); jeśli ρ nadal ≠0 → model tilt+rho zostaje.
- toys: patrz toys_pf w JSON (pokrycie, pull, σ_μ mediana)
- syst (U) po PF: 7 plików; mixing (S5): mix_control.json

Szablon do docu T: kolumny powyżej + wiersz „wrzesień/½ lumi (stare weto)” dla kontekstu. Odsłonięcie: decyzja po tej tabeli.
