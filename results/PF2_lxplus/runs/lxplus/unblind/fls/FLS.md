# f_LS — „rozbieżność” 0,33 ± 0,14 vs 0,22 / 0,34 / 0,49 (19.11.2026, projekt A)

**Wniosek: rozbieżności nie ma. Wpływ na μ (428770 i 427929) = 0 z konstrukcji. Nic nie zmieniać w nominale.**

## Skąd liczby

| liczba | co to | źródło |
|---|---|---|
| 0,22 ± 0,27 / 0,34 ± 0,18 / 0,49 ± 0,35 | f_LS zmierzone w **danych** 428770 per kanał (D⁰D⁰ / D⁺D⁰ / D⁺D⁺), test LS/OS (PD v1.2) | `PD_pairs_diff/runs/lxplus/tfin_diff.json` (`like_sign`), tabela `PD_fls_truth/fls_truth.md` §5 |
| 0,33 ± 0,14 | **średnia tych samych trzech pomiarów** | jw., `NOTA_WYNIK/nota.md` l. 154, 194 |
| 0,27 ± 0,03 (0,329 / 0,270 / 0,184) | oczekiwanie Pythii A3: prawda ważona ε + dylucja znaku D⁰ (w = 0,13) | `PD_fls_truth/fls_truth.md` §5, `UNB_WYNIK.md` l. 153 |
| 0,20 | prawda generatora bez ε i dylucji (slajd 6 prezentacji) | `PREZ_2026-11-19/deck/build_deck.py` l. 171 |

Sprawdzenie średniej (lokalnie): ważona odwrotnością wariancji 0,332 ± 0,138; ważona udziałem par danych (0,27/0,58/0,15) 0,330.
Kanały zgodne między sobą: χ² = 0,38 / 2 ndf. Uśrednienie poprawne.

## Czy f_LS wchodzi do μ?

Nie. Nominalny fit par (`AA_merge/src/charmtag/pairs/fit2d.py`, `model.py`) nie rozdziela par po znaku: szablon sygnału to pary Pythii
z f_LS Pythii wbudowanym. f_LS (`like_sign`, `f_mc`) używa tylko test różnicowy LS/OS (`pairs/diff.py`, μ_LS/OS = 1,46 [1,08; 1,84]).
`budzet_427929.py` i budżet 428770 nie zawierają f_LS. Δμ od zamiany „wspólne → per kanał” = **0** (nie ma czego zamieniać).

## Kontrola na danych (z PD, już policzona)

Dane vs oczekiwanie Pythii (prawda + dylucja): pulle −0,40 / +0,39 / +0,87σ, łącznie 0,33 vs 0,27 → +0,43σ (z profili fit2d 0,21σ).
Czysty DPS (0,5) −1,2σ, LHCb-podobne (0,09–0,12) +1,6σ. Zgodne z Pythią.

Nie zrobione: f_LS w czystym runie 427929. Lokalnie są tylko zbiorcze json (`r427929/pairs/`), bez par z ładunkami, więc trzeba by
puścić `diff.py` na lxplus. Przy ~4× mniejszej statystyce par błąd wyniósłby ~±0,3, więc test nie ma mocy.

## Rekomendacja

- W nocie pisać precyzyjnie: „f_LS (dane 428770) = 0,33 ± 0,14 (kanały 0,22 / 0,34 / 0,49), oczekiwanie Pythii 0,27 ± 0,03 (+0,4σ)”.
- Na slajdzie 6: 0,20 to prawda bez dylucji. Dopisać „(prawda generatora; po dylucji znaku D⁰ obserwowalne 0,27)”, żeby nie zestawiać tego z 0,33.
- Bez biegu na lxplus.
