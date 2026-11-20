# docs — mapa

| plik | dla kogo | co zawiera |
|---|---|---|
| [double_charm.md](double_charm.md) | projekt A | problem → metoda (klasyfikator D z OOF, WP, fit masy, SF z D*, pary D–D z wetem wspólnych śladów, fit2d, normalizacja truth, systematyki) → jak uruchomić → stan (ślepo) → otwarte |
| [signal_like_background.md](signal_like_background.md) | projekt B — **kopia historyczna**; aktualna w repo charmbkg | tło pod pikiem / obiekty „podobne do sygnału”: SB+weto, SALAD, FETA, flow/CATHODE, candrot, pseudo-WS — co działa, co nie, liczby closure; zastosowanie jako systematyka |
| [data_format.md](data_format.md) | wszyscy | format h5 v1.1/v1.2: 60 kolumn, jednostki (GeV vs MeV), decayType, prawda, hipotezy K/π, duplikaty, czego brak |
| [running_on_lxplus.md](running_on_lxplus.md) | operator | HTCondor/EOS, DAG-i, orkiestrator; bez ścieżek prywatnych |

Kod: pakiet `charmtag` (projekt A, double charm). Projekt B (tło podobne do sygnału) od 15.11.2026 w osobnym repo **charmbkg** (pakiet `charmbkg`);
`charmtag.bkggen`, `charmtag.hyp` itd. to aliasy do `charmbkg`. Wejścia CLI: `charmtag …` (oraz `charmtag-bkggen …` = `charmbkg-bkggen …`).

Przegląd projektu (wynik, znane problemy, mapa repo): [`../PRZEGLAD.md`](../PRZEGLAD.md). Wyniki z liczbami: [`../results/README.md`](../results/README.md). Skrypty analizy i odtworzenie wyniku głównego: [`../analysis/README.md`](../analysis/README.md).
