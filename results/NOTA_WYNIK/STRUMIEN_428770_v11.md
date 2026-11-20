# Strumień danych 428770 w v1.1 (Indico) — weryfikacja 20.11.2026 (tylko odczyt)

**Werdykt: `user.malewick.data.minbias.428770...v1.1.Aug7_EXT0` (PanDA 51907090) to strumień physics_Main — te same zdarzenia co v2.0 data_428770.**
„minbias” w nazwie to nazwa analizy/konfiguracji (`charm_minbias_dpt3to150`), nie strumień.

## Dowody
1. **Liczba zdarzeń w ntuplach:** suma `n_events_read` w h5 v1.2 (39 części, źródło v1.1 EXT0 51907090) = **7 989 496**.
   v2.0 data_428770 (97 części, EXT0 53347072, physics_Main): **7 989 496**, a więc identyczna co do zdarzenia.
2. **Nakładanie event_number:** 4 części v1.2 (000002, 000004, 000011, 000012; 581 430 zdarzeń) mają **100,00 %** zdarzeń
   w tabeli zdarzeń v2.0.
   - v2.0 to 7 989 496 unikalnych event_number, LB 4–323.
   - HLT_noalg_L1J12 odpalił w 56,1 % zdarzeń.
3. **Kontenery AOD runu 428770 (rucio):**
   - physics_Main: 11,53 M zdarzeń (r13928_p5279, r14755_p5604, r15869_p6304) oraz 11,50 M (r16487_p5949);
   - physics_MinBias: 16,52 M (r14956_p5604) oraz 16,50 M (r16546_p6832).

   7,99 M zdarzeń w ntuplach mieści się w physics_Main. Strumień MinBias (16,5 M) nie był użyty ani w v1.1, ani w v2.0.
4. **Ograniczenia sprawdzenia:**
   - Kontener v1.1 w rucio jest dziś pusty (0 plików), a `list-parent-dids` nic nie zwraca. Wejścia zadania PanDA nie da się odczytać bez proxy/SSO, więc dowodem jest punkt 2.
   - Lokalne `FEAT_charmnet/data/data_428770/*.h5` to schemat v1.2 z v1.1 (job 000002 i 000025), czyli również physics_Main.

Skrypt: `~/private/charm_q/logs_strm/ov2.py`, log `ov2.log` (lxplus).

## Konsekwencje
- **Dawny wynik 1,50 (v1.1/v1.2) i 1,27 (v2.0)** powstały na tych samych zdarzeniach physics_Main z dominacją L1_J12.
  - Różnica 1,50 → 1,27 to normalizacja (direct → prawda mb_nd), a nie dane.
  - Obciążenie triggerem J12 (test zerowy Z = 3,25) dotyczy w równym stopniu 1,50. Wniosek noty §7 („428770 = obserwacja obciążona triggerem”) bez zmian.
- **Plan „428770 strumień MinBias jako dodatkowe dane”** (INNE_RUNY_MB.md) pozostaje aktualny.
  - To naprawdę nowe zdarzenia: 16,5 M, ten sam run, ta sama kalibracja i pile-up μ≈0,05.
  - Warunek: deduplikacja po event_number względem physics_Main. Zdarzenie może trafić do obu strumieni, gdy przejdzie triggery z obu.
  - Menu/preskale MBTS w tym strumieniu trzeba sprawdzić (trigPS w produkcji v2.0) przed włączeniem do μ.
  - Nadal potrzebne: produkcja ntupli (w ramach produkcji v2.0 lub własne zadania grid).
