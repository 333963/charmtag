"""charmtag.pairs — pary D–D (D0D0, D⁺D0, D⁺D⁺) i fit 2D (etap O, 2026-11-02).

Port modelu M7 v3 + plan B (odsloniecie_wynik_v20260927_01/code/kcat.py) z rozszerzeniami:
systematyki kształtu (wynik_stabilnosc_v20260928_01/code/kcat_syst.py), więz κ partnera
(t12_kappa_partner_v20260927_01), toys/Asimov/pull. Moduł addytywny: nie zmienia reszty charmtag.
"""
SPECIES_DT = {"D0": 10, "Dplus": 4}
CHANNELS = {"D0D0": (10, 10), "DpD0": (4, 10), "DpDp": (4, 4)}
