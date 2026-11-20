"""charmtag.syst — systematyki fitu par D–D liczone ŚLEPO (etap U, 2026-11-03).

Metoda: dane Asimov z modelu w punkcie nuisance'ów dopasowanych na obszarze ślepym (μ=1, bez centrum 2D),
prawda zmieniona o wariant źródła → fit modelem nominalnym → Δμ = μ_wariant − μ_nominal (oba na Asimov, ten sam start).
Więzy (d_ij) mają środki w wartościach dopasowanych (globalne obserwable Asimov), więc nominalny Asimov daje μ = 1.
Wpływ na Z: Z_Asimov(wariant) − Z_Asimov(nominal) oraz Z z systematykami profilowanymi / z kwadraturą σ.
μ z danych nie jest liczone (centrum 2D nie istnieje w histogramach ślepych).

Moduły: core (kształty, SF, κ, kategorie), norm (tune, n_trk, dyfrakcja, stat. MC), files (jackknife po plikach,
podziały — ślepo: nadmierna dyspersja, czułość), wp (warianty WP/kategorii: przebudowa par), report (tabela).
Addytywnie: nie zmienia charmtag.pairs (używa jego API).
"""
SOURCES_DOC = {
    "width_D0": "szerokość piku D0 ±6%", "width_Dp": "szerokość piku D⁺ ±3%",
    "mscale_D0": "skala masy D0 ±1 MeV", "mscale_Dp": "skala masy D⁺ ±1 MeV",
    "bkg_poly3": "tło fałszywych nóg exp(poly2)→exp(poly3)", "swap_D0": "udział swapów K↔π D0 ±20%",
    "refl_Ds": "refleksja Ds→KKπ w D⁺ (5%, płaska)", "refl_Lc": "refleksja Λc→pKπ w D⁺ (1%, płaska)",
    "sf_D0": "SF D0 (score) ±1σ", "sf_Dp": "SF D⁺ (score) ±1σ", "sf_tag": "SF tagu D* ±14%",
    "kappa_channel": "κ: plan B ↔ jedno κ na kanał", "kappa_partner03": "κ: więz partnera σ=0,3",
    "kappa_partner05": "κ: więz partnera σ=0,5", "cat_Y": "wydajności kategorii score ±5% (niezależnie)",
}
