"""Ręcznie napisane wagi (zero etykiet, zero treningu).

Sieć: P (13 cech przygotowanych) --W1,b1--> tanh --W2,b2--> sigmoid.
Każdy neuron ukryty = jedna "miękka decyzja fizyczna": tanh(a·(x - x0))  (rzadkie W1:
1 niezerowa waga, poza bramką AND z 2 wagami).  a = nachylenie (1/szerokość przejścia),
x0 = próg w przestrzeni cechy przygotowanej, w2 = waga głosu w logicie wyjściowym.

Konwencja słownika:  (nazwa, {cecha: a}, b1, w2, uzasadnienie)   gdzie b1 = -a·x0.
Wagi różnią się między D0 (cτ=123 µm, 2 tory) a D+ (cτ=312 µm, 3 tory).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .features import PREPARED, MASS_SAFE

@dataclass
class Unit:
    name: str
    w1: dict            # cecha przygotowana -> nachylenie a
    x0: float           # próg (dla bramki AND: próg sumy ważonej)
    w2: float           # waga wyjściowa
    why: str
    @property
    def b1(self) -> float:
        # tanh(a·(x−x0)) ⇒ b1 = −a·x0 dla neuronu 1-cechowego; dla bramki AND x0 jest progiem sumy ważonej ⇒ b1 = −x0.
        # (rev0 miało b1 = −x0 zawsze → neurony ze |a|≠1 miały przesunięty próg; vertex_chi2 i helicity były martwe. Poprawione 2026-11-04.)
        if len(self.w1) == 1:
            return -float(next(iter(self.w1.values()))) * self.x0
        return -self.x0

def _ln1p(x):  # pomocnik: próg w jednostkach surowych -> przestrzeń log1p
    return float(np.log1p(x))

# ---------------------------------------------------------------- D0 -> K pi
D0_REV0 = [
    Unit("lifetime", {"lxy_sig_t": 1.6}, float(np.arcsinh(1.5 / 3)), 1.5,
         "Lxy/σ: D0 leci ~cτγ≈0.1–0.5 mm, Lxy/σ sygnału 1–6; tło prompt symetryczne wokół 0 → próg 1.5σ, najsilniejszy głos."),
    Unit("pointing", {"point_t": 1.5}, 1.5, 1.2,
         "cosθ_pointing: dla prawdziwego D pęd || wektor lotu, 1-cosθ ~ 1e-3…1e-4 (point_t 3–4); kombinatoryka płaska (point_t 0–1.5) → próg 1.5."),
    Unit("flight_AND_point", {"lxy_sig_t": 1.0, "point_t": 0.6}, float(np.arcsinh(2.5 / 3) + 0.6 * 2.5), 0.6,
         "Bramka AND: duże Lxy/σ ma sens tylko gdy jednocześnie pointing dobry (duże Lxy + złe pointing = V0/konwersja/fałszywy wierzchołek) → premia za koniunkcję."),
    Unit("vertex_chi2", {"chi2_t": -1.4}, _ln1p(3.0), 0.6,
         "χ² wierzchołka 2 torów: sygnał ~1/ndf, tło przypadkowe szerokie → kara powyżej 3 (znak ujemny w a: duże χ² → -1)."),
    Unit("ip_min", {"ipsig_min_t": 1.7}, _ln1p(0.8), 0.8,
         "min IPsig córek: oba tory z D0 są przesunięte (IPsig ~1–3); prompt tło ma przynajmniej jeden tor z IP≈0 → próg 0.8σ."),
    Unit("ip_max", {"ipsig_max_t": 1.2}, _ln1p(2.0), 0.5,
         "max IPsig: przynajmniej jeden tor wyraźnie przesunięty (>2σ) – uzupełnia ip_min gdy drugi tor ma słabą rozdzielczość."),
    Unit("d0_abs", {"d0_min_t": 1.0}, float(np.arcsinh(0.03 / 0.05)), 0.3,
         "min|d0| w mm: ~cτ·sinθ ≈ 30–100 µm dla córek D; niezależne od σ(d0) → mała waga, odporność na błędną skalę błędów w MC."),
    Unit("helicity", {"cos_star": -4.0}, 0.7, 0.8,
         "|cosθ*|: D0 ma spin 0 → rozkład płaski; tło kombinatoryczne (miękki tor sparowany z twardym) piętrzy się przy |cosθ*|→1 → kara powyżej 0.7. UWAGA: skorelowane z m (rzeźbi) → tylko wariant A."),
    Unit("dedx_K", {"kness_K": 0.7}, 0.0, 0.7,
         "dE/dx toru K vs Bethe-Bloch: dla p<1.3 GeV K jonizuje ~1.5–3× MIP, π ~1× → dodatni kness = K potwierdzony; bramka pędowa zeruje tam, gdzie BB nie separuje."),
    Unit("dedx_pi", {"kness_pi1": 0.7}, 0.0, 0.4,
         "dE/dx toru π: odrzuca przypisania K↔π zamienione (odbicie masy) i protony; słabszy głos niż K bo π to hipoteza domyślna."),
    Unit("daughter_p", {"logp_K": 1.2}, float(np.log(2.0)), 0.3,
         "ln p_K: córki prawdziwych D są twarde (p_K>2 GeV przy pT(D)>3), kombinatoryka z miękkich torów → łagodna premia. Korelacja z m przez pT bin → tylko A."),
    Unit("eta", {"abs_eta": -1.5}, 1.8, 0.2,
         "|η|: przy |η|>1.8 rozdzielczość d0/masy spada, tło rośnie → mała kara; zmienia szerokość sygnału → tylko A."),
]
D0_BIAS = -1.2  # próg startowy: bez żadnego pozytywnego dowodu kandydat ma P≈0.23 (tło dominuje a priori)

# ---------------------------------------------------------------- D+ -> K pi pi
DP_REV0 = [
    Unit("lifetime", {"lxy_sig_t": 1.5}, float(np.arcsinh(2.5 / 3)), 2.0,
         "D+ żyje 2.5× dłużej niż D0 (cτ=312 µm) → próg Lxy/σ 2.5, największa waga; rozdzielczość 3-torowego wierzchołka lepsza."),
    Unit("pointing", {"point_t": 1.5}, 1.8, 1.2,
         "Dłuższy lot → kąt pointing jeszcze ostrzejszy; próg point_t 1.8 (1-cosθ < 1.6e-2)."),
    Unit("flight_AND_point", {"lxy_sig_t": 1.0, "point_t": 0.6}, float(np.arcsinh(3.5 / 3) + 0.6 * 2.8), 0.6,
         "AND lot+pointing jak dla D0, progi przesunięte w górę zgodnie z dłuższym cτ."),
    Unit("vertex_chi2", {"chi2_t": -1.4}, _ln1p(3.0), 0.8,
         "χ² 3-torowego wierzchołka: więcej ndf → silniejsza separacja przypadkowych trójek → waga 0.8."),
    Unit("ip_min", {"ipsig_min_t": 1.7}, _ln1p(0.8), 1.0,
         "min IPsig z 3 torów: wymaganie 'wszystkie 3 przesunięte' jest bardzo czyste; kombinatoryka prawie zawsze ma tor prompt → waga 1.0."),
    Unit("ip_max", {"ipsig_max_t": 1.2}, _ln1p(2.5), 0.4,
         "max IPsig: słabszy, bo przy 3 torach łatwo o jeden przesunięty przypadkiem → próg 2.5σ, waga 0.4."),
    Unit("d0_abs", {"d0_min_t": 1.0}, float(np.arcsinh(0.04 / 0.05)), 0.3,
         "min|d0| ≈ 40 µm+ dla córek D+ (dłuższy lot) – skalo-niezależna kontrola IPsig."),
    Unit("helicity", {"cos_star": -4.0}, 0.75, 0.5,
         "|cosθ*| dla 3 ciał mniej ostre (zależy od definicji w charmpy) → mniejsza waga 0.5; rzeźbi m → tylko A."),
    Unit("dedx_K", {"kness_K": 0.7}, 0.0, 0.7,
         "dE/dx K jak dla D0; w D+ ładunek K ≠ ładunek π więc błędne przypisanie rzadsze, ale protony (Λc) nadal odrzucane."),
    Unit("dedx_pi1", {"kness_pi1": 0.7}, 0.0, 0.3,
         "dE/dx π1: odrzuca p/K udające π; po 0.3 na każdy pion (łącznie 0.6 ≈ jeden K)."),
    Unit("dedx_pi2", {"kness_pi2": 0.7}, 0.0, 0.3,
         "dE/dx π2: j.w."),
    Unit("daughter_p", {"logp_pi_min": 1.2}, float(np.log(1.0)), 0.3,
         "ln min p_π: najmiększy pion w kombinatoryce D+ pochodzi z underlying event (p<1 GeV) → premia powyżej 1 GeV; tylko A."),
    Unit("eta", {"abs_eta": -1.5}, 1.8, 0.2,
         "|η| jak dla D0."),
]
DP_BIAS = -1.6  # 3-ciałowa kombinatoryka większa → surowszy prior

# ---------------------------------------------------------------- rewizja 1 (2026-11-04)
# Jedna runda po OBEJRZENIU kwantyli cech przygotowanych w MC (fd truth 0/1 vs mb truth≥8 w SR, podpróbki S6/S) – bez fitowania.
# Zapis zmian i liczb: revision_log.md. Źródło pseudo-istotności: adapt_v11.py (Lxy/0,20 mm; |d0|/σ_d0(pT)).
def _rev1(units, bias, changes, extra):
    out = []
    for u in units:
        c = changes.get(u.name, {})
        out.append(Unit(u.name, dict(u.w1), c.get("x0", u.x0), c.get("w2", u.w2), u.why + c.get("why", "")))
    return out + extra, bias

_D0_CH = {
    "lifetime": {"x0": 0.75, "why": " [rev1] preselekcja Lxy>0,3 mm: mediana tła 0,67, sygnału 0,82 (asinh) → próg 0,75 zamiast 0,48."},
    "pointing": {"x0": 2.1, "why": " [rev1] point_t mediany 2,53 (sig) / 1,78 (bkg) → próg 2,1 zamiast 1,5."},
    "flight_AND_point": {"x0": 0.75 + 0.6 * 2.1, "why": " [rev1] próg sumy przesunięty zgodnie z lifetime/pointing."},
    "vertex_chi2": {"x0": _ln1p(1.5), "why": " [rev1] χ² 2 torów: mediany 0,35/0,44 (log1p) – słaba separacja, próg obniżony do 1,5."},
    "ip_min": {"x0": 0.75, "why": " [rev1] pseudo-IPsig (|d0|/σ_d0(pT)): mediany 1,02/0,48 → próg 0,75."},
    "ip_max": {"x0": 1.25, "why": " [rev1] mediany 1,44/1,02 → próg 1,25."},
    "d0_abs": {"x0": 0.75, "why": " [rev1] min|d0| mediany 1,04/0,44 (asinh/0,05 mm) → próg 0,75."},
    "helicity": {"w2": 1.0, "why": " [rev1] |cosθ*| 0,37 vs 0,88 – najsilniejsza pojedyncza cecha w MC → w2 1,0 (nadal tylko A)."},
}
_DP_CH = {
    "lifetime": {"x0": 1.15, "why": " [rev1] preselekcja Lxy>0,5 mm: mediany 1,37/0,96 (asinh) → próg 1,15 zamiast 0,75."},
    "pointing": {"x0": 2.4, "why": " [rev1] mediany 2,96/1,82 → próg 2,4 zamiast 1,8."},
    "flight_AND_point": {"x0": 1.15 + 0.6 * 2.4, "why": " [rev1] próg sumy przesunięty zgodnie z lifetime/pointing."},
    "ip_max": {"x0": 1.6, "why": " [rev1] mediany 1,78/1,40 → próg 1,6."},
    "helicity": {"w2": 0.8, "why": " [rev1] |cosθ*| 0,48 vs 0,81 → w2 0,8."},
}
_D0_EXTRA = [Unit("D_ip", {"d_ipsig_t": -2.0}, 1.05, 0.5,
    "[rev1] IPsig kandydata D względem PV: prawdziwe D wskazuje na PV (log1p mediana 0,94), kombinatoryka 1,16 → kara powyżej 1,05. Niezależne od m (B)."),]
_DP_EXTRA = [Unit("D_ip", {"d_ipsig_t": -2.0}, 1.2, 0.6,
    "[rev1] jak D0; D⁺: 1,05 (sig) vs 1,39 (bkg) → próg 1,2, w2 0,6."),]
D0_REV1, D0_BIAS_REV1 = _rev1(D0_REV0, D0_BIAS, _D0_CH, _D0_EXTRA)
DP_REV1, DP_BIAS_REV1 = _rev1(DP_REV0, DP_BIAS, _DP_CH, _DP_EXTRA)

REGISTRY = {
    ("D0", 0): (D0_REV0, D0_BIAS),
    ("D0", 1): (D0_REV1, D0_BIAS_REV1),
    ("Dp", 0): (DP_REV0, DP_BIAS),
    ("Dp", 1): (DP_REV1, DP_BIAS_REV1),
}

def build(species: str = "D0", variant: str = "A", rev: int = 1):
    """-> (W1 [H,F], b1 [H], W2 [1,H], b2 [1], units)  float32; wariant B zeruje neurony z cechami rzeźbiącymi."""
    units, bias = REGISTRY[(species, rev)]
    if variant.upper() == "B":
        units = [u for u in units if all(f in MASS_SAFE for f in u.w1)]
    F, H = len(PREPARED), len(units)
    W1 = np.zeros((H, F), np.float32); b1 = np.zeros(H, np.float32); W2 = np.zeros((1, H), np.float32)
    for i, u in enumerate(units):
        for f, a in u.w1.items():
            W1[i, PREPARED.index(f)] = a
        b1[i] = u.b1
        W2[0, i] = u.w2
    b2 = np.array([bias], np.float32)
    return W1, b1, W2, b2, units

def n_nonzero(species="D0", variant="A", rev=1) -> int:
    W1, b1, W2, b2, _ = build(species, variant, rev)
    return int((W1 != 0).sum() + (b1 != 0).sum() + (W2 != 0).sum() + (b2 != 0).sum())
