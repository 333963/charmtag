"""mult_common.py — MULT_charm (16.11.2026): wspólne definicje analizy „uzysk D i par vs krotność, dane vs MC”.

Wersja definicji: MULT-v1 (zamrożona w PROCEDURA.md; zmiana = odstępstwo opisane w raporcie).
Format wejścia: ntuple ROOT CharmAnalysis v2.0 (charm-prod-v2.0, Oct28) — drzewa `CharmAnalysis` + `Tracks` (zsynchronizowane wpisami).
Tryb zgodności v1.1 (Aug7): brak nSelBL/leadTrackPtBL/xi_cal_A/C/GRL/trigPS → zastępstwa opisane w `resolve_event_branches`.
Jednostki w ROOT: MeV (pędy, masy), mm (Lxy) — tu wszystko przeliczane na GeV.
"""
from __future__ import annotations

import numpy as np

DEF_VERSION = "MULT-v1"

# ----------------------------------------------------------------------------- masy [GeV]
M_D0, M_DP, M_DS, DM_DSTAR, M_PI = 1.86484, 1.86966, 2.01026, 0.145426, 0.13957

# ----------------------------------------------------------------------------- akceptacja i cięcia (MULT-v1, zamrożone)
PT_MIN, PT_MAX, ETA_MAX = 3.0, 30.0, 2.5
# wstępne (ekstrakcja; luźniejsze — pozwalają zaostrzyć w analizie bez ponownego czytania ROOT)
PRE = {
    "Dplus": dict(lxy_min=0.5, chi2_max=10.0, cosp_min=0.95, m_lo=1.70, m_hi=2.05),
    "D0": dict(lxy_min=0.1, chi2_max=10.0, m_lo=1.70, m_hi=2.05),
    "Dstar": dict(md0_win=0.060, dm_lo=0.1395, dm_hi=0.170),
}
# nominalne (analiza)
NOM = {
    "Dplus": dict(lxy_min=1.0, chi2_max=8.0, impsig_max=3.0, cosp_min=0.98, dstar_veto=0.160, m_lo=1.73, m_hi=2.00),
    "D0": dict(lxy_min=0.5, chi2_max=5.0, impsig_max=3.0, m_lo=1.75, m_hi=2.00),   # 1,75–2,00 = wspólne okno v1.1 i v2.0
    "Dstar": dict(md0_win=0.040, dm_lo=0.1395, dm_hi=0.165),
}
SPECIES = ("D0", "Dplus", "Dstar")
SP_CODE = {"D0": 10, "Dplus": 4, "Dstar": 1}           # = DMesons_decayType
PDG = {"D0": 421, "Dplus": 411, "Dstar": 413}

# ----------------------------------------------------------------------------- selekcja zdarzeń (MULT-v1)
EVSEL = dict(npv_eq=1, npu_eq=0, pvz_max=150.0, nbl_min=1)
# Schematy triggera (dane). Jeden łańcuch na schemat → waga zdarzenia = trigPS (LB), ta sama w liczniku i mianowniku.
# lead_min: cięcie offline na leadTrackPtBL [GeV] stosowane IDENTYCZNIE w danych i MC (plateau sptrk).
TRIG_SCHEMES = {
    "mbts": dict(chain="HLT_mb_mbts_L1MBTS_1", lead_min=0.0, role="nominał 427929 (physics_MinBias, PS≈1); w 428770 PS≈5·10⁴ → tylko informacyjnie"),
    "mbts2": dict(chain="HLT_noalg_mb_L1MBTS_2", lead_min=0.0, role="alternatywa MBTS (bez wymogu toru)"),
    "sptrk4": dict(chain="HLT_mb_sptrk_pt4_L1MBTS_2", lead_min=5.0, role="nominał 428770 (PS 25/50) + emulacja offline w MC"),
    "sptrk2": dict(chain="HLT_mb_sptrk_pt2_L1MBTS_2", lead_min=2.5, role="wariant 428770 (PS 625/1250)"),
    "j12": dict(chain="HLT_noalg_L1J12", lead_min=0.0, role="TYLKO diagnostyka obciążenia (strumień Main, dominujący trigger kandydatów)"),
    "all": dict(chain=None, lead_min=0.0, role="TYLKO diagnostyka: brak wymogu triggera = próbka ISO_sf/w2 (obciążona)"),
}
TRIG_STORE = sorted({s["chain"] for s in TRIG_SCHEMES.values() if s["chain"]} |
                    {"HLT_noalg_mb_L1MBTS_1", "HLT_noalg_L1AFP_A_OR_C", "HLT_j20_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1",
                     "HLT_noalg_L1RD0_FILLED", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"})
L1_TBP_STORE = ["L1_MBTS_1", "L1_MBTS_2", "L1_J12", "L1_AFP_A_OR_C", "L1_RD0_FILLED"]

# ----------------------------------------------------------------------------- estymatory krotności
# E_bl    : nSelBL (MinBias_BL: pT>0,5, |η|<2,5, |d0|<1,5 mm, …) — z córkami kandydata (jak ALICE „incl.”)
# E_blx   : nSelBL − (córki TEGO kandydata z pT>0,5, |η|<2,5) — licznik; mianownik zdarzeniowy = nSelBL (zdarzenia bez D)
# E_opp   : ślady drzewa Tracks pT>0,5, 0,5<|η|<2,5 w półkuli PRZECIWNEJ do D (minus córki kandydata w tym obszarze);
#           mianownik: każde zdarzenie wchodzi z wagą ½ do bin(n_A) i ½ do bin(n_C)
# E_trk500: EventInfo_nTracksPt500Eta25 (porównanie z gen_norm/w2 i ISO_sf)
ESTIMATORS = ("bl", "blx", "opp", "trk500")
ETA_GAP = 0.5
TRK_PT_MIN = 0.5
# biny względne x = n/⟨n⟩ (⟨n⟩ z wybranej próbki zdarzeń, osobno dane i MC — jak ALICE) i bezwzględne
X_EDGES = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 1e9])
N_EDGES = np.array([0, 5, 10, 15, 20, 30, 40, 60, 1e9])          # bezwzględne n (nSelBL), identyczne dane/MC
X_FIT = (0.5, 4.0)       # zakres dopasowania Ŷ ∝ x^α
# biny w2 (= gen_norm / etap1 w2d.py)
W2_NE = [-1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 1e4]
W2_PE = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 1e4]
W2_MINMC = 20


# ----------------------------------------------------------------------------- nazwy gałęzi
def res(keys, *names, required=True):
    """Pierwsza istniejąca gałąź spośród names, każda z sufiksem _NOSYS albo bez (v1.1 i v2.0 mają _NOSYS przy %SYS%)."""
    for name in names:
        for n in (name + "_NOSYS", name):
            if n in keys:
                return n
    if required:
        raise KeyError(f"brak gałęzi {names}")
    return None


def detect_version(keys) -> str:
    if res(keys, "EventInfo_nSelBL", required=False):
        return "v2.0"
    if res(keys, "EventInfo_xi_cal_plus", required=False):
        return "v1.1"
    return "unknown"


def resolve_event_branches(keys):
    """Mapa nazwa_logiczna → gałąź. Brak w v1.1: nselbl, leadbl, xiA/xiC, grl (→ None, zastępstwa w ekstrakcji:
    nselbl := nTracksPt500Eta25, leadbl := lead z Tracks |η|<2,5, grl := 1). Flaga `compat` w meta npz."""
    m = {
        "evnum": res(keys, "EventInfo_eventNumber"), "run": res(keys, "EventInfo_runNumber"),
        "lb": res(keys, "EventInfo_lumiBlock", required=False), "bcid": res(keys, "EventInfo_bcid", required=False),
        "npv": res(keys, "EventInfo_nPrimaryVertices"), "npu": res(keys, "EventInfo_nPileupVertices"),
        "pvz": res(keys, "EventInfo_pv_z"), "n500": res(keys, "EventInfo_nTracksPt500Eta25"),
        "nselbl": res(keys, "EventInfo_nSelBL", required=False), "leadbl": res(keys, "EventInfo_leadTrackPtBL", required=False),
        "xiA": res(keys, "EventInfo_xi_cal_A", "EventInfo_xi_cal_plus", required=False),
        "xiC": res(keys, "EventInfo_xi_cal_C", "EventInfo_xi_cal_minus", required=False),
        "grl": res(keys, "EventInfo_passGRL_MinBiasAFP", required=False),
        "mcch": res(keys, "EventInfo_mcChannelNumber", required=False),
    }
    return m


def resolve_trigger_branches(keys):
    """trigPassed_<chain> (v1.1 i v2.0: sufiks _NOSYS), trigPS_<chain> (tylko v2.0), trigTBP_<L1> (tylko v2.0)."""
    out = {}
    for c in TRIG_STORE:
        b = res(keys, f"EventInfo_trigPassed_{c}", required=False)
        if b:
            out[f"tp_{c}"] = b
        b = res(keys, f"EventInfo_trigPS_{c}", required=False)
        if b:
            out[f"ps_{c}"] = b
    for c in L1_TBP_STORE:
        b = res(keys, f"EventInfo_trigTBP_{c}", required=False)
        if b:
            out[f"tbp_{c}"] = b
    return out


CAND_BR = ["decayType", "pdgId", "pt", "eta", "phi", "m", "DeltaMass", "D0Index", "mKpi1", "mKpi2",
           "fitOutput__Lxy", "fitOutput__Chi2", "fitOutput__ImpactSignificance", "cosThetaP"]
CAND_DAU = ["pt", "eta", "trackId"]
CAND_MC = ["TruthMatchCategory", "truthPdgId"]
CAND_DAU_MC = ["pdgId", "truthPdgId"]


def resolve_cand_branches(keys, is_mc):
    m = {b: res(keys, f"DMesons_{b}") for b in CAND_BR}
    m.update({f"dau_{b}": res(keys, f"DMesons_daughterInfo__{b}") for b in CAND_DAU})
    if is_mc:
        m.update({b: res(keys, f"DMesons_{b}") for b in CAND_MC})
        m.update({f"dau_{b}": res(keys, f"DMesons_daughterInfo__{b}") for b in CAND_DAU_MC})
    return m


TRUTH_BR = {"pdg": "TruthDmesons_pdgId", "pt": "TruthDmesons_pt", "eta": "TruthDmesons_eta", "fromB": "TruthDmesons_fromBottom"}


# ----------------------------------------------------------------------------- binowanie / wagi
def x_bin(x):
    return np.clip(np.digitize(x, X_EDGES) - 1, 0, len(X_EDGES) - 2)


def n_bin(n):
    return np.clip(np.digitize(n, N_EDGES) - 1, 0, len(N_EDGES) - 2)


def w2_map(d_ntrk, d_lead, m_ntrk, m_lead, d_w=None, m_w=None):
    """w2 = p_dane/p_MC w binach gen_norm (komórki MC < 20 → 0). d_w: wagi zdarzeń danych (trigPS)."""
    hd, _, _ = np.histogram2d(d_ntrk, d_lead, [W2_NE, W2_PE], weights=d_w)
    hm_raw, _, _ = np.histogram2d(m_ntrk, m_lead, [W2_NE, W2_PE])
    hm, _, _ = np.histogram2d(m_ntrk, m_lead, [W2_NE, W2_PE], weights=m_w)
    pd_, pm = hd / max(hd.sum(), 1e-300), hm / max(hm.sum(), 1e-300)
    return np.divide(pd_, pm, out=np.zeros_like(pd_), where=hm_raw >= W2_MINMC)


def w2_lookup(W, ntrk, lead):
    i = np.clip(np.digitize(ntrk, W2_NE) - 1, 0, len(W2_NE) - 2)
    j = np.clip(np.digitize(lead, W2_PE) - 1, 0, len(W2_PE) - 2)
    return W[i, j]


def powerlaw_fit(x, y, ey, lo=X_FIT[0], hi=X_FIT[1]):
    """Ważone LS w log: ln Ŷ = ln a + α ln x na binach x∈[lo, hi] z y>0. Zwraca (α, σ_α, χ², ndf)."""
    x, y, ey = map(np.asarray, (x, y, ey))
    s = (x >= lo) & (x <= hi) & (y > 0) & (ey > 0) & np.isfinite(y)
    if s.sum() < 3:
        return np.nan, np.nan, np.nan, 0
    X = np.vstack([np.ones(s.sum()), np.log(x[s])]).T
    Y = np.log(y[s]); W = 1.0 / (ey[s] / y[s]) ** 2
    A = X.T @ (W[:, None] * X); b = X.T @ (W * Y)
    cov = np.linalg.inv(A); p = cov @ b
    chi2 = float(((Y - X @ p) ** 2 * W).sum())
    return float(p[1]), float(np.sqrt(cov[1, 1])), chi2, int(s.sum() - 2)
