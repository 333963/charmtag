#!/usr/bin/env python3
"""PAIR_PLANES (20.11.2026, projekt A): płaszczyzny mas par D–D z h5 v2.0 (lub v1.2) → histogramy npz.

Ulepszona wersja makra look_for_doubleD.C / plots_doubleD *_margins.png z kodu bazowego charmpp:
  * pary tylko w obrębie zdarzenia, BEZ par dzielących tor (Dmeson_daughter{0,1,2}_trackId, dokładnie);
  * orientacja nóg: różne gatunki — stała (kanał X×Y: x = X); te same — x = wyższe pT;
  * podział OS / LS wg znaku pdgId (zapach powabu: cc̄ vs cc);
  * wariant `vetoDs`: noga D⁺ odrzucona, gdy m(Kππ) − m(Kπ)_i < 160 MeV dla któregoś Kπ (weto D*, jak W+D `MindMDP`;
    w trybie minbias CharmAnalysis go nie ma → refleksja D*⁺ przy m(Kππ) ≈ 2010 MeV);
  * poziomy: `reco` (wszystkie), `matched` (MC: obie nogi Dmeson_truth_category ∈ {0 MatchedC, 1 MatchedB});
  * dla par w pudełku pik×pik i w pudełku pasmo×pasmo: rozkłady |Δφ| i Δy (zamiast bezużytecznej „prawdy” z masą PDG).
Bez score klasyfikatora (h5 go nie ma): to diagnostyka na poziomie findera, jak oryginał. Dane: tylko passGRL_MinBiasAFP == 1
(gdy kolumna jest). MC: waga event_weight (fd: BEZ wag BR wymuszonych rozpadów — tylko kształty).

Kandydaci: D⁺ (decayType 4) m ∈ [1,70; 2,20] GeV; D⁰ (10) m ∈ [1,70; 2,05], L_xy ≥ 0,25 mm (jak analiza A); D* (1) Δm ∈ [0,138; 0,170].

Użycie: pp_extract.py --sample S --files F1 [F2 ...] --out OUT.npz   |   --h5dir DIR --part K --nparts N --out OUT.npz [--wd WD]
"""
import argparse, glob, json, os, subprocess, sys, time
import numpy as np
import h5py

SP = {"Dp": 4, "Ds": 1, "D0": 10}
# gatunek: (kolumna, skala→MeV, zakres MeV, liczba binów, okno piku MeV, PDG MeV)
AXD = {"Dp": ("Dmeson_mass", 1e3, (1700.0, 2200.0), 100, (1840.0, 1900.0), 1869.66),
       "D0": ("Dmeson_mass", 1e3, (1700.0, 2050.0), 70, (1835.0, 1895.0), 1864.84),
       "Ds": ("Dmeson_delta_mass", 1e3, (138.0, 170.0), 64, (143.5, 147.5), 145.43)}
CH = [("Dp", "Dp"), ("D0", "Dp"), ("D0", "D0"), ("Dp", "Ds"), ("D0", "Ds"), ("Ds", "Ds")]
RELS, VARS = ("OS", "LS"), ("nom", "vetoDs")
DPHI = np.linspace(0, np.pi, 33)
DY = np.linspace(-4, 4, 41)
PDGM = {4: 1.86966, 1: 2.01026, 10: 1.86484}
VETO_MEV = 160.0
COLS = ["Dmeson_decayType", "Dmeson_mass", "Dmeson_delta_mass", "Dmeson_pt", "Dmeson_eta", "Dmeson_phi", "Dmeson_charge",
        "Dmeson_Lxy", "Dmeson_mKpi1", "Dmeson_mKpi2", "Dmeson_daughter0_trackId", "Dmeson_daughter1_trackId", "Dmeson_daughter2_trackId"]
OPT = ["event_uid", "event_number", "Dmeson_truth_category", "event_weight", "passGRL_MinBiasAFP"]
XRD = os.environ.get("XRD", "root://eosuser.cern.ch")
MAXPAIRS = 5_000_000


def keys(level):
    return [f"{x}_{y}|{r}|{v}|{level}" for x, y in CH for r in RELS for v in VARS]


def empty(is_mc):
    H = {}
    for lev in (("reco", "matched") if is_mc else ("reco",)):
        for k in keys(lev):
            x, y = k.split("|")[0].split("_")
            H[k + "|h2"] = np.zeros((AXD[x][3], AXD[y][3]))
            for reg in ("box", "sb"):
                H[k + f"|dphi_{reg}"] = np.zeros(len(DPHI) - 1)
                H[k + f"|dy_{reg}"] = np.zeros(len(DY) - 1)
            H[k + "|n_shared"] = np.zeros(1)
    return H


def local_copy(path, wd):
    if wd and path.startswith("/eos/"):
        dst = os.path.join(wd, os.path.basename(path))
        try:
            if subprocess.run(["xrdcp", "-s", "-f", f"{XRD}/{path}", dst], timeout=900).returncode == 0:
                return dst, True
        except Exception:
            pass
        if os.path.exists(dst):
            os.remove(dst)
    return path, False


def load(fn):
    with h5py.File(fn, "r") as f:
        d = f["events"]
        names = d.dtype.names
        miss = [c for c in COLS if c not in names]
        if miss:
            sys.exit(f"{fn}: brak kolumn {miss}")
        use = COLS + [c for c in OPT if c in names]
        out = {c: np.empty(d.shape[0], d.dtype[c]) for c in use}
        step = 500_000
        for i in range(0, d.shape[0], step):
            blk = d.fields(use)[i:i + step]
            for c in use:
                out[c][i:i + step] = blk[c]
    return out


def candidates(c, is_mc):
    dt, m, dm, lxy = c["Dmeson_decayType"], c["Dmeson_mass"], c["Dmeson_delta_mass"], c["Dmeson_Lxy"]
    sel = ((dt == 4) & (m >= 1.70) & (m <= 2.20)) | ((dt == 10) & (m >= 1.70) & (m <= 2.05) & (lxy >= 0.25)) \
        | ((dt == 1) & (dm >= 0.138) & (dm <= 0.170))
    if "passGRL_MinBiasAFP" in c and not is_mc:
        sel &= c["passGRL_MinBiasAFP"] == 1
    ev = c["event_uid"] if "event_uid" in c else c["event_number"]
    idx = np.nonzero(sel)[0]
    idx = idx[np.argsort(ev[idx], kind="stable")]
    g = {k: v[idx] for k, v in c.items()}
    g["ev"] = ev[idx]
    mev = np.where(g["Dmeson_decayType"] == 1, g["Dmeson_delta_mass"], g["Dmeson_mass"]) * 1e3
    g["mev"] = mev.astype(np.float64)
    m_kpipi = g["Dmeson_mass"] * 1e3
    k1, k2 = g["Dmeson_mKpi1"], g["Dmeson_mKpi2"]
    bad = ((k1 > 0) & (m_kpipi - k1 < VETO_MEV)) | ((k2 > 0) & (m_kpipi - k2 < VETO_MEV))
    g["vetoOK"] = ~((g["Dmeson_decayType"] == 4) & bad)
    g["trk"] = np.stack([g[f"Dmeson_daughter{i}_trackId"].astype(np.int64) for i in range(3)], 1)
    pt, eta = g["Dmeson_pt"].astype(np.float64), g["Dmeson_eta"].astype(np.float64)
    mass = np.vectorize(PDGM.get)(g["Dmeson_decayType"]).astype(np.float64)
    pz = pt * np.sinh(eta)
    e = np.sqrt(mass ** 2 + pt ** 2 + pz ** 2)
    g["y"] = 0.5 * np.log((e + pz) / (e - pz))
    g["matched"] = np.isin(g["Dmeson_truth_category"], (0, 1)) if is_mc and "Dmeson_truth_category" in g else np.zeros(len(idx), bool)
    g["w"] = g["event_weight"].astype(np.float64) if is_mc and "event_weight" in g else np.ones(len(idx))
    return g


def pair_chunks(ev):
    """Zwraca kolejne (I, J) indeksów par w zdarzeniu, porcjami ≤ MAXPAIRS."""
    if len(ev) < 2:
        return
    starts = np.r_[0, np.nonzero(ev[1:] != ev[:-1])[0] + 1]
    ends = np.r_[starts[1:], len(ev)]
    cnt = ends - starts
    npair = cnt * (cnt - 1) // 2
    keep = cnt >= 2
    starts, ends, npair = starts[keep], ends[keep], npair[keep]
    cum = np.cumsum(npair)
    lo = 0
    while lo < len(starts):
        hi = int(np.searchsorted(cum, (cum[lo - 1] if lo else 0) + MAXPAIRS, side="right"))
        hi = max(hi, lo + 1)
        n_ev = ends[lo:hi] - starts[lo:hi]
        # kandydaci zdarzeń z ≥ 2 kandydatami (zdarzenia z 1 pominięte) i koniec ich zdarzenia
        i = np.repeat(starts[lo:hi], n_ev) + (np.arange(int(n_ev.sum())) - np.repeat(np.cumsum(n_ev) - n_ev, n_ev))
        evend = np.repeat(ends[lo:hi], n_ev)
        c = evend - i - 1
        tot = int(c.sum())
        if tot:
            I = np.repeat(i, c)
            off = np.arange(tot) - np.repeat(np.cumsum(c) - c, c)
            J = I + 1 + off
            yield I, J
        lo = hi


def fill(H, g, is_mc):
    s = g["Dmeson_decayType"]
    nshared_tot = 0
    for I, J in pair_chunks(g["ev"]):
        ti, tj = g["trk"][I], g["trk"][J]
        shared = np.zeros(len(I), bool)
        for a in range(3):
            for b in range(3):
                shared |= (ti[:, a] == tj[:, b]) & (ti[:, a] >= 0)
        for X, Y in CH:
            if X == Y:
                msk = (s[I] == SP[X]) & (s[J] == SP[X])
                lead = g["Dmeson_pt"][I] >= g["Dmeson_pt"][J]
                xi, yi = np.where(lead, I, J)[msk], np.where(lead, J, I)[msk]
            else:
                m1 = (s[I] == SP[X]) & (s[J] == SP[Y])
                m2 = (s[I] == SP[Y]) & (s[J] == SP[X])
                msk = m1 | m2
                xi, yi = np.where(m1, I, J)[msk], np.where(m1, J, I)[msk]
            if not len(xi):
                continue
            sh = shared[msk]
            rel = np.where(g["Dmeson_charge"][xi] * g["Dmeson_charge"][yi] < 0, "OS", "LS")
            vok = g["vetoOK"][xi] & g["vetoOK"][yi]
            w = g["w"][xi]
            mx, my = g["mev"][xi], g["mev"][yi]
            ax, ay = AXD[X], AXD[Y]
            inx = (mx >= ax[4][0]) & (mx <= ax[4][1]); iny = (my >= ay[4][0]) & (my <= ay[4][1])
            box, sb = inx & iny, ~inx & ~iny
            dphi = np.abs(np.angle(np.exp(1j * (g["Dmeson_phi"][xi].astype(np.float64) - g["Dmeson_phi"][yi]))))
            dy = g["y"][xi] - g["y"][yi]
            levels = [("reco", np.ones(len(xi), bool))]
            if is_mc:
                levels.append(("matched", g["matched"][xi] & g["matched"][yi]))
            for lev, lm in levels:
                for r in RELS:
                    for v in VARS:
                        k = f"{X}_{Y}|{r}|{v}|{lev}"
                        base = lm & (rel == r) & (vok if v == "vetoDs" else True)
                        H[k + "|n_shared"][0] += float((base & sh).sum())
                        ok = base & ~sh
                        H[k + "|h2"] += np.histogram2d(mx[ok], my[ok], bins=(ax[3], ay[3]), range=(ax[2], ay[2]), weights=w[ok])[0]
                        for reg, rm in (("box", box), ("sb", sb)):
                            q = ok & rm
                            H[k + f"|dphi_{reg}"] += np.histogram(dphi[q], bins=DPHI, weights=w[q])[0]
                            H[k + f"|dy_{reg}"] += np.histogram(dy[q], bins=DY, weights=w[q])[0]
        nshared_tot += int(shared.sum())
    return nshared_tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--files", nargs="*")
    ap.add_argument("--h5dir"); ap.add_argument("--part", type=int, default=0); ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--max-files", type=int, default=0); ap.add_argument("--wd", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.files:
        mine = a.files
    else:
        files = sorted(glob.glob(os.path.join(a.h5dir, f"{a.sample}_*_none.h5")))
        if not files:
            sys.exit(f"brak plików {a.h5dir}/{a.sample}_*_none.h5")
        mine = files[a.part::a.nparts]
    if a.max_files:
        mine = mine[:a.max_files]
    is_mc = not a.sample.startswith("data")
    H = empty(is_mc)
    meta = {"sample": a.sample, "part": a.part, "nparts": a.nparts, "is_mc": is_mc, "files": [], "n_cand_sel": 0,
            "n_pairs_shared_any": 0, "axes": {k: [v[2][0], v[2][1], v[3], v[4][0], v[4][1], v[5]] for k, v in AXD.items()},
            "veto_MeV": VETO_MEV, "event_key": None}
    t0 = time.time()
    for fn in mine:
        loc, tmp = local_copy(fn, a.wd)
        try:
            c = load(loc)
        finally:
            if tmp:
                os.remove(loc)
        meta["event_key"] = "event_uid" if "event_uid" in c else "event_number"
        g = candidates(c, is_mc)
        meta["n_cand_sel"] += int(len(g["ev"]))
        meta["n_pairs_shared_any"] += fill(H, g, is_mc)
        meta["files"].append(os.path.basename(fn))
        print(f"{os.path.basename(fn)}: kandydaci {len(g['ev'])} ({time.time() - t0:.0f} s)", flush=True)
    meta["seconds"] = round(time.time() - t0, 1)
    np.savez_compressed(a.out, **H)
    json.dump(meta, open(a.out.replace(".npz", ".meta.json"), "w"), indent=1)
    print(f"OK {a.sample} część {a.part}/{a.nparts}: {len(mine)} plików, {meta['n_cand_sel']} kandydatów, {meta['seconds']} s")


if __name__ == "__main__":
    main()
