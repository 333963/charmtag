#!/usr/bin/env python3
"""UB_TOPO (19.11.2026, projekt A): hipotezy kształtu |Δη| (i Δφ) par D–D z Pythii 8 A3 na poziomie generatora (lokalnie; wynik → ref/topo_gen.json).

Pary (D0/D⁺ ze stanem końcowym kanału analizy, pT 3–30 GeV, |η| < 2,5; jak MULT/mult_truth_gen):
  pythia  — pary z TEGO SAMEGO zdarzenia (SPS + DPS w proporcji Pythii; = oczekiwanie MC nominału),
  dps     — pary NIEZALEŻNE: nogi z różnych zdarzeń (losowe parowanie, ten sam skład kanałów) = granica czystego DPS,
  ls / os — pary like/opposite-sign z tego samego zdarzenia (LS w Pythii ≈ DPS, kontrola).
Waga pary = ε(nogi 1)·ε(nogi 2), ε(pT, |η|) = eps_matched z GEN_mcstat effmap_direct.npz (marginalizacja po n_trk wagą mianownika);
stały czynnik ε znosi się w udziałach binów, liczy się tylko zależność od (pT, η). Wynik: udziały binów per kanał i hipoteza + błąd (bootstrap plików).
  topo_gen.py --gen 'AFP_ddpred/gen/A3/gen_*.tsv.gz' --eff effmap_direct.npz --out ref/topo_gen.json [--max-files N]"""
import argparse
import glob
import gzip
import json

import numpy as np

EDGES = {"dy": [0.0, 1.0, 5.0], "dphi": [0.0, np.pi / 2, np.pi], "dphi_u": [0.0, np.pi / 4, 3 * np.pi / 4, np.pi]}
CH = {(421, 421): "D0D0", (411, 421): "DpD0", (411, 411): "DpDp"}
SP = {421: "D0", 411: "Dplus"}


def read_gen(files):
    """→ tablica (plik, ev, pdg, pt, eta, phi) nóg w akceptancji."""
    rows = []
    for fi, f in enumerate(files):
        with gzip.open(f, "rt") as h:
            for line in h:
                if line[0] != "D":
                    continue
                _, ev, pdg, pt, eta, phi, _fb, mode = line.rstrip("\n").split("\t")
                pdg = int(pdg); pt = float(pt); eta = float(eta)
                if abs(pdg) in SP and int(mode) == 1 and 3.0 <= pt < 30.0 and abs(eta) < 2.5:
                    rows.append((fi, int(ev), pdg, pt, eta, float(phi)))
    return np.array(rows, float).reshape(-1, 6)


def eff_fn(path):
    z = np.load(path)
    pe, ee = z["pt_edges"], z["eta_edges"]
    E = {}
    for sp in ("D0", "Dplus"):
        num, den = z[f"{sp}__num_matched"].sum(2), z[f"{sp}__den"].sum(2)
        E[sp] = np.where(den > 0, num / np.maximum(den, 1), 0.0)

    def f(sp, pt, eta):
        i = np.clip(np.digitize(pt, pe) - 1, 0, len(pe) - 2); j = np.clip(np.digitize(np.abs(eta), ee) - 1, 0, len(ee) - 2)
        return E[sp][i, j]
    return f


def leg_eff(a, eff):
    sp = np.where(np.abs(a[:, 2]) == 421, "D0", "Dplus")
    out = np.empty(len(a))
    for s in ("D0", "Dplus"):
        m = sp == s
        out[m] = eff(s, a[m, 3], a[m, 4]) if m.any() else 0
    return out


def pair_vars(a, i, j):
    dphi = np.abs((a[i, 5] - a[j, 5] + np.pi) % (2 * np.pi) - np.pi)
    dy = np.abs(a[i, 4] - a[j, 4])
    ls = np.sign(a[i, 2]) == np.sign(a[j, 2])
    lo = np.minimum(np.abs(a[i, 2]), np.abs(a[j, 2])).astype(int); hi = np.maximum(np.abs(a[i, 2]), np.abs(a[j, 2])).astype(int)
    ch = np.array([CH[(x, y)] for x, y in zip(np.where(lo == 411, 411, 421), np.where(hi == 421, 421, 411))]) if len(i) else np.array([])
    return dphi, dy, ls, ch


def same_event_pairs(a):
    key = a[:, 0] * 1e10 + a[:, 1]
    o = np.argsort(key, kind="stable"); a = a[o]; key = key[o]
    st = np.r_[0, np.flatnonzero(np.diff(key)) + 1, len(a)]
    I, J = [], []
    for s, e in zip(st[:-1], st[1:]):
        for i in range(s, e):
            for j in range(i + 1, e):
                I.append(i); J.append(j)
    return a, np.array(I, int), np.array(J, int)


def indep_pairs(a, Ich, n_per, rng):
    """Pary niezależne o tym samym składzie kanałów co pary jednego zdarzenia: dla kanału losuje nogi odpowiednich gatunków z RÓŻNYCH zdarzeń."""
    d0 = np.flatnonzero(np.abs(a[:, 2]) == 421); dp = np.flatnonzero(np.abs(a[:, 2]) == 411)
    I, J = [], []
    for ch, n in n_per.items():
        A, B = {"D0D0": (d0, d0), "DpD0": (dp, d0), "DpDp": (dp, dp)}[ch]
        i = rng.choice(A, n); j = rng.choice(B, n)
        ok = (a[i, 0] != a[j, 0]) | (a[i, 1] != a[j, 1])
        I.append(i[ok]); J.append(j[ok])
    return np.concatenate(I), np.concatenate(J)


def fractions(vals, w, edges):
    h = np.histogram(vals, edges, weights=w)[0]
    return h / h.sum() if h.sum() > 0 else np.full(len(edges) - 1, np.nan)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", required=True); ap.add_argument("--eff", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--max-files", type=int, default=None); ap.add_argument("--nboot", type=int, default=50); ap.add_argument("--seed", type=int, default=19)
    a_ = ap.parse_args(argv)
    files = sorted(glob.glob(a_.gen))[: a_.max_files]
    rng = np.random.default_rng(a_.seed)
    a = read_gen(files)
    eff = eff_fn(a_.eff)
    a, I, J = same_event_pairs(a)
    w_leg = leg_eff(a, eff)
    dphi, dy, ls, ch = pair_vars(a, I, J)
    n_per = {c: int((ch == c).sum()) * 20 for c in ("D0D0", "DpD0", "DpDp")}
    Ii, Ji = indep_pairs(a, ch, n_per, rng)
    dphi_i, dy_i, ls_i, ch_i = pair_vars(a, Ii, Ji)
    W = {"w": w_leg[I] * w_leg[J], "u": np.ones(len(I))}
    Wi = {"w": w_leg[Ii] * w_leg[Ji], "u": np.ones(len(Ii))}
    fidx = a[I, 0].astype(int); nf = len(files)
    res = {"files": nf, "pairs_same_event": int(len(I)), "pairs_indep": int(len(Ii)), "edges": {k: list(map(float, v)) for k, v in EDGES.items()},
           "eff": a_.eff, "note": "udziały binów; w = ważone ε(pT,|η|) nóg, u = bez wag; błąd = bootstrap plików (pythia/ls/os)"}
    var_of = {"dy": (dy, dy_i), "dphi": (dphi, dphi_i), "dphi_u": (dphi, dphi_i)}
    for c in ("all", "D0D0", "DpD0", "DpDp"):
        m = np.ones(len(I), bool) if c == "all" else ch == c
        mi = np.ones(len(Ii), bool) if c == "all" else ch_i == c
        R = {"n_pairs": int(m.sum()), "f_LS": {k: float(np.average(ls[m], weights=W[k][m])) for k in W}}
        for v, (x, xi) in var_of.items():
            e = EDGES[v]
            for k in ("w", "u"):
                hyp = {"pythia": fractions(x[m], W[k][m], e), "ls": fractions(x[m & ls], W[k][m & ls], e),
                       "os": fractions(x[m & ~ls], W[k][m & ~ls], e), "dps": fractions(xi[mi], Wi[k][mi], e)}
                boot = []
                for _ in range(a_.nboot):
                    cnt = np.bincount(rng.integers(0, nf, nf), minlength=nf)
                    bw = cnt[fidx] * W[k]
                    boot.append(fractions(x[m], bw[m], e))
                R[f"{v}_{k}"] = {h: [float(t) for t in f] for h, f in hyp.items()}
                R[f"{v}_{k}"]["pythia_err"] = [float(t) for t in np.nanstd(np.array(boot), 0)]
        res[c] = R
    json.dump(res, open(a_.out, "w"), indent=1)
    print(json.dumps({c: {"n": res[c]["n_pairs"], "dy_w": res[c]["dy_w"], "dphi_u_w": res[c]["dphi_u_w"]} for c in ("all", "D0D0", "DpD0", "DpDp")}, indent=1))


if __name__ == "__main__":
    main()
