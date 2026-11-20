#!/usr/bin/env python3
"""mult_truth_gen.py — MULT_charm: przewidywanie kształtu Pythii 8 na poziomie generatora (AFP_ddpred/gen/<tune>/gen_*.tsv.gz).

Format (Q_lxplus/truth/pythia_dd.cc): `E ev proc nch lead` (tylko zdarzenia z ≥1 D o pT>2, |η|<3), `D ev pdg pt eta phi fromB mode`;
`<plik>.json`: n_ok, proc, nch_hist (WSZYSTKIE zdarzenia, mieszanka ND+SD+DD+CD; n_ch = naładowane stabilne pT>0,5, |η|<2,5).
Liczy (inelastyczne, mieszanka jak w danych):
  Ŷ_s(x) = [N_s(bin)/N_ev(bin)] / [N_s/N_ev], x = n_ch/⟨n_ch⟩ (biny X_EDGES), s = D0, D⁺, D* (pT 3–30, |η|<2,5; D0 z D* liczony)
  pary D0+D⁺ (wszystkie / LS / OS) na D vs x; α (Ŷ ∝ x^α na [0,5; 4]); ND-only: licznik z proc=101 / mianownik nieznany → tylko udział ND.
Uwaga: n_ch zawiera produkty rozpadu D (odpowiednik estymatora `bl`, nie `blx`/`opp`).
  python3 mult_truth_gen.py --gen AFP_ddpred/gen --tunes A3 syst_MPIup syst_MPIdown syst_Monash --out out/truth_gen [--max-files N]
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mult_common as C  # noqa: E402

SP = {"D0": 421, "Dplus": 411, "Dstar": 413}


def one_tune(d, max_files):
    files = sorted(glob.glob(os.path.join(d, "gen_*.tsv.gz")))
    files = [f for f in files if os.path.exists(f + ".json")][: max_files or None]
    hist = np.zeros(301); n_ok = 0
    nb = len(C.X_EDGES) - 1
    acc = {}
    rows = []
    for f in files:
        js = json.load(open(f + ".json")); hist += np.asarray(js["nch_hist"], float); n_ok += js["n_ok"]
        df = pd.read_csv(f, sep="\t", header=None, comment="#", names=list(range(8)), dtype={0: str}, low_memory=False)
        E = df[df[0] == "E"]; D = df[df[0] == "D"]
        ev = pd.DataFrame({"ev": E[1].astype(np.int64).to_numpy(), "proc": E[2].astype(int).to_numpy(), "nch": E[3].astype(int).to_numpy()})
        d = pd.DataFrame({"ev": D[1].astype(np.int64).to_numpy(), "pdg": D[2].astype(int).to_numpy(),
                          "pt": D[3].astype(float).to_numpy(), "eta": D[4].astype(float).to_numpy()})
        d = d[(d.pt > C.PT_MIN) & (d.pt < C.PT_MAX) & (d.eta.abs() < C.ETA_MAX)]
        ap = d.pdg.abs()
        g = pd.DataFrame({"ev": d.ev, **{sp: (ap == p).astype(int) for sp, p in SP.items()},
                          "c": (((ap == 421) | (ap == 411)) & (d.pdg > 0)).astype(int),
                          "cb": (((ap == 421) | (ap == 411)) & (d.pdg < 0)).astype(int)}).groupby("ev").sum()
        m = ev.set_index("ev").join(g, how="inner")
        m["file"] = len(rows)
        rows.append(m)
    M = pd.concat(rows)
    nmean = float((hist * np.arange(301)).sum() / hist.sum())
    # mianownik: zdarzenia w binach x
    x_all = np.arange(301) / nmean
    b_all = C.x_bin(x_all)
    N_ev = np.bincount(b_all, hist, nb)
    X_ev = np.bincount(b_all, hist * x_all, nb) / np.where(N_ev > 0, N_ev, 1)
    b = C.x_bin(M.nch.to_numpy() / nmean)
    out = {"files": len(files), "n_events": int(n_ok), "nch_mean": nmean, "x": X_ev.tolist(), "N_ev": N_ev.tolist(), "species": {}}
    for sp in SP:
        N = np.bincount(b, M[sp].to_numpy(), nb); tot = N.sum() / N_ev.sum()
        Y = np.divide(N, N_ev, out=np.zeros(nb), where=N_ev > 0) / tot
        eY = np.divide(np.sqrt(N), N_ev, out=np.zeros(nb), where=N_ev > 0) / tot
        out["species"][sp] = {"N": N.tolist(), "Y": Y.tolist(), "eY": eY.tolist(), "alpha": C.powerlaw_fit(X_ev, Y, eY),
                              "frac_ND": float(M[M.proc == 101][sp].sum() / max(M[sp].sum(), 1))}
    c, cb = M.c.to_numpy().astype(float), M.cb.to_numpy().astype(float); n = c + cb
    P = {"N_D": np.bincount(b, n, nb), "pairs": np.bincount(b, n * (n - 1) / 2, nb),
         "LS": np.bincount(b, (c * (c - 1) + cb * (cb - 1)) / 2, nb), "OS": np.bincount(b, c * cb, nb)}
    out["pairs"] = {k: v.tolist() for k, v in P.items()}
    out["pairs_per_D"] = {k: np.divide(P[k], P["N_D"], out=np.full(nb, np.nan), where=P["N_D"] > 0).tolist() for k in ("pairs", "LS", "OS")}
    # jackknife po 10 grupach plików dla α (stat. prawdy)
    if len(files) >= 10:
        jk = {sp: [] for sp in SP}
        grp = M.file.to_numpy() % 10
        for k in range(10):
            s = grp != k
            for sp in SP:
                N = np.bincount(b[s], M[sp].to_numpy()[s], nb); tot = N.sum() / (N_ev.sum() * 0.9)
                Y = np.divide(N, N_ev * 0.9, out=np.zeros(nb), where=N_ev > 0) / tot
                jk[sp].append(C.powerlaw_fit(X_ev, Y, np.maximum(np.divide(np.sqrt(N), N_ev * 0.9, out=np.zeros(nb), where=N_ev > 0) / tot, 1e-9))[0])
        for sp in SP:
            a = np.asarray(jk[sp]); out["species"][sp]["alpha_jk_err"] = float(np.sqrt(9 / 10 * ((a - a.mean()) ** 2).sum()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "AFP_ddpred", "gen"))
    ap.add_argument("--tunes", nargs="+", default=["A3", "syst_MPIup", "syst_MPIdown", "syst_Monash"])
    ap.add_argument("--out", required=True); ap.add_argument("--max-files", type=int, default=None)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    R = {}
    lines = ["# Przewidywanie Pythii 8 (generator): Ŷ_s vs n_ch/⟨n_ch⟩, inelastyczne (ND+SD+DD+CD)", "",
             "n_ch: naładowane pT>0,5, |η|<2,5 (z produktami D); D: pT 3–30, |η|<2,5. α: Ŷ ∝ x^α na x∈[0,5;4] (σ_stat: jackknife 10 grup).", ""]
    for tn in a.tunes:
        d = os.path.join(a.gen, tn)
        if not glob.glob(os.path.join(d, "gen_*.tsv.gz")):
            print(f"pomijam {tn}: brak plików", file=sys.stderr); continue
        r = one_tune(d, a.max_files); R[tn] = r
        lines.append(f"## {tn}: {r['files']} plików, {r['n_events']:.3g} zdarzeń, ⟨n_ch⟩ = {r['nch_mean']:.2f}")
        lines.append("| j | x | " + " | ".join(f"Ŷ {sp}" for sp in SP) + " | pary/D | LS/D | OS/D |")
        lines.append("|" + "---|" * (5 + len(SP)))
        for j in range(len(C.X_EDGES) - 1):
            lines.append(f"| {j} | {r['x'][j]:.2f} | " + " | ".join(f"{r['species'][sp]['Y'][j]:.3f}±{r['species'][sp]['eY'][j]:.3f}" for sp in SP)
                         + " | " + " | ".join(f"{r['pairs_per_D'][k][j]:.4f}" for k in ("pairs", "LS", "OS")) + " |")
        lines.append("α: " + "; ".join(f"{sp} {r['species'][sp]['alpha'][0]:.3f} ± {r['species'][sp].get('alpha_jk_err', r['species'][sp]['alpha'][1]):.3f} "
                                       f"(udział ND {r['species'][sp]['frac_ND']:.3f})" for sp in SP))
        lines.append("")
    json.dump(R, open(os.path.join(a.out, "truth_gen.json"), "w"), indent=1, default=float)
    open(os.path.join(a.out, "TRUTH_GEN.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
