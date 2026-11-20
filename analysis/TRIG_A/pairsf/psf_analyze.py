#!/usr/bin/env python3
"""TRIG_A/pairsf — SF klasyfikatora H dla D0 w KATEGORIACH SCORE używanych przez pary D–D, per klasa triggera i run (projekt A).

Pojedyncze D0 z tagu D* (fit Δm pass/fail, charmtag.dstar_refit.fit_passfail — ten sam kod co B_trig/EV). ŻADNYCH par D–D.
Wejście: wyniki etapu b_trig (bt_<próbka>.npz lub części parts/bt_<próbka>_*.npz; score = nominalny model D0_H, L_xy ≥ 0,25 mm).
Selekcje („pass”) = przedziały score z progów par (pair_inputs.json → meta.thresholds["10"], globalne w pT):
  c0: ε_S∈[0,30) %   c1: [30,45)   c2: [45,65)   c4: [65,90) (tylko z tagiem D*)   ge65: ≥WP65 (kat. 3 z tagiem)   ge90: ≥WP90 (wejście nóg D0)
SF = ε_dane/ε_MC(fd_nd, prawda), mianownik: wszystkie D0←D* w |m(Kπ) − m0| < 55 MeV. Klasy rozłączne: J12only / MBTSbased / other (428770),
MB427929. Wyjście: --out → psf.json, PSF.md.
"""
import argparse, glob, json, os, sys
from multiprocessing import Pool
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import btrig_common as C  # noqa: E402

THR = {"0.3": 0.9839380383491516, "0.45": 0.9589824676513672, "0.65": 0.8569592237472534, "0.9": 0.38729336857795715}
SEL = {"c0": (THR["0.3"], 9.0), "c1": (THR["0.45"], THR["0.3"]), "c2": (THR["0.65"], THR["0.45"]), "c4": (THR["0.9"], THR["0.65"]),
       "ge65": (THR["0.65"], 9.0), "ge90": (THR["0.9"], 9.0)}
BINS = {"3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0), "5-30": (5.0, 30.0), "3-30": (3.0, 30.0)}
CLASSES = {428770: ["all", "J12only", "MBTSbased", "other"], 427929: ["MB427929"]}
DISJ = [(428770, "J12only"), (428770, "MBTSbased"), (428770, "other"), (427929, "MB427929")]
NEED = ("run", "event", "m", "pt", "dm", "score")


def load(paths, mc=False):
    cols = {}
    for p in paths:
        z = np.load(p, allow_pickle=False)
        keys = [k for k in z.files if k in NEED or k.startswith("e_tp_") or k == "e_grl" or (mc and k in ("truth", "truth_ds"))]
        for k in keys:
            cols.setdefault(k, []).append(z[k])
    n = {k: len(v) for k, v in cols.items()}
    if len(set(n.values())) != 1:
        raise SystemExit(f"STOP: niespójne kolumny w częściach: {n}")
    return {k: np.concatenate(v) for k, v in cols.items()}


def _fit(job):
    key, dp, df_, shape = job
    return key, C.fit_eps(dp, df_, shape, interval=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", action="append", required=True, help="GLOB@RUN (powtarzalne)")
    ap.add_argument("--mc", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--nproc", type=int, default=8); ap.add_argument("--max-rows", type=int, default=0, help="tylko testy")
    ap.add_argument("--bins", default="", help="tylko testy: lista binów"); ap.add_argument("--sels", default="", help="tylko testy: lista selekcji")
    a = ap.parse_args()
    if a.bins:
        for k in [k for k in BINS if k not in a.bins.split(",")]:
            del BINS[k]
    if a.sels:
        for k in [k for k in SEL if k not in a.sels.split(",")]:
            del SEL[k]
    from charmtag import dstar_refit as DR
    M = load(sorted(glob.glob(a.mc)), mc=True)
    f_sig = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < C.NOM_WIN)
    shapes = {b: DR.shape_from_mc(M["dm"][f_sig & (M["pt"] >= lo) & (M["pt"] < hi)] * 1e3) for b, (lo, hi) in BINS.items()}
    eps_mc = {}
    for b, (lo, hi) in BINS.items():
        s = M["score"][f_sig & (M["pt"] >= lo) & (M["pt"] < hi)]
        for sn, (t0, t1) in SEL.items():
            p = (s >= t0) & (s < t1)
            eps_mc[(b, sn)] = (float(p.mean()), float(np.sqrt(p.mean() * (1 - p.mean()) / max(len(s), 1))), int(len(s)))
    jobs, meta = [], {}
    for spec in a.data:
        g, run = spec.rsplit("@", 1); run = int(run)
        fl = sorted(glob.glob(g))
        if not fl:
            raise SystemExit(f"STOP: brak plików {g}")
        D = load(fl)
        if a.max_rows:
            D = {k: v[:a.max_rows] for k, v in D.items()}
        bits = {k[2:]: D[k] for k in D if k.startswith("e_tp_") or k == "e_grl"}
        cls = C.classify(bits, run)
        nom = np.abs(D["m"] - C.M0_D0) < C.NOM_WIN
        meta[str(run)] = {"files": len(fl), "rows": int(len(D["pt"])), "missing_chains": cls.get("_missing", []),
                          "frac": {c: float(cls[c].mean()) for c in CLASSES[run]}}
        for c in CLASSES[run]:
            for b, (lo, hi) in BINS.items():
                md = cls[c] & nom & (D["pt"] >= lo) & (D["pt"] < hi)
                for sn, (t0, t1) in SEL.items():
                    pas = (D["score"] >= t0) & (D["score"] < t1)
                    jobs.append(((run, c, b, sn), D["dm"][md & pas], D["dm"][md & ~pas], shapes[b]))
    with Pool(a.nproc) as P:
        fits = dict(P.map(_fit, jobs, chunksize=1))
    rows = []
    for (run, c, b, sn), f in fits.items():
        e, ee, n = eps_mc[(b, sn)]
        err = C.eps_sym_err(f)
        sf = f["eps"] / e if e > 0 and np.isfinite(f["eps"]) else np.nan
        rows.append({"run": run, "cls": c, "bin": b, "sel": sn, "eps_data": f["eps"], "eps_data_err": err, "determined": f["determined"],
                     "valid": f["valid"], "N_sig": f["N"], "N_sig_err": f["N_err"], "n_pass": f["n_pass"], "n_fail": f["n_fail"],
                     "eps_mc": e, "eps_mc_err": ee, "n_mc": n, "SF": sf,
                     "SF_err": float(np.hypot(err / e, sf * ee / e)) if e > 0 and np.isfinite(sf) else np.nan})
    R = {(r["run"], r["cls"], r["bin"], r["sel"]): r for r in rows}
    ok = lambda r: r is not None and r["determined"] and np.isfinite(r["SF"]) and np.isfinite(r["SF_err"]) and r["SF_err"] > 0  # noqa: E731
    tests = []
    for b in BINS:
        for sn in SEL:
            rs = [R.get((run, c, b, sn)) for run, c in DISJ]
            use = [r for r in rs if ok(r)]
            t = {"bin": b, "sel": sn, "n_classes": len(use)}
            if len(use) >= 2:
                x = C.chi2_consistency([r["SF"] for r in use], [r["SF_err"] for r in use])
                t.update({"chi2": float(x["chi2"]), "ndf": int(x["ndf"]), "p": float(x["p"]), "mean": float(x["mean"]), "mean_err": float(x["err"])})
            j, mb, al = R.get((428770, "J12only", b, sn)), R.get((427929, "MB427929", b, sn)), R.get((428770, "all", b, sn))
            if ok(j) and ok(mb):
                t["J12_minus_MB427929_sigma"] = float((j["SF"] - mb["SF"]) / np.hypot(j["SF_err"], mb["SF_err"]))
                t["J12_over_MB427929"] = float(j["SF"] / mb["SF"])
            # spójność wewnętrzna: Σ N_sig klas rozłącznych 428770 vs „all”
            sub = [R.get((428770, c, b, sn)) for c in ("J12only", "MBTSbased", "other")]
            if al is not None and all(s is not None for s in sub):
                t["Nsig_sum_over_all"] = float(sum(s["N_sig"] for s in sub) / al["N_sig"]) if al["N_sig"] else float("nan")
                if all(ok(s) for s in sub) and ok(al):
                    w = np.array([s["N_sig"] for s in sub]); t["SF_weighted_minus_all"] = float((w * [s["SF"] for s in sub]).sum() / w.sum() - al["SF"])
            tests.append(t)
    os.makedirs(a.out, exist_ok=True)
    json.dump({"thresholds": THR, "sel": SEL, "meta": meta, "eps_mc": {f"{b}:{s}": v for (b, s), v in eps_mc.items()}, "rows": rows, "tests": tests},
              open(os.path.join(a.out, "psf.json"), "w"), indent=1, default=float)
    L = ["# TRIG_A/pairsf — SF klasyfikatora H (D⁰ z tagu D*) w kategoriach score par, per klasa triggera", "",
         "Selekcje: c0/c1/c2 = przedziały ε_S [0,30)/[30,45)/[45,65) %; c4 = [65,90) %; ge65, ge90 = skumulowane. Progi globalne z `pair_inputs.json`.",
         "SF = ε_dane/ε_MC (fd_nd). „(n)” = fit nieokreślony. Pojedyncze D, bez par.", ""]
    for b in BINS:
        L += [f"## pT {b} GeV", "", "| run | klasa | N_sig | " + " | ".join(SEL) + " |", "|---|---|---|" + "---|" * len(SEL)]
        for run, cl in CLASSES.items():
            for c in cl:
                cells = []
                for sn in SEL:
                    r = R[(run, c, b, sn)]
                    cells.append(f"{r['SF']:.2f} ± {r['SF_err']:.2f}" + ("" if r["determined"] else " (n)"))
                L.append(f"| {run} | {c} | {R[(run, c, b, list(SEL)[-1])]['N_sig']:.0f} | " + " | ".join(cells) + " |")
        L += ["", "| selekcja | ε_MC | χ²/ndf klas rozłącznych (p) | J12only − MB427929 [σ] | J12only/MB427929 | ΣN_sig klas / all | SF ważone klas − SF all |", "|---|---|---|---|---|---|---|"]
        for t in tests:
            if t["bin"] != b:
                continue
            g = lambda k, f="{:.2f}": f.format(t[k]) if k in t and np.isfinite(t[k]) else "–"  # noqa: E731
            chi = f"{t['chi2']:.1f}/{t['ndf']} ({t['p']:.3f})" if "chi2" in t else "–"
            L.append(f"| {t['sel']} | {eps_mc[(b, t['sel'])][0]:.3f} | {chi} | {g('J12_minus_MB427929_sigma', '{:+.2f}')} | {g('J12_over_MB427929')} | {g('Nsig_sum_over_all')} | {g('SF_weighted_minus_all', '{:+.3f}')} |")
        L.append("")
    open(os.path.join(a.out, "PSF.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:30]))


if __name__ == "__main__":
    main()
