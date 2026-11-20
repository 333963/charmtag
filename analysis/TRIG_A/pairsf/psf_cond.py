#!/usr/bin/env python3
"""TRIG_A/pairsf/psf_cond (18.11.2026, projekt A) — WARUNKOWA wydajność score H dla D0 z tagu D*: ε(≥ WP_x | ≥ WP90), dane vs MC fd_nd.
Kontrola wyniku pary hybrydowej (PAIRNET hybsf: SF nogi D0 przy pT 3–5 GeV ≈ 0,6 w 428770). Mianownik = D0←D* z score ≥ WP90 (jak nogi par),
więc próbka „fail” (WP90–WP_x) jest czysta — inaczej niż w psf_analyze (mianownik: wszystkie D0←D*, fit nieokreślony przy 3–5 GeV).
Fit Δm pass/fail (btrig_common.fit_eps). Pojedyncze D, bez par. Wejście: B_trig/runs/lxplus/bt_*.npz.
Użycie: psf_cond.py --in DIR --out DIR [--nproc 6]"""
import argparse, json, os, sys
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "B_trig")); sys.path.insert(0, os.path.join(HERE, "..", "..", "AA_merge", "src"))
import btrig_common as C  # noqa: E402

THR = {"0.3": 0.9839380383491516, "0.45": 0.9589824676513672, "0.65": 0.8569592237472534, "0.9": 0.38729336857795715}
SEL = {"ge65": "0.65", "ge45": "0.45", "ge30": "0.3"}
BINS = {"3-4": (3.0, 4.0), "4-5": (4.0, 5.0), "3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0)}
NEED = ("m", "pt", "dm", "score")


def load(p, mc=False):
    z = np.load(p, allow_pickle=False)
    return {k: z[k] for k in z.files if k in NEED or (mc and k in ("truth", "truth_ds"))}


def _fit(job):
    key, dp, df_, shape = job
    return key, C.fit_eps(dp, df_, shape, interval=True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--in", dest="inp", required=True); ap.add_argument("--out", required=True); ap.add_argument("--nproc", type=int, default=6)
    a = ap.parse_args()
    from charmtag import dstar_refit as DR
    M = load(os.path.join(a.inp, "bt_fd_nd.npz"), mc=True)
    sig = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < C.NOM_WIN) & (M["score"] >= THR["0.9"])
    D = {r: load(os.path.join(a.inp, f"bt_data_{r}.npz")) for r in ("428770", "427929")}
    D["oba"] = {k: np.concatenate([D["428770"][k], D["427929"][k]]) for k in NEED}
    jobs, emc = [], {}
    for b, (lo, hi) in BINS.items():
        ms = sig & (M["pt"] >= lo) & (M["pt"] < hi)
        shape = DR.shape_from_mc(M["dm"][ms] * 1e3)
        for sn, t in SEL.items():
            p = M["score"][ms] >= THR[t]
            emc[(b, sn)] = (float(p.mean()), float(np.sqrt(p.mean() * (1 - p.mean()) / max(len(p), 1))), int(len(p)))
            for r, d in D.items():
                base = (np.abs(d["m"] - C.M0_D0) < C.NOM_WIN) & (d["pt"] >= lo) & (d["pt"] < hi) & (d["score"] >= THR["0.9"])
                pas = d["score"] >= THR[t]
                jobs.append(((r, b, sn), d["dm"][base & pas], d["dm"][base & ~pas], shape))
    with Pool(a.nproc) as P:
        F = dict(P.map(_fit, jobs, chunksize=1))
    rows = []
    L = ["# psf_cond — ε(≥ WP_x | ≥ WP90) dla D⁰ z tagu D*: dane / MC fd_nd (fit Δm pass/fail; pojedyncze D)", "",
         "| pT [GeV] | selekcja | ε MC | run | N_sig (≥WP90) | ε dane | SF | fit |", "|---|---|---|---|---|---|---|---|"]
    for (r, b, sn), f in F.items():
        e, ee, n = emc[(b, sn)]; err = C.eps_sym_err(f); sf = f["eps"] / e
        rows.append({"run": r, "bin": b, "sel": sn, "eps_mc": e, "eps_mc_err": ee, "n_mc": n, "eps": f["eps"], "eps_err": err, "N": f["N"], "N_err": f["N_err"],
                     "determined": f["determined"], "valid": f["valid"], "n_pass": f["n_pass"], "n_fail": f["n_fail"], "SF": sf, "SF_err": float(np.hypot(err / e, sf * ee / e))})
        L.append(f"| {b} | {sn} | {e:.3f} | {r} | {f['N']:.0f} ± {f['N_err']:.0f} | {f['eps']:.3f} ± {err:.3f} | **{sf:.3f} ± {rows[-1]['SF_err']:.3f}** | {'ok' if f['determined'] else 'nieokreślony'} |")
    os.makedirs(a.out, exist_ok=True)
    json.dump({"thr": THR, "rows": rows}, open(os.path.join(a.out, "psf_cond.json"), "w"), indent=1, default=float)
    open(os.path.join(a.out, "PSF_COND.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
