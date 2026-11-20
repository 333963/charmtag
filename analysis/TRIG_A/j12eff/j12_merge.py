#!/usr/bin/env python3
"""J12EFF merge: suma histogramów per próbka → wydajność L1_J12 w funkcji wiodącego dżetu torowego / okna torowego / toru wiodącego
względem zbiorów odniesienia MBTS; dopasowanie erf (próg, szerokość, plateau). Wyjście: j12eff_<próbka>.npz, j12eff.json, J12EFF.md."""
import argparse, glob, json, os, re
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
from j12_extract import E_CJ, E_CW, E_LD, E_NS, REFS

EDG = {"cj": E_CJ, "cw": E_CW, "ld": E_LD}
LAB = {"cj": "pT wiodącego dżetu torowego [GeV]", "cw": "Σ pT torów w oknie 0,8×0,8 [GeV]", "ld": "pT toru wiodącego [GeV]"}


def turn(x, p, t, s):
    return 0.5 * p * (1 + erf((x - t) / (np.sqrt(2) * s)))


def eff(k, n):
    with np.errstate(all="ignore"):
        e = np.where(n > 0, k / n, np.nan)
        z = 1.0
        lo = (k + z * z / 2 - z * np.sqrt(k * (1 - e) + z * z / 4)) / (n + z * z)      # Wilson
        hi = (k + z * z / 2 + z * np.sqrt(k * (1 - e) + z * z / 4)) / (n + z * z)
    return e, np.where(n > 0, (hi - lo) / 2, np.nan)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--parts", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    samples = sorted({re.match(r"j12_(.*)_p\d+\.npz", os.path.basename(f)).group(1) for f in glob.glob(os.path.join(a.parts, "j12_*_p*.npz"))})
    J, md = {}, ["# J12EFF — wydajność L1_J12 z danych w zmiennych torowych (projekt A; zdarzenia, bez par D–D)", "",
                 "ε = N(odniesienie ∧ HLT_noalg_L1J12) / N(odniesienie); zdarzenia z ≥ 1 kandydatem (tak jak drzewo). „po preskali” = trigPassed, „TBP” = bit L1 przed preskalą (gdy gałąź jest).", ""]
    for s in samples:
        fs = sorted(glob.glob(os.path.join(a.parts, f"j12_{s}_p*.npz")))
        H, meta = {}, []
        for f in fs:
            z = np.load(f); meta.append(json.loads(str(z["meta"])))
            for k in z.files:
                if k.startswith("H_"):
                    H[k] = H.get(k, 0) + z[k]
        np.savez_compressed(os.path.join(a.out, f"j12eff_{s}.npz"), **H)
        ev = sum(m["events"] for m in meta)
        psj = sum(m["ps_j12_sum"] for m in meta) / max(1, sum(m["ps_j12_n"] for m in meta))
        psm = sum(m["ps_mbts1_sum"] for m in meta) / max(1, sum(m["ps_mbts1_n"] for m in meta))
        S = J[s] = {"parts": len(fs), "events": ev, "tbp_branch": meta[0]["tbp_branch"], "missing": sorted({c for m in meta for c in m["missing"]}),
                    "mean_ps_J12": psj, "mean_ps_mbts1": psm, "refs": {}}
        md += [f"## {s}", "", f"części {len(fs)}, zdarzeń {ev}, gałąź TBP: {S['tbp_branch']}, brak łańcuchów: {S['missing']}, średnia preskala J12 {psj:.2f}, mbts1 {psm:.1f}", ""]
        md += ["| odniesienie | D* | N_ref | N_ref∧J12 | ε całk. | ε TBP całk. | próg (dżet) | szer. | plateau | próg (okno) | szer. | plateau |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for ri, r in enumerate(REFS):
            for di, dn in enumerate(("wszystkie", "z D*")):
                R = S["refs"][f"{r}|{dn}"] = {}
                row = []
                for v in ("cj", "cw", "ld"):
                    h = H["H_" + v][ri, :, di].sum(axis=-1)              # [j, bins]
                    n, k, kt = h[0], h[1], h[2]
                    e, de = eff(k, n); et, det = eff(kt, n)
                    x = 0.5 * (EDG[v][:-1] + np.minimum(EDG[v][1:], EDG[v][:-1] + 20))
                    R[v] = {"x": x.tolist(), "n": n.tolist(), "k": k.tolist(), "k_tbp": kt.tolist(), "eff": np.nan_to_num(e, nan=-1).tolist(), "err": np.nan_to_num(de, nan=-1).tolist()}
                    fit = None
                    ok = (n >= 20) & np.isfinite(e)
                    kk = kt if kt.sum() > k.sum() else k
                    ee, dee = eff(kk, n)
                    if ok.sum() >= 5 and kk.sum() >= 30:
                        try:
                            p, c = curve_fit(turn, x[ok], ee[ok], sigma=np.maximum(dee[ok], 1e-4), p0=[max(ee[ok].max(), 1e-3), 15, 5], bounds=([0, 0, 0.3], [1.5, 80, 40]))
                            fit = {"plateau": p[0], "thr": p[1], "sigma": p[2], "err": np.sqrt(np.diag(c)).tolist(), "chi2": float((((turn(x[ok], *p) - ee[ok]) / np.maximum(dee[ok], 1e-4)) ** 2).sum()), "ndf": int(ok.sum() - 3)}
                        except Exception as ex:  # noqa: BLE001
                            fit = {"error": str(ex)}
                    R[v]["fit"] = fit
                    if v in ("cj", "cw"):
                        row += [f"{fit['thr']:.1f}", f"{fit['sigma']:.1f}", f"{fit['plateau']:.3f}"] if fit and "thr" in fit else ["–", "–", "–"]
                n0 = H["H_cj"][ri, 0, di].sum(); k0 = H["H_cj"][ri, 1, di].sum(); kt0 = H["H_cj"][ri, 2, di].sum()
                R["N_ref"], R["N_j12"], R["N_j12_tbp"] = float(n0), float(k0), float(kt0)
                md.append(f"| {r} | {dn} | {n0:.0f} | {k0:.0f} | {k0 / n0 if n0 else 0:.5f} | {kt0 / n0 if n0 else 0:.5f} | " + " | ".join(row) + " |")
        md.append("")
        for r in ("mbts1", "mbany"):
            for dn in ("wszystkie", "z D*"):
                for v in ("cj", "cw"):
                    R = S["refs"][f"{r}|{dn}"][v]
                    md += [f"### {s}: ε(J12) vs {LAB[v]} — odniesienie {r}, {dn}", "", "| bin | N_ref | N∧J12 | ε | N∧J12 TBP |", "|---|---|---|---|---|"]
                    for i in range(len(R["x"])):
                        if R["n"][i] > 0:
                            md.append(f"| {EDG[v][i]:g}–{EDG[v][i + 1]:g} | {R['n'][i]:.0f} | {R['k'][i]:.0f} | {R['eff'][i]:.4f} ± {R['err'][i]:.4f} | {R['k_tbp'][i]:.0f} |")
                    md.append("")
    json.dump(J, open(os.path.join(a.out, "j12eff.json"), "w"), indent=1)
    open(os.path.join(a.out, "J12EFF.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md[:40]))


if __name__ == "__main__":
    main()
