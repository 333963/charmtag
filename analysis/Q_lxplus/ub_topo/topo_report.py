#!/usr/bin/env python3
"""UB_TOPO raport (19.11.2026, projekt A): kształt |Δη| (lub Δφ) sygnału par D–D vs hipotezy z generatora Pythii (ref/topo_gen.json).

Wejście: fit2d_diff.json + pair_inputs_diff.json z `charmtag pairs/fit2d --unblind --diff <var>` (profile ΔNLL μ_i per kanał × bin).
Hipotezy kształtu h (udziały binów sygnału; ważone ε nóg): pythia (pary z jednego zdarzenia = oczekiwanie MC), dps (pary niezależne).
Rodzina h(λ) = λ·h_dps + (1 − λ)·h_pythia: λ = udział sygnału o kształcie par niezależnych PONAD skład Pythii (λ = 0: kształt Pythii,
λ = 1: czysty DPS). Dla każdego λ wspólne μ (kanały × biny) przy kształcie h(λ): q(λ) = 2·min_μ Σ ΔNLL_i(μ·h_i/f_i) (diff.combine_bins);
profil Δq(λ) = q(λ) − min q → λ̂, 68 % (Δq ≤ 1), Z(λ = 0), Z(λ = 1). Odniesienie: „nadwyżka = DPS” ⇒ λ = (μ − 1)/μ.
Prognoza mocy (σ_λ z biegu, przed liczeniem w PROGNOZA.md) — wynik interpretować tylko przez nią.
  topo_report.py --pairs DIR --gen ref/topo_gen.json --out DIR [--weights w|u] [--mu-ref 1.50] [--selftest]"""
import argparse
import json
import os
import sys

import numpy as np

LAM = np.round(np.arange(-1.5, 2.5001, 0.01), 3)
CH = ("D0D0", "DpD0", "DpDp")


def _diff_mod():
    from charmtag.pairs import diff as DF
    return DF


def gen_var_key(var, nb):
    return {"dy": "dy", "dphi": "dphi"}.get(var, var) if not (var == "dphi" and nb == 3) else "dphi_u"


def lam_scan(prof, ffit, hyp, chans):
    """prof {(ch, k): scan}, ffit {(ch, k): f_fit}, hyp {ch: {"pythia": [..], "dps": [..]}} → Δq(λ), μ(λ)."""
    DF = _diff_mod()
    q, mu = [], []
    for lam in LAM:
        sc = {}
        for (ch, k) in prof:
            h = lam * np.asarray(hyp[ch]["dps"]) + (1 - lam) * np.asarray(hyp[ch]["pythia"])
            sc[(ch, k)] = max(float(h[k]), 1e-6) / ffit[(ch, k)]
        r = DF.combine_bins(prof, sc)
        q.append(r["q_shape"]); mu.append(r["mu"])
    q = np.asarray(q); dq = q - q.min()
    return dq, np.asarray(mu)


def summarize(dq, mu, lam_ref):
    i = int(np.argmin(dq)); ok = LAM[dq <= 1.0]
    at = lambda x: float(np.interp(x, LAM, dq))  # noqa: E731
    return {"lam_hat": float(LAM[i]), "lam_68": [float(ok.min()), float(ok.max())], "sigma_lam": float((ok.max() - ok.min()) / 2),
            "mu_at_hat": float(mu[i]), "Z_lam0_pythia": float(np.sqrt(at(0.0))), "Z_lam1_dps": float(np.sqrt(at(1.0))),
            "lam_ref_excess_dps": lam_ref, "Z_lam_ref": float(np.sqrt(at(lam_ref))) if lam_ref is not None else None,
            "at_scan_edge": bool(i in (0, len(LAM) - 1))}


def load_inputs(pairs_dir):
    rep = json.load(open(os.path.join(pairs_dir, "fit2d_diff.json")))
    D = json.load(open(os.path.join(pairs_dir, "pair_inputs_diff.json")))
    return rep, D


def analyze(rep, D, gen, weights="w", mu_ref=1.5):
    spec = rep["spec"]; labels = rep["labels"]; nb = len(labels)
    gk = f"{gen_var_key(spec['var'], nb)}_{weights}"
    if spec.get("edges") is not None:
        ge = gen["edges"][gen_var_key(spec["var"], nb)]
        if not np.allclose(ge, spec["edges"]):
            raise SystemExit(f"STOP: biny danych {spec['edges']} ≠ biny hipotez z generatora {ge}")
    prof, ffit, hyp = {}, {}, {}
    for ch, rc in rep["channels"].items():
        for k, lab in enumerate(labels):
            b = rc["bins"].get(lab, {})
            if "profile_dNLL" in b:
                prof[(ch, k)] = b["profile_dNLL"]; ffit[(ch, k)] = float(b["f_fit"])
        hyp[ch] = {"pythia": gen[ch][gk]["pythia"], "dps": gen[ch][gk]["dps"], "pythia_err": gen[ch][gk]["pythia_err"]}
    lam_ref = (mu_ref - 1.0) / mu_ref if mu_ref else None
    out = {"var": spec["var"], "edges": spec.get("edges"), "labels": labels, "weights": weights, "lam_ref_mu": mu_ref, "channels": {}, "hyp": hyp}
    for ch in sorted({c for c, _ in prof}):
        pc = {key: v for key, v in prof.items() if key[0] == ch}
        if len(pc) == nb:
            dq, mu = lam_scan(pc, ffit, hyp, [ch])
            out["channels"][ch] = summarize(dq, mu, lam_ref)
    dq, mu = lam_scan(prof, ffit, hyp, list(hyp))
    out["combined"] = summarize(dq, mu, lam_ref)
    out["combined"]["channels"] = sorted({c for c, _ in prof})
    out["scan_combined"] = {"lam": LAM[::10].tolist(), "dq": dq[::10].tolist()}
    return out


def to_md(o, forecast=None):
    L = [f"# UB_TOPO — kształt `{o['var']}` sygnału par D–D vs Pythia/DPS (odsłonięte, 428770)", "",
         f"Biny {o['labels']}; hipotezy z generatora Pythii A3 (ważone ε nóg: `{o['weights']}`), rodzina h(λ) = λ·h_DPS + (1−λ)·h_Pythia.",
         "λ = 0 — kształt Pythii (nominał MC), λ = 1 — czysty DPS (pary niezależne); „nadwyżka = DPS” ⇒ λ = (μ−1)/μ "
         f"= {o['combined']['lam_ref_excess_dps']:.2f} (μ = {o['lam_ref_mu']}).", "",
         "| kanał | h_Pythia | h_DPS | λ̂ | 68 % | Z(λ=0, Pythia) | Z(λ=1, DPS) | Z(λ = nadwyżka-DPS) |", "|---|---|---|---|---|---|---|---|"]
    rows = list(o["channels"].items()) + [("łącznie", o["combined"])]
    for ch, r in rows:
        hp = o["hyp"].get(ch, {}).get("pythia"); hd = o["hyp"].get(ch, {}).get("dps")
        L.append(f"| {ch} | {'' if hp is None else '/'.join(f'{x:.3f}' for x in hp)} | {'' if hd is None else '/'.join(f'{x:.3f}' for x in hd)} | "
                 f"{r['lam_hat']:.2f}{' (krawędź)' if r['at_scan_edge'] else ''} | [{r['lam_68'][0]:.2f}; {r['lam_68'][1]:.2f}] | "
                 f"{r['Z_lam0_pythia']:.2f} | {r['Z_lam1_dps']:.2f} | {r['Z_lam_ref']:.2f} |")
    L += ["", "Interpretacja (ustalona przed wynikiem, PROGNOZA.md): σ_λ ≈ 0,8 → test rozróżnia najwyżej czysty DPS od Pythii (~1,2σ);",
          "pytania „nadwyżka DPS czy dżetowa” ten test nie rozstrzyga (Δλ ≈ 0,3 ≪ σ_λ; obciążenie dżetowe daje kształt SPS-podobny, λ ≤ 0)."]
    if forecast:
        L += ["", f"Prognoza (PROGNOZA.md): {forecast}"]
    return "\n".join(L) + "\n"


def _shift_scan(scan, d):
    """Profil przesunięty w μ o d (min w nowym miejscu), na tej samej siatce kluczy."""
    DF = _diff_mod()
    g = np.array(sorted(float(k) for k in scan))
    v = DF.profile_at(scan, g - d)
    return {f"{x:g}": float(y) for x, y in zip(g, v)}


def asimov_rep(pd_rep, gen, lam_true, mu_true=1.5, scale_sigma=1.0):
    """Profile z biegu PD (sign; like ↔ bin |Δη|≥1, opposite ↔ |Δη|<1 — podobne udziały) przesunięte do wartości Asimova
    μ_i = μ·h_i(λ_true)/f_i; scale_sigma > 1 poszerza profil (np. klasa triggera o mniejszej statystyce)."""
    rep = {"spec": {"var": "dy", "edges": [0.0, 1.0, 5.0]}, "labels": ["[0,1)", "[1,5]"], "channels": {}}
    src = {0: "opposite", 1: "like"}
    for ch, rc0 in pd_rep["channels"].items():
        hp, hd = np.asarray(gen[ch]["dy_w"]["pythia"]), np.asarray(gen[ch]["dy_w"]["dps"])
        h = lam_true * hd + (1 - lam_true) * hp
        rc = {"bins": {}}
        for k, lab in enumerate(rep["labels"]):
            b = rc0["bins"].get(src[k], {})
            if "profile_dNLL" not in b:
                continue
            g = np.array(sorted(float(x) for x in b["profile_dNLL"]))
            v = np.array([b["profile_dNLL"][f"{x:g}"] for x in g]); m0 = g[np.argmin(v)]
            sc = {f"{x:g}": float(y) / scale_sigma ** 2 for x, y in zip(g, v)}
            target = mu_true * h[k] / hp[k]
            rc["bins"][lab] = {"f_fit": float(hp[k]), "profile_dNLL": _shift_scan(sc, target - m0)}
        rep["channels"][ch] = rc
    return rep


def forecast(pd_pairs, gen, out):
    """Prognoza mocy (Asimov) z kształtów profili biegu PD: σ_λ, Z(Pythia | DPS prawdziwy), Z(nadwyżka-DPS | Pythia prawdziwa),
    klasa triggera (σ × 1/√r_c) i 427929 (μ 0,53, σ_μ jak 428770)."""
    pd_rep = json.load(open(os.path.join(pd_pairs, "fit2d_diff.json")))
    F = {}
    for nm, lam_true, mu, ss in (("dps_true_428770", 1.0, 1.5, 1.0), ("pythia_true_428770", 0.0, 1.5, 1.0),
                                 ("J12only_r0.77", 0.0, 1.5, 1 / np.sqrt(0.77)), ("rest_r0.23", 0.0, 1.5, 1 / np.sqrt(0.23)),
                                 ("r427929_mu0.53", 0.0, 0.53, 1.0)):
        o = analyze(asimov_rep(pd_rep, gen, lam_true, mu, ss), {}, gen, "w", 1.5)["combined"]
        F[nm] = {"lam_true": lam_true, "mu": mu, "sigma_scale": ss, "sigma_lam": o["sigma_lam"], "lam_68": o["lam_68"],
                 "Z_pythia": o["Z_lam0_pythia"], "Z_dps": o["Z_lam1_dps"], "Z_excess_dps": o["Z_lam_ref"]}
    s1, s2 = F["J12only_r0.77"]["sigma_lam"], F["rest_r0.23"]["sigma_lam"]
    F["J12_vs_rest_sigma_dlam"] = float(np.hypot(s1, s2))
    json.dump(F, open(os.path.join(out, "forecast.json"), "w"), indent=1)
    for k, v in F.items():
        print(k, v)
    return F


def selftest(gen, out):
    """Domknięcie na syntetycznych profilach parabolicznych: dane z kształtem h(λ_true), σ_i jak w PD (sign) → odzysk λ (bez obciążenia)."""
    rng = np.random.default_rng(3)
    grid = np.round(np.arange(-1.5, 10.01, 0.25), 2)
    res = {}
    for lam_true in (0.0, 1.0):
        hats = []
        for t in range(40):
            rep = {"spec": {"var": "dy", "edges": [0.0, 1.0, 5.0]}, "labels": ["[0,1)", "[1,5]"], "channels": {}}
            for ch, sig in (("D0D0", (0.91, 1.42)), ("DpD0", (0.62, 1.25)), ("DpDp", (1.15, 4.47))):
                hp, hd = np.asarray(gen[ch]["dy_w"]["pythia"]), np.asarray(gen[ch]["dy_w"]["dps"])
                h = lam_true * hd + (1 - lam_true) * hp
                rc = {"bins": {}}
                for k, lab in enumerate(rep["labels"]):
                    ff = hp[k]; mu_true = 1.5 * h[k] / ff; s = sig[k] * 0.25      # σ/4: test samej metody (odzysk λ), nie prognoza
                    m = mu_true + rng.normal(0, s)
                    rc["bins"][lab] = {"f_fit": ff, "profile_dNLL": {f"{g:g}": float(0.5 * ((g - m) / s) ** 2) for g in grid}}
                rep["channels"][ch] = rc
            hats.append(analyze(rep, {}, gen, "w", 1.5)["combined"]["lam_hat"])
        res[lam_true] = (float(np.mean(hats)), float(np.std(hats)))
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(out, "selftest.json"), "w"), indent=1)
    for lt, (m, s) in res.items():
        print(f"selftest λ_true={lt}: λ̂ = {m:.3f} ± {s:.3f} (rozrzut 40 prób)")
        if abs(m - lt) > 3 * s / np.sqrt(40) + 0.05:
            raise SystemExit(f"SELFTEST FAIL: λ_true={lt}, średnie λ̂={m:.3f}")
    print("SELFTEST OK")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs"); ap.add_argument("--gen", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--weights", default="w", choices=["w", "u"]); ap.add_argument("--mu-ref", type=float, default=1.5)
    ap.add_argument("--forecast", default=""); ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--forecast-from", help="katalog pairs biegu PD (sign) → prognoza mocy (Asimov), forecast.json")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    gen = json.load(open(a.gen))
    if a.selftest:
        return selftest(gen, a.out)
    if a.forecast_from:
        return forecast(a.forecast_from, gen, a.out)
    if os.environ.get("CHARM_UNBLIND") != "1":
        raise SystemExit("UB_TOPO: wymaga CHARM_UNBLIND=1")
    rep, D = load_inputs(a.pairs)
    o = analyze(rep, D, gen, a.weights, a.mu_ref)
    o["weights_alt"] = analyze(rep, D, gen, "u" if a.weights == "w" else "w", a.mu_ref)["combined"]   # syst. wagi ε (w ↔ u)
    json.dump(o, open(os.path.join(a.out, "ubtopo.json"), "w"), indent=1, default=float)
    md = to_md(o, a.forecast) + f"\nSyst. wag ε (hipotezy bez wag ε): λ̂ = {o['weights_alt']['lam_hat']:.2f} {o['weights_alt']['lam_68']}.\n"
    open(os.path.join(a.out, "UBTOPO.md"), "w").write(md)
    print(md)


if __name__ == "__main__":
    sys.exit(main())
