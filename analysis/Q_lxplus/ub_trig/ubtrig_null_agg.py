#!/usr/bin/env python3
"""UBTRIG — agregacja testu zerowego (19.11.2026, projekt A).

Klasy zdarzeń przetasowane losowo w warstwach (nogi × D*), liczności klas zachowane → rozkład μ_c i różnicy J12only − reszta
przy BRAKU związku z triggerem. Porównanie z biegiem prawdziwym: jeśli małe klasy przy losowym podziale też dają μ ≪ μ_nom,
to obciążenie metody (model tła w małych klasach), a nie trigger.
  ubtrig_null_agg.py --real ubtrig.json --null 'null/s*/ubtrig.json' --out DIR"""
import argparse
import glob
import json
import os

import numpy as np

CLASSES = ["J12only", "MBTS_j20", "MBTS_noj20", "other"]


def mu_of(rep, c):
    x = rep.get("fits", {}).get(c, {})
    cm = x.get("combined") if isinstance(x, dict) else None
    return float(cm["mu"]) if cm and cm.get("mu") is not None else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True); ap.add_argument("--null", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    R = json.load(open(a.real))
    N = [json.load(open(p)) for p in sorted(glob.glob(a.null))]
    mu_nom = R["fits"]["nominal"]["combined"]["mu"]
    rows, out = [], {"n_null": len(N), "mu_nom": mu_nom, "classes": {}}
    for c in CLASSES:
        v = np.array([mu_of(n, c) for n in N]); v = v[np.isfinite(v)]
        real = mu_of(R, c)
        m, s = (float(v.mean()), float(v.std(ddof=1))) if len(v) > 1 else (np.nan, np.nan)
        z = (real - m) / s if s > 0 else np.nan
        p_low = float((v <= real).mean()) if len(v) else np.nan
        out["classes"][c] = {"real": real, "null_mean": m, "null_std": s, "z_real_vs_null": z, "frac_null_le_real": p_low, "n": int(len(v)),
                             "null_values": v.tolist()}
        rows.append(f"| {c} | {real:.2f} | {m:.2f} ± {s:.2f} | {m - mu_nom:+.2f} | {z:+.2f} | {p_low:.3f} | {len(v)} |")
    d = np.array([n.get("J12_vs_rest", {}).get("diff_sigma", np.nan) for n in N]); d = d[np.isfinite(d)]
    mr = np.array([n.get("J12_vs_rest", {}).get("mu_rest", np.nan) for n in N]); mr = mr[np.isfinite(mr)]
    dr = R.get("J12_vs_rest", {}).get("diff_sigma", np.nan); mrr = R.get("J12_vs_rest", {}).get("mu_rest", np.nan)
    out["J12_vs_rest"] = {"real_diff_sigma": dr, "null_diff_mean": float(d.mean()) if len(d) else None, "null_diff_std": float(d.std(ddof=1)) if len(d) > 1 else None,
                          "frac_null_ge_real": float((d >= dr).mean()) if len(d) else None, "real_mu_rest": mrr,
                          "null_mu_rest_mean": float(mr.mean()) if len(mr) else None, "null_mu_rest_std": float(mr.std(ddof=1)) if len(mr) > 1 else None,
                          "null_diff_values": d.tolist()}
    j = out["J12_vs_rest"]
    zc = (dr - j["null_diff_mean"]) / j["null_diff_std"] if j["null_diff_std"] else np.nan
    j["z_real_vs_null"] = zc
    biased = abs(j["null_mu_rest_mean"] - mu_nom) > 2 * j["null_mu_rest_std"] / np.sqrt(len(mr)) if len(mr) > 1 else None
    out["verdict"] = ("metoda obciążona w małych klasach (μ_reszta w teście zerowym ≠ μ_nom)" if biased else
                      "metoda nieobciążona w małych klasach (μ_reszta w teście zerowym ≈ μ_nom)") if biased is not None else "za mało ziaren"
    md = ["# UBTRIG — test zerowy (klasy losowe o tych samych licznościach)", "",
          f"{len(N)} ziaren; klasy przetasowane w warstwach (zdarzenie z nogami × z kandydatem D*), liczności klas w warstwach jak w biegu prawdziwym. "
          f"μ_nom (próbka 05.11) = {mu_nom:.3f}. Oczekiwanie przy nieobciążonej metodzie: μ_c ≈ μ_nom w każdej klasie.", "",
          "| klasa | μ prawdziwe | μ test zerowy (średnia ± rozrzut) | średnia − μ_nom | Z (prawdziwe vs zerowy) | udział ziaren ≤ prawdziwe | n |",
          "|---|---|---|---|---|---|---|", *rows, "",
          f"J12only − reszta: prawdziwe {dr:+.2f}σ; test zerowy {j['null_diff_mean']:+.2f} ± {j['null_diff_std']:.2f}σ → Z = {zc:+.2f}; "
          f"udział ziaren ≥ prawdziwego: {j['frac_null_ge_real']:.3f}.",
          f"μ reszty: prawdziwe {mrr:.2f}; test zerowy {j['null_mu_rest_mean']:.2f} ± {j['null_mu_rest_std']:.2f} (rozrzut ziaren).", "",
          f"**Werdykt:** {out['verdict']}.", ""]
    json.dump(out, open(os.path.join(a.out, "ubtrig_null.json"), "w"), indent=1, default=float)
    open(os.path.join(a.out, "UBTRIG_NULL.md"), "w").write("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
