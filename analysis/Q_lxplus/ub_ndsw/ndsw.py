#!/usr/bin/env python3
"""UB_NDSW (19.11.2026, projekt A): czy szeroki pik D* w 427929 (skala σ 1,40 × MC) jest prawdziwy?
Ten sam fit Δm co UB427_NDS (2G z MC × skala σ, tło próg×pol2, okno D0 ±40 MeV, 139,6–160 MeV) na D* z 428770 (h5 v1.1, próbka nominału 05.11)
w podzbiorach klas triggera (tabela trg_*.npz paczki ub_trig, lookup po (run, event)) + 427929 (h5 po filtrze UB427) dla porównania.
Skala σ z profilu NLL (siatka 0,60–2,00, ΔNLL = 0,5) — całość i biny pT 3–5, 5–8, 8–30, 3–8.
  ndsw.py --h5-428770 'glob' --h5-427929 'glob' --trg 'trg_*.npz' --ref ndstar_full_428770.json --out DIR
  ndsw.py --selftest --ref ndstar_full_428770.json --out DIR"""
import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "ub_trig"))
import nds_var as V  # noqa: E402
import dstar_refit_nds as R  # noqa: E402

EDGES = np.round(np.arange(139.6, 170.0 + 1e-9, 0.05), 6)
M0_D0 = 1.86484
WIN = 0.040
PT_BINS = {"all": (3, 30), "3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30), "3-8": (3, 8)}
TCLS = ["J12only", "MBTS_j20", "MBTS_noj20", "other"]          # = ub_trig CLASSES (0..3); 4 = brak klucza w tabeli trg
GROUPS = {"all_428770": [0, 1, 2, 3, 4], "MBTS": [1, 2], "J12only": [0], "other": [3], "MBTS_j20": [1], "MBTS_noj20": [2], "nomatch": [4]}
GRID = np.round(np.arange(0.60, 2.00 + 1e-9, 0.05), 4)
REF_428770 = {"n_w40": 482979, "N": 4436.4, "scale": 1.024}


def read_h5(fs, chunk=1_000_000):
    """→ dm [MeV], pt, run, event dla D* (decayType 1, 0<Δm<165 MeV, |m − Δm − m_D0| < 40 MeV) — jak u4_ndstar / n8_ndstar."""
    import h5py
    out = {k: [] for k in ("dm", "pt", "run", "ev")}
    for f in fs:
        with h5py.File(f, "r") as g:
            e = g["events"]; n = e.shape[0]
            for i in range(0, n, chunk):
                dt = e.fields("Dmeson_decayType")[i:i + chunk]
                ix = np.flatnonzero(dt == 1)
                if not len(ix):
                    continue
                x = e.fields(["Dmeson_delta_mass", "Dmeson_mass", "Dmeson_pt", "run_number", "event_number"])[i:i + chunk][ix]
                s = (x["Dmeson_delta_mass"] > 0) & (x["Dmeson_delta_mass"] < 0.165) & \
                    (np.abs(x["Dmeson_mass"] - x["Dmeson_delta_mass"] - M0_D0) < WIN)
                out["dm"].append(x["Dmeson_delta_mass"][s] * 1e3); out["pt"].append(x["Dmeson_pt"][s])
                out["run"].append(x["run_number"][s].astype(np.int64)); out["ev"].append(x["event_number"][s].astype(np.int64))
        print(f"  {os.path.basename(f)}: D* w oknie {sum(len(v) for v in out['dm'])} (narastająco)", flush=True)
    return {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in out.items()}


def hists(d, cls=None):
    """{grupa: {bin: hist}}; cls = klasa 0..4 per kandydat (None → jedna grupa 'r427929')."""
    H = {}
    groups = {"r427929": None} if cls is None else GROUPS
    for g, ks in groups.items():
        sel = np.ones(len(d["dm"]), bool) if ks is None else np.isin(cls, ks)
        H[g] = {b: np.histogram(d["dm"][sel & (d["pt"] >= lo) & (d["pt"] < hi)], EDGES)[0].astype(float) for b, (lo, hi) in PT_BINS.items()}
    return H


def profile_scale(h, shape):
    """Skala σ z profilu NLL (tło próg×pol2, skala stała w siatce) → s_hat, s_lo, s_hi (ΔNLL = 0,5), lub NaN."""
    nl = []
    for s in GRID:
        try:
            r = V.fit_generic(h, EDGES, shape, bkg="B", s_free=False, s0=float(s))
            nl.append(r["nll"] if r["valid"] else np.nan)
        except Exception:  # noqa: BLE001
            nl.append(np.nan)
    nl = np.array(nl); ok = np.isfinite(nl)
    if ok.sum() < 5:
        return {"s_hat": np.nan, "s_lo": np.nan, "s_hi": np.nan, "grid": GRID.tolist(), "nll": nl.tolist()}
    g, y = GRID[ok], nl[ok]; i = int(np.argmin(y))
    j = slice(max(0, i - 2), min(len(g), i + 3)); c = np.polyfit(g[j], y[j], 2)
    s_hat = float(-c[1] / (2 * c[0])) if c[0] > 0 else float(g[i])
    s_hat = float(np.clip(s_hat, g[0], g[-1])); ymin = min(float(np.polyval(c, s_hat)) if c[0] > 0 else y[i], y[i]); d = y - ymin
    fine = np.linspace(g[0], g[-1], 2001); df = np.interp(fine, g, d)
    left, right = fine[fine <= s_hat], fine[fine >= s_hat]
    lo = float(left[np.flatnonzero(df[fine <= s_hat] >= 0.5)[-1]]) if (df[fine <= s_hat] >= 0.5).any() else np.nan
    hi = float(right[np.flatnonzero(df[fine >= s_hat] >= 0.5)[0]]) if (df[fine >= s_hat] >= 0.5).any() else np.nan
    return {"s_hat": s_hat, "s_lo": lo, "s_hi": hi, "grid": GRID.tolist(), "nll": nl.tolist(),
            "dnll_at_1p02": float(np.interp(1.02, g, d)), "dnll_at_1p40": float(np.interp(1.40, g, d))}


def fit_group(hb, shape):
    out = {}
    for b, h in hb.items():
        r = {"n_entries": int(h.sum())}
        try:
            f = V.fit_generic(h, EDGES, shape, bkg="B"); f.pop("attempts", None); r["free"] = f
        except Exception as ex:  # noqa: BLE001
            r["free"] = {"error": repr(ex)}
        try:
            f = V.fit_generic(h, EDGES, shape, bkg="B", s_free=False, s0=1.0); f.pop("attempts", None); r["fixed_mc"] = f
        except Exception as ex:  # noqa: BLE001
            r["fixed_mc"] = {"error": repr(ex)}
        r["profile"] = profile_scale(h, shape)
        p = r["profile"]
        print(f"    {b:5s} n={r['n_entries']:7d} N={r['free'].get('N', np.nan):6.0f}±{r['free'].get('N_err', np.nan):4.0f} "
              f"s={p['s_hat']:.2f} [{p['s_lo']:.2f}; {p['s_hi']:.2f}] N(s=1)={r['fixed_mc'].get('N', np.nan):6.0f}", flush=True)
        out[b] = r
    return out


def verdict(res):
    """Kryterium zapisane w README przed wynikiem (porównanie MBTS(428770) z 427929; binem rozstrzygającym jest 3–8 GeV i całość)."""
    def sc(g, b):
        p = res["groups"].get(g, {}).get(b, {}).get("profile", {}); return p.get("s_hat", np.nan), p.get("s_lo", np.nan), p.get("s_hi", np.nan)
    v = {}
    for b in ("all", "3-8"):
        s, lo, hi = sc("MBTS", b)
        sig = 0.5 * (hi - lo) if np.isfinite(lo) and np.isfinite(hi) else np.nan
        near_mc = np.isfinite(s) and np.isfinite(sig) and abs(s - 1.02) <= 2 * sig
        near_427 = np.isfinite(s) and np.isfinite(sig) and abs(s - 1.40) <= 2 * sig
        v[b] = {"s": s, "sigma": sig, "zgodne_z_1.02": bool(near_mc), "zgodne_z_1.40": bool(near_427)}
    a = v["all"]
    if a["zgodne_z_1.02"] and not a["zgodne_z_1.40"]:
        w = "SZEROKOŚĆ NIEPRAWDZIWA: MBTS(428770) ≈ MC → w 427929 tło wciągnięte w sygnał; N_D* ze stałą skalą (μ427929 ≈ 0,77)"
    elif a["zgodne_z_1.40"] and not a["zgodne_z_1.02"]:
        w = "SZEROKOŚĆ PRAWDZIWA: MBTS(428770) ≈ 1,4 → zdarzenia MB mają szerszy pik (np. rozdzielczość miękkiego π / wielokrotne rozpraszanie); μ427929 = 0,53 zostaje"
    else:
        w = "NIEROZSTRZYGNIĘTE: MBTS(428770) zgodne z oboma (lub z żadnym) w 2σ — syst N_D* 427929 zostaje (±0,22)"
    return v, w


def synth(shape, seed=5):
    """selftest (mechanika, S/B zawyżone): 428770 — J12only (s=1,0, N 6000), MBTS_j20/noj20 (s=1,4, N 1500+1500), other (s=1,0, N 200); 427929 s=1,4 N 2000."""
    rng = np.random.default_rng(seed)
    c = lambda s: np.diff((lambda q: (q - q[0]) / (q[-1] - q[0]))(R._sig_cdf(EDGES, shape["mu"], shape["sc"] * s, shape["k"], shape["f"])))  # noqa: E731
    u = np.clip(0.5 * (EDGES[1:] + EDGES[:-1]) - R.MPI, 1e-6, None); pb = u ** 0.5 * (1 + 0.3 * u / 20); pb /= pb.sum()
    ctr = 0.5 * (EDGES[1:] + EDGES[:-1])

    def draw(ns, nb, s):
        k = rng.poisson(ns * c(s) + nb * pb); dm = np.repeat(ctr, k.astype(int)) + rng.uniform(-0.025, 0.025, int(k.sum()))
        pt = 3 + rng.exponential(3.0, len(dm)); pt = np.where(pt < 30, pt, 29.9); return dm, pt
    parts, cls = [], []
    for k, (ns, nb, s) in {0: (6000, 200000, 1.0), 1: (1500, 25000, 1.4), 2: (1500, 25000, 1.4), 3: (200, 20000, 1.0)}.items():
        dm, pt = draw(ns, nb, s); parts.append((dm, pt)); cls.append(np.full(len(dm), k, np.int8))
    d8 = {"dm": np.concatenate([p[0] for p in parts]), "pt": np.concatenate([p[1] for p in parts])}
    dm, pt = draw(2000, 110000, 1.4)
    return d8, np.concatenate(cls), {"dm": dm, "pt": pt}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5-428770"); ap.add_argument("--h5-427929"); ap.add_argument("--trg"); ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    shape = json.load(open(a.ref))["shape_2g_mc"]
    info = {}
    if a.selftest:
        d8, cls, d7 = synth(shape); info["selftest"] = True
    else:
        from ubtrig_run import load_trg, lookup_cls
        import btrig_common as C
        f8 = sorted(glob.glob(a.h5_428770)); f7 = sorted(glob.glob(a.h5_427929))
        print(f"h5 428770: {len(f8)}; h5 427929: {len(f7)}", flush=True)
        tkey, tcls = load_trg(sorted(glob.glob(a.trg)))
        d8 = read_h5(f8); d7 = read_h5(f7)
        cls, nmiss = lookup_cls(tkey, tcls, C.ev_key(d8["run"], d8["ev"]), 4)
        info.update(files_428770=len(f8), files_427929=len(f7), n_trg_events=int(len(tkey)), dstar_nomatch=int(nmiss), n_dstar_428770=int(len(cls)))
    info["n_w40_428770"] = int(len(d8["dm"][(d8["dm"] > 139.6) & (d8["dm"] < 170)]))
    info["class_counts"] = {g: int(np.isin(cls, ks).sum()) for g, ks in GROUPS.items()}
    print("klasy kandydatów:", info["class_counts"], flush=True)
    H8 = hists(d8, cls); H7 = hists(d7)
    np.savez(os.path.join(a.out, "ndsw_hists.npz"), edges=EDGES, **{f"{g}__{b}": h for G in (H8, H7) for g, hb in G.items() for b, h in hb.items()})
    res = {"info": info, "ref_428770_nominal": REF_428770, "shape_2g_mc": shape, "groups": {}}
    for g, hb in list(H8.items()) + list(H7.items()):
        print(f"[{g}]", flush=True)
        res["groups"][g] = fit_group(hb, shape)
    res["verdict_detail"], res["verdict"] = verdict(res)
    json.dump(res, open(os.path.join(a.out, "ndsw.json"), "w"), indent=1, default=float)
    md = ["# UB_NDSW — szerokość piku D* (skala σ × MC) w klasach triggera 428770 vs 427929", "",
          f"Fit Δm jak UB427_NDS (2G z MC × skala σ, tło próg×pol2, okno D0 ±40 MeV). Skala σ z profilu NLL (ΔNLL = 0,5). "
          f"Kandydaci 428770 w oknie: {info['n_w40_428770']} (nominał 05.11: {REF_428770['n_w40']}); bez klucza w tabeli trg: {info.get('dstar_nomatch', '–')}.", "",
          "| grupa | bin pT | wpisy | N_D* (skala wolna) | skala σ [68 %] | N_D* (skala = 1) | ΔNLL(1,02) | ΔNLL(1,40) |", "|---|---|---|---|---|---|---|---|"]
    for g, gb in res["groups"].items():
        for b, r in gb.items():
            p, f, x = r["profile"], r["free"], r["fixed_mc"]
            md.append(f"| {g} | {b} | {r['n_entries']} | {f.get('N', np.nan):.0f} ± {f.get('N_err', np.nan):.0f} | {p['s_hat']:.2f} [{p['s_lo']:.2f}; {p['s_hi']:.2f}] | "
                      f"{x.get('N', np.nan):.0f} ± {x.get('N_err', np.nan):.0f} | {p.get('dnll_at_1p02', np.nan):.1f} | {p.get('dnll_at_1p40', np.nan):.1f} |")
    md += ["", "Kryterium (README, zapisane przed wynikiem): skala MBTS(428770), całość 3–30 GeV, w 2σ od 1,02 i nie od 1,40 → szerokość w 427929 nieprawdziwa; "
           "w 2σ od 1,40 i nie od 1,02 → prawdziwa; inaczej nierozstrzygnięte. Bin 3–8 GeV pomocniczo.",
           f"Szczegóły: {json.dumps({k: {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in res['verdict_detail'].items()}, ensure_ascii=False)}",
           "", f"**Werdykt:** {res['verdict']}.", ""]
    open(os.path.join(a.out, "NDSW.md"), "w").write("\n".join(md))
    print("\n".join(md[-3:]), flush=True)
    if a.selftest:
        s_m = res["groups"]["MBTS"]["all"]["profile"]["s_hat"]; s_j = res["groups"]["J12only"]["all"]["profile"]["s_hat"]
        ok = abs(s_m - 1.4) < 0.2 and abs(s_j - 1.0) < 0.15
        print(f"SELFTEST: s(MBTS) {s_m:.2f} (prawda 1,4), s(J12only) {s_j:.2f} (prawda 1,0) → {'OK' if ok else 'FAIL'}")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
