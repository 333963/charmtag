#!/usr/bin/env python3
"""mult_gates.py — MULT_charm: bramki G0, G1a, G1b, (G1c część zdarzeniowa) z npz mult_extract.py (PROCEDURA §4). Bez fitów, bez par.
  python3 mult_gates.py --mc DIR --data427929 DIR --data428770 DIR --out gates.json
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mult_common as C  # noqa: E402


def load_ev(d):
    fs = sorted(glob.glob(os.path.join(d, "*.npz")))
    E, metas = {}, []
    for f in fs:
        z = np.load(f); m = json.loads(str(z["meta"])); metas.append(m)
        for k in z.files:
            if k.startswith("ev_"):
                E.setdefault(k[3:], []).append(z[k])
    lens = {k: sum(len(x) for x in v) for k, v in E.items()}
    E = {k: np.concatenate(v) for k, v in E.items() if lens[k] == max(lens.values())}
    return E, metas


def g0(metas, exp_files=None, exp_events=None):
    n = sum(m["n_events"] for m in metas)
    r = {"files": len(metas), "events": n, "formats": sorted({m["format"] for m in metas}),
         "compat_missing": sorted({c for m in metas for c in m["compat_missing"]})}
    ok = r["formats"] == ["v2.0"] and not r["compat_missing"]
    if exp_files:
        ok &= r["files"] == exp_files
    if exp_events:
        ok &= n == exp_events
    r["PASS"] = bool(ok)
    return r


def basesel(E, data=True):
    s = (E["npv"] == 1) & (np.abs(E["pvz"]) < C.EVSEL["pvz_max"]) & (E["nbl"] >= 1)
    if data:
        s &= E["grl"] == 1
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc"); ap.add_argument("--data427929"); ap.add_argument("--data428770"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    R = {}
    if a.mc:
        E, M = load_ev(a.mc); R["G0_mb_nd"] = g0(M, 80, 19714175)
        R["G1c_mc_npu_frac"] = float((E["npu"][basesel(E, False)] > 0).mean())
    for run, d, nf, ne in (("427929", a.data427929, None, None), ("428770", a.data428770, 97, 7989496)):
        if not d:
            continue
        E, M = load_ev(d); R[f"G0_{run}"] = g0(M, nf, ne)
        s = basesel(E)
        R[f"G1c_{run}_npu_frac"] = float((E["npu"][s] > 0).mean())
        ch = "HLT_mb_mbts_L1MBTS_1"
        if f"ps_{ch}" in E:
            lbs = np.unique(E["lb"][s]); ps1 = np.array([np.all(E[f"ps_{ch}"][s & (E["lb"] == lb)] == 1) for lb in lbs])
            R[f"G1a_{run}"] = {"frac_LB_PS1": float(ps1.mean()), "n_LB": int(len(lbs)), "PASS": bool(ps1.mean() >= 0.95)}
            ref = s & (E[f"tp_{ch}"] == 1) & (E["npu"] == 0)
            for sc in ("sptrk4", "sptrk2"):
                c2 = C.TRIG_SCHEMES[sc]["chain"]; lm = C.TRIG_SCHEMES[sc]["lead_min"]
                if f"tp_{c2}" not in E:
                    R[f"G1b_{run}_{sc}"] = {"note": f"brak {c2} w próbce"}; continue
                out = {}
                for thr in (lm, 6.0, 8.0):
                    r = ref & (E["leadbl"] > thr)
                    n = int(r.sum())
                    eff = float((E[f"tp_{c2}"][r] * E[f"ps_{c2}"][r]).sum() / n) if n else float("nan")   # HLT po prescale × PS
                    out[f"lead>{thr}"] = {"n_ref": n, "eff": eff}
                out["PASS"] = bool(out[f"lead>{lm}"]["eff"] >= 0.95) if out[f"lead>{lm}"]["n_ref"] >= 50 else None
                R[f"G1b_{run}_{sc}"] = out
        for sc, info in C.TRIG_SCHEMES.items():
            if info["chain"] and f"tp_{info['chain']}" in E:
                m = s & (E[f"tp_{info['chain']}"] == 1) & (E["leadbl"] > info["lead_min"]) & (E["npu"] == 0)
                R.setdefault(f"counts_{run}", {})[sc] = {"events": int(m.sum()), "sum_PS": float(E[f"ps_{info['chain']}"][m].sum()) if f"ps_{info['chain']}" in E else None}
    json.dump(R, open(a.out, "w"), indent=1)
    print(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
