#!/usr/bin/env python3
"""A3_KOREKTA (20.11.2026): korekta prawdy A3 (mode = ekskluzywny stan końcowy, pythia_dd.cc modeOf) do definicji mc23
(TruthDmesons_decayType topologiczny). Czynniki excl/topo z projektu B (CVER_2026-11-20/RAPORT.md): D0 0,7342±0,0024, D+ 0,9387±0,0013,
D* 0,72–0,734 (centr. 0,727). topo = excl / f. [pary/D*]: D0D0 × f*/f0², DpD0 × f*/(f0 fp), DpDp × f*/fp². Bez danych."""
import json, os
H = os.path.dirname(os.path.abspath(__file__))
T = json.load(open(f"{H}/out/truth_npair.json")); M = json.load(open(f"{H}/out/truth_mbnd.json"))
f0, fp = 0.7342, 0.9387
out = {"f0": f0, "fp": fp, "variants": {}}
for tag, fs in (("central", 0.727), ("fstar_0.720", 0.720), ("fstar_0.734", 0.734)):
    k = {"D0D0": fs / f0**2, "DpD0": fs / (f0 * fp), "DpDp": fs / fp**2}
    fA3 = {ch: T["truth_over_direct"][ch] * k[ch] for ch in k}
    b7 = T["B7"]["mb_nd_dstar_per_event"] / (T["B7"]["A3_dstar_per_ND_event"] / fs)
    out["variants"][tag] = {"fstar": fs, "k_pair_per_dstar": k, "f_A3_kor": fA3,
                            "f_mbnd": M["truth_self_over_direct"], "A3kor_over_mbnd": {ch: fA3[ch] / M["truth_self_over_direct"][ch] for ch in k},
                            "B7_ratio_mb_over_A3kor": b7, "B7_pass_15pct": abs(b7 - 1) < 0.15,
                            "Dp_over_Dstar_A3kor": 5.15 * fs / fp, "Dp_over_Dstar_mbnd": 3.63,
                            "D0_over_Dstar_A3kor": 3.80 * fs / f0, "D0_over_Dstar_mbnd": 3.77}
json.dump(out, open(f"{H}/out/a3_korekta.json", "w"), indent=1)
for t, v in out["variants"].items():
    print(t, {c: round(x, 4) for c, x in v["f_A3_kor"].items()}, "A3kor/mb", {c: round(x, 3) for c, x in v["A3kor_over_mbnd"].items()},
          "B7", round(v["B7_ratio_mb_over_A3kor"], 3), v["B7_pass_15pct"], "D+/D*", round(v["Dp_over_Dstar_A3kor"], 2), "D0/D*", round(v["D0_over_Dstar_A3kor"], 2))
