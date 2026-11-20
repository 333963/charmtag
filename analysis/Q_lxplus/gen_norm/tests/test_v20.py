#!/usr/bin/env python3
"""Test v2.0 (17.11.2026) bez lxplus: (1) ref_by_anchors — referencja w2 po kotwicach (h5 v1.1) w plikach o innym podziale,
(2) gn_root --kind data na atrapie v2.0 (MULT_charm/tests/local/synth) → kolumny nselbl/leadbl/npu/grl/run/lb/tp_/ps_,
(3) gn_root --kind mc: okno D0 1,75–2,00 nie zmienia wyniku na v1.1 (lokalny plik mb_nd v1.1, porcja zdarzeń).
Uruchom: .venv_h5/bin/python tests/test_v20.py [--skip-root]"""
import glob
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

G = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, G)
import gn_combine as C  # noqa: E402

CH = os.path.abspath(os.path.join(G, "..", "..", ".."))


def test_anchors():
    rng = np.random.default_rng(5)
    # 3 „AOD” po 1000 zdarzeń; referencja = AOD 1 (cały) ; pliki v2.0: [AOD0 | pierwsza połowa AOD1], [druga połowa AOD1 | AOD2]
    aod = [np.arange(k * 10**6, k * 10**6 + 1000, dtype=np.uint64) for k in range(3)]
    has_cand = [rng.random(1000) < 0.7 for _ in range(3)]
    A = aod[1][has_cand[1]]
    f1 = np.r_[aod[0], aod[1][:500]]; f2 = np.r_[aod[1][500:], aod[2]]
    DA = [("data_428770_000001", {}, {"ev_evnum": f1}), ("data_428770_000002", {}, {"ev_evnum": f2})]
    masks, info = C.ref_by_anchors(DA, A, gap=50)
    got = np.r_[f1[masks[0]], f2[masks[1]]]
    assert info["anchors_found_unique"] == len(A), info
    # brzegi bloku wewnątrz pliku mogą zgubić zdarzenia bez kandydata (≤ kilka) — dokładnie jak w skanie 17.11
    lost = len(np.setdiff1d(aod[1], got)); extra = len(np.setdiff1d(got, aod[1]))
    assert extra == 0 and lost <= 6, (lost, extra)
    print(f"ref_by_anchors OK: {info['ref_events']} / 1000 zdarzeń AOD ref (zgubione na brzegach {lost}, obce {extra})")


def test_root():
    py = sys.executable
    syn = sorted(glob.glob(os.path.join(CH, "charm_work/MULT_charm/tests/local/synth/data_*.root")))
    v11 = sorted(glob.glob(os.path.join(CH, "mc23_minbiasND_Dpt3to150GeV_v9/*.root")))
    with tempfile.TemporaryDirectory() as d:
        if syn:
            o = os.path.join(d, "da.npz")
            subprocess.run([py, os.path.join(G, "gn_root.py"), "--kind", "data", "--max-events", "3000", syn[0], o], check=True, capture_output=True)
            z = np.load(o); meta = json.loads(str(z["meta"]))
            need = {"ev_evnum", "ev_n_trk", "ev_lead_pt", "ev_nselbl", "ev_leadbl", "ev_npu", "ev_grl", "ev_run", "ev_lb"}
            assert need <= set(z.files), need - set(z.files)
            tp = [k for k in z.files if k.startswith("ev_tp_")]; ps = [k for k in z.files if k.startswith("ev_ps_")]
            assert tp and ps and meta["fmt"] == "v2.0", (tp, ps, meta["fmt"])
            assert float(z["ev_leadbl"].max()) < 1000, "leadbl nie w GeV"
            print(f"gn_root data v2.0 OK: {len(z['ev_evnum'])} zdarzeń, {len(tp)} tp_, {len(ps)} ps_, leadbl max {z['ev_leadbl'].max():.1f} GeV")
        else:
            print("POMINIĘTE: brak atrapy v2.0 (MULT_charm/tests/local/synth)")
        if v11:
            o = os.path.join(d, "mc.npz")
            subprocess.run([py, os.path.join(G, "gn_root.py"), "--kind", "mc", "--max-events", "20000", v11[0], o], check=True, capture_output=True)
            z = np.load(o); meta = json.loads(str(z["meta"]))
            # v1.1: wszystkie poprawne D0 mają m w 1,75–2,00 → okno nic nie zmienia
            m = z["c_m"][z["c_decay_type"] == 10]
            assert meta["fmt"] == "v1.1" and ((m > 1750) & (m < 2000)).all()
            assert int(z["ev_n_corr_d0"].sum()) == int((z["c_decay_type"] == 10).sum())
            print(f"gn_root mc v1.1 OK: {meta['n_events']} zdarzeń, poprawne D0 {int(z['ev_n_corr_d0'].sum())} (okno bez wpływu), extra {meta['extra_columns']}")
        else:
            print("POMINIĘTE: brak lokalnego pliku v1.1 mb_nd")


if __name__ == "__main__":
    test_anchors()
    if "--skip-root" not in sys.argv:
        test_root()
    print("TEST v2.0 OK")
