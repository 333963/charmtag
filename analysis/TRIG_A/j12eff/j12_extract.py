#!/usr/bin/env python3
"""TRIG_A/j12eff — wydajność triggera L1_J12 z DANYCH w funkcji zmiennych torowych (projekt A; tylko zdarzenia, bez par D–D).

Dla każdego zdarzenia ROOT v2.0 (drzewa CharmAnalysis + Tracks, wpis i ↔ wpis i): wiodący dżet torowy anty-kT R=0,4 (tory pT>0,5 GeV, |η|<2,5;
jądro tj_kernel.antikt z paczki trkjet), maksimum Σ pT torów w oknie 0,8×0,8 (odpowiednik naładowany okna L1), nSelBL, pT toru wiodącego,
bity triggera (po preskali) + trigTBP L1_J12 (przed preskalą, gdy gałąź istnieje), flaga „zdarzenie z kandydatem D* w oknie sygnału”.
Wyjście: HISTOGRAMY (nie tabela zdarzeń) per zbiór odniesienia × {wszystkie, J12} × {wszystkie zdarzenia, z D*}: npz.
Zbiory odniesienia (maski zdarzeń, nieobciążone przez J12): mbts1 (HLT_mb_mbts_L1MBTS_1 ∧ GRL), sptrk4 (HLT_mb_sptrk_pt4_L1MBTS_2 ∧ GRL),
sptrk2, mbany (dowolny łańcuch oparty na MBTS bez j20), all (wszystko, ∧ GRL).
Użycie: j12_extract.py --out OUT.npz [--max-events N] PLIK.root [...]
"""
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tj_kernel as K  # noqa: E402
from trg_root_extract import resolve  # noqa: E402
from btrig_common import find_col  # noqa: E402

njit = K.njit
E_CJ = np.concatenate([np.arange(0, 20, 1.0), np.arange(20, 40, 2.5), [40, 50, 70, 1e9]])      # GeV
E_CW = E_CJ.copy()
E_LD = np.concatenate([np.arange(0, 10, 0.5), np.arange(10, 20, 2.0), [20, 30, 1e9]])
E_NS = np.array([0, 10, 20, 30, 40, 50, 60, 80, 100, 150, 1e9])
REFS = ["mbts1", "sptrk4", "sptrk2", "mbany", "all"]
NETA, NPHI = 25, 32


@njit(cache=False)
def per_event(tpt, teta, tphi, off, R, ptmin, etamax):
    n = len(off) - 1
    cj = np.zeros(n); cw = np.zeros(n); ld = np.zeros(n); nt = np.zeros(n, np.int64)
    grid = np.zeros((NETA, NPHI))
    for e in range(n):
        a, b = off[e], off[e + 1]
        m = 0
        for i in range(a, b):
            if tpt[i] > ptmin and abs(teta[i]) < etamax:
                m += 1
        nt[e] = m
        if m == 0:
            continue
        pt = np.empty(m); eta = np.empty(m); phi = np.empty(m); mass = np.full(m, 0.13957); fl = np.zeros(m, np.int64)
        grid[:, :] = 0.0
        k = 0
        for i in range(a, b):
            if tpt[i] > ptmin and abs(teta[i]) < etamax:
                pt[k] = tpt[i]; eta[k] = teta[i]; phi[k] = tphi[i]; k += 1
                if tpt[i] > ld[e]:
                    ld[e] = tpt[i]
                ie = int((teta[i] + 2.5) / 0.2)
                if ie >= NETA:
                    ie = NETA - 1
                ip = int((tphi[i] + np.pi) / (2 * np.pi / NPHI))
                if ip >= NPHI:
                    ip = NPHI - 1
                if ip < 0:
                    ip = 0
                grid[ie, ip] += tpt[i]
        jpt, _je, _jp, _jn, _jf = K.antikt(pt, eta, phi, mass, fl, R)
        cj[e] = jpt[0]
        best = 0.0
        for i in range(NETA - 3):
            for j in range(NPHI):
                s = 0.0
                for x in range(4):
                    for y in range(4):
                        s += grid[i + x, (j + y) % NPHI]
                if s > best:
                    best = s
        cw[e] = best
    return cj, cw, ld, nt


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--max-events", type=int, default=0); ap.add_argument("--chunk", type=int, default=200000)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)
    import awkward as ak
    import uproot
    shp = {"cj": (len(E_CJ) - 1, len(E_NS) - 1), "cw": (len(E_CW) - 1, len(E_NS) - 1), "ld": (len(E_LD) - 1, len(E_NS) - 1), "cjld": (len(E_CJ) - 1, len(E_LD) - 1)}
    # H[var][ref, j (0 wszystkie, 1 J12 po preskali, 2 J12 TBP), d (0 wszystkie, 1 z D*)]
    H = {v: np.zeros((len(REFS), 3, 2) + s) for v, s in shp.items()}
    stat = {"files": 0, "events": 0, "tbp_branch": None, "missing": [], "ps_j12_sum": 0.0, "ps_j12_n": 0, "ps_mbts1_sum": 0.0, "ps_mbts1_n": 0}
    t0 = time.time()
    for path in a.files:
        with uproot.open(path) as R:
            C, T = R["CharmAnalysis"], R["Tracks"]
            if C.num_entries != T.num_entries:
                raise RuntimeError(f"{path}: Tracks {T.num_entries} ≠ CharmAnalysis {C.num_entries}")
            keys = C.keys()
            m = resolve(keys)
            tbp = [k for k in keys if "trigTBP" in k and k.replace("_NOSYS", "").endswith("L1_J12")]
            stat["tbp_branch"] = tbp[0] if tbp else None
            dsb = {k: find_col(keys, [f"DMesons_{k}"]) for k in ("decayType", "pt", "m", "DeltaMass", "D0Index")}
            want = sorted(set(list(m.values()) + tbp[:1] + [b for b in dsb.values() if b]))
            nmax = C.num_entries if not a.max_events else min(C.num_entries, a.max_events)
            for s in range(0, nmax, a.chunk):
                e = min(s + a.chunk, nmax)
                arr = C.arrays(want, entry_start=s, entry_stop=e, library="ak")
                n = len(arr)
                g = lambda k, d=0: ak.to_numpy(arr[m[k]]) if k in m else np.full(n, d)  # noqa: E731
                tp = lambda c: (g(f"tp_{c}", 0) == 1)  # noqa: E731
                grl = g("grl", 1) == 1
                j12 = tp("HLT_noalg_L1J12")
                j12tbp = (ak.to_numpy(arr[tbp[0]]) == 1) if tbp else np.zeros(n, bool)
                refs = {"mbts1": tp("HLT_mb_mbts_L1MBTS_1") & grl, "sptrk4": tp("HLT_mb_sptrk_pt4_L1MBTS_2") & grl,
                        "sptrk2": tp("HLT_mb_sptrk_pt2_L1MBTS_2") & grl,
                        "mbany": (tp("HLT_mb_mbts_L1MBTS_1") | tp("HLT_mb_sptrk_pt4_L1MBTS_2") | tp("HLT_mb_sptrk_pt2_L1MBTS_2")
                                  | tp("HLT_noalg_mb_L1MBTS_1") | tp("HLT_noalg_mb_L1MBTS_2") | tp("HLT_mb_sptrk_L1MBTS_1")) & grl,
                        "all": grl}
                if "ps_HLT_noalg_L1J12" in m:
                    p = g("ps_HLT_noalg_L1J12"); ok = p > 0; stat["ps_j12_sum"] += float(p[ok].sum()); stat["ps_j12_n"] += int(ok.sum())
                if "ps_HLT_mb_mbts_L1MBTS_1" in m:
                    p = g("ps_HLT_mb_mbts_L1MBTS_1"); ok = p > 0; stat["ps_mbts1_sum"] += float(p[ok].sum()); stat["ps_mbts1_n"] += int(ok.sum())
                hasd = np.zeros(n, bool)
                if all(dsb.values()):
                    dm = arr[dsb["DeltaMass"]] / 1e3
                    md0 = (arr[dsb["m"]] - arr[dsb["DeltaMass"]]) / 1e3
                    okd = (arr[dsb["decayType"]] == 1) & (arr[dsb["D0Index"]] >= 0) & (abs(md0 - 1.86484) < 0.055) & (arr[dsb["pt"]] / 1e3 >= 3) \
                        & (dm > 0.1440) & (dm < 0.1470)
                    hasd = ak.to_numpy(ak.sum(okd, axis=1)) > 0
                tr = T.arrays(["Tracks_pt_NOSYS", "Tracks_eta_NOSYS", "Tracks_phi_NOSYS"], entry_start=s, entry_stop=e)
                cnt = ak.to_numpy(ak.num(tr["Tracks_pt_NOSYS"]))
                off = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
                tpt = ak.to_numpy(ak.flatten(tr["Tracks_pt_NOSYS"])).astype(np.float64) * 1e-3
                teta = ak.to_numpy(ak.flatten(tr["Tracks_eta_NOSYS"])).astype(np.float64)
                tphi = ak.to_numpy(ak.flatten(tr["Tracks_phi_NOSYS"])).astype(np.float64)
                cj, cw, ld, _nt = per_event(tpt, teta, tphi, off, 0.4, 0.5, 2.5)
                ns = g("nSelBL", 0).astype(float)
                X = {"cj": (cj, ns, E_CJ, E_NS), "cw": (cw, ns, E_CW, E_NS), "ld": (ld, ns, E_LD, E_NS), "cjld": (cj, ld, E_CJ, E_LD)}
                for ri, rn in enumerate(REFS):
                    for ji, jm in enumerate((np.ones(n, bool), j12, j12tbp)):
                        for di, dmk in enumerate((np.ones(n, bool), hasd)):
                            sel = refs[rn] & jm & dmk
                            if not sel.any():
                                continue
                            for v, (x, y, ex, ey) in X.items():
                                H[v][ri, ji, di] += np.histogram2d(x[sel], y[sel], bins=(ex, ey))[0]
                stat["events"] += n
                stat["missing"] = sorted(set(stat["missing"]) | {c for c in ("HLT_noalg_L1J12", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_pt4_L1MBTS_2") if f"tp_{c}" not in m})
        stat["files"] += 1
        print(f"[{stat['files']}/{len(a.files)}] {os.path.basename(path)} razem {stat['events']} zd. {time.time() - t0:.0f}s", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, meta=np.array(json.dumps(stat)), **{"H_" + v: h for v, h in H.items()})
    print("zapisano", a.out, json.dumps(stat))


if __name__ == "__main__":
    sys.exit(main())
