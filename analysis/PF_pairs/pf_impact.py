"""PF etap 3 (ślepo): wpływ weta wspólnych śladów na pary. Nogi z S9 (legs_<sample>.npz: po WP jak pairs/build.py),
ID śladów z podpróbek npz S/S6 trzema kluczami: old (bity pT,η), tol (1 %), d0z0 (bity d0, z0sinθ_PV).
Liczy per kanał i region (SBxSB, arms=SRxSB, centre tylko MC) liczbę par bez/z wetem, corr(m1,m2), a w MC: pary prawdziwe (truth 0/1 obie nogi).
Dodatkowo tag D*: dopasowanie D* (dt 1, Δm w oknie) do nóg D0 kluczem old vs d0z0. SS w danych: tylko zliczenie — NIE zapisujemy niczego o kształcie."""
import os, sys, glob, json, numpy as np
W = os.path.expanduser("~/mnt/charm/charm_work")
S9 = os.path.join(W, "S9_pairs_ad/out"); OUT = os.path.join(W, "PF_pairs/out")
SUB = {10: {"data_none": "S_bkggen/sub/data_p*.npz", "incl_mc": "S_bkggen/sub/mb_p*.npz", "signal_mc": "S_bkggen/sub/fd_p*.npz"},
       4: {"data_none": "S6_kand/sub/dp_data_p*.npz", "incl_mc": "S6_kand/sub/dp_mb_p*.npz", "signal_mc": "S6_kand/sub/dp_fd_p*.npz"}}
M0 = {10: 1864.84, 4: 1869.66}; SBW = {10: 75.0, 4: 70.0}; WIN = {10: 55.0, 4: 45.0}
def bits(a): return np.asarray(a, np.float32).view(np.uint32).astype(np.int64)
def ids(s, i):
    pt = np.asarray(s[f"Dmeson_daughter{i}_pt"]); ok = pt > 0
    old = (bits(pt) << 31) ^ bits(s[f"Dmeson_daughter{i}_eta"])
    new = (bits(s[f"Dmeson_daughter{i}_d0"]) << 31) ^ bits(s[f"Dmeson_daughter{i}_z0SinThetaPV"])
    tol = np.round(np.log(np.clip(pt, 1e-3, None)) * 100).astype(np.int64) * 100000 + np.round(np.asarray(s[f"Dmeson_daughter{i}_eta"]) * 100).astype(np.int64) + 50000
    return [np.where(ok, x, -1) for x in (old, tol, new)]
def attach(sample):
    z = np.load(os.path.join(S9, f"legs_{sample}.npz")); L = {k: z[k] for k in z.files}
    T = {key: {i: np.full(len(L["row"]), -2, np.int64) for i in range(3)} for key in ("old", "tol", "d0z0")}
    dstar = []
    for dt in (10, 4):
        sel = np.nonzero(L["dt"] == dt)[0]
        for f in sorted(glob.glob(os.path.join(W, SUB[dt][sample]))):
            s = np.load(f); r = s["_row"].astype(np.int64); pos = np.clip(np.searchsorted(r, L["row"][sel]), 0, len(r) - 1); hit = r[pos] == L["row"][sel]
            for i in range(3):
                o, t, n = ids(s, i)
                T["old"][i][sel[hit]] = o[pos[hit]]; T["tol"][i][sel[hit]] = t[pos[hit]]; T["d0z0"][i][sel[hit]] = n[pos[hit]]
            if dt == 10:  # kandydaci D* z tej samej podpróbki (tag)
                dtc = np.asarray(s["Dmeson_decayType"]); dm = np.asarray(s["Dmeson_delta_mass"]).astype(float)
                k = np.nonzero((dtc == 1) & (np.abs(dm - 0.14543) < 0.0010))[0]
                if len(k):
                    ev = (np.asarray(s["run_number"]).astype(np.int64)[k] << 40) ^ np.asarray(s["event_number"]).astype(np.int64)[k]
                    o0, t0, n0 = ids(s, 0); o1, t1, n1 = ids(s, 1)
                    dstar.append({"ev": ev, "old": np.stack([o0[k], o1[k]], 1), "d0z0": np.stack([n0[k], n1[k]], 1)})
    keep = T["old"][0] != -2
    L = {k: v[keep] for k, v in L.items()}; T = {key: {i: T[key][i][keep] for i in range(3)} for key in T}
    D = {k: np.concatenate([d[k] for d in dstar]) for k in ("ev", "old", "d0z0")} if dstar else None
    return L, T, D
def region(m, dt):
    d = np.abs(m - M0[dt]); return np.where(d < WIN[dt], 2, np.where(d > SBW[dt], 0, 1))  # 0 SB, 1 gap, 2 SR
def pairs_idx(ev, dt, a, b):
    o = np.argsort(ev, kind="stable"); ev_s = ev[o]; dt_s = dt[o]
    starts = np.r_[0, np.nonzero(np.diff(ev_s))[0] + 1, len(ev_s)]; ia = []; ib = []
    for s, e in zip(starts[:-1], starts[1:]):
        A = o[s:e][dt_s[s:e] == a]; B = o[s:e][dt_s[s:e] == b] if b != a else A
        if a == b:
            if len(A) < 2: continue
            i, j = np.triu_indices(len(A), 1); ia.append(A[i]); ib.append(A[j])
        else:
            if not len(A) or not len(B): continue
            ia.append(np.repeat(A, len(B))); ib.append(np.tile(B, len(A)))
    return (np.concatenate(ia), np.concatenate(ib)) if ia else (np.array([], int), np.array([], int))
def shared(T, key, ia, ib):
    sh = np.zeros(len(ia), bool)
    for i in range(3):
        for j in range(3):
            sh |= (T[key][i][ia] == T[key][j][ib]) & (T[key][i][ia] >= 0)
    return sh
if __name__ == "__main__":
    sample = sys.argv[1]; L, T, D = attach(sample)
    ev = (L["run"] << 40) ^ L["event"]; res = {"sample": sample, "n_legs": {int(d): int((L["dt"] == d).sum()) for d in (10, 4)}}
    for ch, (a, b) in {"D0D0": (10, 10), "DpD0": (4, 10), "DpDp": (4, 4)}.items():
        ia, ib = pairs_idx(ev, L["dt"], a, b)
        if not len(ia): res[ch] = {"n_pairs": 0}; continue
        ra, rb = region(L["m"][ia], a), region(L["m"][ib], b); m1, m2 = L["m"][ia], L["m"][ib]
        regs = {"SBxSB": (ra == 0) & (rb == 0), "arms": ((ra == 2) & (rb == 0)) | ((ra == 0) & (rb == 2))}
        if sample != "data_none": regs["centre"] = (ra == 2) & (rb == 2)
        r = {"n_pairs": int(len(ia))}
        for key in ("old", "tol", "d0z0"):
            sh = shared(T, key, ia, ib); rk = {"frac_shared_all": float(sh.mean())}
            for rn, msk in regs.items():
                n = int(msk.sum()); k = msk & ~sh
                rk[rn] = {"n": n, "n_after_veto": int(k.sum()), "frac_removed": float(sh[msk].mean()) if n else None,
                          "corr_m1m2_before": float(np.corrcoef(m1[msk], m2[msk])[0, 1]) if n > 10 else None,
                          "corr_m1m2_after": float(np.corrcoef(m1[k], m2[k])[0, 1]) if k.sum() > 10 else None}
            if sample != "data_none":
                tr = np.isin(L["truth"][ia], [0, 1]) & np.isin(L["truth"][ib], [0, 1]) & regs.get("centre", False)
                rk["true_pairs_centre"] = int(tr.sum()); rk["true_pairs_flagged_shared"] = int((tr & sh).sum())
                bk = regs.get("centre", False) & ~tr
                rk["centre_bkg_frac_shared"] = float(sh[bk].mean()) if bk.any() else None
            r[key] = rk
        res[ch] = r
    if D is not None:  # tag D*: udział nóg D0 (SR) z dopasowanym D*
        sel = np.nonzero((L["dt"] == 10) & (region(L["m"], 10) == 2))[0]; tag = {}
        for key in ("old", "d0z0"):
            lo = np.minimum(T[key][0][sel], T[key][1][sel]); hi = np.maximum(T[key][0][sel], T[key][1][sel])
            kl = set(zip(ev[sel].tolist(), lo.tolist(), hi.tolist()))
            dl = np.minimum(D[key][:, 0], D[key][:, 1]); dh = np.maximum(D[key][:, 0], D[key][:, 1])
            kd = set(zip(D["ev"].tolist(), dl.tolist(), dh.tolist()))
            tag[key] = {"n_dstar_cand": int(len(D["ev"])), "n_D0_legs_SR": int(len(sel)), "n_tagged": int(len(kl & kd))}
        res["dstar_tag"] = tag
    os.makedirs(OUT, exist_ok=True); json.dump(res, open(os.path.join(OUT, f"impact_{sample}.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))
