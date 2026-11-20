"""PF etap 2: udział par same-event ze wspólnym śladem per kanał, z kluczem śladu = bity (d0, z0SinThetaPV) [niezależne od refitu]
oraz kluczem bitowym (pt,eta) [stary build._trk_id] i tolerancją 1%. Kanały: D0D0 (10,10), DpDp (4,4), DpD0 (4,10) - łączy dwa npz po (run,event).
Regiony: SR (|m-m0|<window) / SB. Wagi: brak (zliczenia)."""
import sys, json, numpy as np
M0 = {10: 1.86484, 4: 1.86966}; WIN = {10: 0.055, 4: 0.045}; SBMIN = {10: 0.075, 4: 0.070}
def load(f, dtsel, maxrows):
    z = np.load(f, mmap_mode="r"); n = min(maxrows, len(z["event_number"]))
    dt = np.asarray(z["Dmeson_decayType"][:n]); sel = np.nonzero(dt == dtsel)[0]
    g = lambda k: np.asarray(z[k][:n])[sel]
    key = (g("run_number").astype(np.int64) << 40) ^ g("event_number").astype(np.int64)
    out = {"key": key, "m": g("Dmeson_mass").astype(float), "dt": np.full(len(sel), dtsel)}
    for i in range(3):
        for v in ("pt","eta","phi","d0","z0SinThetaPV"): out[(i,v)] = g(f"Dmeson_daughter{i}_{v}")
    return out
def pairs(a, b, same):
    # indeksy par po wspólnym kluczu zdarzenia
    ka, kb = a["key"], b["key"]; oa = np.argsort(ka, kind="stable"); ob = np.argsort(kb, kind="stable")
    ka, kb = ka[oa], kb[ob]; ia=[]; ib=[]
    ua, sa, ca = np.unique(ka, return_index=True, return_counts=True)
    ub, sb_, cb = np.unique(kb, return_index=True, return_counts=True)
    common, pa, pb = np.intersect1d(ua, ub, return_indices=True)
    for x, y in zip(pa, pb):
        na, nb = ca[x], cb[y]
        if na*nb > 4000: continue
        A = oa[sa[x]:sa[x]+na]; B = ob[sb_[y]:sb_[y]+nb]
        if same:
            if na < 2: continue
            i, j = np.triu_indices(na, 1); ia.append(A[i]); ib.append(A[j])
        else:
            ia.append(np.repeat(A, nb)); ib.append(np.tile(B, na))
    return np.concatenate(ia), np.concatenate(ib)
def shared(a, b, ia, ib):
    ex_d0 = np.zeros(len(ia), bool); ex_pt = np.zeros(len(ia), bool); tol = np.zeros(len(ia), bool); nsh = np.zeros(len(ia), int)
    for i in range(3):
        pa = a[(i,"pt")][ia]
        for j in range(3):
            pb = b[(j,"pt")][ib]; ok = (pa>0)&(pb>0)
            e = ok & (a[(i,"d0")][ia]==b[(j,"d0")][ib]) & (a[(i,"z0SinThetaPV")][ia]==b[(j,"z0SinThetaPV")][ib])
            ex_d0 |= e; nsh += e
            ex_pt |= ok & (pa==pb) & (a[(i,"eta")][ia]==b[(j,"eta")][ib])
            dp = np.abs(a[(i,"phi")][ia]-b[(j,"phi")][ib]); dp = np.minimum(dp, 2*np.pi-dp)
            tol |= ok & (np.abs(pa-pb)/np.maximum(pa,1e-6) < 0.01) & (np.abs(a[(i,"eta")][ia]-b[(j,"eta")][ib])<1e-3) & (dp<1e-3)
    return ex_d0, ex_pt, tol, nsh
def region(m, dt):
    d = np.abs(m - M0[dt]); return np.where(d < WIN[dt], "SR", np.where(d > SBMIN[dt], "SB", "gap"))
if __name__ == "__main__":
    fD0, fDp, tag, out, maxrows = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
    res = {"tag": tag}
    d0 = load(fD0, 10, maxrows) if fD0 != "-" else None
    dp = load(fDp, 4, maxrows) if fDp != "-" else None
    for ch, (a, b, same) in {"D0D0": (d0, d0, True), "DpDp": (dp, dp, True), "DpD0": (dp, d0, False)}.items():
        if a is None or b is None: continue
        ia, ib = pairs(a, b, same); ex_d0, ex_pt, tol, nsh = shared(a, b, ia, ib)
        ra, rb = region(a["m"][ia], a["dt"][0]), region(b["m"][ib], b["dt"][0])
        r = {"n_pairs": int(len(ia)), "frac_shared_d0z0": float(ex_d0.mean()), "frac_shared_old_bitwise": float(ex_pt.mean()), "frac_tol1pct": float(tol.mean()),
             "frac_share_2tracks": float((nsh>=2).mean()), "by_region": {}}
        for reg in ("SRxSR", "SRxSB", "SBxSB"):
            x, y = reg.split("x"); msk = ((ra==x)&(rb==y)) | ((ra==y)&(rb==x))
            r["by_region"][reg] = {"n": int(msk.sum()), "frac_shared_d0z0": float(ex_d0[msk].mean()) if msk.any() else None}
        res[ch] = r
    json.dump(res, open(out, "w"), indent=1); print(json.dumps(res, indent=1))
