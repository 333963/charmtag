"""N2 etap 1: D0 z tagu D* (½ lumi ∩ podpróbki S) → histogramy Δm pass/fail per pT × WP, dane i MC.
Tag: D* (decayType 1) ↔ D0 (decayType 10) w tym samym zdarzeniu po bitach (d0, z0sinθ_PV) obu córek (PF; D0 w D* nie jest refitowany).
Score D0 z L_charmtag/runs/v11_full/scores/D0_H/<sample>.h5 (row = indeks w pliku źródłowym = _row w npz)."""
import os, sys, glob, json, numpy as np, h5py
W = os.path.expanduser("~/mnt/charm/charm_work"); OUT = os.path.join(W, "N2_sf/out"); os.makedirs(OUT, exist_ok=True)
SUB = {"data_none": "S_bkggen/sub/data_p*.npz", "signal_mc": "S_bkggen/sub/fd_p*.npz", "incl_mc": "S_bkggen/sub/mb_p*.npz"}
SC = os.path.join(W, "L_charmtag/runs/v11_full/scores/D0_H")
META = json.load(open(os.path.join(W, "L_charmtag/runs/v11_full/models/D0_H/meta.json")))["thresholds"]
PAIR_THR = json.load(open(os.path.join(W, "L_charmtag/runs/v11_full/pairs/pair_inputs.json")))["meta"]["thresholds"]["10"]
PT = [(3, 5), (5, 8), (8, 30)]; M0 = 1.86484; WIN = 0.055
EDGES = np.arange(139.57, 165.0001, 0.1)
def bits(a): return np.asarray(a, np.float32).view(np.uint32).astype(np.int64)
def key(s, i):
    pt = np.asarray(s[f"Dmeson_daughter{i}_pt"]); k = (bits(s[f"Dmeson_daughter{i}_d0"]) << 31) ^ bits(s[f"Dmeson_daughter{i}_z0SinThetaPV"])
    return np.where(pt > 0, k, -1)
def thr_for(wp, pt):
    """próg score: WP 80/70/50 per bin (meta), 65/90 globalny (pary)."""
    if wp in ("0.65", "0.9"): return np.full(len(pt), PAIR_THR[wp])
    t = np.zeros(len(pt))
    for (lo, hi) in PT:
        sel = (pt >= lo) & (pt < hi); t[sel] = META[f"{lo}-{hi}"][f"{float(wp):.2f}"]
    return t
WPS = ["0.80", "0.70", "0.50", "0.65", "0.9"]
def run(sample):
    f = h5py.File(os.path.join(SC, f"{sample}.h5"), "r"); srow = f["row"][:]; sscore = f["score"][:]; o = np.argsort(srow); srow, sscore = srow[o], sscore[o]
    H = {}; n_tag = 0; n_ds = 0; truthinfo = []
    for fn in sorted(glob.glob(os.path.join(W, SUB[sample]))):
        s = np.load(fn); dt = np.asarray(s["Dmeson_decayType"]); ev = (np.asarray(s["run_number"]).astype(np.int64) << 40) ^ np.asarray(s["event_number"]).astype(np.int64)
        k0, k1 = key(s, 0), key(s, 1); m = np.asarray(s["Dmeson_mass"]).astype(float); dm = np.asarray(s["Dmeson_delta_mass"]).astype(float)
        pt = np.asarray(s["Dmeson_pt"]).astype(float); row = np.asarray(s["_row"]).astype(np.int64); tr = np.asarray(s["Dmeson_truth_category"]).astype(float)
        ds = np.nonzero((dt == 1) & (dm > 0.139) & (dm < 0.165))[0]; d0 = np.nonzero(dt == 10)[0]; n_ds += len(ds)
        # słownik D0: (ev, min(k0,k1), max(k0,k1)) -> indeks
        lo = np.minimum(k0[d0], k1[d0]); hi = np.maximum(k0[d0], k1[d0])
        d0key = ev[d0] * 1000003 ^ (lo * 31 + hi)  # hash; kolizje pomijalne
        order = np.argsort(d0key); d0key_s = d0key[order]; d0_s = d0[order]
        lo2 = np.minimum(k0[ds], k1[ds]); hi2 = np.maximum(k0[ds], k1[ds]); dskey = ev[ds] * 1000003 ^ (lo2 * 31 + hi2)
        pos = np.searchsorted(d0key_s, dskey); pos = np.clip(pos, 0, len(d0key_s) - 1); hit = d0key_s[pos] == dskey
        ids, i0 = ds[hit], d0_s[pos[hit]]
        # masa D0 kandydata w oknie + spójność m(D*)−Δm
        ok = (np.abs(m[i0] - M0) < WIN) & (np.abs((m[ids] - dm[ids]) - m[i0]) < 0.0015)
        ids, i0 = ids[ok], i0[ok]
        # score po row
        p = np.searchsorted(srow, row[i0]); p = np.clip(p, 0, len(srow) - 1); hs = srow[p] == row[i0]
        ids, i0, sc = ids[hs], i0[hs], sscore[p[hs]]
        n_tag += len(ids)
        dmm = dm[ids] * 1000.0; ptd = pt[i0]; trd = tr[i0]
        for (lo_, hi_) in PT:
            b = (ptd >= lo_) & (ptd < hi_); name = f"{lo_}-{hi_}"
            for wp in WPS:
                pas = sc >= thr_for(wp, ptd)
                for lab, msk in (("pass", b & pas), ("fail", b & ~pas)):
                    h, _ = np.histogram(dmm[msk], bins=EDGES); H[f"h_{name}_{wp}_{lab}"] = H.get(f"h_{name}_{wp}_{lab}", 0) + h
                    if sample != "data_none":  # prawda: D0 truth 0/1
                        ht, _ = np.histogram(dmm[msk & np.isin(trd, [0, 1])], bins=EDGES); H[f"t_{name}_{wp}_{lab}"] = H.get(f"t_{name}_{wp}_{lab}", 0) + ht
        if sample != "data_none":
            truthinfo.append(np.isin(trd, [0, 1]).mean())
        print(sample, os.path.basename(fn), "D*", len(ds), "tag", len(ids))
    H["edges"] = EDGES; H["n_tag"] = np.array([n_tag]); H["n_dstar"] = np.array([n_ds])
    np.savez_compressed(os.path.join(OUT, f"dstar_hist_{sample}.npz"), **H)
    print(sample, "razem tagów", n_tag, "z D*", n_ds, "truth frac", np.mean(truthinfo) if truthinfo else None)
if __name__ == "__main__":
    run(sys.argv[1])
