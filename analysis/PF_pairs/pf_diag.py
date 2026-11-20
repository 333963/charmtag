"""PF etap 1: diagnoza wspólnych śladów między kandydatami tego samego zdarzenia.
Dla każdego zdarzenia: wszystkie pary (kandydat a, kandydat b), wszystkie pary córek (i,j):
odległość w (pt,eta,phi) oraz zgodność bitowa d0/z0/dEdx. Wyjście: histogramy + udziały."""
import sys, json, numpy as np
f = sys.argv[1]; out = sys.argv[2]; maxrows = int(sys.argv[3]) if len(sys.argv) > 3 else 400000
z = np.load(f, mmap_mode="r")
n = min(maxrows, len(z["event_number"]))
ev = np.asarray(z["event_number"][:n]).astype(np.int64); run = np.asarray(z["run_number"][:n]).astype(np.int64)
key = (run << 40) ^ ev
dt = np.asarray(z["Dmeson_decayType"][:n]); m = np.asarray(z["Dmeson_mass"][:n]).astype(float)
D = {}
for i in range(3):
    for v in ("pt","eta","phi","d0","z0SinThetaPV","dEdx","npts_dEdx","z0SinTheta"):
        D[(i,v)] = np.asarray(z[f"Dmeson_daughter{i}_{v}"][:n])
# sortuj po kluczu
o = np.argsort(key, kind="stable"); key = key[o]; dt = dt[o]; m = m[o]
for k in D: D[k] = D[k][o]
starts = np.r_[0, np.nonzero(np.diff(key))[0] + 1, len(key)]
rel = []; deta=[]; dphi=[]; exact_d0=[]; exact_z0=[]; exact_pt=[]; pairs_tot=0; cand_in_ev=[]
A=[];B=[]
for s,e in zip(starts[:-1], starts[1:]):
    k = e-s; cand_in_ev.append(k)
    if k < 2 or k > 60: continue
    ia, ib = np.triu_indices(k, 1); A.append(ia+s); B.append(ib+s)
A = np.concatenate(A); B = np.concatenate(B); pairs_tot = len(A)
# min po córkach (i,j) z ważnymi pt>0
best_rel = np.full(pairs_tot, np.inf); best_deta = np.full(pairs_tot, np.inf); best_dphi = np.full(pairs_tot, np.inf)
any_exact_d0 = np.zeros(pairs_tot, bool); any_exact_pt = np.zeros(pairs_tot, bool); any_exact_z0 = np.zeros(pairs_tot, bool)
any_tol = np.zeros(pairs_tot, bool)
for i in range(3):
    pa = D[(i,"pt")][A]; 
    for j in range(3):
        pb = D[(j,"pt")][B]; ok = (pa > 0) & (pb > 0)
        r = np.abs(pa-pb)/np.maximum(pa,1e-6); de = np.abs(D[(i,"eta")][A]-D[(j,"eta")][B]); dp = np.abs(D[(i,"phi")][A]-D[(j,"phi")][B]); dp = np.minimum(dp, 2*np.pi-dp)
        d = np.where(ok, r + de + dp, np.inf)
        better = d < (best_rel + best_deta + best_dphi)
        best_rel = np.where(better, r, best_rel); best_deta = np.where(better, de, best_deta); best_dphi = np.where(better, dp, best_dphi)
        any_exact_d0 |= ok & (D[(i,"d0")][A] == D[(j,"d0")][B]) & (D[(i,"z0SinThetaPV")][A] == D[(j,"z0SinThetaPV")][B])
        any_exact_z0 |= ok & (D[(i,"z0SinTheta")][A] == D[(j,"z0SinTheta")][B])
        any_exact_pt |= ok & (pa == pb) & (D[(i,"eta")][A] == D[(j,"eta")][B])
        any_tol |= ok & (r < 0.01) & (de < 1e-3) & (dp < 1e-3)
res = {"file": f, "rows": int(n), "events": int(len(starts)-1), "pairs_same_event": int(pairs_tot),
       "frac_exact_pt_eta": float(any_exact_pt.mean()), "frac_exact_d0_z0PV": float(any_exact_d0.mean()), "frac_exact_z0": float(any_exact_z0.mean()),
       "frac_tol_1pct_1e-3": float(any_tol.mean()),
       "frac_rel_lt": {str(t): float((best_rel < t).mean()) for t in (1e-4,1e-3,3e-3,1e-2,3e-2,0.1)},
       "frac_deta_lt": {str(t): float((best_deta < t).mean()) for t in (1e-5,1e-4,1e-3,1e-2)},
       "frac_dphi_lt": {str(t): float((best_dphi < t).mean()) for t in (1e-5,1e-4,1e-3,1e-2)},
       "hist_log10_rel": np.histogram(np.log10(np.clip(best_rel[np.isfinite(best_rel)],1e-7,10)), bins=np.arange(-7,1.01,0.25))[0].tolist(),
       "hist_log10_deta": np.histogram(np.log10(np.clip(best_deta[np.isfinite(best_deta)],1e-8,10)), bins=np.arange(-8,1.01,0.25))[0].tolist(),
       "n_cand_per_event_pct": np.percentile(cand_in_ev,[50,90,99,100]).tolist(),
       # zgodność: tolerancja vs exact d0
       "tol_and_exactd0": float((any_tol & any_exact_d0).mean()), "tol_not_exactd0": float((any_tol & ~any_exact_d0).mean()), "exactd0_not_tol": float((~any_tol & any_exact_d0).mean())}
json.dump(res, open(out,"w"), indent=1); print(json.dumps({k:v for k,v in res.items() if not k.startswith("hist")}, indent=1))
