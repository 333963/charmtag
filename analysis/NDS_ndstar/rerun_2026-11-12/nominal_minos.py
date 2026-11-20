"""rerun 2026-11-12: nominał N_D* (threshold_pol2, 2G kształt V MC, bin 0,05, [m_π+0,9; 160]) — HESSE, MINOS i skan profilu ΔNLL=0,5.
Ten sam model i starty co dstar_refit.fit_single (charmtag z 10.11, NLL względem modelu nasyconego). Kontrola: wynik N == fit_single.
Użycie: PYTHONPATH=REPO/CHARM_git python nominal_minos.py dstar_nom.npz dstar_post.json OUT.json"""
import sys, json, numpy as np
from iminuit import Minuit
from charmtag import dstar_refit as R

npz, post, out = sys.argv[1:4]
z = np.load(npz); e0 = z["edges"]; shape = json.load(open(post))["norm"]["shape_mc"]
h, e = R._rebin(z["norm_all"].astype(float), e0)
hs = np.where(h > 0, h, 1.0)
def cdf_s(xe, mu, s):
    c = R._sig_cdf(xe, mu, shape["sc"] * s, shape["k"], shape["f"]); return (c - c[0]) / (c[-1] - c[0])
def nll(N, mu, s, nb, a, b, c1, c2):
    lam = np.clip(np.diff(N * cdf_s(e, mu, s) + nb * R._bkg_cdf_t(e, a, b, c1, c2, False)), 1e-9, None)
    return float(np.sum(lam - h - h * np.log(lam / hs)))
best = None
for a0 in (0.5, 0.6, 0.4):
    m = Minuit(nll, N=max(0.01 * h.sum(), 50), mu=shape["mu"], s=1.0, nb=h.sum(), a=a0, b=0.0, c1=0.0, c2=0.0)
    m.errordef = 0.5; m.strategy = 2
    m.limits["N"] = (0, None); m.limits["s"] = (0.8, 1.6); m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
    m.limits["a"] = (0.01, 3); m.limits["b"] = (-0.5, 1.5); m.limits["c1"] = (-5, 20); m.limits["c2"] = (-10, 40); m.limits["nb"] = (0, None)
    m.fixed["b"] = True
    m.simplex(); m.migrad(ncall=50000); m.hesse()
    if best is None or (m.valid and (not best.valid or m.fval < best.fval)): best = m
    if best.valid: break
m = best
res = {"N": m.values["N"], "hesse": m.errors["N"], "valid": m.valid, "edm": m.fmin.edm, "fval": m.fval,
       "params": {k: m.values[k] for k in m.parameters}, "params_at_limit": [p.name for p in m.params if p.has_limits and ((p.lower_limit is not None and abs(p.value - p.lower_limit) < 1e-3*max(1,abs(p.error))) or (p.upper_limit is not None and abs(p.value - p.upper_limit) < 1e-3*max(1,abs(p.error))))]}
print("fit", res, flush=True)
m.minos("N")
me = m.merrors["N"]
res["minos"] = {"lower": me.lower, "upper": me.upper, "is_valid": me.is_valid, "lower_valid": me.lower_valid, "upper_valid": me.upper_valid,
                "at_lower_limit": me.at_lower_limit, "at_upper_limit": me.at_upper_limit}
print("minos", res["minos"], flush=True)
json.dump(res, open(out, "w"), indent=1, default=float)
# niezależny skan profilu: N stałe, reszta minimalizowana od minimum; interpolacja ΔNLL=0,5
N0, f0 = res["N"], res["fval"]; v0 = dict(res["params"])
grid = N0 + np.linspace(-1.6, 1.6, 33) * res["hesse"]
prof = []
for Nf in grid:
    mm = Minuit(nll, **dict(v0, N=Nf)); mm.errordef = 0.5; mm.strategy = 1
    for k in ("s", "mu", "a", "b", "c1", "c2", "nb"): mm.limits[k] = m.limits[k]
    mm.fixed["N"] = True; mm.fixed["b"] = True
    mm.migrad(ncall=50000)
    prof.append((float(Nf), float(mm.fval - f0), bool(mm.valid)))
    print(f"N={Nf:7.1f} dNLL={mm.fval - f0:7.4f} valid={mm.valid}", flush=True)
P = np.array([(a, b) for a, b, _ in prof]); i0 = int(np.argmin(P[:, 1]))
lo = np.interp(0.5, P[:i0 + 1, 1][::-1], P[:i0 + 1, 0][::-1]); hi = np.interp(0.5, P[i0:, 1], P[i0:, 0])
res["profile_scan"] = {"points": prof, "lower": float(lo - N0), "upper": float(hi - N0), "min_dnll": float(P[:, 1].min()),
                       "all_valid": all(v for *_, v in prof)}
print("profile", res["profile_scan"]["lower"], res["profile_scan"]["upper"], flush=True)
json.dump(res, open(out, "w"), indent=1, default=float)
