"""NDS etap 1: dlaczego fit_single nie zbiega na pełnej lumi (norm_all z V dstar_nom.npz)."""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])  # AA_merge/src
from charmtag import dstar_refit as R
from iminuit import Minuit
z = np.load(sys.argv[2]); e = z["edges"]; h = z["norm_all"].astype(float)
shape = json.load(open(sys.argv[3]))["norm"]["shape_mc"]
print("shape", shape, "entries [LO,HI]:", R._rebin(h, e)[0].sum())
hh, ee = R._rebin(h, e)
def cdf_s(xe, mu, s):
    c = R._sig_cdf(xe, mu, shape["sc"] * s, shape["k"], shape["f"]); return (c - c[0]) / (c[-1] - c[0])
def nll(N, mu, s, nb, a, b, c1, c2):
    lam = np.clip(np.diff(N * cdf_s(ee, mu, s) + nb * R._bkg_cdf(ee, a, b, c1, c2)), 1e-9, None)
    return float(np.sum(lam - hh * np.log(lam)))
for b0 in (0.05, 0.15, 0.4):
    m = Minuit(nll, N=max(0.01 * hh.sum(), 50), mu=shape["mu"], s=1.0, nb=hh.sum(), a=0.5, b=b0, c1=0.0, c2=0.0); m.errordef = 0.5
    m.limits["N"] = (0, None); m.limits["s"] = (0.8, 1.6); m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
    m.limits["a"] = (0.01, 3); m.limits["b"] = (-1, 3); m.limits["c1"] = (-0.5, 2); m.limits["c2"] = (-0.05, 0.2); m.limits["nb"] = (0, None)
    m.migrad(ncall=20000); m.hesse()
    print(f"\n== b0={b0}: valid={m.valid} fval={m.fval:.1f} edm={m.fmin.edm:.2e} has_accurate_covar={m.fmin.has_accurate_covar} has_posdef_covar={m.fmin.has_posdef_covar} at_limit={m.fmin.has_parameters_at_limit} reached_call_limit={m.fmin.has_reached_call_limit} above_max_edm={m.fmin.is_above_max_edm}")
    for p in m.parameters: print(f"   {p:3s} = {m.values[p]:12.5g} ± {m.errors[p]:10.3g}  lim {m.limits[p]}  fixed={m.fixed[p]}")
