import numpy as np, json
rng=np.random.default_rng(1); n=1_000_000
g1,g2,g3=rng.standard_normal((3,n))
def mu_up(mu,st,syst,nd,cl):
    x=(mu+st*g1)*np.exp(np.log1p(syst)*g2)*np.exp(np.log1p(nd)*g3)
    return float(np.quantile(x,cl))
def seff(mu_u,f): return 28*f/(mu_u-1+f) if mu_u>1-f else float('inf')
out={}
for lab,mu,syst in [("0.50 syst sym (stare 0,387)",0.50,0.387),("0.50 syst +0,551 (final)",0.50,0.551),("0.72 skala z MC (final)",0.72,0.551),("0.70 frag ALICE (final)",0.50/0.71,0.551)]:
    st=0.41*mu/0.50
    r={}
    for cl in (0.95,0.68):
        u=mu_up(mu,st,syst,0.414,cl); r[cl]=(u,[seff(u,f) for f in (0.40,0.54,0.20)])
    out[lab]=r; print(lab, {cl:(round(v[0],2),[round(s,1) for s in v[1]]) for cl,v in r.items()})
