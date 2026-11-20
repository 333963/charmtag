import numpy as np
rng=np.random.default_rng(1)
# ALICE pp 5.02 TeV (arXiv 2105.06335 Tab.1), stat (+syst_up, -syst_dn) [%]
A={"f0":(39.1,1.7,2.5,3.7),"fp":(17.3,1.8,1.7,2.1),"fs":(15.5,1.2,4.1,1.9)}
# Pythia (mc23 mb_nd): ratios z prawdy (TRUTH_NPAIR: D0/D* 3.77, D+/D* 3.63 po BR) + bariony/Ds jak e+e- -> e+e- (Lisovyi 2016)
BR={"D0":0.03947,"Dp":0.0938,"Ds":0.677*0.03947}
r0=3.77*BR["Ds"]/BR["D0"]; rp=3.63*BR["Ds"]/BR["Dp"]
print(f"Pythia mb_nd truth: f(D0)/f(D*)={r0:.3f}, f(D+)/f(D*)={rp:.3f}; e+e-: {0.6086/0.2429:.3f}, {0.2404/0.2429:.3f}")
P={"f0":0.6086,"fp":0.2404,"fs":0.2429}
def k(a,p):
    return {"D0D0":(a["f0"]**2/a["fs"])/(p["f0"]**2/p["fs"]),"DpD0":(a["f0"]*a["fp"]/a["fs"])/(p["f0"]*p["fp"]/p["fs"]),"DpDp":(a["fp"]**2/a["fs"])/(p["fp"]**2/p["fs"])}
W={"D0D0":0.334,"DpD0":0.478,"DpDp":0.188}
def keff(kk): return 1/sum(W[c]/kk[c] for c in W)
c={x:A[x][0]/100 for x in A}; k0=k(c,P); print({x:round(v,3) for x,v in k0.items()},"k_eff",round(keff(k0),3))
def draw(m,s,u,d,n):
    z=rng.standard_normal(n); sy=np.where(z>0,u,d); return (m+s*rng.standard_normal(n)+sy*np.abs(z)*np.sign(z))/100
n=200000; T={x:draw(*A[x],n) for x in A}
kk=k(T,P); ke=1/sum(W[c]/kk[c] for c in W)
for name,v in list(kk.items())+[("eff",ke)]:
    q=np.percentile(v,[16,50,84]); print(name, f"{q[1]:.3f} [{q[0]:.3f}; {q[2]:.3f}]")
for mu in (0.50,0.53,0.72,1.27,1.50): print(mu, "->", round(mu/keff(k0),2), "(CRqcd -11%:", round(mu/0.89,2),")")
