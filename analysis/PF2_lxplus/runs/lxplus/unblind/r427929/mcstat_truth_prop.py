"""Propozycja 19.11: wiersz 'statystyka MC par' w budżecie 427929 z normalizacji z prawdy mb_nd (GEN_mcstat/out/truth_mbnd.json)
zamiast metody bezpośredniej (n_eff). Nie zmienia budzet_427929.json ani nominału."""
import json, math, os
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.abspath(os.path.join(H, "../../../../.."))
b = json.load(open(os.path.join(H, "budzet_427929.json")))
t = json.load(open(os.path.join(R, "GEN_mcstat/out/truth_mbnd.json")))
CH = ["D0D0", "DpD0", "DpDp"]; w = b["w"]
mu = 0.49894733917589273  # truthnorm nom427
rows = b["rows"]; mc = [r for r in rows if r["source"].startswith("statystyka MC par")][0]
var = {
 "A_w2_jk+boot (zachowawczo, baza 428770 dziedziczona)": [t["stat_rel"]["total"][c] for c in CH],
 "B_unw_poisson+boot (wagi 427929 ≈ 1, dolne oszacowanie)": [math.sqrt(1 / t["n_true_pairs"][c] + t["stat_rel"]["bootstrap_maps"][c] ** 2) for c in CH],
}
def totals(new_ch):
    new_c = sum(wi * x for wi, x in zip(w, new_ch))
    g = {}
    for r in rows:
        c = new_c if r is mc else r["common"]
        g.setdefault(r["group"], []).append(c)
    gs = {k: math.sqrt(sum(x * x for x in v)) for k, v in g.items()}
    tot_dn = math.sqrt(sum(x * x for x in gs.values()))
    tr = [r for r in rows if r["source"].startswith("tune Pythii")][0]
    up_c = sum(wi * x for wi, x in zip(w, b["tune_up"]))
    nup = math.sqrt(gs["Normalizacja"] ** 2 - tr["common"] ** 2 + up_c ** 2)
    tot_up = math.sqrt(tot_dn ** 2 - gs["Normalizacja"] ** 2 + nup ** 2)
    return new_c, tot_up, tot_dn
out = {"mu": mu, "old": {"mc_common": mc["common"], "rel_up": b["syst_rel_up"], "rel_dn": b["syst_rel_dn"]}}
print(f"stary: mc_stat {mc['common']:.3f}; syst +{b['syst_rel_up']:.3f}/-{b['syst_rel_dn']:.3f} -> +{mu*b['syst_rel_up']:.3f}/-{mu*b['syst_rel_dn']:.3f}")
for k, v in var.items():
    c, up, dn = totals(v)
    out[k] = {"ch": v, "common": c, "rel_up": up, "rel_dn": dn, "abs_up": mu * up, "abs_dn": mu * dn}
    print(f"{k}: ch {[round(x,3) for x in v]} wspólne {c:.3f}; syst +{up:.3f}/-{dn:.3f} -> +{mu*up:.3f}/-{mu*dn:.3f}")
json.dump(out, open(os.path.join(H, "mcstat_truth_prop.json"), "w"), indent=1)
