"""N2 v4 — SF D0 (species 10) zintegrowane po pT z binów B2b (EV pełna lumi), D⁺ (4) bez zmian (v3).
Wagi binów = N_D* per bin z B2b (pf95, ±55): 3–5 735, 5–8 1375, 8–30 2110 (proxy widma pT D0 w parach; nadpisywalne --w).
Bin 8–30: SF abs B2b; 5–8: SF95 = 0,91±0,10 (syst pokrywający 0,82–1,00, natura składnika nierozstrzygnięta) × SF_rel sW;
3–5: niemierzalne → 0,91±0,15 (jak 5–8). WP<0,5: skalowanie v3 czynnikiem v4/v3 @0,5. Wynik: sf_for_pairs_v4_proposed.json."""
import json, sys, math
w = dict(zip(["3-5","5-8","8-30"], [735.,1375.,2110.]))
if "--w" in sys.argv:
    w = dict(zip(["3-5","5-8","8-30"], map(float, sys.argv[sys.argv.index("--w")+1].split(","))))
tot = sum(w.values()); w = {k: v/tot for k, v in w.items()}
v3 = json.load(open("sf_for_pairs_v3.json"))
b2b = json.load(open("sf_for_pairs_v3_proposed.json"))["B2b_EV"]
sf830, u830 = b2b["SF"]["8-30"], b2b["SF_unc"]["8-30"]
rel58 = {"0.95":1.0, "0.90":1.01, "0.80":1.00, "0.70":1.01, "0.50":1.06}   # sW, B2b §0
sf58 = {k: 0.91*rel58[k] for k in rel58}; u58 = {k: math.hypot(0.10, 0.91*0.05) for k in rel58}
sf35 = {k: 0.91 for k in rel58}; u35 = {k: 0.15 for k in rel58}
out, unc = {}, {}
for wp in ["0.95","0.90","0.80","0.70","0.50"]:
    out[wp] = w["3-5"]*sf35[wp] + w["5-8"]*sf58[wp] + w["8-30"]*sf830[wp]
    unc[wp] = math.sqrt((w["3-5"]*u35[wp])**2 + (w["5-8"]*u58[wp])**2 + (w["8-30"]*u830[wp])**2)
def interp(x, a, b, fa, fb): return fa + (fb-fa)*(x-a)/(b-a)
d0 = {"0.9": out["0.90"], "0.8": out["0.80"], "0.7": out["0.70"], "0.65": interp(.65,.5,.7,out["0.50"],out["0.70"]),
      "0.6": interp(.6,.5,.7,out["0.50"],out["0.70"]), "0.5": out["0.50"]}
r = out["0.50"]/v3["10"]["0.5"]
d0["0.4"] = v3["10"]["0.4"]*r; d0["0.3"] = v3["10"]["0.3"]*r
d0u = {"0.9": unc["0.90"], "0.8": unc["0.80"], "0.7": unc["0.70"], "0.65": max(unc["0.50"],unc["0.70"]), "0.5": unc["0.50"]}
v4 = {"_note": "N2 v4 (2026-11-05, B2b EV pelna lumi): D0 = srednia wazona N_D* binow 3-5/5-8/8-30 (%s); 8-30 SF abs B2b ±0,03-0,04; "
      "5-8 SF95 0,91±0,10 (syst 0,82-1,00) x SF_rel sW; 3-5 0,91±0,15. D+ = v3 (brak nowej probki). SF95(D0) = %.3f (v3 zakladalo 1). "
      "Nominal unblind nadal v3 do decyzji uzytkownika; v4 = kandydat (zgodny z v3 w 1 sigma, blad 2x mniejszy)." % (", ".join("%s %.2f"%kv for kv in w.items()), out["0.95"]),
      "10": {k: round(v,4) for k,v in d0.items()}, "4": v3["4"],
      "unc": {"10": {k: round(v,3) for k,v in d0u.items()}, "4": v3["unc"]["4"]},
      "SF95_D0": {"val": round(out["0.95"],4), "unc": round(unc["0.95"],3)},
      "per_bin_D0": {"8-30": {"SF": sf830, "unc": u830}, "5-8": {"SF": {k: round(v,3) for k,v in sf58.items()}, "unc": round(u58["0.95"],3)},
                     "3-5": {"SF": 0.91, "unc": 0.15}, "weights": w},
      "source": {"v3": "sf_for_pairs_v3.json", "B2b": "charm_work/B2b_dstar/ (charm-B2b-dstar-anomalia.md §7)"}}
json.dump(v4, open("sf_for_pairs_v4_proposed.json","w"), indent=1, ensure_ascii=False)
print("D0 v4:", {k: round(v,3) for k,v in d0.items()}); print("unc:", {k: round(v,3) for k,v in d0u.items()})
print("v3 :", v3["10"]); print("SF95:", round(out["0.95"],3), "±", round(unc["0.95"],3))
