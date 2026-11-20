"""TUNE — systematyka normalizacji μ od tune Pythii: [pary/D*]_reco wariant vs A3, per kanał i grupa (MPI, CR, Monash).
Wejście: ../runs/lxplus/truth/combined_A3 + combined_syst_<wariant>/truth_norm.json (etap ORCH `tune`). Logika = charmtag.syst.norm.collect.
Uruchom z charm_work: ../.venv_h5/bin/python TUNE_pythia/tune_syst.py
"""
import json, math, os, sys

sys.path.insert(0, "Q_lxplus/syst")
from charmtag_syst import norm as N  # noqa: E402

TD = "../runs/lxplus/truth"
CH = ("D0D0", "DpD0", "DpDp")
nom = json.load(open(f"{TD}/combined_A3/truth_norm.json"))
rows, miss = {}, []
for grp, names in N.TUNE_GROUPS.items():
    for nm in names:
        fn = f"{TD}/combined_syst_{nm}/truth_norm.json"
        if not os.path.exists(fn):
            miss.append(nm); continue
        rv = json.load(open(fn))
        rows[nm] = {"group": grp, "n_events": (rv.get("generator") or {}).get("n_events"),
                    "delta": {ch: (N._ppd(rv, ch) / N._ppd(nom, ch) - 1.0) if (N._ppd(rv, ch) and N._ppd(nom, ch)) else None for ch in CH}}
out = {"variants": rows, "missing": miss, "groups": {}}
for grp, names in N.TUNE_GROUPS.items():
    have = [n for n in names if n in rows]
    out["groups"][grp] = {ch: max([abs(rows[n]["delta"][ch] or 0) for n in have] or [0.0]) for ch in CH} if have else None
tot = {ch: math.sqrt(sum((g[ch] ** 2) for g in out["groups"].values() if g)) for ch in CH}
out["total_rel"] = tot
os.makedirs("TUNE_pythia/out", exist_ok=True)
json.dump(out, open("TUNE_pythia/out/tune_syst.json", "w"), indent=1)
md = ["# Tune Pythii — δ([pary/D*]_reco) względem A3", "", "| wariant | grupa | N zd. | D⁰D⁰ | D⁺D⁰ | D⁺D⁺ |", "|---|---|---|---|---|---|"]
f = lambda v: "—" if v is None else f"{100*v:+.1f} %"  # noqa: E731
for nm, r in rows.items():
    md.append(f"| {nm} | {r['group']} | {r['n_events']} | " + " | ".join(f(r["delta"][ch]) for ch in CH) + " |")
md += ["", "| grupa (max \\|δ\\|) | D⁰D⁰ | D⁺D⁰ | D⁺D⁺ |", "|---|---|---|---|"]
for g, v in out["groups"].items():
    md.append(f"| {g} | " + (" | ".join(f"{100*v[ch]:.1f} %" for ch in CH) if v else "CZEKA | CZEKA | CZEKA") + " |")
md.append("| **suma kwadratowa** | " + " | ".join(f"**{100*tot[ch]:.1f} %**" for ch in CH) + " |")
md += ["", "δμ/μ = −δ/(1+δ) (μ ∝ 1/oczekiwanie). Błąd statystyczny wariantu (10⁸ zd. vs 5·10⁸ A3) nie jest tu odjęty — przy |δ| rzędu kilku % sprawdzić rozrzut między seedami.",
       f"Brakujące warianty: {', '.join(miss) if miss else 'brak'}."]
open("TUNE_pythia/tune_syst.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
