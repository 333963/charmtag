#!/usr/bin/env python3
"""T-PREP: wypełnia tfin_template.md liczbami z biegu pf (compare_pf.json, toys, syst, mix) i — gdy jest — z odsłonięcia (fit2d_unblind.json).
Użycie: tfin_fill.py --pf <PF_OUT> [--unblind <UB_OUT>] --template tfin_template.md --out tfin  → tfin.md + tfin.json
Brak pliku/klucza → "—" (nigdy nie zgaduje). Progi bramek G1–G4, G6, G11 jak w unblinding_checklist.md."""
import argparse, glob, json, os, time
CH = ("D0D0", "DpD0", "DpDp")
NAME = {"D0D0": "D0D0", "DpD0": "D⁺D0", "DpDp": "D⁺D⁺"}
SHARED_BAND = {"D0D0": (0.30, 0.65), "DpD0": (0.21, 0.50), "DpDp": (0.72, 0.93)}   # PF3: [dane½ − 10 pp; MC + 2 pp], SB×SB po progach par
def jl(p):
    try: return json.load(open(p))
    except Exception: return None
def g(d, *ks):
    for k in ks:
        if not isinstance(d, dict) or k not in d: return None
        d = d[k]
    return d
def first(*vals):  # pierwsza wartość nie-None (0.0 jest wartością!)
    for v in vals:
        if v is not None: return v
    return None
def f(x, nd=2):
    try: return "—" if x is None else f"{float(x):.{nd}f}"
    except Exception: return str(x)
def pct(a, b):
    try: return f"{100*(float(a)/float(b)-1):+.0f} %"
    except Exception: return "—"
def ok(cond):
    return "—" if cond is None else ("✔" if cond else "✖")
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pf", required=True); ap.add_argument("--unblind"); ap.add_argument("--template", required=True); ap.add_argument("--out", required=True); ap.add_argument("--ctrl", help="compare_ctrl.json (PF3) — izolowany udział usuniętych par per kanał (G6)")
    a = ap.parse_args(); P = os.path.join(a.pf, "pairs")
    cmp_ = jl(os.path.join(a.pf, "compare_pf.json")) or {}
    counts = cmp_.get("counts_pf") or jl(os.path.join(P, "pairs_counts.json")) or {}
    toys = cmp_.get("toys_pf") or jl(os.path.join(P, "toys", "fit2d_toys.json")) or jl(os.path.join(P, "fit2d_toys.json")) or {}
    fpf, fno = cmp_.get("fit2d_pf") or {}, cmp_.get("fit2d_nominal") or {}
    sb_pf, sb_no = cmp_.get("sbsb_pf") or {}, cmp_.get("sbsb_nominal") or {}
    syst_files = cmp_.get("syst_pf") or {os.path.basename(x): jl(x) for x in glob.glob(os.path.join(P, "syst", "*.json"))}
    sblind = jl(os.path.join(P, "syst", "syst_blind.json")) or (syst_files.get("syst_blind.json") if isinstance(syst_files, dict) else None) or {}
    ysyst = jl(os.path.join(P, "syst", "syst_bkg_model_Dplus.json")) or {}
    salad_files = {os.path.basename(x): jl(x) for x in glob.glob(os.path.join(P, "syst", "bkg_shape_salad_*.json"))}
    mix = cmp_.get("mix_pf") or jl(os.path.join(a.pf, "mix", "pairs", "mix_control.json")) or {}
    ub = jl(os.path.join(a.unblind, "pairs", "fit2d_unblind.json")) if a.unblind else None
    out = {"generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()), "pf": a.pf, "mode": "ODSŁONIĘTE" if ub else "ŚLEPO",
           "unblind_src": f" i `{a.unblind}`" if ub else "", "norm_mode": cmp_.get("norm_mode_pf") or "—"}
    b, s = counts.get("pairs_before_veto"), counts.get("pairs_shared_track")
    out.update(track_id=counts.get("track_id", "—"), pairs_before_veto=b if b is not None else "—", pairs_shared_track=s if s is not None else "—",
               shared_pct=f(100*s/b, 1) if (b and s is not None) else "—", pairs_kept=counts.get("pairs_kept", "—"), region=counts.get("region", "—"))
    rows, gates = [], {}
    def sbv(d, ch):
        if isinstance(d, dict) and ch not in d and any(k.startswith(ch + "_") for k in d):
            return sum(v for k, v in d.items() if k.startswith(ch + "_") and isinstance(v, (int, float)))
        v = d.get(ch) if isinstance(d, dict) else None
        if isinstance(v, dict): v = v.get("n") or v.get("total") or sum(x for x in v.values() if isinstance(x, (int, float)))
        return v
    for ch in CH:
        cp, cn = dict(fpf.get(ch) or {}), fno.get(ch) or {}
        if cp.get("rho") is None: cp["rho"] = g(mix, "channels", ch, "corr_m1m2_same_sbsb")
        rows.append(f"| {NAME[ch]} | {f(sbv(sb_no, ch),0)} | {f(sbv(sb_pf, ch),0)} | {pct(sbv(sb_pf, ch), sbv(sb_no, ch))} | {f(cn.get('Z_asimov_arms'))} | {f(cp.get('Z_asimov_arms'))} | "
                    f"{f(cn.get('sigma_mu'))} | {f(cp.get('sigma_mu'))} | {f(cn.get('dev_ndf_blind'))} | {f(cp.get('dev_ndf_blind'))} | {f(cp.get('rho'),3)} | {f(cp.get('tilt'),3)} |")
        gch = {}
        dn = cp.get("dev_ndf_blind"); gch["G4"] = None if dn is None else float(dn) < 1.5
        za, zn = cp.get("Z_asimov_arms"), cn.get("Z_asimov_arms")
        gch["G11"] = None if (za is None or zn is None) else float(za) >= float(zn) - 0.1
        gates[ch] = gch
    out["blind_rows"] = "\n".join(rows)
    out["Zq_nominal"], out["Zq_pf"] = f(fno.get("Z_quadrature")), f(fpf.get("Z_quadrature"))
    trows = []
    for ch in CH:
        t = (toys.get("channels") or toys).get(ch) if isinstance(toys, dict) else None
        t = t or {}
        sp, bias, cov, n = first(t.get("pull_sigma"), t.get("sigma_pull"), t.get("pull_std")), t.get("bias"), first(t.get("coverage68"), t.get("coverage")), first(t.get("n"), t.get("ntoy"))
        smu = (fpf.get(ch) or {}).get("sigma_mu")
        bs = (float(bias)/float(smu)) if (bias is not None and smu) else None
        g1 = None if sp is None else 0.9 <= float(sp) <= 1.1
        g2 = None if bs is None else abs(bs) < 0.10
        g3 = None if cov is None else 0.62 <= float(cov) <= 0.74
        gates[ch].update(G1=g1, G2=g2, G3=g3)
        trows.append(f"| {NAME[ch]} | {f(sp)} | {f(bs)} | {f(cov)} | {n or '—'} | {ok(g1)} | {ok(g2)} | {ok(g3)} |")
    out["toys_rows"] = "\n".join(trows)
    srows = []
    rep = None
    for k, v in syst_files.items():
        if isinstance(v, dict) and ("per_channel" in v or "total" in v): rep = v
    for ch in CH:
        stat = first(g(sblind, "sigma_stat", ch), (fpf.get(ch) or {}).get("sigma_mu"))
        pc = (g(rep, "per_channel", ch) or {}) if rep else {}
        T = sblind.get("totals") or {}
        def tq(*names):
            v = [float(T[n][ch]) for n in names if isinstance(T.get(n), dict) and T[n].get(ch) is not None]
            return (sum(x * x for x in v) ** 0.5) if v else None
        shape = first(pc.get("shape_bkg_kappa"), tq("Kształty", "κ, kategorie", "Tło"))
        sf, norm = first(pc.get("sf"), tq("SF")), first(pc.get("norm"), tq("Normalizacja"))
        tot = first(pc.get("total"), g(sblind, "total_all", ch))
        y = g(ysyst, "data_rel", "wp70", ch) or g(ysyst, ch)
        sal = None
        for fn, v in salad_files.items():
            if ch.startswith("Dp") and "Dplus" in fn: sal = g(v, "rel", "wp70") or g(v, "syst_rel")
            if ch == "D0D0" and "D0" in fn: sal = g(v, "rel", "wp70") or g(v, "syst_rel")
        ratio = (float(tot)/float(stat)) if (tot is not None and stat) else None
        g7 = None if ratio is None else ratio < 1.0
        gates[ch]["G7"] = g7
        srows.append(f"| {NAME[ch]} | {f(stat)} | {f(shape)} | {f(sf)} | {f(norm)} | {f(y,3)} | {f(sal,3)} | {f(tot)} | {f(ratio)} | {ok(g7)} |")
    out["syst_rows"] = "\n".join(srows); out["n_syst_json"] = len(syst_files)
    if mix:
        out["mix_line"] = " · ".join(f"{NAME[ch]}: χ²/ndf SB×SB same vs mix {f(first(g(mix, 'channels', ch, 'sbsb_chi2_ndf'), g(mix, ch, 'chi2ndf')))}, corr same {f(first(g(mix, 'channels', ch, 'corr_m1m2_same_sbsb'), g(mix, ch, 'corr_same')),3)} / mix {f(first(g(mix, 'channels', ch, 'corr_m1m2_mix_sbsb'), g(mix, ch, 'corr_mix')),3)}, ramiona same/mix {f(g(mix, 'channels', ch, 'arms_ratio_same_over_mix'),3)}" for ch in CH)
    else: out["mix_line"] = "— (brak mix_control.json)"
    grows = []
    sh_ok = {}
    ctrl = jl(a.ctrl) if a.ctrl else None
    g6src, g6fr = "brak", {}
    for ch in CH:
        lo, hi = SHARED_BAND[ch]
        fr = first(g(ctrl, "channels", ch, "SBxSB", "frac_removed"), g(ctrl, ch, "sbsb", "frac_removed"), g(ctrl, ch, "frac_removed_sbsb")) if ctrl else None
        if fr is not None: g6src = "compare_ctrl (izolowane)"
        else:
            no, pf_ = sbv(sb_no, ch), sbv(sb_pf, ch)
            fr = (1 - pf_ / no) if (no and pf_ is not None) else None
            if fr is not None and g6src == "brak": g6src = "SB×SB nominal→PF (NIEizolowane: miesza weto ze zmianą kodu)"
        g6fr[ch] = fr
        sh_ok[ch] = None if fr is None else lo <= fr <= hi
        gates[ch]["G6"] = sh_ok[ch]
    out["g6_line"] = f"G6 źródło: {g6src}; udział usuniętych SB×SB: " + ", ".join(f"{NAME[ch]} {f(100*g6fr[ch],0) if g6fr[ch] is not None else '—'} % (pas {int(100*SHARED_BAND[ch][0])}–{int(100*SHARED_BAND[ch][1])} %)" for ch in CH)
    for gid, lab in (("G1", "toys σ(pull) 0,9–1,1"), ("G2", "toys |bias| < 0,1σ"), ("G3", "pokrycie 68 % 0,62–0,74"), ("G4", "dev/ndf ślepy < 1,5"),
                     ("G6", "udział usuniętych par SB×SB w pasie (per kanał)"), ("G7", "syst/stat < 1"), ("G11", "Z_A PF ≥ Z_A nom − 0,1")):
        grows.append(f"| {gid} | {lab} | " + " | ".join(ok(gates[ch].get(gid)) for ch in CH) + " |")
    out["gate_rows"] = "\n".join(grows)
    core = [gates[ch].get(k) for ch in CH for k in ("G1", "G2", "G3", "G4", "G6", "G7", "G11")]
    out["verdict"] = "NO-GO" if any(v is False for v in core) else ("GO (automatyczne)" if all(v is True for v in core) else "NIEKOMPLETNE")
    if ub:
        L = ["| kanał | μ̂ | 68 % | Z lok. | dev/ndf | μ bez tagu | κ=channel μ (Z) |", "|---|---|---|---|---|---|---|"]
        for ch in CH:
            r = g(ub, "channels", ch) or {}
            p68 = r.get("mu_profile_68"); p68s = f"[{f(p68[0])}, {f(p68[1])}]" if isinstance(p68, (list, tuple)) and len(p68) == 2 else "—"
            nt = r.get("mu_no_tag"); nts = f"{f(nt[0])} ± {f(nt[1])}" if isinstance(nt, (list, tuple)) else "—"
            kc = r.get("crosscheck_kchannel") or {}
            L.append(f"| {NAME[ch]} | {f(r.get('mu_hat'))} ± {f(r.get('mu_err_hesse'))} | {p68s} | {f(r.get('Z'))} (exp {f((fpf.get(ch) or {}).get('Z_asimov_arms'))}) | {f(r.get('dev_ndf'))} | {nts} | {f(kc.get('mu_hat'))} ({f(kc.get('Z'))}) |")
        c = ub.get("combined") or {}; c68 = c.get("mu_common_68"); c68s = f"[{f(c68[0])}, {f(c68[1])}]" if isinstance(c68, (list, tuple)) and len(c68) == 2 else "—"
        sp = g(ub, "syst_profiled", "combined") or {}
        L += ["", f"**Wspólne μ = {f(c.get('mu_common'))} {c68s}, Z = {f(c.get('Z'))}** (stat); z systematykami profilowanymi: μ = {f(sp.get('mu_common'))}, Z = {f(sp.get('Z'))}.",
              "Wynik vs ślepy Asimov: Z_A PF = " + out["Zq_pf"] + "; dane/MC μ i przedział jak wyżej (normalizacja: " + str(out.get("norm_mode", "—")) + ")."]
        out["unblind_section"] = "\n".join(L); out["Z_obs"] = f(c.get("Z")); out["Z_obs_note"] = "stat; z syst: " + f(sp.get("Z"))
    else:
        out["unblind_section"] = "— (nie odsłonięte; `submit_unblind.sh submit` po decyzji)"; out["Z_obs"] = "—"; out["Z_obs_note"] = "—"
    tpl = open(a.template, encoding="utf-8").read()
    for k, v in out.items(): tpl = tpl.replace("{{" + k + "}}", str(v))
    open(a.out + ".md", "w", encoding="utf-8").write(tpl)
    json.dump({"fields": out, "gates": gates}, open(a.out + ".json", "w"), indent=1, default=str)
    print(tpl)
if __name__ == "__main__":
    main()
