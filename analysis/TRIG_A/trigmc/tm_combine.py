#!/usr/bin/env python3
"""TRIG_A/trigmc — test zamknięcia wag w2(n_trk, lead_pt) wobec emulacji triggera L1_J12 na poziomie generatora (projekt A).

Pytanie: jeżeli dane to zdarzenia wybrane triggerem kalorymetrycznym (E_T okna 0,8×0,8 > próg), a MC minbias jest przeważane
tylko w (n_ch, lead_pt) do rozkładu tych zdarzeń (jak w2 w nominale), to o ile [pary/D*] przeważonego MC różni się od prawdziwego
[pary/D*] w zdarzeniach wyzwolonych?  kappa = R_true / R_w2  → μ pozorne = kappa przy danych ≡ Pythia.
Tylko poziom generatora (bez wydajności rekonstrukcji); warianty wag nóg pT^k (k = 0, 1, 2) pokazują czułość na kształt ε(pT).
Wejście: --gen KATALOG (gen_*.tsv.gz, .json, .h3em, .h3all z pythia_dd_jet). Wyjście: --out KATALOG → trigmc.json, TRIGMC.md.
"""
import argparse, glob, gzip, io, json, math, os, sys
import numpy as np
import pandas as pd

NN, NL, NE = 200, 60, 80
ECOL = ["t", "ev", "proc", "nch", "lead", "pthat", "etem", "etall", "weta", "wphi", "cjlead"]
DCOL = ["t", "ev", "pdg", "pt", "eta", "phi", "fromB", "mode", "cjpt"]
CH = ["D0D0", "DpD0", "DpDp"]
NJK = 20


def erf_eff(et, thr, sig):
    from scipy.special import erf
    return 0.5 * (1.0 + erf((np.asarray(et, float) - thr) / (math.sqrt(2.0) * sig)))


def read_h3(fn):
    h = np.zeros(NN * NL * NE)
    if os.path.exists(fn):
        a = np.loadtxt(fn, dtype=np.int64, ndmin=2)
        if a.size:
            np.add.at(h, a[:, 0], a[:, 1])
    return h.reshape(NN, NL, NE)


def read_gen(fn):
    txt = gzip.open(fn, "rt").read()
    E, D = [], []
    for l in txt.split("\n"):
        if not l or l[0] == "#":
            continue
        (E if l[0] == "E" else D).append(l)
    e = pd.read_csv(io.StringIO("\n".join(E)), sep="\t", header=None, names=ECOL)
    d = pd.read_csv(io.StringIO("\n".join(D)), sep="\t", header=None, names=DCOL)
    return e, d


def per_file(fn, fidx):
    """Zwraca (zdarzenia z D* lub parą: tabela zdarzeń, tabela D*, tabela par) — wszystko z kluczem pliku."""
    e, d = read_gen(fn)
    e = e[e.proc == 101].copy()
    d = d[(d.pt > 3) & (d.pt < 30) & (d.eta.abs() < 2.5)].merge(e, on="ev", suffixes=("", "_e"))
    d["jk"] = fidx % NJK
    ds = d[d.pdg.abs() == 413].copy()
    g = d[d.pdg.abs().isin([421, 411])].reset_index(drop=True)
    g["i"] = g.index
    keep = ["ev", "pdg", "pt", "eta", "phi", "fromB", "mode", "cjpt", "i"]
    m = g.merge(g[keep], on="ev", suffixes=("", "2"))
    m = m[m.i < m.i2].copy()
    a, b = m.pdg.abs(), m.pdg2.abs()
    m["ch"] = np.where((a == 421) & (b == 421), 0, np.where((a == 411) & (b == 411), 2, 1))
    m["ls"] = np.sign(m.pdg) == np.sign(m.pdg2)
    dphi = np.abs(np.angle(np.exp(1j * (m.phi - m.phi2))))
    m["dphi"] = dphi
    m["dr"] = np.hypot(dphi, m.eta - m.eta2)
    m["anyB"] = (m.fromB + m.fromB2) > 0
    m["mode11"] = (m["mode"] == 1) & (m.mode2 == 1)
    ecols = ["nch", "lead", "etem", "etall", "cjlead", "jk"]
    d0 = g[g.pdg.abs() == 421]
    return (ds[["pt", "mode", "cjpt"] + ecols], m[["ch", "ls", "dr", "dphi", "anyB", "mode11", "pt", "pt2"] + ecols],
            g[["pdg", "pt", "mode", "cjpt"] + ecols], len(e))


def cellw(h3, eff_e, nedges, ledges):
    """Wagi typu w2: w(komórka n_ch × lead) = N_trig/N_all; h3 (NN, NL, NE); eff_e (NE,) wydajność w binie E_T."""
    nb = np.clip(np.digitize(np.arange(NN) + 0.5, nedges) - 1, 0, len(nedges) - 2)
    lb = np.clip(np.digitize((np.arange(NL) + 0.5) * 0.5, ledges) - 1, 0, len(ledges) - 2)
    allc = np.zeros((len(nedges) - 1, len(ledges) - 1)); trg = np.zeros_like(allc)
    h2a = h3.sum(2); h2t = (h3 * eff_e[None, None, :]).sum(2)
    np.add.at(allc, (nb[:, None].repeat(NL, 1), lb[None, :].repeat(NN, 0)), h2a)
    np.add.at(trg, (nb[:, None].repeat(NL, 1), lb[None, :].repeat(NN, 0)), h2t)
    w = np.divide(trg, allc, out=np.zeros_like(trg), where=allc > 0)
    return w, float(h2t.sum() / max(h2a.sum(), 1.0)), allc


def look(w, nch, lead, nedges, ledges):
    i = np.clip(np.digitize(np.minimum(nch, NN - 1) + 0.5, nedges) - 1, 0, len(nedges) - 2)
    j = np.clip(np.digitize(np.minimum(lead, 29.99), ledges) - 1, 0, len(ledges) - 2)
    return w[i, j]


def jk_ratio(num_w, num_jk, den_w, den_jk):
    """R = Σnum/Σden z błędem jackknife po NJK grupach plików."""
    N = np.bincount(num_jk, weights=num_w, minlength=NJK); Dn = np.bincount(den_jk, weights=den_w, minlength=NJK)
    R = N.sum() / Dn.sum() if Dn.sum() > 0 else np.nan
    r = (N.sum() - N) / np.maximum(Dn.sum() - Dn, 1e-300)
    return R, float(np.sqrt((NJK - 1) / NJK * ((r - r.mean()) ** 2).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", required=True); ap.add_argument("--out", required=True); ap.add_argument("--max-files", type=int, default=0)
    ap.add_argument("--calib", help="j12eff.json: kalibracja emulacji na zmierzonej ε(J12) w zdarzeniach z kandydatem D* (427929) → siatka drobna + TRIGCAL.md")
    a = ap.parse_args()
    fs = sorted(glob.glob(os.path.join(a.gen, "gen_*.tsv.gz")))
    if a.max_files:
        fs = fs[:a.max_files]
    if not fs:
        sys.exit("STOP: brak gen_*.tsv.gz w " + a.gen)
    H = {"em": np.zeros((NN, NL, NE)), "all": np.zeros((NN, NL, NE))}
    DS, PR, SG, nev, nok = [], [], [], 0, 0
    for k, fn in enumerate(fs):
        try:
            ds, m, g, _ = per_file(fn, k)
            j = json.load(open(fn + ".json"))
        except Exception as ex:  # noqa: BLE001
            print("pomijam", fn, ex); continue
        H["em"] += read_h3(fn[:-3] + ".h3em") if os.path.exists(fn[:-3] + ".h3em") else read_h3(fn + ".h3em")
        H["all"] += read_h3(fn[:-3] + ".h3all") if os.path.exists(fn[:-3] + ".h3all") else read_h3(fn + ".h3all")
        DS.append(ds); PR.append(m); SG.append(g); nev += j["proc"].get("101", {}).get("n", 0); nok += 1
    DS, PR, SG = pd.concat(DS), pd.concat(PR), pd.concat(SG)
    ecen = (np.arange(NE) + 0.5) * 0.5
    BINN = {"nom": (np.array([0, 10, 20, 30, 40, 50, 60, 80, 100, 200.0]), np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 30.0])),
            "fine": (np.arange(0, 201, 4.0), np.arange(0, 30.01, 0.5))}
    EMU = [("em", t, s) for t in (6, 8, 10, 12, 15, 20) for s in (1.5, 3.0)] + [("all", t, s) for t in (15, 20, 25, 30) for s in (2.0, 4.0)]
    CAL = None
    if a.calib:
        EMU = [("em", t, sg) for t in (10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 35, 40) for sg in (2.0, 4.0, 6.0)] \
            + [("all", t, sg) for t in (20, 25, 30, 35, 40, 45, 50, 60) for sg in (3.0, 6.0, 9.0)]
        BINN = {"nom": BINN["nom"]}
        cj = json.load(open(a.calib))["data_427929"]["refs"]["mbts1|z D*"]["cj"]
        x, n, k = np.array(cj["x"]), np.array(cj["n"]), np.array(cj["k"])
        use = (x > 5) & (x < 20) & (n > 0)                           # zakres wiarygodny (powyżej ~20 GeV fałszywe tory w danych)
        edges = np.arange(0, 21.0, 1.0)
        CAL = {"x": x[use], "n": n[use], "k": k[use], "tot_n": float(n.sum()), "tot_k": float(k.sum()),
               "bin": np.clip(np.digitize(DS.cjlead.values, edges) - 1, 0, 19), "inr": DS.cjlead.values < 20, "ib": (x[use] - 0.5).astype(int)}
    KP = (0, 1, 2)
    res = {"n_files": nok, "n_events_nd": int(nev), "n_dstar": int(len(DS)), "n_pairs": int(len(PR)), "emulations": []}

    def block(wD, wP, wS, tag):
        """Stosunki dla zadanych wag zdarzeń (na wierszach D*, par, singli)."""
        o = {}
        for kp in KP:
            den = wD * DS.pt.values ** kp
            lw = wP * (PR.pt.values * PR.pt2.values) ** kp
            r = {}
            for c, nm in enumerate(CH + ["all"]):
                sel = (PR.ch.values == c) if c < 3 else np.ones(len(PR), bool)
                r[nm] = jk_ratio(lw[sel], PR.jk.values[sel], den, DS.jk.values)
            for nm, sel in (("LS", PR.ls.values), ("OS", ~PR.ls.values), ("OS_dr08", ~PR.ls.values & (PR.dr.values < 0.8)),
                            ("anyB", PR.anyB.values), ("mode11", PR.mode11.values)):
                r[nm] = jk_ratio(lw[sel], PR.jk.values[sel], den, DS.jk.values)
            r["f_LS"] = jk_ratio(lw[PR.ls.values], PR.jk.values[PR.ls.values], lw, PR.jk.values)
            for sp, pid in (("D0", 421), ("Dp", 411)):
                s = SG.pdg.abs().values == pid
                r[sp + "_per_Dstar"] = jk_ratio((wS * SG.pt.values ** kp)[s], SG.jk.values[s], den, DS.jk.values)
            o[f"k{kp}"] = r
        d0 = (SG.pdg.abs().values == 421)
        o["frac_D0_in_cj10"] = float((wS[d0] * (SG.cjpt.values[d0] > 10)).sum() / max(wS[d0].sum(), 1e-300))
        o["mean_lead_Dstar"] = float((wD * DS.lead.values).sum() / max(wD.sum(), 1e-300))
        o["mean_nch_Dstar"] = float((wD * DS.nch.values).sum() / max(wD.sum(), 1e-300))
        return o

    res["minbias"] = block(np.ones(len(DS)), np.ones(len(PR)), np.ones(len(SG)), "mb")
    for scale, thr, sig in EMU:
        col = "etem" if scale == "em" else "etall"
        eff_e = erf_eff(ecen, thr, sig)
        one = {"scale": scale, "thr": thr, "sigma": sig}
        tD, tP, tS = (erf_eff(X[col].values, thr, sig) for X in (DS, PR, SG))
        one["true"] = block(tD, tP, tS, "true")
        for bn, (ne, le) in BINN.items():
            w, frac, _ = cellw(H[scale], eff_e, ne, le)
            one["frac_events"] = frac
            one["w2_" + bn] = block(look(w, DS.nch.values, DS.lead.values, ne, le), look(w, PR.nch.values, PR.lead.values, ne, le),
                                    look(w, SG.nch.values, SG.lead.values, ne, le), "w2")
        if CAL is not None:                                          # przewidywana ε(J12) w zdarzeniach z D* vs dżet naładowany, χ² względem danych
            num = np.bincount(CAL["bin"][CAL["inr"]], weights=tD[CAL["inr"]], minlength=20); den = np.bincount(CAL["bin"][CAL["inr"]], minlength=20)
            pred = (num / np.maximum(den, 1))[CAL["ib"]]
            obs = CAL["k"] / CAL["n"]; err = np.sqrt(np.maximum(pred * (1 - pred), 1e-6) / CAL["n"])
            one["calib"] = {"chi2": float((((obs - pred) / err) ** 2).sum()), "ndf": int(len(obs)), "eps_Dstar_pred": float(tD.mean()),
                            "eps_Dstar_data": CAL["tot_k"] / CAL["tot_n"], "pred": pred.tolist(), "obs": obs.tolist(), "x": CAL["x"].tolist()}
        res["emulations"].append(one)
    os.makedirs(a.out, exist_ok=True)
    json.dump(res, open(os.path.join(a.out, "trigmc.json"), "w"), indent=1, default=float)
    if CAL is not None:
        val = lambda v: v[0] if isinstance(v, (list, tuple)) else v  # noqa: E731
        rows = []
        for o in res["emulations"]:
            kap = [val(o["true"][f"k{kp}"]["all"]) / val(o["w2_nom"][f"k{kp}"]["all"]) for kp in KP]
            kls = val(o["true"]["k0"]["LS"]) / val(o["w2_nom"]["k0"]["LS"]); kos = val(o["true"]["k0"]["OS"]) / val(o["w2_nom"]["k0"]["OS"])
            cw = val(o["w2_nom"]["k0"]["all"]) / val(res["minbias"]["k0"]["all"])
            rows.append((o["calib"]["chi2"], o, kap, kls, kos, cw))
        rows.sort(key=lambda r: r[0])
        T = ["# TRIG_A/trigcal — emulacja L1_J12 skalibrowana na danych (ε w zdarzeniach z kandydatem D*, run 427929, dżet torowy 5–20 GeV)", "",
             "χ² = porównanie przewidywanej ε(J12 | zdarzenie z D*, pT wiodącego dżetu naładowanego) z danymi w binach 1 GeV. Uwaga: w danych to kandydaci D*",
             "w oknie sygnału (z tłem kombinatorycznym) i dżety z torów zrekonstruowanych; w MC prawdziwe D* i dżety z cząstek naładowanych.",
             f"ε(J12 | D*) w danych = {rows[0][1]['calib']['eps_Dstar_data']:.4f}; ułamek zdarzeń MBTS z J12 w danych = 0,00256.", "",
             "| skala | próg | σ | χ²/ndf | ε(D*) MC | ułamek zd. MC | κ k=0 | κ k=1 | κ k=2 | κ LS | κ OS | c_w |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for chi2, o, kap, kls, kos, cw in rows:
            T.append(f"| {o['scale']} | {o['thr']} | {o['sigma']} | {chi2:.0f}/{o['calib']['ndf']} | {o['calib']['eps_Dstar_pred']:.4f} | {o['frac_events']:.2e} | "
                     + " | ".join(f"{v:.3f}" for v in kap) + f" | {kls:.3f} | {kos:.3f} | {cw:.3f} |")
        open(os.path.join(a.out, "TRIGCAL.md"), "w").write("\n".join(T) + "\n")
        print("OK TRIGCAL"); return

    # ---- raport
    L = ["# TRIG_A/trigmc — zamknięcie wag w2 wobec emulacji triggera dżetowego (generator, Pythia A3 ND)", "",
         f"Pliki: {nok}, zdarzenia ND: {nev:.3g}, D* (3–30 GeV, |η|<2,5): {len(DS)}, pary: {len(PR)}. Błędy: jackknife po {NJK} grupach plików.",
         "κ = [pary/D*] w zdarzeniach wyzwolonych / [pary/D*] w minbias przeważonym w (n_ch, lead_pt) do tych zdarzeń. κ ≠ 1 = obciążenie μ.",
         "Poziom generatora, bez wydajności rekonstrukcji; k = wykładnik wagi pT^k na nogę (proxy kształtu ε).", "",
         "## Minimum bias (odniesienie)", "",
         "| k | pary/D* | f_LS | udział par z b | D⁰ w dżecie nał. >10 GeV | ⟨lead⟩ D* |", "|---|---|---|---|---|---|"]
    mb = res["minbias"]
    for kp in KP:
        r = mb[f"k{kp}"]
        L.append(f"| {kp} | {r['all'][0]:.4f} | {r['f_LS'][0]:.3f} | {r['anyB'][0] / r['all'][0]:.3f} | {mb['frac_D0_in_cj10']:.3f} | {mb['mean_lead_Dstar']:.2f} |")
    for kp in KP:
        L += ["", f"## Emulacje, k = {kp}", "",
              "| skala | próg | σ | ułamek zd. | D⁰ w dżecie >10 | pary/D* prawda | κ (biny nom.) | κ (biny drobne) | κ D⁰D⁰ | κ D⁺D⁰ | κ D⁺D⁺ | κ LS | κ OS | f_LS prawda / w2 | c_w = w2/minbias |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for o in res["emulations"]:
            t, wn, wf = o["true"][f"k{kp}"], o["w2_nom"][f"k{kp}"], o["w2_fine"][f"k{kp}"]

            def kap(nm, w=wn):
                R = t[nm][0] / w[nm][0] if w[nm][0] else float("nan")
                e = R * math.hypot(t[nm][1] / max(t[nm][0], 1e-300), w[nm][1] / max(w[nm][0], 1e-300))
                return f"{R:.3f} ± {e:.3f}"
            L.append(f"| {o['scale']} | {o['thr']} | {o['sigma']} | {o['frac_events']:.2e} | {o['true']['frac_D0_in_cj10']:.3f} | {t['all'][0]:.4f} | "
                     f"{kap('all')} | {kap('all', wf)} | {kap('D0D0')} | {kap('DpD0')} | {kap('DpDp')} | {kap('LS')} | {kap('OS')} | "
                     f"{t['f_LS'][0]:.3f} / {wn['f_LS'][0]:.3f} | {wn['all'][0] / mb[f'k{kp}']['all'][0]:.3f} |")
    L += ["", "Punkty odniesienia z danych do wyboru emulacji: ułamek zdarzeń MBTS z L1_J12 w 427929 ≈ 0,0026; D⁰ (z D*) w dżecie torowym > 10 GeV:",
          "0,61 ± 0,02 (428770, J12-only), 0,16 ± 0,02 (427929); c_w nominału (pary/D*, ważone/nieważone): 1,37 / 1,20 / 0,98 (D⁰D⁰ / D⁺D⁰ / D⁺D⁺)."]
    open(os.path.join(a.out, "TRIGMC.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:12]))


if __name__ == "__main__":
    main()
