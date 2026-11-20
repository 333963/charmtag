#!/usr/bin/env python3
"""Prezentacja projektu A (19.11.2026) wg kartki → deck.pptx + deck.pdf (ta sama specyfikacja slajdów).

Uruchomienie (z dowolnego miejsca):
  PYTHONPATH=<katalog z python-pptx> ~/Documents/charm/.venv_h5/bin/python build_deck.py
Bez python-pptx powstaje tylko PDF (matplotlib).

Podmiana treści:
  - liczby w toku: deck/numbers.json (klucze {{...}}; brak klucza → na slajdzie „{{KLUCZ}} (in progress)”),
  - wykresy punktu 3 (BDT vs NN): gdy w ../p3/runs/ są roc_*.png / scores_*.png / importance_*.png / sculpt_*.png,
    slajdy z nimi tworzą się automatycznie (inaczej slajd-placeholder z opisem modeli).
"""
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PREZ = os.path.dirname(HERE)
P1, P3, P6 = (os.path.join(PREZ, x) for x in ("p1", "p3", "p6"))
FX = os.path.join(HERE, "figs_extra")                       # wykresy dodatkowe (make_extra.py) + p1sim (make_p1_sim.py: MC → "")
P1S = os.path.join(FX, "p1sim")
KRT = os.path.join(FX, "kartka")                            # uzupełnienia wg kartki (make_kartka.py, 20.11)
FD = os.path.join(HERE, "figs_data")                        # fity na danych (fork b; slajdy tworzą się tylko gdy pliki istnieją)
NUM = json.load(open(os.path.join(HERE, "numbers.json"))) if os.path.exists(os.path.join(HERE, "numbers.json")) else {}
TAG = ""
DATE = "19 Nov 2026"


def fill(s):
    """{{KLUCZ}} → wartość z numbers.json albo widoczny placeholder."""
    return re.sub(r"\{\{(\w+)\}\}", lambda m: str(NUM[m.group(1)]) if m.group(1) in NUM else f"[{m.group(1)}: in progress]", s)


def img(*names, d=None):
    out = []
    for n in names:
        nn = n if n.endswith(".png") else n + ".png"
        cands = [os.path.join(d, nn)] if d else [os.path.join(P1S, nn), os.path.join(P1, nn)]   # p1sim ma pierwszeństwo
        p = next((c for c in cands if os.path.exists(c)), cands[0])
        if os.path.exists(p):
            out.append(p)
        else:
            print("BRAK wykresu:", p, file=sys.stderr)
    return out


# ------------------------------------------------------------------------------------------------ specyfikacja slajdów
def data_figs():
    """figs_data/manifest.json (fork b): lista lub słownik wpisów {file|png|path, title, caption|description, slot?}."""
    out = []
    if not os.path.isdir(FD):
        return out
    man = os.path.join(FD, "manifest.json")
    ent = []
    if os.path.exists(man):
        m = json.load(open(man))
        if isinstance(m, dict):
            m = m.get("figures", m.get("plots", [dict(v, file=v.get("file", k)) if isinstance(v, dict) else {"file": k, "title": str(v)} for k, v in m.items()]))
        for e in m:
            f = e.get("file") or e.get("png") or e.get("path")
            if not f:
                continue
            f = f if os.path.isabs(f) else os.path.join(FD, f)
            if not f.endswith(".png") and os.path.exists(os.path.splitext(f)[0] + ".png"):
                f = os.path.splitext(f)[0] + ".png"
            if os.path.exists(f) and f.endswith(".png"):
                ent.append({"file": f, "title": e.get("title", os.path.basename(f)), "caption": e.get("caption") or e.get("description") or "",
                            "slot": e.get("slot", "")})
    else:
        ent = [{"file": f, "title": os.path.splitext(os.path.basename(f))[0], "caption": "", "slot": ""} for f in sorted(glob.glob(os.path.join(FD, "*.png")))]
    return ent


def slot_of(e):
    k = (e["slot"] + " " + os.path.basename(e["file"]) + " " + e["title"]).lower()
    if "datamc" in k or "data_mc" in k or "data/mc" in k or "goodbad" in k:
        return "datamc"
    if "excess" in k or "margin" in k:
        return "excess"
    if "dstar" in k or "d*" in k or "dm" in k.split("_") or "delta" in k:
        return "dstar"
    return "pairs"


def unc_table():
    """Tabela rozbicia niepewności 427929 (figs_extra/unc_breakdown.json z make_extra.py)."""
    f = os.path.join(FX, "unc_breakdown.json")
    if not os.path.exists(f):
        return None
    u = json.load(open(f))
    pc = lambda a, b: f"+{100 * a:.0f} / −{100 * b:.0f}" if abs(a - b) > 0.005 else f"±{100 * a:.0f}"  # noqa: E731
    rows = [["source (run 427929, μ = 0.50)", "Δμ/μ [%]", "", "source", "Δμ/μ [%]"]]
    R = [(a, pc(x, y)) for a, x, y in u["unc_rows"]]
    half = (len(R) + 1) // 2
    for i in range(half):
        l = R[i]; r = R[half + i] if half + i < len(R) else ("", "")
        rows.append([l[0], l[1], "", r[0], r[1]])
    rows.append(["Total systematic (w/o N(D*))", pc(u["syst_up"], u["syst_dn"]), "", "Total N(D*)", pc(u["nds_total"], u["nds_total"])])
    rows.append(["", "", "", "TOTAL (stat ⊕ syst ⊕ N(D*))", pc(*u["total"])])
    return rows


def slides():
    S = []
    DF = data_figs()
    add = lambda title, bullets=(), images=(), table=None, section="", note="": S.append(  # noqa: E731
        {"title": title, "bullets": list(bullets), "images": list(images), "table": table, "section": section, "note": note})

    add("Double-charm production (D–D pairs) in 2022 minimum-bias data",
        ["Result: no D–D pair excess over Pythia in clean minimum bias; the excess of run 428770 is a jet-trigger (J12) bias",
         "[author, affiliation — to be filled]",
         "pp at √s = 13.6 TeV, 2022: run 427929 (clean MBTS minimum bias, ≈ 29.0 M events in v2.0 ntuples) and run 428770 (J12-dominated, ≈ 8.0 M events)",
         f"{DATE}", "Outline: (1) signal & background definitions  (2) invariant-mass distributions  (3) classifier models: BDT vs NN",
         "(4) MC results: score, GOOD/BAD  (5) data GOOD/BAD  (6) pair network & pair topologies  (7) results: μ, trigger, clean run"],
        section="Title")

    # ---------------- (1)
    add("(1) Signal and background definitions — truth matching",
        ["Labels from Dmeson_truth_category (CharmAnalysis ntuples), grouped as in charmpy categorize_dmeson",
         "Background in training and score plots = everything except signal and K↔π reflections",
         "D*: signal = true D* with a true D⁰ daughter"],
        table=[["group", "codes", "meaning"],
               ["Signal", "0, 1", "all daughters matched to a true D of this species/mode (0 prompt c, 1 from b)"],
               ["Reflection K↔π (D⁰)", "100, 101", "same tracks as signal, K and π hypotheses swapped"],
               ["Other charm / other mode", "2, 3, 6, 7", "tracks from a charm hadron of another species or decay mode"],
               ["Partially reconstructed", "8–11", "some daughters from a c/b hadron, others random"],
               ["Combinatorial", "4, 5, 12–15", "no relation to a charm hadron (incl. b and light)"]],
        section="1")
    add("(1) Composition of the signal window: before → after the BDT cut (WP70), MC minimum bias",
        ["Signal window |m − m₀| < 55 MeV (D⁰), 45 MeV (D⁺); MC mb_nd (Pythia 8 A3, non-diffractive)",
         "After WP70 the window is still background-dominated in minimum bias → background modelling from data sidebands",
         "Forced-decay MC after WP70: D⁰ 64 % signal (20 % reflections), D⁺ 76 % signal (21 % partially reco)"],
        table=[["", "signal", "partially reco", "combinatorial", "reflections"],
               ["D⁰", "0.6 % → 9.7 %", "5.3 % → 10.3 %", "93.6 % → 74.8 %", "0.15 % → 2.9 %"],
               ["D⁺", "0.4 % → 20.9 %", "8.7 % → 34.8 %", "90.7 % → 42.0 %", "–"]],
        section="1")

    # ---------------- (2)
    add("(2) D⁰ → K⁻π⁺ invariant mass: signal and background components",
        ["Left: minimum-bias MC — signal invisible in the stack (S/B = 0.6 % in window), overlaid ×10",
         "Right: forced-decay MC — 24 % signal in window, large partially-reconstructed component (22 %)"],
        images=img("p2_mass_D0_incl_mc", "p2_mass_D0_signal_mc"), section="2")
    add("(2) D⁺ → K⁻π⁺π⁺ invariant mass: signal and background components",
        ["Left: minimum-bias MC — S/B = 0.4 % in window", "Right: forced-decay MC — 14 % signal, 32 % partially reconstructed"],
        images=img("p2_mass_Dplus_incl_mc", "p2_mass_Dplus_signal_mc"), section="2")
    add("(2) D*⁺ → D⁰π⁺: Δm = m(Kππ) − m(Kπ)",
        ["Minimum-bias MC: 510 true D* among 325 k candidates (24 files) — small peak at 145.4 MeV",
         "Forced-decay MC: 55.7 k signal among 148 k candidates", "Data 428770 (next slides): peak visible only for GOOD D⁰ legs"],
        images=img("p2_dstar_dm_mb_nd", "p2_dstar_dm_fd_nd"), section="2")
    add("(2) Cutflow — minimum-bias MC (full statistics, nominal BDT \"H\" scores)",
        ["Finder preselection (CharmAnalysis) already applied in the ntuples: mass/pT ranges remove almost nothing",
         "S/B grows mainly with the signal window and the BDT cut"],
        table=[["step", "D⁰ S", "D⁰ B", "D⁰ S/B", "D⁺ S", "D⁺ B", "D⁺ S/B"],
               ["mass range, 3 < pT < 30 GeV", "9,534", "3,170,411", "0.0030", "8,727", "11,078,038", "0.00079"],
               ["D* veto (D⁺)", "–", "–", "–", "8,665", "10,846,291", "0.00080"],
               ["signal window", "8,825", "1,404,246", "0.0063", "7,912", "1,910,643", "0.0041"],
               ["score ≥ WP80", "7,098", "122,334", "0.058", "6,391", "47,687", "0.134"],
               ["score ≥ WP70", "6,187", "55,933", "0.111", "5,600", "21,196", "0.264"],
               ["score ≥ WP50", "4,432", "11,360", "0.39", "3,972", "5,358", "0.74"]],
        section="2")

    # ---------------- (3)
    add("(3) Classifier inputs: low-level vs high-level",
        ["Task: true D vs combinatorial background, candidate by candidate (D⁰: 2 tracks, D⁺: 3 tracks)",
         "II. High-level = 18 features of the nominal BDT \"H\": log pT, |η|, vertex χ², L_xy, IP significance, |z0|, pointing, |cos θ*|,",
         "    daughter min/max |d0|, max |z0|, max |η|, PID: dE/dx, log p, n(dE/dx hits) of K; max dE/dx, min log p of π",
         "I. Low-level = daughter tracks (11 numbers each: log pT, η, Δφ/Δη to D, log p, dE/dx (+missing flag), n hits, d0, z0, K/π role) + L_xy",
         "    vertex-fit quantities (χ², cos θ*, pointing) NOT given — the model must infer them",
         "Same training sample for all variants: signal = truth-matched D in forced-decay MC, background = combinatorial in minimum-bias MC;",
         "training/evaluation inside the signal window (no mass learning); 5-fold by event, out-of-fold metrics"],
        section="3")
    add("(3) Inputs requested in the plan vs inputs used",
        ["max ΔR between daughters is NOT among the 18 BDT \"H\" features (distribution below: separating, partly redundant with pT and L_xy)",
         "Daughter charge q enters only through the K/π role (opposite-sign Kπ, D⁺: Kππ charge pattern fixed by the finder)",
         "Rapidity y of tracks: η used (y ≈ η for pions/kaons at these momenta); daughter pT enters low-level directly, high-level only via log p"],
        table=[["plan (sheet)", "low-level (I) — used as", "high-level (II) — used as"],
               ["y, pT, p, φ", "η, log pT, log p, Δφ/Δη to D", "|η| of D, log pT of D, min log p of π, log p of K"],
               ["q (charge)", "K/π role flag", "– (fixed by finder)"],
               ["dE/dx", "dE/dx + missing flag, n hits", "dE/dx K, max dE/dx π, n(dE/dx hits) K"],
               ["d0, L_xy", "d0, z0 per track; L_xy", "min/max |d0|, max |z0|, L_xy, IP significance"],
               ["max ΔR", "– (implicit in Δφ/Δη)", "missing in H (shown below)"],
               ["χ², cos θ*", "not given (must be inferred)", "vertex χ², |cos θ*|, pointing"],
               ["daughter pT", "log pT", "not explicit (log p, max |η|)"]],
        section="3")
    add("(3) max ΔR between daughters: signal vs background (not used in BDT \"H\")",
        ["Forced-decay MC truth-matched signal vs minimum-bias MC background, signal window, 3 < pT < 30 GeV"],
        images=img("k3_maxdR", d=KRT), section="3")
    add("(3) Models: gradient-boosted trees vs neural networks",
        ["Deep Sets for tracks: unordered set, interchangeable pions in D⁺, same net for 2 and 3 tracks, weight sharing",
         "NN training: rank→Gauss inputs, weighted BCE, AdamW (lr 2·10⁻³), one-cycle, batch 4096, ≤ 30 epochs, early stopping",
         "Reference: production BDT \"H\" (trained vs data sidebands) evaluated on the same candidates"],
        table=[["", "BDT (XGBoost)", "Neural network"],
               ["high-level (18)", "≤ 1500 trees, depth 6, lr 0.05 (nominal)", "MLP 18→128→128→64→1, SiLU, ≈ 27.5 k params"],
               ["low-level (tracks)", "flat vector K, π (pT-ordered) + L_xy", "Deep Sets: φ 11→64→64→64, sum⊕max, ρ 129→128→64→1, ≈ 34 k params"]],
        section="3")
    roc = sorted(glob.glob(os.path.join(P3, "runs", "**", "roc_*.png"), recursive=True))
    other = [sorted(glob.glob(os.path.join(P3, "runs", "**", f"{k}_*.png"), recursive=True)) for k in ("scores", "importance", "sculpt")]
    if roc:
        for i in range(0, len(roc), 2):
            add("(3) BDT vs NN × low/high-level — ROC", [fill("{{P3_SUMMARY}}")], images=roc[i:i + 2], section="3")
        for imgs, t in zip(other, ("score distributions", "feature importance (high-level)", "mass sculpting test")):
            per = 1 if t != "feature importance (high-level)" else 2          # 2×2 panele → jeden obraz na slajd (czytelność)
            for i in range(0, len(imgs), per):
                sp = "" if per == 2 else (" — D⁺" if "Dplus" in imgs[i] else " — D⁰")
                add(f"(3) BDT vs NN — {t}{sp}", [], images=imgs[i:i + per], section="3")
    else:
        add("(3) BDT vs NN × low/high-level — results",
            ["Training running on lxplus (4 variants × D⁰/D⁺; NN up to ~10 h)",
             "Metrics: ROC/AUC, ε_B at ε_S = 90 % (nominal WP), ε_S at fixed background, score shapes, feature importance, mass sculpting",
             fill("{{P3_SUMMARY}}")],
            section="3")

    # ---------------- (4)
    add("(4) BDT score: signal vs background, working point → GOOD / BAD",
        ["Signal from forced-decay MC, background from minimum-bias MC; three pT-dependent WP70 thresholds",
         "At ε_S = 0.70: ε_B = 0.040 for D⁰ (rejection ×25), 0.0094 for D⁺ (×106)",
         "GOOD = score ≥ WP (pT-bin threshold), BAD = below"],
        images=img("p4_score_D0", "p4_score_Dplus"), section="4")
    add("(4) D⁰ invariant mass, GOOD vs BAD — minimum-bias MC",
        ["GOOD: clear peak; BAD: flat background (signal shown ×100)"],
        images=img("p4_mass_goodbad_D0_incl_mc"), section="4")
    add("(4) D⁰ invariant mass, GOOD vs BAD — forced-decay MC",
        ["Forced decay after WP70: window = 64 % signal + 20 % K↔π reflections"],
        images=img("p4_mass_goodbad_D0_signal_mc"), section="4")
    add("(4) D⁺ invariant mass, GOOD vs BAD — minimum-bias MC",
        ["GOOD: peak on a reduced background; BAD: flat background"],
        images=img("p4_mass_goodbad_Dplus_incl_mc"), section="4")
    add("(4) D⁺ invariant mass, GOOD vs BAD — forced-decay MC",
        ["Forced decay after WP70: window = 76 % signal + 21 % partially reconstructed"],
        images=img("p4_mass_goodbad_Dplus_signal_mc"), section="4")
    add("(4) D* Δm, GOOD vs BAD (by the D⁰-leg score, WP90) — minimum-bias MC",
        ["D* has no own classifier: GOOD/BAD defined by the BDT score of the D⁰ daughter (WP90, same as in data and in pairs)",
         "{K4DS_MB}"],
        images=img("k4_dstar_goodbad_mb_nd", d=KRT), section="4")
    add("(4) D* Δm, GOOD vs BAD (by the D⁰-leg score, WP90) — forced-decay MC",
        ["{K4DS_FD}"],
        images=img("k4_dstar_goodbad_fd_nd", d=KRT), section="4")
    add("(4) Mass vs score in 2D — D⁰, signal and background (minimum-bias MC)",
        ["Signal concentrated at m₀ and high score; background at score ≈ 0 and flat in mass → no visible sculpting"],
        images=img("p4_2d_score_mass_D0_mc"), section="4")
    add("(4) Mass vs score in 2D — D⁺, signal and background (minimum-bias MC)",
        ["Same picture for D⁺: no visible mass sculpting by the score"],
        images=img("p4_2d_score_mass_Dplus_mc"), section="4")
    add("(4) Δm vs D⁰-leg score in 2D — D*, signal (forced decay) and background (minimum bias)",
        ["True D* concentrated at Δm ≈ 145.4 MeV and high D⁰-leg score; background spread in Δm and at low score"],
        images=img("k4_2d_dm_score_dstar_mc", d=KRT), section="4")

    # ---------------- (5)
    add("(5) Data 428770: D⁰ invariant mass, GOOD vs BAD",
        ["D⁰: 272 k GOOD out of 8.58 M candidates", "Window/sideband ratio: GOOD 1.37, BAD 1.03 → peak only in GOOD"],
        images=img("p5_data_mass_goodbad_D0"), section="5")
    add("(5) Data 428770: D⁺ invariant mass, GOOD vs BAD",
        ["D⁺: 242 k GOOD out of 32.3 M candidates", "Window/sideband ratio: GOOD 1.46, BAD 0.95 → peak only in GOOD"],
        images=img("p5_data_mass_goodbad_Dplus"), section="5")
    for e in [e for e in DF if slot_of(e) == "datamc"]:
        add(f"(5) {e['title']}", [e["caption"]] if e["caption"] else [], images=[e["file"]], section="5")
    add("(5) Data 428770: D* Δm and mass–score 2D",
        ["D* peak visible only when the D⁰ leg is GOOD (WP90): 144 k GOOD among 906 k candidates",
         "2D: data show the same signal region (m₀, high score) as MC"],
        images=img("p2_dstar_dm_data_428770", "p5_data_2d_score_mass_D0"), section="5")
    add("(5) Data 427929 (main-result run): D⁰ invariant mass, GOOD vs BAD",
        ["{K5_D0}"], images=img("k5_427929_mass_goodbad_D0", d=KRT), section="5")
    add("(5) Data 427929 (main-result run): D⁺ invariant mass, GOOD vs BAD",
        ["{K5_DP}"], images=img("k5_427929_mass_goodbad_Dplus", d=KRT), section="5")
    add("(5) Data 427929: mass vs score in 2D — D⁰ and D⁺",
        ["Same signal region (m₀, high score) as in MC and in run 428770"],
        images=img("k5_427929_2d_score_mass_D0", "k5_427929_2d_score_mass_Dplus", d=KRT), section="5")
    add("(5) Data 427929: D* Δm, GOOD vs BAD by the D⁰-leg score",
        ["{K5_DS}"], images=img("k5_427929_dstar_dm_goodbad", d=KRT), section="5")

    # ---------------- (6)
    add("(6) Production topologies of cc̄ pairs",
        ["Flavour creation: c and c̄ back-to-back (|Δφ| ≈ π)", "Gluon splitting g → cc̄: collinear (small Δφ, Δη) — typical inside jets",
         "Double parton scattering (DPS): two independent scatterings → uncorrelated, also like-sign (cc, c̄c̄) pairs"],
        images=img("p6_schematic_topologies", d=P6), section="6")
    add("(6) Pair topologies in Pythia 8 A3 (generator level, analysis acceptance)",
        ["15,556 D⁰/D⁺ pairs (pT 3–30 GeV, |η| < 2.5): 12,469 OS, 3,087 LS (f_LS = 0.20 generator truth, before wrong-sign D⁰ correction)",
         "OS: two peaks in |Δφ| (0: gluon splitting, π: back-to-back); LS ≈ flat, like independent pairs (DPS proxy)",
         "P(|Δη| < 1): OS 0.60, LS 0.38, independent pairs 0.37",
         "Not comparable directly with data: f_LS(data, 428770) = 0.33 ± 0.14 vs Pythia after sign correction 0.27 ± 0.03 (+0.4σ)"],
        images=img("p6_gen_dphi", "p6_gen_deta", d=P6), section="6")
    add("(6) Pair kinematics: Δφ–Δη and pT₁ vs pT₂ (Pythia 8 A3)",
        ["OS pairs: collinear peak at (0, 0) and away-side ridge", "Mean pT asymmetry of OS pairs: 0.18"],
        images=img("p6_gen_dphi_deta_os", "p6_gen_pt1pt2", d=P6), section="6")
    add("(6) Training on pairs of D mesons: pair network",
        ["Per-leg vector (18 features + species + single-D BDT logit) → shared MLP encoder 20→128→128→64",
         "Symmetric combination (sum ⊕ |difference|) → head → P(both legs true); 5 folds",
         "AUC (DpD0): single-D score sum 0.945 → pair BDT 0.973 → pair NN 0.973; background at ε_S = 50 %: 0.012 → 0.0047"],
        images=img("p6_pnet_arch", "p6_pnet_roc", d=P6), section="6")
    add("(6) Pair network: score and expected gain on μ",
        ["Score distributions: 31.7 k signal and 37.6 k background pairs (v1, MC background)",
         "Expected σ_μ / σ_nominal: single-D sum 0.92, pair BDT 0.81, pair NN (data-sideband background) 0.79",
         fill("Measurement with the pair network on data (ubpnet): {{UBPNET_MU}}")],
        images=img("p6_pnet_score", "p6_pnet_sigma", d=P6), section="6")

    # ---------------- (7) wyniki
    for slot, lab in (("dstar", "D* normalisation fit on data"), ("pairs", "pair fit on data"), ("excess", "correlated pair excess (no fit model)")):
        for e in [e for e in DF if slot_of(e) == slot]:
            add(f"(7) {e['title']}", [e["caption"]] if e["caption"] else [lab], images=[e["file"]], section="7")
    add("(7) μ = (pairs/D*)_data / (pairs/D*)_Pythia — normalisation from truth (mb_nd)",
        ["New normalisation (decision 19.11): all true pairs in mb_nd × direct efficiency maps → 5× more effective MC pairs",
         "Systematics without N_D* (428770) reduced ±0.61 → ±0.36 (MC pair statistics ±0.46 → ±0.18)",
         "Statistical errors include the event-correlation (jackknife) factor: ×1.18 (428770), ×1.13 (427929)",
         "Significance vs μ = 0 (background-only): 428770 Z = 4.05 (≈ 3.4 with event correlations); 427929 Z = 1.60",
         "  428770 vs Pythia: (1.27 − 1)/σ_tot,↓ with σ_tot,↓ = 0.61 from the quoted errors; 427929: BUDZET (syst at μ = 1) / TRUTHNORM (syst at μ̂)"],
        table=[["run", "μ", "stat", "syst", "N(D*)", "vs Pythia (μ = 1)"],
               ["428770 (J12-dominated)", "1.27", "+0.41/−0.40", "±0.36", "±0.28", "+0.4σ (trigger-biased, see next slides)"],
               ["427929 (clean min. bias) — MAIN RESULT", "0.50", "+0.41/−0.37", "+0.28/−0.19", fill("±0.21 {{NDS427_NOTE}}"), "−0.6σ … −0.9σ (consistent)"]],
        section="7")
    add("(7) Summary of μ measurements vs Pythia 8 A3",
        ["Bands: main result 427929 (green stat., yellow stat ⊕ syst ⊕ N(D*)); red line: Pythia (μ = 1)",
         "Clean minimum bias consistent with Pythia; the 428770 excess sits only in J12-triggered events"],
        images=[f for f in [os.path.join(FX, "mu_summary.png")] if os.path.exists(f)], section="7")
    ut = unc_table()
    if ut:
        add("(7) Uncertainty breakdown of the main result (run 427929, μ = 0.50)",
            ["Relative uncertainties Δμ/μ in %, sorted; common μ (channel weights 1/σ²); upper tune row = ALICE fragmentation (k = 0.71)",
             "Largest syst.: N(D*) (Δm fit variants, peak width ×1.40) and MC pair statistics (conservative — still direct-method n_eff)"],
            table=ut, section="7")
    add("(7) Uncertainty ranking — run 427929",
        ["Data statistics dominate; asymmetric: fragmentation measured by ALICE moves μ up only"],
        images=[f for f in [os.path.join(FX, "unc_ranking.png")] if os.path.exists(f)], section="7")
    add("(7) 428770: μ per trigger class — the excess sits in J12-triggered events",
        ["Classes from trigger bits (disjoint); r_c = share of N_D* from Δm fit; background (fake legs, SB×SB) per class",
         "J12only vs rest: +3.6σ (≈ +3.0σ with event correlations); class consistency χ² = 15.7/3, p = 0.001",
         "Negative μ in small classes is unphysical (Pythia = 1) → null test with random classes (next slide)",
         "Note: \"other\" ≈ 40 % AFP + untagged; MBTS_noj20 mostly hard-track sptrk_pt4+ → no clean subset of 428770"],
        table=[["class", "share of N_D*", "μ (truth norm.)", "68 % interval"],
               ["J12only", "0.769 ± 0.019", "1.37", "[1.05; 1.71]"],
               ["MBTS_j20", "0.135 ± 0.015", "−0.51", "[−0.95; −0.05]"],
               ["MBTS_noj20", "0.010 ± 0.004", "1.28", "[0.40; 2.87]"],
               ["other", "0.084 ± 0.012", "−0.59", "[−1.31; 0.11]"]],
        section="7")
    add("(7) Null test: random event classes with the same sizes",
        ["Event classes shuffled within strata (event has legs × has D* candidate); class sizes identical to the trigger classes",
         "If the method is unbiased: μ_c ≈ μ_nominal (1.19 on the 05.11 sample) in every class",
         fill("Random class of J12 size: μ = {{NULL_J12}} (real J12only: 1.63)"),
         fill("Rest (small classes): μ = {{NULL_MU_REST}} (real: −0.53)"),
         fill("J12 − rest: real +3.7σ vs null {{NULL_DIFF}} → {{NULL_VERDICT}}")],
        section="7")
    add("(7) Clean run 427929 (MBTS minimum bias, prescale 1) — main result",
        ["Same selection, models, templates and 2D fit as 428770; N(D*) from Δm fit (2046 ± 485); normalisation from truth (mb_nd)",
         "Event jackknife (200 groups): σ_jk/σ = 1.13 applied to the statistical error",
         "Main result (decision 19.11): μ = 0.50 +0.41/−0.37 (stat) +0.28/−0.19 (syst) ± 0.21 (N(D*)) → consistent with Pythia (−0.6 … −0.9σ)",
         "Upper syst includes charm fragmentation measured in pp (ALICE: fewer D mesons per c quark, k = 0.71 → μ/k ≈ 0.70)",
         "Difference to 428770: ≈ −1.3σ",
         "  cross-check: direct (old) normalisation gives μ = 0.57 ± 0.45 — not the quoted result"],
        section="7")
    add("(7) D* peak width in data vs pT — origin of the N(D*) uncertainty",
        ["427929: Δm peak 1.40× wider than MC overall, but only at low pT (large background); at 8–30 GeV 1.12 [0.99; 1.26]",
         "428770 (all and MBTS-triggered): compatible with MC width in every bin",
         "Pre-registered test (MBTS events of 428770: 0.96 ± 0.24) is inconclusive → width kept free, N(D*) ±41 %",
         "Variant: width fixed to MC → N(D*) = 1417 → μ ≈ 0.72 (not nominal)"],
        images=[f for f in [os.path.join(FX, "width_vs_pt.png")] if os.path.exists(f)], section="7")
    add("(7) Pair topology |Δη| in data (428770): λ = DPS fraction above Pythia",
        ["h(λ) = λ·h_DPS + (1 − λ)·h_Pythia; λ = 0 Pythia, λ = 1 pure DPS, \"excess = DPS\" ⇒ λ = 0.33",
         "Signal more collinear than Pythia (|Δη| < 1 enhanced) — opposite to a DPS excess",
         "Pure DPS disfavoured at 2.8σ, \"excess = DPS\" at 2.1σ; consistent with jet-induced gluon-splitting enhancement"],
        table=[["channel", "λ̂", "68 %", "Z(Pythia)", "Z(DPS)", "Z(excess = DPS)"],
               ["D⁰D⁰", "−0.60", "[−1.50; 0.72]", "0.47", "1.16", "0.74"],
               ["D⁺D⁰", "−1.50 (edge)", "[−1.50; −0.95]", "2.01", "2.71", "2.28"],
               ["D⁺D⁺", "0.18", "[−1.40; 2.50]", "0.12", "0.39", "0.08"],
               ["combined", "−1.11", "[−1.50; −0.51]", "1.71", "2.80", "2.13"]],
        section="7")
    add("(7) Charm yield vs event multiplicity (427929, MBTS)",
        ["Ŷ ∝ x^α with x = n_trk/⟨n_trk⟩; Δα_fit = α_data − α_MC (same mass fit in data and MC)",
         "D*: Δα = 0.03 ± 0.22 (Z = 0.1); D⁺: Δα = −0.00 ± 0.26 → no extra charm–multiplicity correlation beyond Pythia",
         "→ \"event-weight model\" systematic unchanged; D⁰ data fit not converged (open)",
         "Fit closure in MC (G2) fails for D⁺/D*: α_fit − α_true ≈ +0.15 → both data and MC fitted the same way"],
        section="7")
    add("(7) Trigger systematic (pre-registered rules, TRIG_A/SYST_TRIG.md)",
        ["d = (μ_nominal − μ_clean)/μ_nominal from two measurements:",
         "  d₁ (427929 vs 428770, jackknife ×1.13) = +0.65 ± 0.46 (1.4σ);  d₂ (J12only vs rest within 428770) = +1.39 ± 0.44 (3.2σ)",
         "  combined d = +1.04 ± 0.32 (Z = 3.2); J12 emulation predicted only d_emu = 0.16 (κ = 1.24) → data 2.8σ above emulation",
         "Status: effect → main result moved to the clean run 427929 (decision 19.11)",
         fill("Null test outcome: {{NULL_VERDICT}}")],
        section="7")
    add("Summary and next steps",
        ["Full chain for D–D pairs in ATLAS 2022 minimum bias: BDT selection, 2D pair fit, D* normalisation, systematics",
         "428770: μ = 1.27 (Z ≈ 3.4) — but the excess sits in J12 (jet) triggered events; topology prefers collinear (jet-like) pairs",
         "MAIN RESULT — 427929 (clean minimum bias): μ = 0.50 +0.41/−0.37 (stat) +0.28/−0.19 (syst) ± 0.21 (N_D*) → consistent with Pythia; σ_eff > ~6 mb (95 % CL)",
         "Variants: D* width from MC → μ ≈ 0.72; ALICE fragmentation → μ ≈ 0.70; both → ≈ 1.0 — no extra DPS in clean minimum bias",
         "428770 reported as trigger-biased observation (null test: J12 Z = 3.25, J12 − rest Z ≈ 8; collinear topology)",],
        section="End")
    # liczby do slajdów uzupełnień wg kartki (figs_extra/kartka/kartka_numbers.json)
    kp = os.path.join(KRT, "kartka_numbers.json")
    if os.path.exists(kp):
        K = json.load(open(kp)); tok = {}
        for nm, key in (("mb_nd", "K4DS_MB"), ("fd_nd", "K4DS_FD")):
            k = K.get(f"dstar_{nm}", {})
            if k:
                es = f"ε_S(GOOD) = {k['eff_S_good']:.2f}, " if k.get("eff_S_good") is not None else ""
                tok[key] = f"{k['n']:,} D⁰←D* candidates, {k['n_sig']:,} true D*; {es}ε_B(GOOD) = {k['eff_B_good']:.3f}"
        for sp, key, nm in (("D0", "K5_D0", "D⁰"), ("Dplus", "K5_DP", "D⁺")):
            k = K.get(f"data427929_{sp}", {})
            if k:
                tok[key] = (f"{nm}: {k['n_good']:,} GOOD out of {k['n_all']:,} candidates (WP70, same thresholds as 428770); window/sideband ratio: "
                            f"GOOD {k['good_window_over_sb_scaled']:.2f}, BAD {k['bad_window_over_sb_scaled']:.2f}")
        k = K.get("data427929_dstar", {})
        if k:
            tok["K5_DS"] = f"{k['n_good']:,} GOOD among {k['n']:,} D⁰←D* candidates (log scale); D* peak at 145.4 MeV after the D⁰-leg cut"
        for e in S:
            e["bullets"] = [tok.get(b[1:-1], b) if b.startswith("{K") and b.endswith("}") else b for b in e["bullets"]]
    S = [dict(e, bullets=[b for b in e["bullets"] if not (b.startswith("{K") and b.endswith("}"))]) for e in S]
    return S


# ------------------------------------------------------------------------------------------------ PPTX
def build_pptx(S, out):
    from pptx import Presentation
    from pptx.dml.color import RGBColor
    from pptx.util import Emu, Inches, Pt
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]
    W, H = 13.333, 7.5
    for i, s in enumerate(S, 1):
        sl = prs.slides.add_slide(blank)
        bar = sl.shapes.add_shape(1, 0, 0, prs.slide_width, Inches(0.95))
        bar.fill.solid(); bar.fill.fore_color.rgb = RGBColor(0x1F, 0x3A, 0x5F); bar.line.fill.background()
        tb = sl.shapes.add_textbox(Inches(0.3), Inches(0.12), Inches(W - 0.6), Inches(0.75)).text_frame
        tb.word_wrap = True; p = tb.paragraphs[0]; p.text = fill(s["title"]); p.font.size = Pt(26 if i == 1 else 22); p.font.bold = True
        p.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        nb = len(s["bullets"]); has_img = bool(s["images"]); top = 1.1
        nl = sum(max(1, -(-len(fill(b)) // 150)) for b in s["bullets"])     # zawijanie długich punktów (≈150 znaków/wiersz przy 14 pt)
        if nb:
            bh = 0.30 * nl + 0.12 * nb + 0.2 if (has_img or s["table"]) else H - 1.8
            tf = sl.shapes.add_textbox(Inches(0.4), Inches(top), Inches(W - 0.8), Inches(bh)).text_frame
            tf.word_wrap = True
            for j, b in enumerate(s["bullets"]):
                p = tf.paragraphs[0] if j == 0 else tf.add_paragraph()
                txt = fill(b); ind = txt.startswith("    ") or txt.startswith("  ")
                p.text = ("" if ind else "• ") + txt.strip(); p.font.size = Pt(14 if (has_img or s["table"]) else 18); p.level = 1 if ind else 0
            top += bh
        if s["table"]:
            rows, cols = len(s["table"]), len(s["table"][0])
            th = min(0.42 * rows, H - top - 0.6)
            t = sl.shapes.add_table(rows, cols, Inches(0.4), Inches(top + 0.1), Inches(W - 0.8), Inches(th)).table
            for r, row in enumerate(s["table"]):
                for c, v in enumerate(row):
                    cell = t.cell(r, c); cell.text = fill(str(v))
                    for pp in cell.text_frame.paragraphs:
                        pp.font.size = Pt(12); pp.font.bold = r == 0
            top += th + 0.2
        if has_img:
            n = len(s["images"]); avail_h = H - top - 0.5; slot_w = (W - 0.8) / n
            for k, f in enumerate(s["images"]):
                h_, w_ = mpimg.imread(f).shape[:2]; ar = w_ / h_
                w = min(slot_w - 0.2, avail_h * ar); hh = w / ar
                x = 0.4 + k * slot_w + (slot_w - w) / 2
                sl.shapes.add_picture(f, Inches(x), Inches(top + 0.05), Emu(int(Inches(w))), Emu(int(Inches(hh))))
        ft = sl.shapes.add_textbox(Inches(0.3), Inches(H - 0.45), Inches(W - 0.6), Inches(0.35)).text_frame
        p = ft.paragraphs[0]; p.text = f"{TAG}   ·   {DATE}   ·   {i}/{len(S)}"
        p.font.size = Pt(10); p.font.italic = True; p.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    prs.save(out)


# ------------------------------------------------------------------------------------------------ PDF (matplotlib)
def build_pdf(S, out):
    W, H = 13.333, 7.5
    with PdfPages(out) as pdf:
        for i, s in enumerate(S, 1):
            fig = plt.figure(figsize=(W, H))
            fig.patches.append(plt.Rectangle((0, 1 - 0.95 / H), 1, 0.95 / H, transform=fig.transFigure, color="#1F3A5F", zorder=0))
            fig.text(0.025, 1 - 0.5 / H, fill(s["title"]), color="white", fontsize=20 if i == 1 else 17, weight="bold", va="center")
            y = 1 - 1.2 / H; small = bool(s["images"] or s["table"])
            fs = 11 if small else 15; dy = (0.36 if small else 0.5) / H
            for b in s["bullets"]:
                txt = fill(b); ind = txt.startswith("  ")
                fig.text(0.04 if ind else 0.03, y, ("" if ind else "• ") + txt.strip(), fontsize=fs, va="top", wrap=True)
                y -= dy * max(1, -(-len(txt) // (170 if small else 130)))
            if s["table"]:
                rows = s["table"]; th = min(0.4 * len(rows), y * H - 0.7) / H
                ax = fig.add_axes([0.03, y - th - 0.01, 0.94, th]); ax.axis("off")
                L = [max(len(fill(str(r[c]))) for r in rows) + 2 for c in range(len(rows[0]))]
                tab = ax.table(cellText=[[fill(str(v)) for v in r] for r in rows[1:]], colLabels=rows[0], loc="center", cellLoc="center",
                               colWidths=[x / sum(L) for x in L])
                tab.auto_set_font_size(False); tab.set_fontsize(10); tab.scale(1, 1.25)
                for (r, c), cell in tab.get_celld().items():
                    if r == 0:
                        cell.set_text_props(weight="bold"); cell.set_facecolor("#DDE5EE")
                y -= th + 0.03
            if s["images"]:
                n = len(s["images"]); bottom = 0.07; hh = y - bottom - 0.01; ww = 0.94 / n
                for k, f in enumerate(s["images"]):
                    ax = fig.add_axes([0.03 + k * ww, bottom, ww - 0.01, hh]); ax.imshow(mpimg.imread(f)); ax.axis("off")
            fig.text(0.025, 0.02, f"{TAG}   ·   {DATE}   ·   {i}/{len(S)}",
                     fontsize=9, style="italic", color="#555555")
            pdf.savefig(fig); plt.close(fig)


def main():
    S = slides()
    pdf = os.path.join(HERE, "deck.pdf"); build_pdf(S, pdf); print("PDF:", pdf, len(S), "slajdów")
    try:
        build_pptx(S, os.path.join(HERE, "deck.pptx")); print("PPTX:", os.path.join(HERE, "deck.pptx"))
    except ImportError:
        print("brak python-pptx — tylko PDF", file=sys.stderr)
    ph = sorted({m for s in S for t in [s["title"], *s["bullets"], *(str(v) for r in (s["table"] or []) for v in r)]
                 for m in re.findall(r"\{\{(\w+)\}\}", t) + re.findall(r"\[(\w+): in progress\]", t) if m not in NUM})
    print("placeholdery bez wartości:", ph)
    json.dump({"slides": [s["title"] for s in S], "missing_placeholders": ph}, open(os.path.join(HERE, "deck_index.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
