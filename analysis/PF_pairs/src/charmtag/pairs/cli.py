"""CLI: `charmtag pairs ...` i `charmtag fit2d ...` (wywoływane z charmtag/cli.py lub `python -m charmtag.pairs.cli`).

  charmtag pairs  -c config/mac_v11.yaml [--unblind] [--norm direct|frozen|truth]
  charmtag fit2d  -c config/mac_v11.yaml [--unblind] [--ntoys 50] [--asimov-only]
  charmtag fit2d  --closure [~/Documents/charm/unblinding_freeze_v20260927_01]     # test domknięcia z 09-27
Domyślnie ślepo: centrum 2D (obie nogi w oknie) nie trafia do histogramów i μ z danych nie jest liczone.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time

import numpy as np

log = logging.getLogger("charmtag")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(prog="charmtag pairs|fit2d", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["pairs", "fit2d"])
    ap.add_argument("-c", "--config", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--h5-dir", default=None)
    ap.add_argument("--unblind", action="store_true", help="odsłonięcie centrum 2D (pełna płaszczyzna, μ, Z)")
    ap.add_argument("--norm", default="direct", choices=["direct", "frozen", "truth"])
    ap.add_argument("--ntoys", type=int, default=0)
    ap.add_argument("--asimov-only", action="store_true")
    ap.add_argument("--closure", nargs="?", const="~/Documents/charm/unblinding_freeze_v20260927_01", default=None)
    ap.add_argument("--stageb-hists", default="~/Documents/charm/odsloniecie_wynik_v20260927_01/results/stageB_hists.npz")
    a = ap.parse_args(argv)

    from .. import config as C
    from ..cli import _setup_log
    ov = {"paths": {}}
    if a.out_dir:
        ov["paths"]["out_dir"] = a.out_dir
    if a.h5_dir:
        ov["paths"]["h5_dir"] = a.h5_dir
    cfg = C.load(a.config, ov)
    from .build import pdir
    out = pdir(cfg)
    _setup_log(out)
    t0 = time.time()
    log.info("charmtag %s unblind=%s out=%s", a.command, a.unblind, out)
    if a.command == "pairs":
        from .build import cmd_pairs
        cmd_pairs(cfg, unblind=a.unblind, norm_mode=a.norm)
    else:
        from . import fit2d as F
        from .model import Model, load_inputs
        if a.closure:
            r = F.run_closure(os.path.expanduser(a.closure), os.path.join(out, "closure"), os.path.expanduser(a.stageb_hists), ntoy=a.ntoys)
            print(json.dumps({"PASS": r["PASS"], **r["check"]}, indent=1))
            log.info("koniec fit2d --closure (%.0f s)", time.time() - t0)
            return 0 if r["PASS"] else 3
        I = load_inputs(os.path.join(out, "pair_inputs"))
        M = Model(I)
        if a.asimov_only:
            rep = {"asimov_nominal": F.asimov(M)}
            F._js(os.path.join(out, "fit2d_asimov.json"), rep)
            print(json.dumps(rep, indent=1, default=float))
        elif a.unblind:
            p = os.path.join(out, "pairs_full.npz")
            if not os.path.exists(p):
                raise SystemExit(f"brak {p}: najpierw `charmtag pairs --unblind` (odsłonięcie wymaga jawnej flagi na obu etapach)")
            H = dict(np.load(p))
            rep = F.run_unblind(M, H, out, ntoy=a.ntoys)
            write_report(out, I, rep, blind=False)
        else:
            H = dict(np.load(os.path.join(out, "pairs_blind.npz")))
            rep = F.run_blind(M, H, out, ntoy=a.ntoys)
            write_report(out, I, rep, blind=True)
    log.info("koniec %s (%.0f s)", a.command, time.time() - t0)
    return 0


def write_report(out, I, rep, blind=True):
    from . import CHANNELS
    L = [f"# Pary D–D — fit 2D ({'ŚLEPO' if blind else 'ODSŁONIĘTE'})", ""]
    m = I.meta
    L.append(f"- Normalizacja: {m.get('norm_mode')}; N_D* dane = {m.get('dstar_fit_data', {}).get('N', float('nan')):.0f}; "
             f"lumi = {m.get('lumi_fraction', float('nan')):.3f} pełnego 428770 (skraca się w μ).")
    L.append(f"- Oczekiwane pary (pełna płaszczyzna, przed kategoriami): " + ", ".join(f"{k} {v:.1f}" for k, v in I.n_pair.items()))
    L.append("")
    if blind:
        L += ["| kanał | dev/ndf (bez centrum) | ramiona dev/bin | Z Asimov (nuisance z ramion) | σ_μ | S oczek. |", "|---|---|---|---|---|---|"]
        for ch in CHANNELS:
            r = rep["channels"][ch]
            if "skipped" in r:
                L.append(f"| {ch} | {r['skipped']} | | | | |")
                continue
            L.append(f"| {ch} | {r['blind_region_dev_ndf']:.2f} | {r['arms_dev_per_bin']:.2f} | {r['asimov_from_arms']['Z']:.2f} | "
                     f"{r['asimov_from_arms']['sigma_mu']:.2f} | {r['S_expected_full_plane']:.1f} |")
        L.append(f"\nŁącznie (kwadratura) Z_Asimov = {rep['Z_asimov_quadrature']:.2f}σ. Przesunięcia μ (Asimov, prawda zmieniona, fit nominalny):")
        for g, v in rep["asimov_shifts"].items():
            L.append(f"- {g}: " + (", ".join(f"{k} {x:+.3f}" for k, x in v.items()) if isinstance(v, dict) else str(v)))
        for s in ("0.3", "0.5"):
            a = rep.get(f"asimov_partner_sigma{s}")
            if a:
                L.append(f"- więz κ partnera σ={s}: Z_Asimov " + ", ".join(f"{ch} {a[ch]['Z']:.2f}" for ch in CHANNELS) + f", kwadr. {a['Z_quadrature']:.2f}")
    else:
        L += ["| kanał | μ | 68% | Z | dev/ndf | κ-kanał μ (Z) | bez tagu |", "|---|---|---|---|---|---|---|"]
        for ch in CHANNELS:
            r = rep["channels"][ch]
            k = r.get("crosscheck_kchannel", {})
            nt = r.get("mu_no_tag", [np.nan, np.nan])
            L.append(f"| {ch} | {r['mu_hat']:.2f} | {r['mu_profile_68'][0]:.2f}–{r['mu_profile_68'][1]:.2f} | {r['Z']:.2f} | {r['dev_ndf']:.2f} | "
                     f"{k.get('mu_hat', np.nan):.2f} ({k.get('Z', np.nan):.1f}) | {nt[0]:.2f} ± {nt[1]:.2f} |")
        c = rep["combined"]
        L.append(f"\nWspólne μ = {c['mu_common']:.2f} [{c['mu_common_68'][0]:.2f}, {c['mu_common_68'][1]:.2f}], Z = {c['Z']:.2f}σ.")
        if "syst_profiled" in rep:
            s = rep["syst_profiled"]["combined"]
            L.append(f"- z systematykami kształtu (profil): wspólne μ = {s['mu_common']:.2f}, Z = {s['Z']:.2f}σ")
        for sg, row in rep.get("kappa_partner", {}).items():
            L.append(f"- więz κ partnera σ={sg}: wspólne μ = {row['combined']['mu_common']:.2f}, Z = {row['combined']['Z']:.2f}σ")
    with open(os.path.join(out, "fit2d_blind.md" if blind else "fit2d_unblind.md"), "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
