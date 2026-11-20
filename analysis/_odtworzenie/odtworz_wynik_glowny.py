#!/usr/bin/env python3
"""Odtwarza z samego repo fit wyniku głównego (run 427929, normalizacja z prawdy mb_nd) i porównuje z zapisanym wynikiem.

Uruchamia NIEZMIENIONY skrypt `PF2_lxplus/runs/lxplus/unblind/truthnorm/truthnorm_run.py` w drzewie złożonym z repo
(zob. zloz_charm_work.py). Wejścia: `pair_inputs.json` + tablice z `*.npz.json` (results/), czynniki f_ch z
`GEN_mcstat/out/truth_mbnd.json`, minimalizator `Q_lxplus/ub_jk/ubjk_run.py:replica_fit`.

  python analysis/_odtworzenie/odtworz_wynik_glowny.py                      # 427929, kod zamrożony (analysis/PF_pairs/src)
  python analysis/_odtworzenie/odtworz_wynik_glowny.py --charmtag repo      # ten sam fit bieżącym pakietem charmtag z repo
  python analysis/_odtworzenie/odtworz_wynik_glowny.py --what nom428        # 428770 (próbka 05.11, baza N_D* 5584)
  ... --f-one                                                               # kontrola f ≡ 1 (normalizacja „direct”)

Czas: kilka minut na jeden fit (3 kanały + profil wspólny). Wymaga numpy, scipy, pandas, iminuit.
Wynik odniesienia: results/PF2_lxplus/runs/lxplus/unblind/truthnorm/res/tn_<what>[_f1].json.
Uwaga: to odtwarza OSTATNI krok (fit z gotowych histogramów i normalizacji). Wcześniejsze kroki (ntuple → h5 → score →
nogi → histogramy, N_D*, c_w) wymagają danych spoza repo — opis w analysis/README.md.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import zloz_charm_work as Z  # noqa: E402

CH = ("D0D0", "DpD0", "DpDp")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--what", default="nom427", choices=["nom427", "nom428"])
    ap.add_argument("--f-one", action="store_true", help="kontrola: f_ch ≡ 1")
    ap.add_argument("--charmtag", default="frozen", choices=["frozen", "repo"], help="frozen = analysis/PF_pairs/src (kod, który policzył wynik); repo = bieżący pakiet")
    ap.add_argument("--workdir", default=None, help="katalog roboczy (domyślnie tymczasowy)")
    ap.add_argument("--tol", type=float, default=1e-6, help="dopuszczalne |Δμ| względem zapisu")
    a = ap.parse_args()
    work = os.path.abspath(a.workdir) if a.workdir else tempfile.mkdtemp(prefix="charm_odtworzenie_")
    cw, _ = Z.assemble(work)
    src = os.path.join(cw, "PF_pairs", "src") if a.charmtag == "frozen" else REPO
    env = dict(os.environ, TN_SP=os.path.join(work, "scratch"), PYTHONPATH=src, PYTHONDONTWRITEBYTECODE="1")
    out = os.path.join(work, "out")
    cmd = [sys.executable, os.path.join(cw, Z.UNB, "truthnorm", "truthnorm_run.py"), "--what", a.what, "--out", out] + (["--f-one"] if a.f_one else [])
    print("uruchamiam:", " ".join(cmd), f"\n  PYTHONPATH={src}\n  TN_SP={env['TN_SP']}", flush=True)
    subprocess.run(cmd, env=env, check=True)
    name = f"tn_{a.what}{'_f1' if a.f_one else ''}.json"
    new = json.load(open(os.path.join(out, name)))
    ref = json.load(open(os.path.join(REPO, "results", Z.UNB, "truthnorm", "res", name)))
    rows = [("μ wspólne", new["fit"]["combined"]["mu"], ref["fit"]["combined"]["mu"]),
            ("68 % dół", new["fit"]["combined"]["mu_68"][0], ref["fit"]["combined"]["mu_68"][0]),
            ("68 % góra", new["fit"]["combined"]["mu_68"][1], ref["fit"]["combined"]["mu_68"][1]),
            ("Z", new["fit"]["combined"]["Z"], ref["fit"]["combined"]["Z"])]
    rows += [(f"μ {c}", new["fit"][c]["mu_hat"], ref["fit"][c]["mu_hat"]) for c in CH]
    print(f"\nkod: {new['charmtag']}\n{'wielkość':12s} {'odtworzone':>20s} {'zapis w repo':>20s} {'różnica':>10s}")
    worst = 0.0
    for n, x, y in rows:
        worst = max(worst, abs(x - y))
        print(f"{n:12s} {x:20.14f} {y:20.14f} {x - y:10.1e}")
    ok = worst < a.tol
    print(f"\n{'ZGODNE' if ok else 'NIEZGODNE'}: max |różnica| = {worst:.1e} (próg {a.tol:g}); katalog roboczy: {work}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
