#!/usr/bin/env python3
"""UBPNET — nogi z DANYCH 428770 w CAŁYM zakresie mas fitu (okno sygnału otwarte; zatwierdzone 18.11, CHARM_UNBLIND=1).

Dokładnie pn_extract.py z łańcucha PAIRNET v3 (te same: h5 v2.0, cechy H, modele H, preselekcja WP90, arbitraż D0 --arb, uid/ftag/tid),
jedyna zmiana: maska zaślepienia `sb_mask` (tylko pasma boczne) zastąpiona maską zakresu fitu (25 binów gatunku). Dzięki temu nogi z pasm
są identyczne jak w częściach SB×SB, na których trenowano pnetd (te same ek → te same foldy OOF). meta: sb_only = False, full_window = True.
Użycie: upn_extract.py --models DIR --out OUT.npz FILE.h5 [...]   (wymusza --tag data --arb)"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PN = os.environ.get("UPN_PN") or os.path.join(os.path.dirname(HERE), "pairnet")
sys.path.insert(0, PN)


def range_mask(m, isp):
    import pn_extract as E
    bb = np.digitize(np.asarray(m, float), E.SB_EDGES[isp]) - 1
    return (bb >= 0) & (bb < 25)


def main(argv=None):
    if os.environ.get("CHARM_UNBLIND") != "1":
        raise SystemExit("upn_extract: wymaga CHARM_UNBLIND=1 (okno sygnału danych)")
    import pn_extract as E
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--single" in argv or "--tag" in argv:
        raise SystemExit("upn_extract: --tag/--single ustawiane wewnętrznie")
    E.sb_mask = range_mask                                         # jedyna różnica względem pn_extract --tag data
    out = argv[argv.index("--out") + 1]
    E.main(["--tag", "data", "--arb"] + argv)
    z = dict(np.load(out, allow_pickle=False))
    meta = json.loads(str(z["meta"]))
    meta.update({"sb_only": False, "full_window": True, "extractor": "ub_pnet/upn_extract.py (pn_extract + maska zakresu)"})
    z["meta"] = np.array(json.dumps(meta))
    np.savez_compressed(out, **z)
    print(f"[upn_extract] okno otwarte, nóg {len(z['m'])} → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
