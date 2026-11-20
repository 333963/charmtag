"""09.11: bootstrap σ(N) po zdarzeniach w cmd_fit (`pipeline._boot_err`) — wiersze bliźniacze zawyżają σ względem NLL."""
import numpy as np
import pandas as pd

from charmtag import massfit as MF
from charmtag import pipeline as P

LO, HI, M0 = 1.75, 2.0, 1.86484
FC = {"nbins": 100, "signal_model": "gauss", "background_model": "exp", "sigma_mode": "constrained",
      "sigma_constraint_rel": 0.10, "reflection_mode": "off", "bootstrap_reps": 60}


def _sample(rng, n_ev, dup):
    s = rng.normal(M0, 0.02, int(0.15 * n_ev))
    b = LO + rng.exponential(0.3, 3 * n_ev)
    m = np.concatenate([s, b[b < HI][: n_ev - len(s)]])
    ev = np.arange(len(m))
    if dup:   # każde zdarzenie dwa razy (te same masy) — skrajny przypadek „bliźniaków”
        m, ev = np.repeat(m, 2), np.repeat(ev, 2)
    return pd.DataFrame({"m": m, "run": 1, "event": ev})


def _ratio(dup, seed):
    rng = np.random.default_rng(seed)
    d = _sample(rng, 4000, dup)
    shape = MF.fit_mc_shape(rng.normal(M0, 0.02, 20000), np.ones(20000), LO, HI, M0, "gauss", 100)
    f = MF.fit_data(d["m"].to_numpy(), LO, HI, M0, shape, FC)
    b = P._boot_err(d, LO, HI, M0, shape, FC, None, 0.0, seed)
    assert b and "N_sig_err_boot" in b, b
    return b["N_sig_err_boot"] / f["N_sig_err"]


def test_boot_independent_rows_close_to_one():
    assert 0.7 < _ratio(False, 1) < 1.35


def test_boot_duplicated_rows_inflate_sigma():
    assert 1.2 < _ratio(True, 2) < 1.8   # oczekiwane √2


def test_boot_disabled_or_no_event_columns():
    d = pd.DataFrame({"m": np.linspace(LO, HI, 100)})
    assert P._boot_err(d, LO, HI, M0, {}, FC, None, 0.0, 0) is None
    assert P._boot_err(d.assign(run=1, event=0), LO, HI, M0, {}, {**FC, "bootstrap_reps": 0}, None, 0.0, 0) is None
