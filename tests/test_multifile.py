"""Etap P: próbka = lista plików / glob — wynik identyczny jak dla jednego pliku."""
import json
import os

import h5py
import numpy as np
import pandas as pd
import pytest

from charmtag import config as C
from charmtag import samples as SA
from charmtag.io import n_rows, open_events, read_rows
from charmtag.synth import make


@pytest.fixture(scope="module")
def two(tmp_path_factory):
    d = tmp_path_factory.mktemp("mf")
    make(str(d / "h5"), n=3000, seed=5)
    base = C.load(None, {"paths": {"h5_dir": str(d / "h5"), "out_dir": str(d / "out")}})
    src = C.sample_path(base, "data_none")
    with h5py.File(src) as f:
        a = f["events"][:]
    k = len(a) // 3
    os.makedirs(d / "parts", exist_ok=True)
    names = []
    for j, (lo, hi) in enumerate([(0, k), (k, len(a))], start=1):
        q = d / "parts" / f"data_428770_{j:06d}_none.h5"
        with h5py.File(q, "w") as f:
            f.create_dataset("events", data=a[lo:hi])
            f.attrs["n_events_read"] = 100 * j
        names.append(str(q))
    with open(d / "parts" / "data_428770_000003_none.h5.part", "w") as f:
        f.write("x")  # niedokończony — glob musi pominąć
    multi = C.load(None, {"paths": {"h5_dir": str(d / "h5"), "out_dir": str(d / "out2"),
                                    "samples": {"data_none": ["../parts/data_428770_*_none.h5"]},
                                    "sample_info": {"data_none": {"jobs_total": 47}}}})
    return base, multi, src, names, a


def test_paths_glob(two):
    base, multi, src, names, a = two
    assert C.sample_paths(multi, "data_none") == names
    assert C.is_multi(multi, "data_none") and not C.is_multi(base, "data_none")
    assert C.source(base, "data_none") == src
    assert C.sample_path(multi, "data_none") == names[0]


def test_multievents_slices(two):
    *_, names, a = two
    with open_events(names) as ds:
        assert ds.shape[0] == len(a) == n_rows(names)
        for lo, hi in [(0, 10), (len(a) // 3 - 5, len(a) // 3 + 5), (0, len(a)), (len(a) - 3, len(a) + 10)]:
            x = ds.fields(["Dmeson_mass", "event_number"])[lo:hi]
            assert np.array_equal(x["Dmeson_mass"], a["Dmeson_mass"][lo:hi])
        i, loc = ds.locate([0, len(a) // 3, len(a) - 1])
        assert list(i) == [0, 1, 1] and loc[1] == 0


def test_read_rows_list(two):
    *_, names, a = two
    df = read_rows(names, ["Dmeson_mass"], chunk_rows=777, max_rows=len(a) // 3 + 7)
    assert np.array_equal(df["Dmeson_mass"].to_numpy(), a["Dmeson_mass"][: len(a) // 3 + 7])


def test_iter_species_equal(two):
    base, multi, *_ = two
    for cfg in (base, multi):
        cfg["io"]["chunk_rows"] = 500
    A = pd.concat([L for L, _ in SA.iter_species(base, "data_none", "D0")], ignore_index=True)
    B = pd.concat([L for L, _ in SA.iter_species(multi, "data_none", "D0")], ignore_index=True)
    pd.testing.assert_frame_equal(A, B)
    assert SA.count_species(base, "data_none", "D0") == SA.count_species(multi, "data_none", "D0")
    La = SA.load_species(base, "data_none", "D0", cache=True)
    Lb = SA.load_species(multi, "data_none", "D0", cache=True)
    assert len(La) == len(Lb)


def test_meta(two):
    base, multi, *_ = two
    m = C.sample_meta(multi, "data_none")
    assert m["n_files"] == 2 and m["jobs"] == [1, 2] and m["n_events_read"] == 300
    assert abs(m["lumi_fraction_jobs"] - 2 / 47) < 1e-12
    assert json.loads(C.sample_meta_json(base, "data_none"))["n_files"] == 1


def test_apply_multi_rows(two):
    """apply na liście plików: row = wiersz globalny; attrs z metadanymi."""
    from charmtag import pipeline as P
    base, multi, src, names, a = two
    for cfg in (base, multi):
        cfg["training"]["xgb"]["n_estimators"] = 30
    P.cmd_train(base, "D0", "H")
    multi["paths"]["out_dir"] = base["paths"]["out_dir"]  # ten sam model
    done = P.cmd_apply(multi, "D0", "H", samples=["data_none"])
    assert done["data_none"] > 0
    sc = P.load_scores(multi, "D0", "H", "data_none")
    assert np.allclose(a["Dmeson_mass"][sc["row"].to_numpy()], sc["m"].to_numpy())
    p = os.path.join(P.sdir(multi, "D0", "H"), "data_none.h5")
    with h5py.File(p) as f:
        assert int(f.attrs["n_files"]) == 2
        assert json.loads(f.attrs["sample_meta"])["jobs"] == [1, 2]
    assert os.path.exists(os.path.join(P.sdir(multi, "D0", "H"), "samples_meta.json"))
