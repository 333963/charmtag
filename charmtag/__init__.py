"""charmtag — double charm (pary D–D) w danych ATLAS minimum-bias 2022: klasyfikator pojedynczego D (XGBoost, k-fold),
fit masy, pary D–D z fitem 2D, systematyki.

Projekt B (tło podobne do sygnału: bkggen, HYP, SALAD, …) jest od 2.1.0 osobnym pakietem `charmbkg` (repo charmbkg);
stare importy `charmtag.bkggen`, `charmtag.hyp` itd. są aliasami do `charmbkg` (zob. _bkg_alias.py).
"""
__version__ = "2.1.1"

from . import _bkg_alias as _alias  # noqa: E402

_alias.install()
