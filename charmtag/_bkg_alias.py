"""Moduły projektu B przeniesione do osobnego pakietu `charmbkg` (repo charmbkg) — 15.11.2026, charmtag 2.1.0.

Stare ścieżki importu (`charmtag.bkggen…`, `charmtag.hyp`, …) działają dalej: ten finder podstawia pod nie moduły
`charmbkg`, więc configi i snapshoty na lxplus nie wymagają zmian. To ALIAS (ten sam obiekt modułu), nie kopia:
`charmtag.hyp is charmbkg.hyp`.
"""
import importlib
import importlib.abc
import importlib.util
import sys

MOVED = ("bkggen", "hyp", "hyp_extract", "hyp_toy", "dstar_tag", "refldata", "bkgsrc", "s8salad", "val_epsb", "sf_rel")
_PREFIX = "charmtag."


def _target(fullname):
    if not fullname.startswith(_PREFIX):
        return None
    rest = fullname[len(_PREFIX):]
    if rest.split(".")[0] in MOVED:
        return "charmbkg." + rest
    return None


class _Loader(importlib.abc.Loader):
    def __init__(self, target):
        self.target = target

    def create_module(self, spec):
        mod = importlib.import_module(self.target)
        self._orig_spec = mod.__spec__   # importlib nadpisze __spec__ specem aliasu — przywracamy w exec_module
        return mod

    def exec_module(self, module):
        module.__spec__ = self._orig_spec   # moduł zostaje `charmbkg.*` (importy względne, pickle, __package__)


class _Finder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        t = _target(fullname)
        if t is None:
            return None
        if importlib.util.find_spec(t) is None:
            raise ImportError(f"{fullname} jest teraz w pakiecie charmbkg (repo charmbkg): pip install charmbkg")
        mod = importlib.import_module(t)
        spec = importlib.util.spec_from_loader(fullname, _Loader(t), is_package=hasattr(mod, "__path__"))
        if hasattr(mod, "__path__"):
            spec.submodule_search_locations = list(mod.__path__)
        return spec


def install():
    if not any(isinstance(f, _Finder) for f in sys.meta_path):
        sys.meta_path.insert(0, _Finder())
