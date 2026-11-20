"""charmtag.hand – ręczny (bez treningu) dyskryminant D0/D+ jako kontrola rzeźbienia i prior.
Docelowo: skopiować ten pakiet jako charmtag/hand/ (addytywnie) – API: charmtag.hand.score(df)."""
from .model import score, forward, export_onnx, score_onnx
from .features import prepare, PREPARED, MASS_SAFE, COLMAP
from .weights import build, n_nonzero, REGISTRY
__all__ = ["score", "forward", "export_onnx", "score_onnx", "prepare", "PREPARED", "MASS_SAFE", "COLMAP", "build", "n_nonzero", "REGISTRY"]
__version__ = "0.1.0"
