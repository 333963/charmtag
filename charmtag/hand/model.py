"""Forward w numpy, API `score(df)` zgodne stylem z charmtag.scores, eksport/odczyt ONNX."""
from __future__ import annotations
import json, os
import numpy as np
import pandas as pd
from .features import prepare, PREPARED
from .weights import build

SPECIES_ALIASES = {"D0": "D0", "d0": "D0", "D+": "Dp", "Dp": "Dp", "dp": "Dp", "Dplus": "Dp", "D_plus": "Dp"}

def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))

def forward(P: np.ndarray, W1, b1, W2, b2):
    """P [N,F] -> (prob [N], logit [N], hidden [N,H])."""
    P = np.asarray(P, np.float32)
    Hd = np.tanh(P @ W1.T + b1)
    logit = (Hd @ W2.T + b2)[:, 0]
    return _sigmoid(logit), logit, Hd

def score(df: pd.DataFrame, species: str = "D0", variant: str = "A", rev: int = 1,
          colmap: dict | None = None, dedx_scale: float = 1.0, return_hidden: bool = False) -> pd.DataFrame:
    """Kolumny wyjściowe (konwencja charmtag.scores):
       score_hand_<species>_<variant>  – prawdopodobieństwo (0..1)
       logit_hand_<species>_<variant>  – logit
       opcjonalnie h_<unit> – aktywacje neuronów (diagnostyka)."""
    sp = SPECIES_ALIASES[species]
    W1, b1, W2, b2, units = build(sp, variant, rev)
    P = prepare(df, colmap=colmap, dedx_scale=dedx_scale).to_numpy()
    p, z, Hd = forward(P, W1, b1, W2, b2)
    tag = f"hand_{sp}_{variant.upper()}"
    out = pd.DataFrame({f"score_{tag}": p.astype(np.float32), f"logit_{tag}": z.astype(np.float32)}, index=df.index)
    if return_hidden:
        for i, u in enumerate(units):
            out[f"h_{u.name}"] = Hd[:, i].astype(np.float32)
    return out

# ------------------------------------------------------------------ ONNX
def export_onnx(path: str, species="D0", variant="A", rev=1, opset=17) -> str:
    """Graf: input P[N,14] (cechy PRZYGOTOWANE, patrz features.prepare) -> Gemm -> Tanh -> Gemm -> Sigmoid.
    Przygotowanie cech (log/asinh/pulle BB) pozostaje w Pythonie – zapisane w metadanych."""
    import onnx
    from onnx import helper, TensorProto, numpy_helper
    W1, b1, W2, b2, units = build(SPECIES_ALIASES[species], variant, rev)
    F = W1.shape[1]
    X = helper.make_tensor_value_info("P", TensorProto.FLOAT, ["N", F])
    Y = helper.make_tensor_value_info("score", TensorProto.FLOAT, ["N", 1])
    L = helper.make_tensor_value_info("logit", TensorProto.FLOAT, ["N", 1])
    inits = [numpy_helper.from_array(W1, "W1"), numpy_helper.from_array(b1, "b1"),
             numpy_helper.from_array(W2, "W2"), numpy_helper.from_array(b2, "b2")]
    nodes = [
        helper.make_node("Gemm", ["P", "W1", "b1"], ["z1"], transB=1),
        helper.make_node("Tanh", ["z1"], ["h"]),
        helper.make_node("Gemm", ["h", "W2", "b2"], ["logit"], transB=1),
        helper.make_node("Sigmoid", ["logit"], ["score"]),
    ]
    g = helper.make_graph(nodes, f"charmtag_hand_{species}_{variant}_rev{rev}", [X], [Y, L], inits)
    m = helper.make_model(g, producer_name="charmtag.hand", opset_imports=[helper.make_opsetid("", opset)])
    m.ir_version = 9
    meta = {"features": PREPARED, "species": species, "variant": variant, "rev": rev,
            "units": [u.name for u in units], "note": "input = charmtag.hand.features.prepare(df)"}
    for k, v in meta.items():
        e = m.metadata_props.add(); e.key = k; e.value = json.dumps(v)
    onnx.checker.check_model(m)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    onnx.save(m, path)
    return path

def score_onnx(path: str, P: np.ndarray) -> np.ndarray:
    import onnxruntime as ort
    s = ort.InferenceSession(path, providers=["CPUExecutionProvider"])
    return s.run(["score"], {"P": np.asarray(P, np.float32)})[0][:, 0]
