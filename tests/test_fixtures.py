"""expected.json is what MATLAB is checked against, so it must be what Python computes now."""
import json
from pathlib import Path

import numpy as np
import pytest

from specleash import Spectrum, compare, discrete_frechet

FIX = Path(__file__).resolve().parent / "fixtures"


def load(name):
    w, r = np.loadtxt(FIX / name, delimiter=",", unpack=True)
    return w, r


def test_expected_values_are_current():
    e = json.loads((FIX / "expected.json").read_text())
    for p in e["pairs"]:
        P, Q = np.column_stack(load(p["p"])), np.column_stack(load(p["q"]))
        assert discrete_frechet(P, Q) == pytest.approx(p["distance"], abs=1e-12), p
    for c in e["comparisons"]:
        d = compare(Spectrum(*load(c["a"])), Spectrum(*load(c["b"])), tuple(c["wl_range"]),
                    normalise_at=c["normalise_at"], x_scale=c["x_scale"]).distance
        assert d == pytest.approx(c["distance"], abs=1e-12), c
