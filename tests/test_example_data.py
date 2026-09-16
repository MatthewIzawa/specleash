"""The example data load as documented in examples/README.md."""
import hashlib
from pathlib import Path

import numpy as np
import pytest

from specleash import read_spectrum

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "data"


def test_the_gartrelle_fit_is_stored_byte_for_byte():
    path = EXAMPLES / "A17274_Gartrelle2021_fit.txt"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == \
        "9dca499f3039c282549a8cb66d169383c978853ffe1eab6076d92f8df5c2cee8"


def test_the_gartrelle_fit_loads_over_0_70_to_2_45_um():
    s = read_spectrum(EXAMPLES / "A17274_Gartrelle2021_fit.txt", unit="um")
    assert len(s.wavelength) == 388
    assert s.wavelength[0] == pytest.approx(0.70) and s.wavelength[-1] == pytest.approx(2.45)
    assert np.interp(1.47, s.wavelength, s.reflectance) == pytest.approx(1.0, abs=0.001)
