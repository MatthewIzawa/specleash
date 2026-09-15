import os
from pathlib import Path

import numpy as np
import pytest

from specleash import compare, normalisation_sweep, read_spectrum
from specleash.frechet import discrete_frechet
from specleash.spectrum import CoverageError, Spectrum


def red(w, slope=0.4, level=1.0):
    return Spectrum(w, level * (1 + slope * (w - 0.55)))


def test_a_spectrum_against_itself_is_zero():
    s = red(np.linspace(0.4, 2.5, 200))
    assert compare(s, s, (0.7, 2.45), normalise_at=1.5).distance == 0.0


def test_the_range_is_required_and_refused_when_not_covered():
    a, b = red(np.linspace(0.4, 2.5, 200)), red(np.linspace(0.8, 2.5, 200))
    with pytest.raises(TypeError):
        compare(a, b)
    with pytest.raises(CoverageError):
        compare(a, b, (0.7, 2.45))


def test_normalisation_happens_before_cropping():
    """A normalising wavelength outside the range is allowed when both spectra cover it."""
    a, b = red(np.linspace(0.4, 2.5, 200)), red(np.linspace(0.4, 2.5, 150), slope=0.5)
    r = compare(a, b, (0.7, 2.45), normalise_at=0.55)
    assert np.interp(0.7, r.wavelength_a, r.curve_a[:, 1]) == pytest.approx(1 + 0.4 * 0.15, rel=1e-6)


def test_unit_scaling_maps_the_range_onto_zero_to_one():
    a, b = red(np.linspace(0.4, 2.5, 200)), red(np.linspace(0.4, 2.5, 150), slope=0.5)
    r = compare(a, b, (0.7, 2.45), normalise_at=1.5, x_scale="unit")
    assert r.curve_a[0, 0] == 0.0 and r.curve_a[-1, 0] == 1.0
    assert r.curve_b[0, 0] == 0.0 and r.curve_b[-1, 0] == 1.0


def test_the_distance_is_that_of_the_curves_it_reports():
    a, b = red(np.linspace(0.4, 2.5, 120)), red(np.linspace(0.4, 2.5, 90), slope=0.6)
    r = compare(a, b, (0.7, 2.45), normalise_at=1.5, x_scale=2.0)
    assert r.distance == pytest.approx(discrete_frechet(r.curve_a, r.curve_b), abs=1e-12)
    assert np.allclose(r.curve_a[:, 0], 2.0 * r.wavelength_a)


def test_the_bottleneck_names_the_wavelengths_of_the_longest_link():
    a, b = red(np.linspace(0.4, 2.5, 120)), red(np.linspace(0.4, 2.5, 90), slope=0.9)
    r = compare(a, b, (0.7, 2.45), normalise_at=0.7)
    i, j = r.path[r.bottleneck_index]
    assert r.bottleneck == (r.wavelength_a[i], r.wavelength_b[j])
    assert r.link_lengths[r.bottleneck_index] == pytest.approx(r.distance)


def test_a_common_grid_gives_equal_sampling_and_must_fit_the_range():
    a, b = red(np.linspace(0.4, 2.5, 120)), red(np.linspace(0.4, 2.5, 90), slope=0.6)
    r = compare(a, b, (0.70, 2.45), normalise_at=1.5, grid=0.05)
    assert len(r.curve_a) == len(r.curve_b) == 36
    with pytest.raises(ValueError, match="grid"):
        compare(a, b, (0.70, 2.43), grid=0.05)


def test_settings_are_recorded():
    s = red(np.linspace(0.4, 2.5, 120))
    r = compare(s, s, (0.7, 2.45), normalise_at=1.5, x_scale="unit", grid=0.05)
    assert r.settings == dict(wl_range=(0.7, 2.45), normalise_at=1.5, x_scale="unit", grid=0.05)


def test_normalisation_sweep_matches_individual_comparisons():
    a, b = red(np.linspace(0.4, 2.5, 120)), red(np.linspace(0.4, 2.5, 90), slope=0.6)
    anchors = [0.8, 1.2, 2.0]
    sweep = normalisation_sweep(a, b, (0.7, 2.45), anchors)
    assert np.allclose(sweep, [compare(a, b, (0.7, 2.45), normalise_at=k).distance for k in anchors])


GARTRELLE = Path(os.environ.get("SPECLEASH_GARTRELLE_DIR",
                                r"C:\Users\Matthew Izawa\Documents\Gordon_Gartrelle\Frechet"))


@pytest.mark.skipif(not (GARTRELLE / "A17274.txt").exists(), reason="Gartrelle's test files not present")
def test_gartrelle_pair_as_supplied():
    """Gartrelle's two test files, compared as his MATLAB script compares them:
    raw columns, no renormalisation, native sampling, wavelength in µm."""
    a = read_spectrum(GARTRELLE / "A17274.txt", unit="um")
    m = read_spectrum(GARTRELLE / "A17274MM38.txt", unit="um")
    assert compare(a, m, (0.70, 2.45)).distance == pytest.approx(0.142488, abs=5e-7)
