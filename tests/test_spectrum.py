import numpy as np
import pytest

from specleash.spectrum import CoverageError, Spectrum, crop, normalise, read_spectrum


def test_the_unit_is_required(tmp_path):
    f = tmp_path / "s.txt"
    f.write_text("0.5 1.0\n1.0 1.1\n")
    with pytest.raises(TypeError):
        read_spectrum(f)
    with pytest.raises(ValueError, match="unit"):
        read_spectrum(f, unit="angstrom")


def test_nanometres_are_converted_to_micrometres(tmp_path):
    f = tmp_path / "s.csv"
    f.write_text("# comment\nwavelength_nm,reflectance\n500,0.04\n600,0.05\n700,0.06\n")
    s = read_spectrum(f, unit="nm", delimiter=",")
    assert np.allclose(s.wavelength, [0.5, 0.6, 0.7])
    assert np.allclose(s.reflectance, [0.04, 0.05, 0.06])


def test_missing_values_are_dropped_and_rows_sorted(tmp_path):
    f = tmp_path / "s.txt"
    f.write_text("0.6 1.1 0.01 1\n0.5 1.0 0.01 1\n0.7 -1.0 -1.0 0\n0.8 1.3 0.02 1\n")
    s = read_spectrum(f, unit="um", sigma_column=2, missing=-1.0)
    assert np.allclose(s.wavelength, [0.5, 0.6, 0.8])
    assert np.allclose(s.sigma, [0.01, 0.01, 0.02])


def test_wavelength_must_increase():
    with pytest.raises(ValueError, match="increasing"):
        Spectrum(np.array([0.5, 0.5, 0.6]), np.array([1.0, 1.0, 1.1]))


def test_crop_refuses_a_range_the_spectrum_does_not_cover():
    s = Spectrum(np.linspace(0.8, 2.4, 50), np.ones(50), name="short")
    with pytest.raises(CoverageError, match="short"):
        crop(s, 0.7, 2.45)


def test_crop_places_points_exactly_at_both_ends():
    w = np.linspace(0.40, 2.50, 43)
    s = Spectrum(w, 1 + 0.2 * w)
    c = crop(s, 0.7, 2.45)
    assert c.wavelength[0] == 0.7 and c.wavelength[-1] == 2.45
    assert c.reflectance[0] == pytest.approx(1 + 0.2 * 0.7)
    assert c.reflectance[-1] == pytest.approx(1 + 0.2 * 2.45)
    assert np.all(np.diff(c.wavelength) > 0)


def test_normalise_sets_unity_at_the_wavelength():
    w = np.linspace(0.5, 2.5, 81)
    s = normalise(Spectrum(w, 0.05 + 0.01 * w, 0.001 * np.ones(81)), 1.5)
    assert np.interp(1.5, s.wavelength, s.reflectance) == pytest.approx(1.0)
    assert s.sigma[0] == pytest.approx(0.001 / (0.05 + 0.015))


def test_normalise_refuses_outside_the_spectrum_or_nonpositive_values():
    w = np.linspace(0.5, 2.5, 81)
    with pytest.raises(CoverageError):
        normalise(Spectrum(w, np.ones(81)), 2.6)
    with pytest.raises(ValueError, match="positive"):
        normalise(Spectrum(w, np.zeros(81)), 1.5)
