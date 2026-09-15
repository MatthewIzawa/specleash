"""Reflectance spectra: reading, cropping and normalisation.

Wavelength is held in micrometres. A reader is always told the unit of a file;
nothing here infers it from the numbers.
"""
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

_TO_UM = {"um": 1.0, "nm": 1e-3}
_TOL = 1e-9


class CoverageError(ValueError):
    """A spectrum does not cover a requested wavelength or range."""


@dataclass(frozen=True)
class Spectrum:
    wavelength: np.ndarray
    reflectance: np.ndarray
    sigma: Optional[np.ndarray] = None
    name: str = ""
    source: str = ""

    def __post_init__(self):
        w = np.asarray(self.wavelength, dtype=float)
        r = np.asarray(self.reflectance, dtype=float)
        if w.ndim != 1 or r.shape != w.shape:
            raise ValueError("wavelength and reflectance must be 1-D arrays of equal length")
        if len(w) < 2:
            raise ValueError("a spectrum needs at least two points")
        if not (np.all(np.isfinite(w)) and np.all(np.isfinite(r))):
            raise ValueError("wavelength and reflectance must be finite")
        if np.any(np.diff(w) <= 0):
            raise ValueError("wavelength must be strictly increasing")
        object.__setattr__(self, "wavelength", w)
        object.__setattr__(self, "reflectance", r)
        if self.sigma is not None:
            s = np.asarray(self.sigma, dtype=float)
            if s.shape != w.shape:
                raise ValueError("sigma must match wavelength in length")
            object.__setattr__(self, "sigma", s)

    def label(self):
        return self.name or "spectrum"


def read_spectrum(path, unit, *, columns=(0, 1), sigma_column=None, delimiter=None,
                  comments="#", missing=None, name=None) -> Spectrum:
    """Read a text spectrum. `unit` ("um" or "nm") is required.

    Lines that do not parse as numbers are skipped, so a column-name line after
    the comments is harmless. Rows whose reflectance equals `missing` (MITHNEOS
    files use -1) or holds a non-finite value are dropped, and rows are sorted
    by wavelength.
    """
    if unit not in _TO_UM:
        raise ValueError(f"unit must be one of {sorted(_TO_UM)}; got {unit!r}")
    cols = list(columns) + ([sigma_column] if sigma_column is not None else [])
    data = np.genfromtxt(path, delimiter=delimiter, comments=comments, usecols=cols,
                         invalid_raise=False, encoding="utf-8")
    data = np.atleast_2d(data)
    keep = np.all(np.isfinite(data), axis=1)
    if missing is not None:
        keep &= data[:, 1] != missing
    data = data[keep]
    data = data[np.argsort(data[:, 0], kind="stable")]
    return Spectrum(data[:, 0] * _TO_UM[unit], data[:, 1],
                    data[:, 2] if sigma_column is not None else None,
                    name=name if name is not None else Path(path).stem, source=str(path))


def normalise(spec: Spectrum, at: float) -> Spectrum:
    """Divide reflectance (and sigma) by the reflectance interpolated at `at` µm."""
    w = spec.wavelength
    if not (w[0] - _TOL <= at <= w[-1] + _TOL):
        raise CoverageError(f"{spec.label()} covers {w[0]:.4f}-{w[-1]:.4f} um; cannot normalise at {at} um")
    value = float(np.interp(at, w, spec.reflectance))
    if value <= 0:
        raise ValueError(f"{spec.label()} reflectance at {at} um is {value}; it must be positive to normalise")
    sigma = None if spec.sigma is None else spec.sigma / value
    return Spectrum(w, spec.reflectance / value, sigma, spec.name, spec.source)


def crop(spec: Spectrum, lo: float, hi: float) -> Spectrum:
    """The spectrum on [lo, hi], with points placed exactly at lo and hi.

    Both ends are interpolated, so two spectra cropped to one range start and
    end at the same wavelengths whatever their sampling. A range the spectrum
    does not cover is refused.
    """
    if not lo < hi:
        raise ValueError(f"range must have lo < hi; got {lo}, {hi}")
    w, r = spec.wavelength, spec.reflectance
    if w[0] > lo + _TOL or w[-1] < hi - _TOL:
        raise CoverageError(f"{spec.label()} covers {w[0]:.4f}-{w[-1]:.4f} um, not {lo}-{hi} um")
    inside = (w > lo + _TOL) & (w < hi - _TOL)
    wn = np.concatenate([[lo], w[inside], [hi]])
    rn = np.concatenate([[np.interp(lo, w, r)], r[inside], [np.interp(hi, w, r)]])
    sn = None
    if spec.sigma is not None:
        s = spec.sigma
        sn = np.concatenate([[np.interp(lo, w, s)], s[inside], [np.interp(hi, w, s)]])
    return Spectrum(wn, rn, sn, spec.name, spec.source)
