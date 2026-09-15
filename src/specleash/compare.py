"""Compare two spectra by the discrete Fréchet distance between their curves.

Each spectrum is normalised (optionally, on its whole extent), cropped to the
comparison range with points placed exactly at both ends, optionally resampled
onto a common grid, and turned into the polyline (x_scale × wavelength,
reflectance). Every choice that changes the number is an argument and is
recorded on the result.
"""
from dataclasses import dataclass
from typing import Tuple

import numpy as np

from .frechet import frechet_coupling
from .spectrum import Spectrum, crop, normalise


@dataclass(frozen=True)
class Comparison:
    distance: float
    path: np.ndarray
    link_lengths: np.ndarray
    bottleneck_index: int
    bottleneck: Tuple[float, float]   # wavelengths (µm) of the pair that sets the distance
    curve_a: np.ndarray               # (n, 2): scaled wavelength, reflectance
    curve_b: np.ndarray
    wavelength_a: np.ndarray          # µm
    wavelength_b: np.ndarray
    settings: dict


def _grid(lo, hi, step):
    n = int(round((hi - lo) / step)) + 1
    w = lo + step * np.arange(n)
    if abs(w[-1] - hi) > 1e-9:
        raise ValueError(f"range {lo}-{hi} um is not a whole number of {step} um grid steps")
    w[-1] = hi
    return w


def _prepare(spec, lo, hi, normalise_at, grid):
    s = spec if normalise_at is None else normalise(spec, normalise_at)
    s = crop(s, lo, hi)
    if grid is not None:
        w = _grid(lo, hi, grid)
        s = Spectrum(w, np.interp(w, s.wavelength, s.reflectance), None, s.name, s.source)
    return s


def _x(w, lo, hi, x_scale):
    if x_scale == "unit":
        return (w - lo) / (hi - lo)
    if isinstance(x_scale, str) or not x_scale > 0:
        raise ValueError(f"x_scale must be a positive number or 'unit'; got {x_scale!r}")
    return float(x_scale) * w


def compare(a: Spectrum, b: Spectrum, wl_range, *, normalise_at=None, x_scale=1.0, grid=None) -> Comparison:
    """Discrete Fréchet comparison of `a` and `b` over `wl_range` (µm; required)."""
    lo, hi = (float(v) for v in wl_range)
    if not lo < hi:
        raise ValueError(f"wl_range must have lo < hi; got {wl_range}")
    sa = _prepare(a, lo, hi, normalise_at, grid)
    sb = _prepare(b, lo, hi, normalise_at, grid)
    ca = np.column_stack([_x(sa.wavelength, lo, hi, x_scale), sa.reflectance])
    cb = np.column_stack([_x(sb.wavelength, lo, hi, x_scale), sb.reflectance])
    r = frechet_coupling(ca, cb)
    i, j = r.path[r.bottleneck]
    return Comparison(
        distance=r.distance, path=r.path, link_lengths=r.link_lengths, bottleneck_index=r.bottleneck,
        bottleneck=(float(sa.wavelength[i]), float(sb.wavelength[j])),
        curve_a=ca, curve_b=cb, wavelength_a=sa.wavelength, wavelength_b=sb.wavelength,
        settings=dict(wl_range=(lo, hi), normalise_at=normalise_at, x_scale=x_scale, grid=grid))


def normalisation_sweep(a, b, wl_range, anchors, **kwargs) -> np.ndarray:
    """The distance at each normalisation wavelength in `anchors`."""
    return np.array([compare(a, b, wl_range, normalise_at=k, **kwargs).distance for k in anchors])
