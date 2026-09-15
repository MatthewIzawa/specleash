"""Write the fixture spectra and the values both implementations must reproduce.

Run from specleash/ with src on the path:  python tests/fixtures/make_expected.py
The CSVs are headerless (wavelength_um, reflectance) so MATLAB's readmatrix reads them directly.
"""
import json
from pathlib import Path

import numpy as np

from specleash import Spectrum, compare, discrete_frechet

HERE = Path(__file__).resolve().parent


def spectra():
    rng = np.random.default_rng(2026)
    w1 = np.round(np.linspace(0.40, 2.50, 141), 6)
    w2 = np.round(np.linspace(0.45, 2.48, 88), 6)
    s1 = 0.05 * (1 + 0.45 * (w1 - 0.55) - 0.05 * (w1 - 1.5) ** 2) + rng.normal(0, 3e-4, w1.size)
    s2 = 0.04 * (1 + 0.60 * (w2 - 0.55)) + rng.normal(0, 3e-4, w2.size)
    s3 = 0.03 * (1 + 0.20 * (w2 - 0.55) + 0.02 * np.sin(6 * w2))
    return {"fixture_a.csv": (w1, s1), "fixture_b.csv": (w2, s2), "fixture_c.csv": (w2, s3)}


def main():
    for name, (w, r) in spectra().items():
        np.savetxt(HERE / name, np.column_stack([w, r]), delimiter=",", fmt="%.10g")
    # Expected values come from the files as written, which is what MATLAB and
    # test_fixtures.py read, so rounding to 10 digits cannot open a gap.
    data = {name: tuple(np.loadtxt(HERE / name, delimiter=",", unpack=True)) for name in spectra()}
    pairs = []
    for p, q in (("fixture_a.csv", "fixture_b.csv"), ("fixture_b.csv", "fixture_c.csv")):
        P, Q = (np.column_stack(data[p]), np.column_stack(data[q]))
        pairs.append(dict(p=p, q=q, distance=discrete_frechet(P, Q)))
    comparisons = []
    for a, b, rng_, norm, xs in (("fixture_a.csv", "fixture_b.csv", [0.7, 2.45], 1.5, 1.0),
                                 ("fixture_a.csv", "fixture_c.csv", [0.5, 2.4], None, "unit"),
                                 ("fixture_b.csv", "fixture_c.csv", [0.6, 2.3], 0.55, 2.0)):
        A, B = Spectrum(*data[a]), Spectrum(*data[b])
        d = compare(A, B, tuple(rng_), normalise_at=norm, x_scale=xs).distance
        comparisons.append(dict(a=a, b=b, wl_range=rng_, normalise_at=norm, x_scale=xs, distance=d))
    (HERE / "expected.json").write_text(json.dumps(dict(pairs=pairs, comparisons=comparisons), indent=1))
    print("wrote", HERE / "expected.json")


if __name__ == "__main__":
    main()
