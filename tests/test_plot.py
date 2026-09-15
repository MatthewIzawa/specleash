import matplotlib
matplotlib.use("Agg")
import numpy as np

from specleash import Spectrum, compare
from specleash.plot import leash_plot


def test_leash_plot_draws_curves_links_and_the_bottleneck():
    w = np.linspace(0.4, 2.5, 80)
    a, b = Spectrum(w, 1 + 0.4 * w), Spectrum(w, 1 + 0.6 * w)
    r = compare(a, b, (0.7, 2.45), normalise_at=1.5)
    ax = leash_plot(r, labels=("A", "B"), max_links=1000)
    i, j = r.path[r.bottleneck_index]
    xs = [tuple(np.round(line.get_xdata(), 12)) for line in ax.get_lines()]
    assert (round(r.curve_a[i, 0], 12), round(r.curve_b[j, 0], 12)) in xs
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "A" in labels and "B" in labels and any(t.startswith("d = ") for t in labels)
    assert len(ax.get_lines()) == len(r.path) + 3
