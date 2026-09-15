"""The leash plot: two spectra, the matching between them, and the link that sets the distance."""
import numpy as np


def leash_plot(comparison, ax=None, *, labels=("a", "b"), colors=("#1f5fa8", "#b03a2e"),
               link_color="0.78", bottleneck_color="#e08e0b", max_links=400, equal_aspect=True,
               y_label="Normalised reflectance", x_label=None):
    """Draw both curves as compared, every matched pair (thinned to `max_links`),
    and the longest link. Equal aspect keeps link lengths true to the distance."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    ca, cb, path = comparison.curve_a, comparison.curve_b, comparison.path
    step = max(1, int(np.ceil(len(path) / max_links)))
    for k in range(0, len(path), step):
        i, j = path[k]
        ax.plot([ca[i, 0], cb[j, 0]], [ca[i, 1], cb[j, 1]], color=link_color, lw=0.4, zorder=1)
    ax.plot(ca[:, 0], ca[:, 1], color=colors[0], lw=1.2, label=labels[0], zorder=2)
    ax.plot(cb[:, 0], cb[:, 1], color=colors[1], lw=1.2, label=labels[1], zorder=2)
    i, j = path[comparison.bottleneck_index]
    ax.plot([ca[i, 0], cb[j, 0]], [ca[i, 1], cb[j, 1]], color=bottleneck_color, lw=2.2, zorder=3,
            label=f"d = {comparison.distance:.3f}")
    if x_label is None:
        x_label = "Scaled wavelength" if comparison.settings["x_scale"] != 1.0 else "Wavelength (µm)"
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if equal_aspect:
        ax.set_aspect("equal", adjustable="datalim")
    ax.legend(frameon=False)
    return ax
