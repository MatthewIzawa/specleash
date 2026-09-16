# specleash

Compare the shapes of largely featureless reflectance spectra, such as those of
D-type asteroids and carbonaceous chondrites, with the discrete Fréchet distance.

Think of a person and a dog walking along two curves, each always moving
forward. The Fréchet distance is the shortest leash that lets them both reach
the end. `specleash` computes it exactly, returns the matching between the two
spectra, and reports the pair of wavelengths that sets the distance.

```python
from specleash import read_spectrum, compare
from specleash.plot import leash_plot

asteroid = read_spectrum("asteroid.txt", unit="um")
meteorite = read_spectrum("meteorite.txt", unit="nm")
result = compare(asteroid, meteorite, wl_range=(0.70, 2.45), normalise_at=1.5)
print(result.distance, result.bottleneck)
leash_plot(result, labels=("asteroid", "meteorite"))
```

Every choice that changes the number is an argument: the wavelength range (no
default), the normalisation wavelength, the scaling of the wavelength axis, and
the sampling. A comparison over a range that either spectrum does not cover is
refused.

The method follows Gartrelle et al. (2021, Icarus 361, 114349). A MATLAB version
is in `matlab/`. `examples/` holds the polynomial fit to the spectrum of
(17274) 2000 LC16 that Gartrelle et al. (2021) compared with Tagish Lake spectra.

Tests: `python -m pytest`
