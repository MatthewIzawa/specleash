"""specleash: discrete Fréchet distance for the shapes of featureless reflectance spectra."""
from .compare import Comparison, compare, normalisation_sweep
from .frechet import FrechetResult, discrete_frechet, frechet_coupling
from .spectrum import CoverageError, Spectrum, crop, normalise, read_spectrum

__version__ = "0.1.0.dev0"
__all__ = ["Comparison", "compare", "normalisation_sweep", "FrechetResult", "discrete_frechet",
           "frechet_coupling", "CoverageError", "Spectrum", "crop", "normalise", "read_spectrum"]
