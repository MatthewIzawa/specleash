# MATLAB

`specleash_frechet.m` computes the exact discrete Fréchet distance and the matching.
`specleash_compare.m` loads two spectra and compares them with the same conventions
and numbers as the Python package; `'Plot', true` draws the leash plot.

This replaces the calculation and plotting modules Gordon Gartrelle used for
Gartrelle et al. (2021): one function, no recursion limit, and a stated wavelength
range. Example, with his convention P = asteroid and Q = meteorite:

    d = specleash_compare('A17274.txt', 'A17274MM38.txt', [0.70 2.45], 'Plot', true)

Tests: `runtests('matlab/tests')` from the repository root.
Written for R2021a or later. Not yet run on MATLAB as of the first commit; see CI.
