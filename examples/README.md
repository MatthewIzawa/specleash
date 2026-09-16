# Example data

## `data/A17274_Gartrelle2021_fit.txt`

The polynomial fit to the visible and near-infrared spectrum of the Amor asteroid
(17274) 2000 LC16 that G. M. Gartrelle compared with Tagish Lake spectra in
Gartrelle et al. (2021, Icarus 361, 114349).

- **Columns:** wavelength (µm) and normalised reflectance, tab-separated, 388 rows,
  0.700–2.450 µm in 0.005 µm steps. There are no rows between 0.755 and 0.770,
  1.700 and 1.710, or 2.435 and 2.445 µm.
- **What it is:** a fourth-order polynomial fitted over 0.70–2.45 µm to the spectrum
  of 17274 in the Bus-DeMeo taxonomy (DeMeo et al., 2009, Icarus 202, 160–180; visible
  segment from Binzel et al., 2004, Icarus 170, 259–294), after that spectrum was
  normalised to 1 at 1.5 µm. Refitting that spectrum the same way reproduces this
  file with a root-mean-square difference of 0.0006.
- **Where the input spectrum is:** the spectrum of 17274, and the Tagish Lake spectra
  it was compared with, are in the supplementary data of Izawa et al. (2015,
  Icarus 254, 324–332).
- **Integrity:** the file is stored byte for byte as supplied (`.gitattributes`
  turns off line-ending conversion). SHA-256:
  `9dca499f3039c282549a8cb66d169383c978853ffe1eab6076d92f8df5c2cee8`.

```python
from specleash import read_spectrum, compare

fit = read_spectrum("examples/data/A17274_Gartrelle2021_fit.txt", unit="um")
mm38 = read_spectrum("mm38.csv", unit="nm", delimiter=",")   # MM-38 from Izawa et al. (2015)
result = compare(fit, mm38, wl_range=(0.70, 2.45), normalise_at=1.5)
print(result.distance, result.bottleneck)   # 0.152, set at 2.45 µm on both spectra
```
