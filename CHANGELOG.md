# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). Python, R and MATLAB/Octave share one version number.

## [0.1.0] - unreleased

First public release: Python (`geoareaweight` on PyPI), R and MATLAB/Octave implementations.

### Added
- Weights: `cos_lat_weights`, exact spherical `band_weights`, WGS84 (authalic latitude),
  `area_weights`, `weights_2d`, `cell_area`, `infer_bounds`.
- Statistics: `area_mean`, `region_mean`, `area_integral`, `area_fraction`, keeping all
  non-spatial dimensions, with NaN-aware renormalisation.
- Python: xarray support (CF-based coordinate detection, label alignment of masks and
  weights, lazy dask evaluation), `numpy.ma` support, type hints (`py.typed`), doctests.
- `area_integral` accepts `weights=`, `radius=` and `skipna=` (all languages).
- Validation of longitude wrap-around (e.g. 350, 355, 0, 5), monotonic coordinates,
  latitude range (catches swapped lat/lon), method names, mask/weight shapes and
  negative weights, with error messages that say how to fix the input.
- MATLAB/Octave: case-insensitive name-value options with unknown-option errors,
  `'Valid'` option for `gaw_area_fraction`, `Contents.m`; internal helpers in `private/`.
- English documentation, figures, `llms.txt`, CI workflow, citation metadata.

### Fixed (relative to the pre-release code)
- xarray masks/weights were applied by position, not by label: wrong results for data
  with descending latitude, and an error for data stored as (lon, lat).
- Dimensions named `x`/`y` (projected grids) were taken as longitude/latitude.
- 1-D latitude/longitude coordinates that are not dimension coordinates raised an error.
- With 2-D latitudes, unknown methods and the alias `"sphere_band"` silently returned
  equal weights.
- Masked arrays (e.g. from netCDF4) used their fill values as data.
- Numeric masks with NaN (e.g. from `regionmask`) included the NaN cells.
- Longitudes crossing the 0/360 meridian produced a huge cell width.
- `area_integral` returned 0 instead of NaN when no valid cell remained, and returned a
  NumPy array for DataArray input.
- `area_mean` copied the weights to the full data shape (memory); now a matrix product.
