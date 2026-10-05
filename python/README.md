# geoareaweight

**Correct area weights and area-weighted means, totals and area fractions for latitude–longitude grids (NumPy and xarray).**

On a regular lat–lon grid every latitude row has the same number of cells, but a
cell's area shrinks with cos φ. A plain `mean()` over latitude and longitude
therefore over-weights the poles: on NCEP/NCAR Reanalysis 1 it gives a 1980–2020
global mean temperature of 4.8 °C instead of 14.2 °C and overstates the warming
trend by 75 %. The same mistake led to the retraction of a *Nature* paper in
2022 and was one of the reasons for another retraction in 2026.

`geoareaweight` computes exact cell areas (sphere or WGS84 ellipsoid) and uses
them for NaN-aware spatial statistics that keep time and all other dimensions.

![Cell area versus latitude](https://raw.githubusercontent.com/GISWLH/GeoAreaWeight/main/figs/fig1_concept_cell_area_vs_lat.png)

## Install

```bash
pip install geoareaweight            # NumPy only
pip install "geoareaweight[xarray]"  # with xarray support
```

Python ≥ 3.10, NumPy ≥ 1.23. xarray and dask are optional.

## Usage

```python
import xarray as xr
import geoareaweight as gaw

da = xr.open_dataset("air.mon.mean.nc")["air"]       # (time, lat, lon), any order

gaw.area_mean(da)                                    # area-weighted global mean, DataArray(time)
gaw.area_mean(da, mask=land)                         # land-only mean (boolean (lat, lon) mask)
gaw.area_mean(da, method="ellipsoid")                # exact WGS84 cell areas
gaw.area_integral(runoff_m_per_yr)                   # total volume, m3/yr
gaw.area_fraction(spei <= -1, valid=spei.notnull())  # drought area fraction
```

NumPy arrays need the coordinates and, if not `(..., lat, lon)`, the axis positions:

```python
import numpy as np
lat = np.arange(-89.5, 90, 1.0); lon = np.arange(0.5, 360, 1.0)
x = np.random.rand(12, lat.size, lon.size)            # (time, lat, lon)
gaw.area_mean(x, lat, lon)                            # shape (12,)
gaw.area_mean(np.moveaxis(x, 0, -1), lat, lon, lat_axis=0, lon_axis=1)
gaw.cell_area(lat, lon, units="km2").sum()            # 510 065 881 km2 (sphere, R = 6 371 008.8 m)
gaw.area_weights(lat)                                 # 1-D weights summing to 1
```

## API

| Function | Purpose |
|---|---|
| `area_mean(x, lat=None, lon=None, *, method="band", weights=None, mask=None, lat_bounds=None, lon_bounds=None, lat_axis=-2, lon_axis=-1, lat_name=None, lon_name=None, skipna=True)` | area-weighted mean over lat/lon |
| `region_mean(x, lat, lon, mask, **kw)` | `area_mean` with a required mask |
| `area_integral(x, lat=None, lon=None, *, ellipsoid=False, radius=None, weights=None, mask=None, ..., units="m2", skipna=True)` | `sum(x * cell_area)` |
| `area_fraction(cond, lat=None, lon=None, *, valid=None, **kw)` | area fraction where `cond` is true |
| `cell_area(lat, lon=None, lat_bounds=None, lon_bounds=None, ellipsoid=False, radius=None, units="m2")` | cell areas, `(nlat, nlon)` |
| `weights_2d(lat, lon=None, lat_bounds=None, lon_bounds=None, method="band", normalize=True)` | 2-D weight field |
| `area_weights(lat, lat_bounds=None, method="band", normalize=True)` | 1-D latitude weights |
| `band_weights`, `cos_lat_weights`, `infer_bounds`, `authalic_sin_lat`, `authalic_radius` | building blocks |

Methods: `"band"` (default, exact spherical cell area), `"cos"` (cos φ),
`"ellipsoid"` (exact WGS84), `"none"` (plain mean, for comparison).

Key behaviour:

* Coordinates are cell centres in degrees, ascending or descending, regular or
  irregular; longitudes may cross 0/360 or ±180. Bounds are inferred as mid-points
  unless `lat_bounds` / `lon_bounds` are given (`n+1` edges or `(n, 2)`).
* NaN and masked (`numpy.ma`) values are skipped and the weights renormalised per
  time step; NaN is returned where no valid cell remains.
* Masks and weights are `(nlat, nlon)`; for DataArrays they are aligned by
  coordinate labels (e.g. a descending-latitude dataset with an ascending mask).
* DataArray input returns a DataArray; dask arrays stay lazy.
* Curvilinear grids (2-D latitude) fall back to cos φ with a warning; pass the
  model's cell areas (`areacella`) as `weights=` for exact results.

Documentation, experiments, R and MATLAB/Octave versions:
<https://github.com/GISWLH/GeoAreaWeight>

## License

MIT
