"""GeoAreaWeight: correct area weights and area-weighted statistics on latitude-longitude grids.

Quick reference
---------------
area_mean(x, lat, lon)            area-weighted mean (intensive quantities)
region_mean(x, lat, lon, mask)    area-weighted mean inside a boolean mask
area_integral(x, lat, lon)        sum(x * cell_area)  (per-area flux -> total)
area_fraction(cond, lat, lon)     fraction of the area where cond is True
cell_area(lat, lon)               cell areas [m2 or km2], shape (nlat, nlon)
weights_2d(lat, lon)              2-D weight field, shape (nlat, nlon)
area_weights(lat)                 1-D latitude weights
band_weights / cos_lat_weights    the underlying weight formulas

Weighting methods: "band" (default, exact sphere), "cos", "ellipsoid"
(exact WGS84), "none" (plain mean, for comparison). Angles are in degrees.
"""
from .means import area_fraction, area_integral, area_mean, region_mean
from .weights import (
    EARTH_RADIUS,
    METHODS,
    WGS84_A,
    WGS84_E2,
    WGS84_F,
    area_weights,
    authalic_radius,
    authalic_sin_lat,
    band_weights,
    cell_area,
    cos_lat_weights,
    infer_bounds,
    weights_2d,
)

__version__ = "0.1.0"

__all__ = [
    "area_mean",
    "region_mean",
    "area_integral",
    "area_fraction",
    "area_weights",
    "cos_lat_weights",
    "band_weights",
    "cell_area",
    "weights_2d",
    "infer_bounds",
    "authalic_sin_lat",
    "authalic_radius",
    "EARTH_RADIUS",
    "WGS84_A",
    "WGS84_F",
    "WGS84_E2",
    "METHODS",
    "__version__",
]
