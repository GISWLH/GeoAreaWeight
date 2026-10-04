"""GeoAreaWeight: area-weighted means on latitude-longitude grids."""
from .weights import (EARTH_RADIUS, WGS84_A, WGS84_F, area_weights, authalic_radius,
                      authalic_sin_lat, band_weights, cell_area, cos_lat_weights,
                      infer_bounds, weights_2d)
from .means import area_fraction, area_integral, area_mean, region_mean

__version__ = "0.1.0"
__all__ = [
    "area_mean", "region_mean", "area_integral", "area_fraction",
    "area_weights", "cos_lat_weights", "band_weights", "cell_area", "weights_2d",
    "infer_bounds", "authalic_sin_lat", "authalic_radius",
    "EARTH_RADIUS", "WGS84_A", "WGS84_F",
]
