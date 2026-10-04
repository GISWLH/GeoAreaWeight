"""Area-weighted means, regional means, integrals and area fractions."""
from __future__ import annotations

import numpy as np

from .weights import cell_area, weights_2d

__all__ = ["area_mean", "area_integral", "area_fraction", "region_mean"]

_LAT_NAMES = ("lat", "latitude", "Latitude", "y", "nav_lat", "rlat")
_LON_NAMES = ("lon", "longitude", "Longitude", "x", "nav_lon", "rlon")


def _find(da, names, given):
    if given is not None:
        return given
    for n in names:
        if n in da.dims or n in da.coords:
            return n
    raise KeyError(f"could not find a coordinate among {names}; pass lat_name/lon_name")


def _core(x, w2d, lat_axis, lon_axis, mask, skipna):
    """Weighted mean of ``x`` over (lat_axis, lon_axis) with renormalised weights."""
    x = np.asarray(x, dtype=float)
    nd = x.ndim
    la, lo = lat_axis % nd, lon_axis % nd
    if la == lo:
        raise ValueError("lat_axis and lon_axis must differ")
    xm = np.moveaxis(x, (la, lo), (-2, -1))
    w = np.broadcast_to(np.asarray(w2d, dtype=float), xm.shape[-2:])
    if w.shape != xm.shape[-2:]:
        raise ValueError(f"weights {w.shape} do not match data grid {xm.shape[-2:]}")
    wfull = np.broadcast_to(w, xm.shape).copy()
    if mask is not None:
        m = np.broadcast_to(np.asarray(mask, dtype=bool), xm.shape[-2:])
        wfull = wfull * m
    valid = np.isfinite(xm) if skipna else np.ones(xm.shape, dtype=bool)
    wv = np.where(valid, wfull, 0.0)
    num = np.sum(np.where(valid, xm, 0.0) * wv, axis=(-2, -1))
    den = np.sum(wv, axis=(-2, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0, num / den, np.nan)
    if not skipna:  # NaN anywhere -> NaN (plain np.mean behaviour)
        out = np.where(np.any(~np.isfinite(xm) & (wfull > 0), axis=(-2, -1)), np.nan, out)
    return out, den


def area_mean(x, lat=None, lon=None, *, method="band", weights=None, mask=None,
              lat_bounds=None, lon_bounds=None, lat_axis=-2, lon_axis=-1,
              lat_name=None, lon_name=None, skipna=True):
    """Area-weighted spatial mean (works for time series: extra dims are kept).

    Parameters
    ----------
    x : ndarray or xarray.DataArray with a lat and a lon dimension.
    lat, lon : 1-D (or 2-D lat for curvilinear grids; then cos(lat) is used).
        Optional for DataArray (taken from the coordinates).
    method : "band" (default, exact spherical cell area), "cos", "ellipsoid"
        (WGS84), "none" (=plain unweighted mean, for comparison).
    weights : optional 2-D cell-area field (e.g. CMIP ``areacella``); overrides ``method``.
    mask : optional boolean field (True = include) for regional / land-only means.
    skipna : weights are renormalised over valid (finite) cells at every step.

    Returns an array with the lat/lon dimensions removed (DataArray in -> DataArray out).
    """
    try:
        import xarray as xr
    except Exception:  # pragma: no cover
        xr = None
    if xr is not None and isinstance(x, xr.DataArray):
        ln = _find(x, _LAT_NAMES, lat_name)
        on = _find(x, _LON_NAMES, lon_name)
        lat_v = x[ln].values if lat is None else lat
        lon_v = x[on].values if lon is None else lon
        if isinstance(mask, xr.DataArray):
            mask = mask.transpose(*[d for d in x.dims if d in mask.dims]).values
        if isinstance(weights, xr.DataArray):
            weights = weights.transpose(*[d for d in x.dims if d in weights.dims]).values
        dims = list(x.dims)
        la, lo = dims.index(ln) if ln in dims else None, dims.index(on) if on in dims else None
        if la is None or lo is None:   # 2-D coords on curvilinear grid: (y, x) dims
            la, lo = dims.index(x[ln].dims[0]), dims.index(x[ln].dims[1])
        # data layout must be (..., lat, lon) after moveaxis; weights built for (lat, lon) order
        res, _ = _core(x.values, _make_w(lat_v, lon_v, method, weights, lat_bounds, lon_bounds),
                       la, lo, mask, skipna)
        keep = [d for i, d in enumerate(dims) if i not in (la % x.ndim, lo % x.ndim)]
        coords = {k: v for k, v in x.coords.items()
                  if set(v.dims).issubset(keep)}
        return xr.DataArray(res, dims=keep, coords=coords, name=x.name, attrs=x.attrs)
    x = np.asarray(x, dtype=float)
    if lat is None and weights is None:
        raise ValueError("lat (or weights) is required for plain arrays")
    res, _ = _core(x, _make_w(lat, lon, method, weights, lat_bounds, lon_bounds),
                   lat_axis, lon_axis, mask, skipna)
    return res


def _make_w(lat, lon, method, weights, lat_bounds, lon_bounds):
    if weights is not None:
        return np.asarray(weights, dtype=float)
    return weights_2d(lat, lon, lat_bounds, lon_bounds, method, normalize=False)


def region_mean(x, lat=None, lon=None, mask=None, **kw):
    """Area-weighted mean inside ``mask`` (True = in region). Alias of ``area_mean``."""
    if mask is None:
        raise ValueError("mask is required")
    return area_mean(x, lat, lon, mask=mask, **kw)


def area_integral(x, lat, lon, *, ellipsoid=False, mask=None, lat_bounds=None,
                  lon_bounds=None, lat_axis=-2, lon_axis=-1, units="m2"):
    """Sum of x * cell_area over the grid (for fluxes per unit area -> totals).

    e.g. precipitation in mm/yr  ->  pass x*1e-3 (m/yr) and units="m2" to get m3/yr.
    NaN cells are skipped.
    """
    a = cell_area(lat, lon, lat_bounds, lon_bounds, ellipsoid=ellipsoid, units=units)
    x = np.asarray(x, dtype=float)
    xm = np.moveaxis(x, (lat_axis % x.ndim, lon_axis % x.ndim), (-2, -1))
    w = a if mask is None else a * np.asarray(mask, dtype=bool)
    return np.nansum(np.where(np.isfinite(xm), xm, 0.0) * w, axis=(-2, -1))


def area_fraction(cond, lat=None, lon=None, *, valid=None, **kw):
    """Area fraction (0-1) of cells where ``cond`` is True, among valid cells.

    The drought-area style statistic: ``area_fraction(spei < -1, lat, lon)``.
    NaN in a float ``cond`` / cells where ``valid`` is False are excluded.
    """
    try:
        import xarray as xr
    except Exception:  # pragma: no cover
        xr = None
    if xr is not None and isinstance(cond, xr.DataArray):
        c = cond.astype(float)
        if valid is not None:
            c = c.where(valid)
        return area_mean(c, lat, lon, **kw)
    c = np.asarray(cond).astype(float)
    if valid is not None:
        c = np.where(np.asarray(valid), c, np.nan)
    return area_mean(c, lat, lon, **kw)
