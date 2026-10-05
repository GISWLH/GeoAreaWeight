"""Area-weighted means, regional means, area integrals and area fractions.

Every public function accepts either

* a NumPy array (or ``numpy.ma.MaskedArray``) plus ``lat``/``lon`` arrays and
  the positions of the latitude/longitude axes (``lat_axis``/``lon_axis``,
  default ``-2``/``-1``), or
* an ``xarray.DataArray``, whose latitude/longitude coordinates are found
  automatically (CF ``standard_name``/``units`` attributes, then the names
  ``lat``/``latitude``/``nav_lat`` and ``lon``/``longitude``/``nav_lon``,
  then rotated-pole ``rlat``/``rlon``). Dimension order is irrelevant, dask
  arrays stay lazy, and DataArray masks/weights are aligned by coordinate
  labels.

All dimensions other than latitude and longitude (time, level, ensemble ...)
are kept in the result.
"""
from __future__ import annotations

import sys
from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .weights import _canonical_method, cell_area, weights_2d

__all__ = ["area_mean", "region_mean", "area_integral", "area_fraction"]

_STD_NAMES = {"lat": ("latitude",), "lon": ("longitude",)}
_ROT_STD_NAMES = {"lat": ("grid_latitude",), "lon": ("grid_longitude",)}
_UNITS = {
    "lat": ("degrees_north", "degree_north", "degree_n", "degrees_n", "degreen", "degreesn"),
    "lon": ("degrees_east", "degree_east", "degree_e", "degrees_e", "degreee", "degreese"),
}
_NAMES = {"lat": ("lat", "latitude", "nav_lat"), "lon": ("lon", "longitude", "nav_lon")}
_ROT_NAMES = {"lat": ("rlat",), "lon": ("rlon",)}


# --------------------------------------------------------------------------- helpers
def _is_dataarray(x: Any) -> bool:
    xr = sys.modules.get("xarray")
    return xr is not None and isinstance(x, xr.DataArray)


def _is_dataset(x: Any) -> bool:
    xr = sys.modules.get("xarray")
    return xr is not None and isinstance(x, xr.Dataset)


def _to_float(x: Any) -> NDArray[np.float64]:
    """Float array; masked elements (e.g. from netCDF4) become NaN."""
    if np.ma.isMaskedArray(x):
        return np.ma.filled(x.astype(float), np.nan)
    return np.asarray(x, dtype=float)


def _grid_field(f: Any, shape: tuple[int, int], name: str, kind: str) -> NDArray:
    """Validate a mask or weight field against the ``(nlat, nlon)`` grid shape.

    A field given as ``(nlon, nlat)`` is transposed when ``nlat != nlon``.
    Masks: NaN / masked / 0 -> excluded. Weights: NaN / masked -> 0, negative -> error.
    """
    if np.ma.isMaskedArray(f):
        f = np.ma.filled(f.astype(float), np.nan)
    a = np.asarray(f)
    if a.ndim == 2 and a.shape != shape and a.shape == shape[::-1]:
        a = a.T
    try:
        a = np.broadcast_to(a, shape)
    except ValueError:
        raise ValueError(
            f"{name} has shape {np.shape(f)} but the data grid (lat, lon) has shape {shape}"
        ) from None
    if kind == "mask":
        if a.dtype == bool:
            return a
        a = a.astype(float)
        return np.isfinite(a) & (a != 0)
    a = a.astype(float)
    if np.any(a < 0):
        raise ValueError(f"{name} must be non-negative")
    return np.where(np.isfinite(a), a, 0.0)


def _reduce_core(
    xm: NDArray[np.float64],
    w: NDArray[np.float64],
    mask: NDArray[np.bool_] | None,
    skipna: bool,
    how: str,
) -> NDArray[np.float64]:
    """Weighted mean or weighted sum of ``xm`` over its last two axes (lat, lon)."""
    if mask is not None:
        w = w * mask
    wf = w.reshape(-1)
    X = xm.reshape(xm.shape[:-2] + (-1,))
    finite = np.isfinite(X)
    num = np.where(finite, X, 0.0) @ wf
    used = (wf > 0).astype(float)
    n_valid = finite.astype(float) @ used  # number of valid cells that carry weight (bool @ bool is bool)
    if how == "mean":
        den = finite @ wf if skipna else np.full(num.shape, wf.sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            out = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    else:
        out = num
    out = np.where(n_valid > 0, out, np.nan)
    if not skipna:
        out = np.where(n_valid < used.sum(), np.nan, out)
    return out[()] if out.ndim == 0 else out


def _find_coord(da: Any, axis: str) -> str | None:
    best: tuple[int, str] | None = None
    for name, c in da.coords.items():
        if c.ndim not in (1, 2):
            continue
        sn = str(c.attrs.get("standard_name", "")).lower()
        un = str(c.attrs.get("units", "")).lower()
        nm = str(name).lower()
        if sn in _STD_NAMES[axis] or un in _UNITS[axis] or nm in _NAMES[axis]:
            rank = 0 if c.ndim == 1 else 2
        elif sn in _ROT_STD_NAMES[axis] or nm in _ROT_NAMES[axis]:
            rank = 1 if c.ndim == 1 else 3  # rotated grids are regular in rotated coords
        else:
            continue
        if best is None or rank < best[0]:
            best = (rank, str(name))
    return None if best is None else best[1]


def _locate(da: Any, axis: str, name: str | None, values: Any) -> tuple[tuple, Any]:
    """Return (dims, values) of the latitude or longitude coordinate of ``da``."""
    if name is None:
        name = _find_coord(da, axis)
        if name is None:
            raise KeyError(
                f"could not identify the {'latitude' if axis == 'lat' else 'longitude'} "
                f"coordinate among coords {list(da.coords)} / dims {list(da.dims)}; "
                f"pass {axis}_name='...'"
            )
    if name in da.coords:
        c = da[name]
        return c.dims, (c.values if values is None else np.asarray(values, dtype=float))
    if name in da.dims:
        if values is None:
            raise ValueError(f"dimension {name!r} has no coordinate values; pass {axis}= explicitly")
        return (name,), np.asarray(values, dtype=float)
    raise KeyError(f"{name!r} is not a coordinate or dimension of the DataArray ({list(da.dims)})")


def _align_field(f: Any, da: Any, core: tuple[str, str], name: str) -> Any:
    """DataArray mask/weights -> NumPy array in ``core`` order, aligned by labels."""
    if not _is_dataarray(f):
        return f
    if f.ndim != 2 or set(f.dims) != set(core):
        raise ValueError(f"{name} must have exactly the dimensions {core}; got {f.dims}")
    for d in core:
        if d in f.indexes and d in da.indexes and not f.indexes[d].equals(da.indexes[d]):
            try:
                f = f.sel({d: da.indexes[d]})
            except KeyError:
                raise ValueError(
                    f"{name} coordinates along {d!r} do not match the data; "
                    "reindex or regrid it to the data grid first"
                ) from None
    return f.transpose(*core).values


def _dispatch(
    x: Any,
    lat: Any,
    lon: Any,
    *,
    how: str,
    build_w: Callable[[Any, Any], NDArray[np.float64]],
    weights: Any,
    mask: Any,
    lat_axis: int,
    lon_axis: int,
    lat_name: str | None,
    lon_name: str | None,
    skipna: bool,
    keep_attrs: bool,
) -> Any:
    if _is_dataset(x):
        raise TypeError("pass a DataArray (e.g. ds['tas']), not a Dataset")
    if _is_dataarray(x):
        import xarray as xr

        lat_dims, lat_v = _locate(x, "lat", lat_name, lat)
        lon_dims, lon_v = _locate(x, "lon", lon_name, lon)
        if len(lat_dims) == 1 and len(lon_dims) == 1:
            core = (lat_dims[0], lon_dims[0])
        elif len(lat_dims) == 2:
            core = tuple(lat_dims)
        else:
            raise ValueError(f"unsupported coordinate layout: lat dims {lat_dims}, lon dims {lon_dims}")
        if core[0] == core[1]:
            raise ValueError(f"latitude and longitude share the dimension {core[0]!r}")
        shape = (x.sizes[core[0]], x.sizes[core[1]])
        w_in = _align_field(weights, x, core, "weights")
        w = _grid_field(build_w(lat_v, lon_v) if w_in is None else w_in, shape, "weights", "weights")
        m_in = _align_field(mask, x, core, "mask")
        m = None if m_in is None else _grid_field(m_in, shape, "mask", "mask")
        return xr.apply_ufunc(
            lambda a: _reduce_core(np.asarray(a, dtype=float), w, m, skipna, how),
            x,
            input_core_dims=[list(core)],
            dask="parallelized",
            output_dtypes=[float],
            dask_gufunc_kwargs={"allow_rechunk": True},
            keep_attrs=keep_attrs,
        )

    xa = _to_float(x)
    if xa.ndim < 2:
        raise ValueError(f"x must have at least 2 dimensions (lat, lon); got shape {xa.shape}")
    la, lo = lat_axis % xa.ndim, lon_axis % xa.ndim
    if la == lo:
        raise ValueError("lat_axis and lon_axis must differ")
    xm = np.moveaxis(xa, (la, lo), (-2, -1))
    shape = xm.shape[-2:]
    if weights is None:
        if lat is None or lon is None:
            raise ValueError("lat and lon (or weights=) are required for NumPy input")
        w_raw = build_w(lat, lon)
        if w_raw.shape != shape:
            raise ValueError(
                f"lat/lon give a {w_raw.shape} grid but the data has {shape} along "
                f"(lat_axis={lat_axis}, lon_axis={lon_axis}); check the axis arguments"
            )
    else:
        w_raw = weights
    w = _grid_field(w_raw, shape, "weights", "weights")
    m = None if mask is None else _grid_field(mask, shape, "mask", "mask")
    return _reduce_core(xm, w, m, skipna, how)


# --------------------------------------------------------------------------- public API
def area_mean(
    x: Any,
    lat: ArrayLike | None = None,
    lon: ArrayLike | None = None,
    *,
    method: str = "band",
    weights: ArrayLike | None = None,
    mask: ArrayLike | None = None,
    lat_bounds: ArrayLike | None = None,
    lon_bounds: ArrayLike | None = None,
    lat_axis: int = -2,
    lon_axis: int = -1,
    lat_name: str | None = None,
    lon_name: str | None = None,
    skipna: bool = True,
) -> Any:
    """Area-weighted spatial mean over latitude and longitude.

    Use this for **intensive** quantities (temperature, precipitation rate,
    concentration, an index ...). Other dimensions such as time are kept.

    Parameters
    ----------
    x : ndarray, numpy.ma.MaskedArray or xarray.DataArray
        Data with a latitude and a longitude dimension. Masked elements count
        as missing.
    lat, lon : array_like, optional
        Cell-centre coordinates [deg]: 1-D, or 2-D ``lat`` for curvilinear
        grids (then ``cos(lat)`` is used). Required for NumPy input unless
        ``weights`` is given; optional overrides for a DataArray.
    method : {"band", "cos", "ellipsoid", "none"}, default "band"
        ``"band"``: exact spherical cell area. ``"cos"``: cos(lat).
        ``"ellipsoid"``: exact WGS84 cell area. ``"none"``: plain unweighted
        mean (for comparison only).
    weights : array_like or DataArray, shape (nlat, nlon), optional
        Your own cell areas (e.g. CMIP ``areacella``, or ``area * land_fraction``).
        Overrides ``method``. NaN weights count as 0.
    mask : array_like or DataArray, shape (nlat, nlon), optional
        ``True`` = include (regional / land-only means). For numeric masks,
        0 and NaN mean "exclude". A ``(nlon, nlat)`` mask is transposed
        automatically when ``nlat != nlon``.
    lat_bounds, lon_bounds : array_like, optional
        ``(n+1,)`` edges or ``(n, 2)`` bounds; inferred from centres if omitted.
    lat_axis, lon_axis : int, default -2, -1
        Positions of the latitude and longitude axes (NumPy input only).
    lat_name, lon_name : str, optional
        Coordinate names (DataArray only) if auto-detection fails.
    skipna : bool, default True
        If True, NaN cells are skipped and the weights are renormalised over the
        valid cells of every time step. If False, any NaN inside the (masked)
        domain gives NaN.

    Returns
    -------
    ndarray, scalar or DataArray
        ``x`` with the latitude and longitude dimensions removed. NaN where no
        valid cell remains.

    Examples
    --------
    >>> import numpy as np, geoareaweight as gaw
    >>> lat = np.arange(-89.5, 90, 1.0); lon = np.arange(0.5, 360, 1.0)
    >>> x = np.broadcast_to(np.abs(lat)[:, None], (lat.size, lon.size))
    >>> round(float(gaw.area_mean(x, lat, lon)), 2)   # sphere, analytic: 90 - 180/pi = 32.70
    32.71
    >>> round(float(gaw.area_mean(x, lat, lon, method="none")), 2)   # plain mean
    45.0
    """
    m = _canonical_method(method)

    def build(lat_v, lon_v):
        return weights_2d(lat_v, lon_v, lat_bounds, lon_bounds, m, normalize=False)

    return _dispatch(x, lat, lon, how="mean", build_w=build, weights=weights, mask=mask,
                     lat_axis=lat_axis, lon_axis=lon_axis, lat_name=lat_name,
                     lon_name=lon_name, skipna=skipna, keep_attrs=True)


def region_mean(
    x: Any,
    lat: ArrayLike | None = None,
    lon: ArrayLike | None = None,
    mask: ArrayLike | None = None,
    **kwargs: Any,
) -> Any:
    """Area-weighted mean inside ``mask`` (``True`` = in region).

    Same as ``area_mean(x, lat, lon, mask=mask, **kwargs)`` but ``mask`` is required.
    """
    if mask is None:
        raise ValueError("mask is required")
    return area_mean(x, lat, lon, mask=mask, **kwargs)


def area_integral(
    x: Any,
    lat: ArrayLike | None = None,
    lon: ArrayLike | None = None,
    *,
    ellipsoid: bool = False,
    radius: float | None = None,
    weights: ArrayLike | None = None,
    mask: ArrayLike | None = None,
    lat_bounds: ArrayLike | None = None,
    lon_bounds: ArrayLike | None = None,
    lat_axis: int = -2,
    lon_axis: int = -1,
    lat_name: str | None = None,
    lon_name: str | None = None,
    units: str = "m2",
    skipna: bool = True,
) -> Any:
    """Area integral ``sum(x * cell_area)`` over latitude and longitude.

    Use this to turn a quantity **per unit area** into a total, e.g. a water
    flux in m/yr into m3/yr, or an indicator (0/1) into an area in m2.

    Parameters
    ----------
    x : ndarray, numpy.ma.MaskedArray or xarray.DataArray
        Quantity per unit area. Convert units first: precipitation in mm/yr
        times ``1e-3`` gives m/yr, and with ``units="m2"`` the result is m3/yr.
    lat, lon : array_like, optional
        1-D cell-centre coordinates [deg] (required for NumPy input unless
        ``weights`` is given).
    ellipsoid : bool, default False
        Use WGS84 cell areas instead of a sphere.
    radius : float, optional
        Sphere radius [m]; default 6 371 008.8 m.
    weights : array_like or DataArray, shape (nlat, nlon), optional
        Your own cell areas, in the units you want (overrides the computed areas
        and ``units``).
    mask : array_like or DataArray, shape (nlat, nlon), optional
        ``True`` = include.
    lat_bounds, lon_bounds, lat_axis, lon_axis, lat_name, lon_name
        As in :func:`area_mean`.
    units : {"m2", "km2"}, default "m2"
        Area unit of the computed cell areas.
    skipna : bool, default True
        Skip NaN cells (they contribute 0). The result is NaN only when no valid
        cell remains. If False, any NaN inside the (masked) domain gives NaN.

    Returns
    -------
    ndarray, scalar or DataArray

    Examples
    --------
    >>> import numpy as np, geoareaweight as gaw
    >>> lat = np.arange(-89.5, 90, 1.0); lon = np.arange(0.5, 360, 1.0)
    >>> ones = np.ones((lat.size, lon.size))
    >>> round(float(gaw.area_integral(ones, lat, lon, units="km2")) / 1e6, 2)  # Earth, 1e6 km2
    510.07
    """
    def build(lat_v, lon_v):
        lat_a = np.asarray(lat_v, dtype=float)
        if lat_a.ndim != 1:
            raise ValueError("area_integral needs 1-D lat/lon; pass cell areas via weights=")
        return cell_area(lat_a, lon_v, lat_bounds, lon_bounds, ellipsoid=ellipsoid,
                         radius=radius, units=units)

    return _dispatch(x, lat, lon, how="sum", build_w=build, weights=weights, mask=mask,
                     lat_axis=lat_axis, lon_axis=lon_axis, lat_name=lat_name,
                     lon_name=lon_name, skipna=skipna, keep_attrs=False)


def area_fraction(
    cond: Any,
    lat: ArrayLike | None = None,
    lon: ArrayLike | None = None,
    *,
    valid: ArrayLike | None = None,
    **kwargs: Any,
) -> Any:
    """Fraction (0-1) of the valid area where ``cond`` is True.

    The "drought area" / "sea-ice extent fraction" statistic, e.g.
    ``area_fraction(spei <= -1, lat, lon)``. It is the area-weighted mean of
    the 0/1 indicator.

    Parameters
    ----------
    cond : array_like or DataArray
        Boolean condition, or a float array of 0/1 with NaN for invalid cells.
        Note that a comparison such as ``spei <= -1`` is ``False`` (not NaN)
        where ``spei`` is NaN; pass ``valid=np.isfinite(spei)`` (or
        ``spei.notnull()``) so that missing cells are excluded from the
        denominator.
    lat, lon : array_like, optional
        As in :func:`area_mean`.
    valid : array_like or DataArray, optional
        Boolean, broadcastable to ``cond``: cells to include in the denominator.
    **kwargs
        Passed to :func:`area_mean` (``method``, ``mask``, ``lat_axis`` ...).

    Examples
    --------
    >>> import numpy as np, geoareaweight as gaw
    >>> lat = np.arange(-89.5, 90, 1.0); lon = np.arange(0.5, 360, 1.0)
    >>> north_of_30 = np.broadcast_to(lat[:, None] >= 30, (lat.size, lon.size))
    >>> round(float(gaw.area_fraction(north_of_30, lat, lon)), 4)
    0.25
    """
    if _is_dataarray(cond):
        c = cond.astype(float)
        if valid is not None:
            c = c.where(valid)
    else:
        c = _to_float(cond)
        if valid is not None:
            c = np.where(np.asarray(valid, dtype=bool), c, np.nan)
    return area_mean(c, lat, lon, **kwargs)
