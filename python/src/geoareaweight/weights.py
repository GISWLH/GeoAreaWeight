"""Area weights and cell areas for latitude-longitude grids.

All functions take plain NumPy-compatible arrays. Angles are in **degrees**.

Conventions used throughout the package
---------------------------------------
* ``lat`` is a 1-D array of cell-centre latitudes (ascending or descending,
  regular or irregular), or a 2-D array for curvilinear grids.
* ``lat_bounds`` / ``lon_bounds`` are either ``n + 1`` cell edges or an
  ``(n, 2)`` array of ``[lower, upper]`` bounds (the CF ``*_bnds`` layout).
  When omitted they are inferred from the centres (see :func:`infer_bounds`).
* 2-D weight / area fields always have shape ``(nlat, nlon)``.
"""
from __future__ import annotations

import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "WGS84_A",
    "WGS84_F",
    "WGS84_E2",
    "EARTH_RADIUS",
    "METHODS",
    "infer_bounds",
    "cos_lat_weights",
    "band_weights",
    "authalic_sin_lat",
    "authalic_radius",
    "area_weights",
    "cell_area",
    "weights_2d",
]

#: WGS84 semi-major axis [m].
WGS84_A = 6378137.0
#: WGS84 flattening.
WGS84_F = 1.0 / 298.257223563
#: WGS84 first eccentricity squared.
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
#: IUGG mean Earth radius R1 [m]; default sphere radius.
EARTH_RADIUS = 6371008.8

#: Canonical weighting methods. Aliases accepted by ``method=`` arguments:
#: ``"sphere"``/``"sphere_band"`` -> ``"band"``, ``"wgs84"`` -> ``"ellipsoid"``,
#: ``"equal"``/``"unweighted"`` -> ``"none"``.
METHODS = ("band", "cos", "ellipsoid", "none")

_ALIASES = {
    "band": "band", "sphere": "band", "sphere_band": "band",
    "cos": "cos",
    "ellipsoid": "ellipsoid", "wgs84": "ellipsoid",
    "none": "none", "equal": "none", "unweighted": "none",
}

_LAT_TOL = 1e-6  # tolerance for |lat| slightly above 90 from float round-off


def _canonical_method(method: str) -> str:
    try:
        return _ALIASES[str(method).lower()]
    except KeyError:
        raise ValueError(
            f"unknown method {method!r}; expected one of {METHODS} "
            "(aliases: 'sphere', 'wgs84', 'unweighted')"
        ) from None


def _as_1d(x: ArrayLike, name: str) -> NDArray[np.float64]:
    a = np.asarray(x, dtype=float)
    if a.ndim == 0:
        a = a.reshape(1)
    if a.ndim != 1:
        raise ValueError(f"{name} must be 1-D (got shape {a.shape})")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} contains NaN or infinite values")
    return a


def _check_lat_range(a: NDArray[np.float64], name: str) -> NDArray[np.float64]:
    finite = a[np.isfinite(a)]
    if finite.size and np.max(np.abs(finite)) > 90.0 + _LAT_TOL:
        raise ValueError(
            f"{name} must be within [-90, 90] degrees (max |value| = "
            f"{np.max(np.abs(finite)):g}); did you swap lat and lon?"
        )
    return np.clip(a, -90.0, 90.0)


def _check_monotonic(c: NDArray[np.float64], name: str) -> None:
    d = np.diff(c)
    if not (np.all(d > 0) or np.all(d < 0)):
        raise ValueError(
            f"{name} must be strictly increasing or strictly decreasing to infer "
            "cell bounds; sort the grid or pass explicit bounds"
        )


def _unwrap_lon(lon: NDArray[np.float64]) -> NDArray[np.float64]:
    """Remove 360-degree jumps (e.g. 350, 355, 0, 5 -> 350, 355, 360, 365)."""
    if lon.size < 2:
        return lon
    d = np.mod(np.diff(lon) + 180.0, 360.0) - 180.0
    return np.concatenate([lon[:1], lon[0] + np.cumsum(d)])


def infer_bounds(
    centers: ArrayLike, lo: float | None = None, hi: float | None = None
) -> NDArray[np.float64]:
    """Infer ``n + 1`` cell edges from ``n`` cell centres.

    Inner edges are the mid-points between neighbouring centres; the two
    outer edges are extrapolated by half of the adjacent spacing and then
    clipped to ``[lo, hi]``. Works for irregular spacing and for ascending
    or descending centres.

    Parameters
    ----------
    centers : array_like, shape (n,)
        Strictly monotonic cell centres, n >= 2.
    lo, hi : float, optional
        Clip the edges to this range (e.g. ``-90, 90`` for latitude).

    Returns
    -------
    ndarray, shape (n + 1,)

    Examples
    --------
    >>> infer_bounds([-60.0, 0.0, 60.0], -90, 90)
    array([-90., -30.,  30.,  90.])
    """
    c = _as_1d(centers, "centers")
    n = c.size
    if n < 2:
        raise ValueError("need >= 2 grid points to infer cell bounds; pass explicit bounds")
    _check_monotonic(c, "centers")
    b = np.empty(n + 1)
    b[1:-1] = 0.5 * (c[:-1] + c[1:])
    b[0] = c[0] - 0.5 * (c[1] - c[0])
    b[-1] = c[-1] + 0.5 * (c[-1] - c[-2])
    if lo is not None or hi is not None:
        b = np.clip(b, -np.inf if lo is None else lo, np.inf if hi is None else hi)
    return b


def _bounds_to_edges(bounds: ArrayLike, n: int, name: str) -> NDArray[np.float64]:
    """Accept ``(n+1,)`` edges or ``(n, 2)`` bounds and return an ``(n, 2)`` array."""
    b = np.asarray(bounds, dtype=float)
    if b.ndim == 1 and b.size == n + 1:
        out = np.column_stack([b[:-1], b[1:]])
    elif b.ndim == 2 and b.shape == (n, 2):
        out = b
    else:
        raise ValueError(f"{name} must have shape (n+1,) or (n, 2) with n={n}; got {b.shape}")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} contains NaN or infinite values")
    return out


def _lat_edges(lat: NDArray[np.float64], lat_bounds: ArrayLike | None) -> NDArray[np.float64]:
    if lat_bounds is None:
        b = infer_bounds(lat, -90.0, 90.0)
        return np.column_stack([b[:-1], b[1:]])
    return _check_lat_range(_bounds_to_edges(lat_bounds, lat.size, "lat_bounds"), "lat_bounds")


def _lon_widths_deg(
    lon: ArrayLike, lon_bounds: ArrayLike | None, *, allow_single: bool
) -> NDArray[np.float64]:
    """Width [deg] of each longitude cell; handles the 0/360 or +-180 wrap."""
    lon = _as_1d(lon, "lon")
    if lon_bounds is not None:
        ed = _bounds_to_edges(lon_bounds, lon.size, "lon_bounds")
        w = np.abs(ed[:, 1] - ed[:, 0])
        # a bound pair such as [359, 1] or [179, -179] crosses the wrap point
        return np.where(w > 180.0, 360.0 - w, w)
    if lon.size == 1:
        if allow_single:
            return np.ones(1)  # any constant: relative weights are unaffected
        raise ValueError("cannot infer the width of a single longitude; pass lon_bounds")
    return np.abs(np.diff(infer_bounds(_unwrap_lon(lon))))


def cos_lat_weights(lat: ArrayLike) -> NDArray[np.float64]:
    """Classic ``cos(latitude)`` weights.

    Proportional to cell area only for regular grids and only approximately
    (error O(dlat^2)); a cell centred on a pole gets weight 0. Prefer
    :func:`band_weights` unless you need to reproduce the textbook formula.

    Examples
    --------
    >>> cos_lat_weights([0.0, 60.0, 90.0]).round(3)
    array([1. , 0.5, 0. ])
    """
    lat = _check_lat_range(np.asarray(lat, dtype=float), "lat")
    return np.clip(np.cos(np.deg2rad(lat)), 0.0, None)


def authalic_sin_lat(lat_deg: ArrayLike, e2: float = WGS84_E2) -> NDArray[np.float64]:
    """``sin(beta)``, where ``beta`` is the authalic latitude of geodetic latitude(s).

    ``sin(beta) = q(phi) / q(90 deg)`` with
    ``q = (1 - e2) * [sin / (1 - e2 sin^2) - ln((1 - e sin) / (1 + e sin)) / (2 e)]``.
    The area of an ellipsoidal lat-lon cell is exactly
    ``Rq^2 * dlon * (sin(beta_north) - sin(beta_south))``, where ``Rq`` is
    :func:`authalic_radius`.
    """
    s = np.sin(np.deg2rad(np.asarray(lat_deg, dtype=float)))
    e = np.sqrt(e2)

    def q(sn):
        return (1.0 - e2) * (sn / (1.0 - e2 * sn**2)
                             - np.log((1.0 - e * sn) / (1.0 + e * sn)) / (2.0 * e))

    return q(s) / q(1.0)


def authalic_radius(a: float = WGS84_A, e2: float = WGS84_E2) -> float:
    """Radius [m] of the sphere with the same surface area as the ellipsoid.

    Examples
    --------
    >>> round(authalic_radius(), 3)
    6371007.181
    """
    e = np.sqrt(e2)
    qp = 1.0 + (1.0 - e2) / e * np.arctanh(e)
    return float(a * np.sqrt(qp / 2.0))


def band_weights(
    lat: ArrayLike,
    lat_bounds: ArrayLike | None = None,
    ellipsoid: bool = False,
    normalize: bool = True,
) -> NDArray[np.float64]:
    """Exact area weight of each latitude band, ``|sin(phi_north) - sin(phi_south)|``.

    Parameters
    ----------
    lat : array_like, shape (nlat,)
        Cell-centre latitudes [deg]; ascending, descending or irregular.
    lat_bounds : array_like, optional
        ``(nlat+1,)`` edges or ``(nlat, 2)`` bounds. If omitted they are
        inferred with :func:`infer_bounds` and clipped to +-90, so pole points
        such as 90N/90S on NCEP-style grids become half-height polar caps.
    ellipsoid : bool, default False
        Use the WGS84 ellipsoid (sin of the authalic latitude) instead of a sphere.
    normalize : bool, default True
        Scale the weights to sum to 1.

    Returns
    -------
    ndarray, shape (nlat,)

    Examples
    --------
    >>> band_weights([-45.0, 45.0], lat_bounds=[-90, 0, 90])
    array([0.5, 0.5])
    """
    lat = _check_lat_range(_as_1d(lat, "lat"), "lat")
    edges = _lat_edges(lat, lat_bounds)
    s = authalic_sin_lat(edges) if ellipsoid else np.sin(np.deg2rad(edges))
    w = np.abs(s[:, 1] - s[:, 0])
    if normalize:
        w = w / w.sum()
    return w


def area_weights(
    lat: ArrayLike,
    lat_bounds: ArrayLike | None = None,
    method: str = "band",
    normalize: bool = True,
) -> NDArray[np.float64]:
    """1-D latitude weights (one value per latitude row).

    Parameters
    ----------
    lat : array_like, shape (nlat,)
        Cell-centre latitudes [deg].
    lat_bounds : array_like, optional
        Cell bounds; used by ``"band"`` and ``"ellipsoid"`` only.
    method : {"band", "cos", "ellipsoid", "none"}, default "band"
        ``"band"``: exact spherical band area from sin(lat) edges.
        ``"cos"``: cos(lat) of the centres (approximation).
        ``"ellipsoid"``: exact WGS84 band area (authalic latitude).
        ``"none"``: all ones, i.e. what a plain ``np.mean`` does.
    normalize : bool, default True
        Scale the weights to sum to 1.

    Returns
    -------
    ndarray, shape (nlat,)

    Examples
    --------
    >>> area_weights([0.0, 60.0], method="cos", normalize=False).round(3)
    array([1. , 0.5])
    """
    m = _canonical_method(method)
    if m == "band":
        return band_weights(lat, lat_bounds, False, normalize)
    if m == "ellipsoid":
        return band_weights(lat, lat_bounds, True, normalize)
    lat = _check_lat_range(_as_1d(lat, "lat"), "lat")
    w = cos_lat_weights(lat) if m == "cos" else np.ones_like(lat)
    return w / w.sum() if normalize else w


def cell_area(
    lat: ArrayLike,
    lon: ArrayLike | None = None,
    lat_bounds: ArrayLike | None = None,
    lon_bounds: ArrayLike | None = None,
    ellipsoid: bool = False,
    radius: float | None = None,
    units: str = "m2",
) -> NDArray[np.float64]:
    """Area of each grid cell, exact for latitude-longitude (graticule) cells.

    Parameters
    ----------
    lat : array_like, shape (nlat,)
        Cell-centre latitudes [deg].
    lon : array_like, shape (nlon,), optional
        Cell-centre longitudes [deg]. Grids crossing the 0/360 or +-180
        meridian (e.g. ``[350, 355, 0, 5]``) are handled.
    lat_bounds, lon_bounds : array_like, optional
        ``(n+1,)`` edges or ``(n, 2)`` bounds; inferred from centres if omitted.
    ellipsoid : bool, default False
        WGS84 ellipsoid (authalic radius 6 371 007.18 m) instead of a sphere.
    radius : float, optional
        Sphere radius [m]; default :data:`EARTH_RADIUS` (6 371 008.8 m).
        Ignored when ``ellipsoid=True``.
    units : {"m2", "km2"}, default "m2"

    Returns
    -------
    ndarray
        Shape ``(nlat, nlon)`` if ``lon`` is given, otherwise ``(nlat,)``
        holding the area per radian of longitude.

    Examples
    --------
    >>> a = cell_area([0.0], [0.5], lat_bounds=[-0.5, 0.5], lon_bounds=[0, 1], units="km2")
    >>> round(float(a[0, 0]), 1)
    12364.2
    """
    if units not in ("m2", "km2"):
        raise ValueError("units must be 'm2' or 'km2'")
    R = authalic_radius() if ellipsoid else (EARTH_RADIUS if radius is None else float(radius))
    f = band_weights(lat, lat_bounds, ellipsoid, normalize=False)  # d sin(phi)
    if lon is None:
        area = R**2 * f
    else:
        dlon = np.deg2rad(_lon_widths_deg(lon, lon_bounds, allow_single=False))
        area = R**2 * np.outer(f, dlon)
    return area / 1e6 if units == "km2" else area


def weights_2d(
    lat: ArrayLike,
    lon: ArrayLike | None = None,
    lat_bounds: ArrayLike | None = None,
    lon_bounds: ArrayLike | None = None,
    method: str = "band",
    normalize: bool = True,
) -> NDArray[np.float64]:
    """2-D weight field of shape ``(nlat, nlon)``.

    * 1-D ``lat`` and ``lon``: latitude weights (see :func:`area_weights`)
      times longitude cell width (the width only matters for irregular
      longitude spacing). Proportional to :func:`cell_area` for ``"band"``
      and ``"ellipsoid"``.
    * 2-D ``lat`` (curvilinear grids, shape ``(ny, nx)``): only ``cos(lat)``
      is possible, and a :class:`UserWarning` is issued for ``"band"`` or
      ``"ellipsoid"``. For exact results pass the model's own cell areas
      (e.g. CMIP ``areacella``) as ``weights=`` to :func:`area_mean`.

    Parameters
    ----------
    lat, lon : array_like
        Cell centres [deg]. ``lon`` is required for 1-D ``lat``.
    lat_bounds, lon_bounds : array_like, optional
        Cell bounds (1-D grids only).
    method : {"band", "cos", "ellipsoid", "none"}, default "band"
    normalize : bool, default True
        Scale the field to sum to 1.

    Returns
    -------
    ndarray, shape (nlat, nlon) or (ny, nx)
    """
    m = _canonical_method(method)
    lat_a = np.asarray(lat, dtype=float)
    if lat_a.ndim == 2:
        if m == "none":
            w = np.ones_like(lat_a)
        else:
            if m != "cos":
                warnings.warn(
                    f"method={method!r} needs 1-D latitudes; using cos(lat) for the 2-D "
                    "(curvilinear) grid. Pass cell areas via weights= for exact results.",
                    UserWarning, stacklevel=2,
                )
            w = cos_lat_weights(lat_a)
    elif lat_a.ndim <= 1:
        if lon is None:
            raise ValueError("lon is required to build a 2-D weight field")
        lat_w = area_weights(lat_a, lat_bounds, m, normalize=False)
        if m == "none":
            w = np.ones((lat_w.size, _as_1d(lon, "lon").size))
        else:
            w = np.outer(lat_w, _lon_widths_deg(lon, lon_bounds, allow_single=True))
    else:
        raise ValueError(f"lat must be 1-D or 2-D (got shape {lat_a.shape})")
    return w / w.sum() if normalize else w
