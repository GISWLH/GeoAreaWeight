"""Area weights for latitude-longitude grids.

All functions work on plain NumPy arrays; latitudes/longitudes in degrees.
"""
from __future__ import annotations

import numpy as np

# WGS84 ellipsoid
WGS84_A = 6378137.0                      # semi-major axis [m]
WGS84_F = 1.0 / 298.257223563            # flattening
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)     # first eccentricity squared
# IUGG mean radius R1, used as the default sphere radius [m]
EARTH_RADIUS = 6371008.8

__all__ = [
    "WGS84_A", "WGS84_F", "EARTH_RADIUS",
    "infer_bounds", "cos_lat_weights", "band_weights", "authalic_sin_lat",
    "authalic_radius", "area_weights", "cell_area", "weights_2d",
]


def _as_1d(x, name):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError(f"{name} must be 1-D (got shape {x.shape})")
    return x


def infer_bounds(centers, lo=None, hi=None):
    """Cell edges (n+1,) from cell centres (n,), also for irregular spacing.

    Edges are the mid-points between neighbouring centres; the two outer edges
    are extrapolated by half the adjacent spacing and clipped to [lo, hi].
    Works for ascending or descending centres.
    """
    c = _as_1d(centers, "centers")
    n = c.size
    if n < 2:
        raise ValueError("need >= 2 grid points to infer cell bounds; "
                         "pass explicit bounds instead")
    b = np.empty(n + 1)
    b[1:-1] = 0.5 * (c[:-1] + c[1:])
    b[0] = c[0] - 0.5 * (c[1] - c[0])
    b[-1] = c[-1] + 0.5 * (c[-1] - c[-2])
    if lo is not None or hi is not None:
        b = np.clip(b, lo if lo is not None else -np.inf,
                    hi if hi is not None else np.inf)
    return b


def _bounds_to_edges(bounds, n, name):
    """Accept (n+1,) edges or (n,2) [lower, upper] bounds -> (n, 2) array."""
    b = np.asarray(bounds, dtype=float)
    if b.ndim == 1 and b.size == n + 1:
        return np.column_stack([b[:-1], b[1:]])
    if b.ndim == 2 and b.shape == (n, 2):
        return b
    raise ValueError(f"{name} must have shape (n+1,) or (n, 2) with n={n}; got {b.shape}")


def cos_lat_weights(lat):
    """Classic cos(latitude) weights (negative values from rounding clipped to 0)."""
    lat = np.asarray(lat, dtype=float)
    return np.clip(np.cos(np.deg2rad(lat)), 0.0, None)


def authalic_sin_lat(lat_deg, e2=WGS84_E2):
    """sin(beta) where beta is the authalic latitude of geodetic latitude(s).

    sin(beta) = q(phi) / q(pi/2),  q = (1-e2) * [ sin/(1-e2 sin^2) - ln((1-e sin)/(1+e sin))/(2e) ]
    The area of a lat-lon cell on the ellipsoid is exactly Rq^2 * dlon * d sin(beta).
    """
    s = np.sin(np.deg2rad(np.asarray(lat_deg, dtype=float)))
    e = np.sqrt(e2)
    def q(sn):
        return (1.0 - e2) * (sn / (1.0 - e2 * sn ** 2)
                             - np.log((1.0 - e * sn) / (1.0 + e * sn)) / (2.0 * e))
    return q(s) / q(1.0)


def authalic_radius(a=WGS84_A, e2=WGS84_E2):
    """Radius of the sphere with the same surface area as the ellipsoid [m]."""
    e = np.sqrt(e2)
    qp = 1.0 + (1.0 - e2) / e * np.arctanh(e)
    return a * np.sqrt(qp / 2.0)


def band_weights(lat, lat_bounds=None, ellipsoid=False, normalize=True):
    """Exact area weight of each latitude band: |sin(phi_n) - sin(phi_s)|.

    With ``ellipsoid=True`` sin(phi) is replaced by sin(authalic latitude)
    (exact for the WGS84 ellipsoid).  Missing ``lat_bounds`` are inferred from
    the centres (``infer_bounds``, clipped to +-90), so poles points such as
    90N/90S on NCEP-style grids become half-cells (polar caps).
    """
    lat = _as_1d(lat, "lat")
    if lat_bounds is None:
        edges = np.column_stack([infer_bounds(lat, -90.0, 90.0)[:-1],
                                 infer_bounds(lat, -90.0, 90.0)[1:]])
    else:
        edges = _bounds_to_edges(lat_bounds, lat.size, "lat_bounds")
    if ellipsoid:
        s = authalic_sin_lat(edges)
    else:
        s = np.sin(np.deg2rad(edges))
    w = np.abs(s[:, 1] - s[:, 0])
    if normalize:
        w = w / w.sum()
    return w


def area_weights(lat, lat_bounds=None, method="band", normalize=True):
    """1-D latitude weights.

    method: ``"cos"``   cos(lat)  (proportional to area for *centre* latitudes,
                        approximate; ~exact for fine grids)
            ``"band"``  exact spherical band area from sin(lat) edges (default)
            ``"ellipsoid"`` exact WGS84 authalic-area weights
            ``"none"``  all ones (what you get with a plain np.mean)
    """
    lat = _as_1d(lat, "lat")
    m = method.lower()
    if m == "cos":
        w = cos_lat_weights(lat)
    elif m in ("band", "sphere", "sphere_band"):
        return band_weights(lat, lat_bounds, False, normalize)
    elif m in ("ellipsoid", "wgs84"):
        return band_weights(lat, lat_bounds, True, normalize)
    elif m in ("none", "equal", "unweighted"):
        w = np.ones_like(lat)
    else:
        raise ValueError(f"unknown method {method!r}")
    return w / w.sum() if normalize else w


def cell_area(lat, lon=None, lat_bounds=None, lon_bounds=None, ellipsoid=False,
              radius=None, units="m2"):
    """Area of each grid cell.

    Returns shape (nlat, nlon) if ``lon`` is given, else (nlat,) per unit of
    longitude (dlon = 1 rad).  Exact for lat-lon (graticule) cells.
    Irregular spacing is supported through the inferred/explicit bounds.
    ``ellipsoid=True`` uses WGS84 (authalic radius), otherwise ``radius``
    (default 6371008.8 m).  units: ``"m2"`` or ``"km2"``.
    """
    lat = _as_1d(lat, "lat")
    if ellipsoid:
        R = authalic_radius()
    else:
        R = EARTH_RADIUS if radius is None else float(radius)
    f = band_weights(lat, lat_bounds, ellipsoid, normalize=False)  # d sin(phi)
    if lon is None:
        area = R ** 2 * f
    else:
        lon = _as_1d(lon, "lon")
        if lon_bounds is None:
            lb = infer_bounds(lon)
            ed = np.column_stack([lb[:-1], lb[1:]])
        else:
            ed = _bounds_to_edges(lon_bounds, lon.size, "lon_bounds")
        dlon = np.deg2rad(np.abs(ed[:, 1] - ed[:, 0]))
        area = R ** 2 * np.outer(f, dlon)
    if units == "km2":
        area = area / 1e6
    elif units != "m2":
        raise ValueError("units must be 'm2' or 'km2'")
    return area


def weights_2d(lat, lon=None, lat_bounds=None, lon_bounds=None, method="band",
               normalize=True):
    """2-D weight field (nlat, nlon).

    * 1-D ``lat`` (+1-D ``lon``): cos / band / ellipsoid weights x longitude width
      (width only matters for irregular longitude spacing).
    * 2-D ``lat`` (curvilinear/rotated grids): cos(lat) approximation only; pass your
      own cell areas (e.g. CMIP ``areacella``) through ``weights=`` of ``area_mean``
      for exact results.
    """
    lat_a = np.asarray(lat, dtype=float)
    if lat_a.ndim == 2:
        w = cos_lat_weights(lat_a) if method.lower() in ("cos", "band", "sphere", "ellipsoid", "wgs84") \
            else np.ones_like(lat_a)
    else:
        w1 = area_weights(lat_a, lat_bounds, method, normalize=False)
        if method.lower() in ("none", "equal", "unweighted"):
            lon = np.arange(lat_a.size * 0 + (1 if lon is None else np.asarray(lon).size), dtype=float)
            w = np.ones((lat_a.size, lon.size))
            return w / w.sum() if normalize else w
        if lon is None:
            raise ValueError("lon is required to build a 2-D weight field")
        lon = _as_1d(lon, "lon")
        if lon_bounds is None and lon.size >= 2:
            lb = infer_bounds(lon)
            ed = np.column_stack([lb[:-1], lb[1:]])
        elif lon_bounds is not None:
            ed = _bounds_to_edges(lon_bounds, lon.size, "lon_bounds")
        else:
            ed = np.array([[0.0, 1.0]])
        dlon = np.abs(ed[:, 1] - ed[:, 0])
        w = np.outer(w1, dlon)
    if normalize:
        w = w / w.sum()
    return w
