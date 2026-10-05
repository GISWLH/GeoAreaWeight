"""Tests for the xarray interface and independent cross-checks (optional dependencies)."""
import numpy as np
import pytest

import geoareaweight as gaw

from conftest import grid

R = gaw.EARTH_RADIUS

xr = pytest.importorskip("xarray")


def test_xarray_roundtrip_matches_numpy():
    lat, lon = grid(5.0)
    t = np.arange(6)
    data = np.random.default_rng(4).normal(size=(6, lat.size, lon.size))
    da = xr.DataArray(data, dims=("time", "lat", "lon"), coords={"time": t, "lat": lat, "lon": lon}, name="v")
    r = gaw.area_mean(da)
    assert isinstance(r, xr.DataArray) and r.dims == ("time",)
    np.testing.assert_allclose(r.values, gaw.area_mean(data, lat, lon))
    # (lon, lat, time) order and descending lat
    da2 = da.transpose("lon", "lat", "time").isel(lat=slice(None, None, -1))
    np.testing.assert_allclose(gaw.area_mean(da2).values, r.values)
    r_ref = xr.DataArray(data, dims=("time", "lat", "lon"), coords={"lat": lat, "lon": lon}).weighted(
        xr.DataArray(gaw.cell_area(lat, lon), dims=("lat", "lon"), coords={"lat": lat, "lon": lon})
    ).mean(("lat", "lon"))
    np.testing.assert_allclose(r.values, r_ref.values)


def test_wgs84_cell_area_matches_cylindrical_equal_area_projection():
    """Independent check: EPSG:6933 (WGS84 cylindrical equal-area) -> cell area = dx*dy."""
    pyproj = pytest.importorskip("pyproj")
    tr = pyproj.Transformer.from_crs(4326, 6933, always_xy=True)
    for lat0, lat1 in [(0, 1), (40, 41), (40, 50), (-89, -80), (80, 90), (65, 66)]:
        x0, y0 = tr.transform(10.0, lat0)
        x1, y1 = tr.transform(11.0, lat1)
        ref = abs((x1 - x0) * (y1 - y0))
        a = gaw.cell_area(np.array([(lat0 + lat1) / 2]), np.array([10.5]), lat_bounds=[lat0, lat1],
                          lon_bounds=[10, 11], ellipsoid=True)[0, 0]
        assert a == pytest.approx(ref, rel=1e-9), (lat0, lat1)


def test_xarray_mask_follows_dimension_order_and_labels():
    lat, lon = grid(5.0)
    data = np.random.default_rng(5).normal(size=(3, lat.size, lon.size))
    da = xr.DataArray(data, dims=("time", "lat", "lon"), coords={"lat": lat, "lon": lon})
    mask = xr.DataArray(np.broadcast_to(lat[:, None] > 30, (lat.size, lon.size)),
                        dims=("lat", "lon"), coords={"lat": lat, "lon": lon})
    ref = gaw.area_mean(da, mask=mask).values
    np.testing.assert_allclose(ref, gaw.area_mean(data, lat, lon, mask=mask.values))
    # data transposed -> mask must still be read as (lat, lon)
    np.testing.assert_allclose(gaw.area_mean(da.transpose("lon", "lat", "time"), mask=mask).values, ref)
    # data with descending latitude (ERA5 style) -> mask aligned by label, not position
    np.testing.assert_allclose(gaw.area_mean(da.isel(lat=slice(None, None, -1)), mask=mask).values, ref)
    # weights given as a DataArray are aligned the same way
    area = xr.DataArray(gaw.cell_area(lat, lon), dims=("lat", "lon"), coords={"lat": lat, "lon": lon})
    np.testing.assert_allclose(gaw.area_mean(da.isel(lat=slice(None, None, -1)), weights=area).values,
                               gaw.area_mean(da).values)
    with pytest.raises(ValueError, match="do not match"):
        gaw.area_mean(da, mask=mask.assign_coords(lat=mask.lat + 0.1))


def test_projected_xy_dims_are_not_mistaken_for_lat_lon():
    da = xr.DataArray(np.ones((2, 3, 4)), dims=("time", "y", "x"),
                      coords={"y": [0.0, 1e3, 2e3], "x": [0.0, 1e3, 2e3, 3e3]})
    with pytest.raises(KeyError, match="lat_name"):
        gaw.area_mean(da)


def test_cf_attributes_and_non_dimension_coords():
    lat, lon = grid(10.0)
    data = np.random.default_rng(6).normal(size=(lat.size, lon.size))
    ref = gaw.area_mean(data, lat, lon)
    da = xr.DataArray(data, dims=("j", "i"), coords={
        "phi": ("j", lat, {"standard_name": "latitude"}),
        "lam": ("i", lon, {"units": "degrees_east"})})
    assert float(gaw.area_mean(da)) == pytest.approx(ref)
    da2 = xr.DataArray(data, dims=("j", "i"))
    assert float(gaw.area_mean(da2, lat, lon, lat_name="j", lon_name="i")) == pytest.approx(ref)


def test_rotated_pole_grid_uses_rotated_coordinates():
    """A rotated lat-lon grid is regular on the rotated sphere: band weights in rlat are exact."""
    rlat, rlon = grid(10.0)
    data = np.random.default_rng(7).normal(size=(rlat.size, rlon.size))
    lat2d = np.full(data.shape, 45.0)   # dummy geographic latitude (would give equal weights)
    da = xr.DataArray(data, dims=("rlat", "rlon"),
                      coords={"rlat": rlat, "rlon": rlon, "lat": (("rlat", "rlon"), lat2d),
                              "lon": (("rlat", "rlon"), np.zeros(data.shape))})
    assert float(gaw.area_mean(da)) == pytest.approx(gaw.area_mean(data, rlat, rlon))


def test_input_validation():
    with pytest.raises(ValueError, match="swap lat and lon"):
        gaw.area_weights(np.arange(0, 360, 10.0))
    with pytest.raises(ValueError, match="strictly"):
        gaw.band_weights([0.0, 10.0, 5.0])
    with pytest.raises(ValueError, match="single longitude"):
        gaw.cell_area([0.0, 10.0], [5.0])
    with pytest.raises(ValueError, match="lat_axis"):
        gaw.area_mean(np.ones((4, 3)), [0.0, 10.0, 20.0], [0.0, 1.0, 2.0, 3.0])
    with pytest.raises(TypeError, match="DataArray"):
        gaw.area_mean(xr.Dataset({"v": (("lat", "lon"), np.ones((2, 2)))}))
    # lat with tiny round-off beyond the pole is accepted
    assert gaw.band_weights([-90.0000000001, 0.0, 90.0000000001]).sum() == pytest.approx(1.0)


def test_area_integral_xarray_nan_and_skipna():
    lat, lon = grid(10.0)
    x = np.ones((2, lat.size, lon.size))
    x[1] = np.nan
    x[0, 0, 0] = np.nan
    r = gaw.area_integral(x, lat, lon)
    a = gaw.cell_area(lat, lon)
    assert r[0] == pytest.approx(a.sum() - a[0, 0])
    assert np.isnan(r[1])                                     # nothing valid -> NaN, not 0
    assert np.isnan(gaw.area_integral(x, lat, lon, skipna=False)[0])
    da = xr.DataArray(x, dims=("time", "lat", "lon"), coords={"lat": lat, "lon": lon}, attrs={"units": "m"})
    rx = gaw.area_integral(da, units="km2")
    assert isinstance(rx, xr.DataArray) and "units" not in rx.attrs
    np.testing.assert_allclose(rx.values, r / 1e6)
    assert float(gaw.area_integral(da[0], weights=a)) == pytest.approx(r[0])


def test_area_fraction_valid_and_xarray():
    lat, lon = grid(5.0)
    spei = np.random.default_rng(9).normal(size=(lat.size, lon.size))
    spei[lat < -60, :] = np.nan
    f = gaw.area_fraction(spei <= -1, lat, lon, valid=np.isfinite(spei))
    ref = gaw.area_mean(np.where(np.isfinite(spei), spei <= -1, np.nan), lat, lon)
    assert f == pytest.approx(ref)
    da = xr.DataArray(spei, dims=("lat", "lon"), coords={"lat": lat, "lon": lon})
    assert float(gaw.area_fraction(da <= -1, valid=da.notnull())) == pytest.approx(ref)


def test_xarray_keeps_attrs_name_and_coords():
    lat, lon = grid(10.0)
    t = np.arange(3)
    da = xr.DataArray(np.ones((3, lat.size, lon.size)), dims=("time", "lat", "lon"),
                      coords={"time": t, "lat": lat, "lon": lon}, name="tas", attrs={"units": "K"})
    r = gaw.area_mean(da)
    assert r.name == "tas" and r.attrs["units"] == "K" and list(r.time.values) == list(t)


def test_dask_stays_lazy():
    pytest.importorskip("dask")
    lat, lon = grid(5.0)
    data = np.random.default_rng(10).normal(size=(8, lat.size, lon.size))
    da = xr.DataArray(data, dims=("time", "lat", "lon"), coords={"lat": lat, "lon": lon}).chunk(
        {"time": 2, "lat": 12})
    r = gaw.area_mean(da)
    assert r.chunks is not None
    np.testing.assert_allclose(r.compute().values, gaw.area_mean(data, lat, lon))
