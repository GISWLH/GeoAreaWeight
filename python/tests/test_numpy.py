"""Tests for the NumPy interface (no optional dependencies)."""
import numpy as np
import pytest

import geoareaweight as gaw

from conftest import grid

R = gaw.EARTH_RADIUS


@pytest.mark.parametrize("res,pole", [(1.0, False), (2.5, True), (0.25, False), (5.0, False)])
def test_sphere_area_sums_to_total(res, pole):
    lat, lon = grid(res, pole)
    a = gaw.cell_area(lat, lon)
    assert a.sum() == pytest.approx(4 * np.pi * R ** 2, rel=1e-12)


def test_ellipsoid_total_area_wgs84():
    lat, lon = grid(1.0)
    a = gaw.cell_area(lat, lon, ellipsoid=True)
    # WGS84 surface area = 510 065 621.7 km2 (authalic sphere)
    assert a.sum() / 1e6 == pytest.approx(510065621.7, rel=1e-7)
    assert gaw.authalic_radius() == pytest.approx(6371007.1809, abs=1e-3)


def test_known_cell_area_1deg_equator():
    # 1x1 deg cell at the equator on the IUGG sphere ~ 12 364 km2
    a = gaw.cell_area(np.array([0.0]), np.array([0.5]), lat_bounds=[-0.5, 0.5],
                      lon_bounds=[0.0, 1.0], units="km2")
    assert a[0, 0] == pytest.approx(R ** 2 * np.deg2rad(1) * 2 * np.sin(np.deg2rad(0.5)) / 1e6)
    assert 12300 < a[0, 0] < 12400


def test_cos_weights_close_to_band_for_fine_grids():
    lat, _ = grid(0.25)
    wc = gaw.area_weights(lat, method="cos")
    wb = gaw.area_weights(lat, method="band")
    assert np.max(np.abs(wc - wb) / wb) < 1e-3 * 5   # excl. poles tiny differences


def test_descending_and_irregular_lat():
    lat = np.array([80., 60., 30., 10., -5., -40., -75.])   # irregular, descending
    w = gaw.band_weights(lat)
    assert w.sum() == pytest.approx(1.0)
    w_asc = gaw.band_weights(lat[::-1])
    np.testing.assert_allclose(w, w_asc[::-1])
    a = gaw.cell_area(lat)          # per radian of longitude
    assert 2 * np.pi * a.sum() == pytest.approx(4 * np.pi * R ** 2, rel=1e-12)


def test_explicit_bounds_forms():
    lat = np.array([-45., 45.])
    e1 = gaw.band_weights(lat, lat_bounds=[-90, 0, 90])
    e2 = gaw.band_weights(lat, lat_bounds=[[-90, 0], [0, 90]])
    np.testing.assert_allclose(e1, [0.5, 0.5])
    np.testing.assert_allclose(e1, e2)
    with pytest.raises(ValueError):
        gaw.band_weights(lat, lat_bounds=[0, 1])


def test_constant_field_any_weights():
    lat, lon = grid(2.0)
    x = np.full((lat.size, lon.size), 3.7)
    for m in ("cos", "band", "ellipsoid", "none"):
        assert gaw.area_mean(x, lat, lon, method=m) == pytest.approx(3.7)


def test_mean_of_sin_lat_is_zero_and_unweighted_is_not():
    lat, lon = grid(1.0)
    x = np.broadcast_to(np.sin(np.deg2rad(lat))[:, None] ** 3, (lat.size, lon.size))
    assert abs(gaw.area_mean(x, lat, lon)) < 1e-12
    # field = |lat| : area-weighted mean over sphere = 90*(pi/2-1)*... analytic: (pi/2 - 1) rad
    y = np.broadcast_to(np.abs(lat)[:, None], (lat.size, lon.size))
    exact = np.rad2deg(np.pi / 2 - 1)
    assert gaw.area_mean(y, lat, lon) == pytest.approx(exact, abs=3e-3)   # O(dlat^2) discretisation
    assert gaw.area_mean(y, lat, lon, method="none") == pytest.approx(45.0, abs=1e-9)


def test_cos_vs_band_on_coarse_grid_differs_slightly():
    lat, lon = grid(10.0)
    y = np.broadcast_to(np.abs(lat)[:, None], (lat.size, lon.size))
    mc, mb = gaw.area_mean(y, lat, lon, method="cos"), gaw.area_mean(y, lat, lon, method="band")
    assert mc != mb and abs(mc - mb) < 0.5


def test_nan_renormalisation_and_mask():
    lat, lon = grid(5.0)
    rng = np.random.default_rng(0)
    x = rng.normal(size=(12, lat.size, lon.size))
    x[:, :5, :] = np.nan          # missing southern rows
    x[3, 10:20, 10:30] = np.nan   # time-dependent hole
    m = gaw.area_mean(x, lat, lon)
    assert m.shape == (12,)
    w = gaw.weights_2d(lat, lon, normalize=False)
    for t in range(12):
        v = np.isfinite(x[t])
        ref = np.sum(x[t][v] * w[v]) / np.sum(w[v])
        assert m[t] == pytest.approx(ref, rel=1e-12)
    # skipna=False propagates NaN
    assert np.isnan(gaw.area_mean(x, lat, lon, skipna=False)).all()
    # mask
    mask = np.zeros((lat.size, lon.size), bool)
    mask[20:30, 5:40] = True
    mm = gaw.area_mean(x, lat, lon, mask=mask)
    ref = np.sum(x[0][mask] * w[mask]) / np.sum(w[mask])
    assert mm[0] == pytest.approx(ref)
    # all-masked -> NaN
    assert np.isnan(gaw.area_mean(x, lat, lon, mask=np.zeros_like(mask))).all()


def test_axis_handling():
    lat, lon = grid(10.0)
    x = np.random.default_rng(1).normal(size=(lon.size, 4, lat.size))   # (lon, time, lat)
    r = gaw.area_mean(x, lat, lon, lat_axis=2, lon_axis=0)
    r2 = gaw.area_mean(np.transpose(x, (1, 2, 0)), lat, lon)
    np.testing.assert_allclose(r, r2)


def test_user_supplied_weights_areacella():
    lat, lon = grid(10.0)
    x = np.random.default_rng(2).normal(size=(lat.size, lon.size))
    a = gaw.cell_area(lat, lon)
    assert gaw.area_mean(x, weights=a) == pytest.approx(gaw.area_mean(x, lat, lon))


def test_area_integral_and_fraction():
    lat, lon = grid(1.0)
    x = np.ones((lat.size, lon.size))
    assert gaw.area_integral(x, lat, lon) == pytest.approx(4 * np.pi * R ** 2)
    cond = np.broadcast_to(lat[:, None] > 0, x.shape)   # northern hemisphere
    assert gaw.area_fraction(cond, lat, lon) == pytest.approx(0.5)
    # a 30-90N cap covers 25 % of the sphere
    cond2 = np.broadcast_to(lat[:, None] >= 30, x.shape)
    assert gaw.area_fraction(cond2, lat, lon) == pytest.approx(0.25, abs=1e-12)
    # unweighted fraction is wrong
    assert gaw.area_mean(cond2.astype(float), lat, lon, method="none") == pytest.approx(60 / 180)


def test_2d_lat_cos_fallback():
    lat1, lon1 = grid(10.0)
    lon2, lat2 = np.meshgrid(lon1, lat1)
    x = np.random.default_rng(3).normal(size=lat2.shape)
    r2 = gaw.area_mean(x, lat2, lon2, method="cos")
    r1 = gaw.area_mean(x, lat1, lon1, method="cos")
    assert r2 == pytest.approx(r1)


def test_irregular_longitude():
    lat = np.array([-60., 0., 60.])
    lon = np.array([0., 10., 30., 90.])
    a = gaw.cell_area(lat, lon)
    # wider lon cells -> bigger area in same latitude row
    assert np.all(np.diff(a, axis=1)[:, 1] > 0)


def test_unknown_method_is_rejected_everywhere():
    lat, lon = grid(10.0)
    lat2, _ = np.meshgrid(lat, lon, indexing="ij")
    for f in (lambda: gaw.area_weights(lat, method="nonsense"),
              lambda: gaw.weights_2d(lat, lon, method="nonsense"),
              lambda: gaw.weights_2d(lat2, method="nonsense"),
              lambda: gaw.area_mean(np.ones(lat2.shape), lat, lon, method="nonsense")):
        with pytest.raises(ValueError, match="unknown method"):
            f()


def test_2d_lat_band_warns_and_aliases_work():
    lat, lon = grid(10.0)
    lat2, _ = np.meshgrid(lat, lon, indexing="ij")
    with pytest.warns(UserWarning, match="cos"):
        w = gaw.weights_2d(lat2, method="sphere_band", normalize=False)
    np.testing.assert_allclose(w, np.cos(np.deg2rad(lat2)))
    np.testing.assert_allclose(gaw.area_weights(lat, method="WGS84"),
                               gaw.area_weights(lat, method="ellipsoid"))
    np.testing.assert_allclose(gaw.weights_2d(lat, lon, method="unweighted"), 1 / lat2.size)


def test_masked_arrays_are_treated_as_missing():
    lat, lon = grid(10.0)
    x = np.ma.masked_array(np.ones((lat.size, lon.size)), mask=False)
    x.data[0, 0] = 9.96921e36   # netCDF fill value
    x[0, 0] = np.ma.masked
    assert gaw.area_mean(x, lat, lon) == pytest.approx(1.0)


def test_numeric_masks_nan_and_zero_exclude():
    lat, lon = grid(10.0)
    x = np.random.default_rng(8).normal(size=(lat.size, lon.size))
    sel = np.broadcast_to(lat[:, None] > 30, x.shape)
    ref = gaw.area_mean(x, lat, lon, mask=sel)
    assert gaw.area_mean(x, lat, lon, mask=np.where(sel, 1.0, np.nan)) == pytest.approx(ref)
    assert gaw.area_mean(x, lat, lon, mask=sel.astype(int)) == pytest.approx(ref)
    # (nlon, nlat) masks / weights are transposed automatically
    assert gaw.area_mean(x, lat, lon, mask=sel.T) == pytest.approx(ref)
    a = gaw.cell_area(lat, lon)
    assert gaw.area_mean(x, weights=a.T) == pytest.approx(gaw.area_mean(x, lat, lon))
    with pytest.raises(ValueError, match="non-negative"):
        gaw.area_mean(x, weights=-a)
    with pytest.raises(ValueError, match="shape"):
        gaw.area_mean(x, lat, lon, mask=np.ones((3, 3), bool))


def test_longitudes_crossing_the_wrap_point():
    lat = np.array([-45.0, 45.0])
    lon_wrapped = np.r_[np.arange(182.5, 360, 5.0), np.arange(2.5, 180, 5.0)]
    a = gaw.cell_area(lat, lon_wrapped)
    np.testing.assert_allclose(a, a[:, :1] * np.ones_like(a))
    assert a.sum() == pytest.approx(4 * np.pi * R ** 2 * (gaw.band_weights(lat, normalize=False).sum() / 2))
    np.testing.assert_allclose(gaw.cell_area(lat, [359.5, 0.5], lon_bounds=[[359, 0], [0, 1]]),
                               gaw.cell_area(lat, [0.5, 1.5], lon_bounds=[0, 1, 2]))
    np.testing.assert_allclose(gaw.cell_area(lat, np.arange(-177.5, 180, 5.0)), a)


def test_skipna_false_on_complete_data():
    lat, lon = grid(10.0)
    x = np.random.default_rng(11).normal(size=(3, lat.size, lon.size))
    np.testing.assert_allclose(gaw.area_mean(x, lat, lon, skipna=False), gaw.area_mean(x, lat, lon))
    np.testing.assert_allclose(gaw.area_integral(x, lat, lon, skipna=False), gaw.area_integral(x, lat, lon))
    x[1, 0, 0] = np.nan
    r = gaw.area_mean(x, lat, lon, skipna=False)
    assert np.isnan(r[1]) and np.isfinite(r[[0, 2]]).all()
    # a NaN outside the mask does not matter
    mask = np.ones(x.shape[1:], bool)
    mask[0, 0] = False
    assert np.isfinite(gaw.area_mean(x, lat, lon, mask=mask, skipna=False)).all()
