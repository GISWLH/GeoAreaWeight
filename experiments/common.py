import warnings, numpy as np, xarray as xr
warnings.filterwarnings("ignore")
import geoareaweight as gaw
D = "/workspace/geoareaweight/data"
SEC = 86400.0

def op(fn, **kw):
    return xr.open_dataset(f"{D}/{fn}", drop_variables=["time_bnds", "time_bounds", "lat_bnds", "lon_bnds"],
                           decode_timedelta=False, **kw)

def annual(da):
    """Annual mean of a monthly series (day-weighted), complete years only."""
    da = da.sel(time=slice(None, None))
    yrs = da.time.dt.year
    cnt = da.groupby(yrs).count("time")  # not used
    n = da.time.groupby(yrs).count()
    full = n.year[n == 12]
    da = da.sel(time=yrs.isin(full))
    w = da.time.dt.days_in_month
    out = (da * w).groupby(da.time.dt.year).sum("time", skipna=False) / w.groupby(da.time.dt.year).sum()
    return out.rename(year="year")

def trend_per_decade(y, years):
    m = np.isfinite(y)
    return np.polyfit(np.asarray(years)[m], np.asarray(y)[m], 1)[0] * 10

def wmean(da, method="band", **kw):
    return gaw.area_mean(da, method=method, **kw)

def land_mask_for(lat, lon, lon360=True):
    """Boolean (nlat, nlon) land mask at cell centres from Natural Earth 110m land (any lat/lon order)."""
    import regionmask
    land = regionmask.defined_regions.natural_earth_v5_1_2.land_110
    lat = np.asarray(lat, float); lon = np.asarray(lon, float)
    lon_ = np.where(lon > 180, lon - 360, lon)
    io, jo = np.argsort(lat), np.argsort(lon_)
    m = land.mask(lon_[jo], lat[io]).notnull().values       # (nlat_sorted, nlon_sorted)
    out = np.empty_like(m)
    out[np.ix_(io, jo)] = m
    return out
