"""Prepare compact annual fields from the big NetCDF files (run once)."""
import os, warnings, numpy as np, xarray as xr
warnings.filterwarnings("ignore")
D = "/workspace/geoareaweight/data"; O = D + "/derived"

def cru_annual(var, fn, how):
    out = O + f"/cru_{var}_annual.nc"
    if os.path.exists(out): return
    ds = xr.open_dataset(f"{D}/{fn}", decode_times=True)[var]
    yrs = np.unique(ds.time.dt.year.values)
    res = []
    for y in yrs:
        a = ds.sel(time=slice(f"{y}-01-01", f"{y}-12-31")).load()
        if var == "pre":   # monthly totals -> annual total (mm/yr)
            r = a.sum("time", skipna=False)
        else:
            r = a.mean("time", skipna=False)
        res.append(r.expand_dims(year=[y]))
        print(var, y, flush=True) if y % 20 == 0 else None
    xr.concat(res, "year").astype("float32").to_dataset(name=var).to_netcdf(out)

cru_annual("tmp", "cru_ts4.09.1901.2024.tmp.dat.nc", "mean")
if os.path.exists(f"{D}/cru_ts4.09.1901.2024.pre.dat.nc"):
    cru_annual("pre", "cru_ts4.09.1901.2024.pre.dat.nc", "sum")

# SPEI-12: December value = SPEI of calendar year
out = O + "/spei12_dec.nc"
if not os.path.exists(out):
    s = xr.open_dataset(f"{D}/spei12.nc")["spei"]
    d = s.sel(time=s.time.dt.month == 12)
    d = d.assign_coords(year=("time", d.time.dt.year.values)).swap_dims(time="year").drop_vars("time")
    d.astype("float32").to_dataset(name="spei12").to_netcdf(out)
print("prep done")
