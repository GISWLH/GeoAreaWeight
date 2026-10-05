"""Prepare compact annual fields from the large CRU TS and SPEIbase NetCDF files (run once)."""
import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import D  # noqa: E402

O = D + "/derived"
os.makedirs(O, exist_ok=True)

def cru_annual(var, fn):
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
        if y % 20 == 0:
            print(var, y, flush=True)
    xr.concat(res, "year").astype("float32").to_dataset(name=var).to_netcdf(out)

cru_annual("tmp", "cru_ts4.09.1901.2024.tmp.dat.nc")
if os.path.exists(f"{D}/cru_ts4.09.1901.2024.pre.dat.nc"):
    cru_annual("pre", "cru_ts4.09.1901.2024.pre.dat.nc")

# SPEI-12: the December value summarises the calendar year
out = O + "/spei12_dec.nc"
if not os.path.exists(out):
    s = xr.open_dataset(f"{D}/spei12.nc")["spei"]
    d = s.sel(time=s.time.dt.month == 12)
    d = d.assign_coords(year=("time", d.time.dt.year.values)).swap_dims(time="year").drop_vars("time")
    d.astype("float32").to_dataset(name="spei12").to_netcdf(out)
print("prep done")
