"""Unweighted vs area-weighted global mean of NCEP/NCAR R1 near-surface air temperature.

Data: air.mon.mean.nc from
https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.derived.surface.html
Usage: python example_python.py [path/to/air.mon.mean.nc]
"""
import sys

import numpy as np
import xarray as xr

import geoareaweight as gaw

path = sys.argv[1] if len(sys.argv) > 1 else "air.mon.mean.nc"
da = xr.open_dataset(path)["air"].sel(time=slice("1980", "2020")).groupby("time.year").mean("time")

wrong = da.mean(("lat", "lon"))               # plain mean: the mistake
right = gaw.area_mean(da)                     # exact spherical cell areas (default method="band")
cosw = gaw.area_mean(da, method="cos")        # textbook cos(lat) weights


def trend(s):
    return np.polyfit(s.year, s, 1)[0] * 10   # per decade


print(f"mean : unweighted {float(wrong.mean()):.2f}  area-weighted {float(right.mean()):.2f}  "
      f"cos {float(cosw.mean()):.2f} degC")
print(f"trend: unweighted {trend(wrong):.3f}  area-weighted {trend(right):.3f} degC/decade")

# share of the globe warmer than 20 degC in the last year (area fraction)
print(f"area > 20 degC: {float(gaw.area_fraction(da.isel(year=-1) > 20)):.1%}")
