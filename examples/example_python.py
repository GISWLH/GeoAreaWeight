"""Unweighted vs area-weighted global mean of NCEP/NCAR R1 near-surface air temperature."""
import sys
import numpy as np, xarray as xr
import geoareaweight as gaw

path = sys.argv[1] if len(sys.argv) > 1 else "air.mon.mean.nc"   # https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.derived.surface.html
da = xr.open_dataset(path).air.sel(time=slice("1980", "2020")).groupby("time.year").mean("time")
wrong = da.mean(("lat", "lon"))                                   # plain mean  <-- the mistake
right = gaw.area_mean(da)                                         # area weighted (exact band areas)
cosw = gaw.area_mean(da, method="cos")
tr = lambda s: np.polyfit(s.year, s, 1)[0] * 10
print(f"mean  : unweighted {wrong.mean():.2f}  area-weighted {right.mean():.2f}  cos {cosw.mean():.2f} degC")
print(f"trend : unweighted {tr(wrong):.3f}  area-weighted {tr(right):.3f} degC/decade")
