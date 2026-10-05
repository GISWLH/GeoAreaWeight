"""Extra analyses: latitude-window, resolution, band contribution, synthetic checks."""
import json
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (  # noqa: E402
    D, DOCS, annual, gaw, land_mask_for, np, op, trend_per_decade, xr,
)
E = {}
Dd = D + "/derived"

# ---- A. latitude-window dependence ------------------------------------------------------------
a = annual(op("air.mon.mean.nc").air).sel(year=slice(1980, 2020)); y = a.year.values
T = (a + 273.15)                     # Kelvin
lat = a.lat.values
Tc = T.mean("year")                  # climatology (lat, lon)
gp = annual(op("gpcp.precip.mon.mean.nc").precip).sel(year=slice(1980, 2020)) * 365.25
Pc = gp.mean("year")
lat_p = Pc.lat.values
Ta = a - a.mean("year")              # anomalies
rows = []
for L in range(10, 91, 10):
    r = dict(L=L)
    # symmetric window |lat| <= L
    for nm, da in [("T", Tc), ("P", Pc)]:
        la = da.lat.values; m = np.zeros(da.shape, bool); m[np.abs(la) <= L, :] = True
        u = float(gaw.area_mean(da.values, la, da.lon.values, method="none", mask=m))
        w = float(gaw.area_mean(da.values, la, da.lon.values, method="band", mask=m))
        r[f"{nm}_sym_unw"], r[f"{nm}_sym_w"] = u, w
    m = np.zeros(Ta.shape[1:], bool); m[np.abs(lat) <= L, :] = True
    su = gaw.area_mean(Ta.values, lat, Ta.lon.values, method="none", mask=m); sw = gaw.area_mean(Ta.values, lat, Ta.lon.values, method="band", mask=m)
    r["Ttrend_sym_unw"], r["Ttrend_sym_w"] = trend_per_decade(su, y), trend_per_decade(sw, y)
    rows.append(r)
# poleward windows lat >= L  (L = 0,30,60 ...)
rows2 = []
for L in (0, 30, 45, 60, 70, 80):
    r = dict(L=L)
    m = np.zeros(Ta.shape[1:], bool); m[lat >= L, :] = True
    su = gaw.area_mean(Ta.values, lat, Ta.lon.values, method="none", mask=m); sw = gaw.area_mean(Ta.values, lat, Ta.lon.values, method="band", mask=m)
    r["Ttrend_unw"], r["Ttrend_w"] = trend_per_decade(su, y), trend_per_decade(sw, y)
    mt = np.zeros(Tc.shape, bool); mt[lat >= L, :] = True
    r["T_unw"] = float(gaw.area_mean(Tc.values, lat, Tc.lon.values, method="none", mask=mt))
    r["T_w"] = float(gaw.area_mean(Tc.values, lat, Tc.lon.values, method="band", mask=mt))
    rows2.append(r)
E["lat_window_sym"] = rows; E["lat_window_pole"] = rows2

# ---- B. band contribution (10-degree bands), NCEP T climatology ------------------------------
edges = np.arange(-90, 91, 10)
w_unw = np.ones(Tc.shape) / Tc.size
w_area = gaw.weights_2d(lat, Tc.lon.values)
band = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (lat >= lo) & (lat < hi) if hi < 90 else (lat >= lo)
    mm = np.zeros(Tc.shape, bool); mm[m, :] = True
    band.append(dict(lo=int(lo), hi=int(hi), cell_share=float(w_unw[mm].sum()), area_share=float(w_area[mm].sum()),
                     T_mean_K=float(np.sum(Tc.values[mm] * w_area[mm]) / w_area[mm].sum()),
                     contrib_diff_K=float(np.sum(Tc.values[mm] * (w_unw[mm] - w_area[mm])))))
E["bands"] = band
E["bands_total_diff_K"] = float(sum(b["contrib_diff_K"] for b in band))

# ---- C. resolution dependence (CRU land T, GPCC land P) --------------------------------------
res_rows = []
def coarsen(da, f):
    if f == 1: return da
    return da.coarsen(lat=f, lon=f, boundary="trim").reduce(np.nanmean)
if os.path.exists(Dd + "/cru_tmp_annual.nc"):
    ct = xr.open_dataset(Dd + "/cru_tmp_annual.nc").tmp.sel(year=slice(1991, 2020)).mean("year")
    ctr = xr.open_dataset(Dd + "/cru_tmp_annual.nc").tmp.sel(year=slice(1980, 2020))
    for f in (1, 2, 4, 10, 20, 40):
        c = coarsen(ct, f); r = dict(res=0.5 * f)
        r["T_unw"] = float(gaw.area_mean(c.values, c.lat.values, c.lon.values, method="none"))
        r["T_w"] = float(gaw.area_mean(c.values, c.lat.values, c.lon.values, method="band"))
        cs = coarsen(ctr, f)
        su = gaw.area_mean(cs.values, cs.lat.values, cs.lon.values, method="none"); sw = gaw.area_mean(cs.values, cs.lat.values, cs.lon.values, method="band")
        r["trend_unw"], r["trend_w"] = trend_per_decade(su, cs.year.values), trend_per_decade(sw, cs.year.values)
        res_rows.append(r)
E["resolution_cru_T"] = res_rows
pres = []
gpc = op("precip.mon.total.1x1.v2020.nc").precip.sel(time=slice("1991-01-01", "2019-12-31"))
gpa = gpc.groupby("time.year").sum("time", skipna=False).mean("year")
for f in (1, 2, 3, 5, 10, 20):
    c = coarsen(gpa, f)
    pres.append(dict(res=1.0 * f, P_unw=float(gaw.area_mean(c.values, c.lat.values, c.lon.values, method="none")),
                     P_w=float(gaw.area_mean(c.values, c.lat.values, c.lon.values, method="band"))))
E["resolution_gpcc_P"] = pres

# NCEP global T climatology interpolated (bilinear) to finer / coarser grids: resolution dependence of the gap
Tcs = Tc.sortby("lat")
lonp = np.concatenate([Tcs.lon.values, [360.0]])
vals = np.concatenate([Tcs.values, Tcs.values[:, :1]], axis=1)
from scipy.interpolate import RegularGridInterpolator
rgi = RegularGridInterpolator((Tcs.lat.values, lonp), vals)
nres = []
for res in (10.0, 5.0, 2.5, 1.0, 0.5, 0.25):
    la = np.arange(-90 + res / 2, 90, res); lo = np.arange(res / 2, 360, res)
    LA, LO = np.meshgrid(la, lo, indexing="ij")
    x = rgi(np.column_stack([LA.ravel(), LO.ravel()])).reshape(LA.shape)
    nres.append(dict(res=res, T_unw=float(gaw.area_mean(x, la, lo, method="none")) - 273.15, T_w=float(gaw.area_mean(x, la, lo)) - 273.15))
E["resolution_ncep_interp_T"] = nres

# ---- D. synthetic fields ----------------------------------------------------------------------
lat1 = np.arange(-89.5, 90, 1.0); lon1 = np.arange(0.5, 360, 1.0)
syn = {}
f = lambda g: np.broadcast_to(g(lat1)[:, None], (lat1.size, lon1.size))
for nm, g in [("constant field", lambda l: np.ones_like(l) * 5.0), ("cos(lat)", lambda l: np.cos(np.deg2rad(l))),
              ("|lat| (linear meridional gradient)", lambda l: np.abs(l)), ("temperature-like profile 30cos²φ-10", lambda l: 30 * np.cos(np.deg2rad(l)) ** 2 - 10)]:
    x = f(g)
    syn[nm] = dict(unw=float(gaw.area_mean(x, lat1, lon1, method="none")), w=float(gaw.area_mean(x, lat1, lon1)))
# random (no latitudinal structure) fields: unweighted vs weighted mean difference ~ 0 on average
rng = np.random.default_rng(0)
d = []
for k in range(500):
    x = rng.normal(size=(lat1.size, lon1.size))
    d.append(gaw.area_mean(x, lat1, lon1, method="none") - gaw.area_mean(x, lat1, lon1))
syn["white-noise fields: unweighted - weighted, mean/std"] = dict(mean=float(np.mean(d)), std=float(np.std(d)))
# polar amplification: trend s(lat)=a+b*|lat|/90 degC/decade
tr = {}
for b in (0.0, 0.1, 0.2, 0.4, 0.8):
    s = 0.15 + b * np.abs(lat1) / 90
    x = np.broadcast_to(s[:, None], (lat1.size, lon1.size))
    tr[str(b)] = dict(unw=float(gaw.area_mean(x, lat1, lon1, method="none")), w=float(gaw.area_mean(x, lat1, lon1)))
syn["polar-amplified trend s=0.15+b|φ|/90"] = tr
E["synthetic"] = syn

# ---- F. normalised gap vs share of variance carried by latitude --------------------------------
def gapstat(x, la, lo, name):
    x = np.asarray(x, float)
    w2 = gaw.weights_2d(la, lo, normalize=False)
    v = np.isfinite(x)
    wv = np.where(v, w2, 0.0); x0 = np.where(v, x, 0.0)
    mw = (x0 * wv).sum() / wv.sum(); mu = x0[v].mean()
    var = ((x0 - mw) ** 2 * wv).sum() / wv.sum()
    zw = wv.sum(1); zm = np.where(zw > 0, (x0 * wv).sum(1) / np.where(zw > 0, zw, 1), np.nan)
    varz = np.nansum(zw * (zm - mw) ** 2) / zw.sum()
    return dict(name=name, mean_unw=float(mu), mean_w=float(mw), sigma=float(np.sqrt(var)), lat_var_frac=float(varz / var),
                gap_over_sigma=float((mu - mw) / np.sqrt(var)))
gs = []
gs.append(gapstat(Tc.values, lat, Tc.lon.values, "NCEP R1 air temperature (global)"))
ct_ = xr.open_dataset(Dd + "/cru_tmp_annual.nc").tmp.sel(year=slice(1980, 2020)).mean("year")
gs.append(gapstat(ct_.values, ct_.lat.values, ct_.lon.values, "CRU air temperature (land)"))
gs.append(gapstat(gpa.values if False else gpc.groupby("time.year").sum("time", skipna=False).mean("year").values, gpc.lat.values, gpc.lon.values, "GPCC precipitation (land)"))
gs.append(gapstat(Pc.values, Pc.lat.values, Pc.lon.values, "GPCP precipitation (global)"))
sp_ = xr.open_dataset(Dd + "/spei12_dec.nc").spei12.sel(year=slice(1981, 2022))
gs.append(gapstat(sp_.mean("year").values, sp_.lat.values, sp_.lon.values, "SPEI-12 time mean"))
gs.append(gapstat((sp_ <= -1).where(sp_.notnull()).astype(float).mean("year").values, sp_.lat.values, sp_.lon.values, "SPEI <= -1 frequency"))
E["gap_vs_latvar"] = gs

# ---- E. land share ----------------------------------------------------------------------------
lm = land_mask_for(lat1, lon1)
E["land_share_1deg"] = dict(cell_share=float(lm.mean()), area_share=float(gaw.area_mean(lm.astype(float), lat1, lon1)))
json.dump(E, open(DOCS / "extras.json", "w"), ensure_ascii=False, indent=1)
for k in ("lat_window_sym", "lat_window_pole", "resolution_cru_T", "resolution_gpcc_P", "land_share_1deg", "bands_total_diff_K"):
    print(k); print(json.dumps(E[k], indent=0)[:2500])
print(json.dumps(E["synthetic"], ensure_ascii=False, indent=0))
