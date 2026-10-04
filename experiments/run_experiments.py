"""Unweighted vs area-weighted means/trends on common open datasets (GeoAreaWeight)."""
import sys, json
sys.path.insert(0, "/workspace/geoareaweight/experiments")
from common import *

R = {"rows": [], "extra": {}}
S = {}   # series for figures

def add(name, kind, unit, period, grid, u, w, extra=""):
    """u, w: annual series (pandas-like arrays) of unweighted / weighted mean."""
    yrs = np.asarray(period)
    u = np.asarray(u, float); w = np.asarray(w, float)
    row = dict(name=name, kind=kind, unit=unit, period=f"{yrs[0]}-{yrs[-1]}", grid=grid,
               mean_unw=float(np.nanmean(u)), mean_w=float(np.nanmean(w)),
               trend_unw=float(trend_per_decade(u, yrs)), trend_w=float(trend_per_decade(w, yrs)), note=extra)
    row["mean_diff"] = row["mean_unw"] - row["mean_w"]
    row["mean_rel_pct"] = 100 * row["mean_diff"] / abs(row["mean_w"]) if abs(row["mean_w"]) > 1e-9 else np.nan
    row["trend_diff"] = row["trend_unw"] - row["trend_w"]
    row["trend_rel_pct"] = 100 * row["trend_diff"] / abs(row["trend_w"]) if abs(row["trend_w"]) > 1e-12 else np.nan
    R["rows"].append(row)
    S[name] = dict(years=yrs, unw=u, w=w)
    print(f"{name:55s} mean {row['mean_unw']:9.3f} vs {row['mean_w']:9.3f} ({row['mean_rel_pct']:+7.1f}%) | trend {row['trend_unw']:8.4f} vs {row['trend_w']:8.4f} ({row['trend_rel_pct']:+7.1f}%)", flush=True)

def pair(da, years, mask=None, **kw):
    u = wmean(da, "none", mask=mask, **kw).values; w = wmean(da, "band", mask=mask, **kw).values
    return u, w

# ---------------------------------------------------------------- 1. NCEP/NCAR R1 near-surface temperature (2.5 deg)
a = annual(op("air.mon.mean.nc").air).sel(year=slice(1980, 2020)); y = a.year.values
lm = op("land.nc").land.isel(time=0).values == 1
u, w = pair(a, y); add("NCEP R1 气温 全球(含海洋)", "T", "°C", y, "2.5°", u, w)
u, w = pair(a, y, mask=lm); add("NCEP R1 气温 仅陆地(含南极)", "T", "°C", y, "2.5°", u, w, "land.nc mask")
lat = a.lat.values
lm2 = lm & (lat[:, None] > -60)
u, w = pair(a, y, mask=lm2); add("NCEP R1 气温 仅陆地(不含南极)", "T", "°C", y, "2.5°", u, w, "land.nc, lat>-60")
S["ncep_air_global"] = S["NCEP R1 气温 全球(含海洋)"]
# Arctic / tropics windows
for nm, sel in [("60-90N", lat >= 60), ("30S-30N", np.abs(lat) <= 30), ("0-90N", lat >= 0)]:
    m = np.zeros(lm.shape, bool); m[sel, :] = True
    u, w = pair(a, y, mask=m); add(f"NCEP R1 气温 {nm}", "T", "°C", y, "2.5°", u, w)
# cos vs band vs ellipsoid for global mean (method comparison)
ab = {}
for m in ("none", "cos", "band", "ellipsoid"):
    s = wmean(a, m).values; ab[m] = dict(mean=float(s.mean()), trend=float(trend_per_decade(s, y)))
R["extra"]["method_compare_ncep_T"] = ab
print("method compare", ab)

# ---------------------------------------------------------------- 2. GISTEMP v4 anomalies (2 deg)
g = annual(op("gistemp1200_GHCNv4_ERSSTv5.nc").tempanomaly).sel(year=slice(1980, 2020)); y = g.year.values
u, w = pair(g, y); add("GISTEMP v4 温度异常 全球", "Tanom", "K", y, "2°", u, w, "NaN 处权重重新归一")
# 1880-2020 too
g2 = annual(op("gistemp1200_GHCNv4_ERSSTv5.nc").tempanomaly).sel(year=slice(1880, 2020)); y2 = g2.year.values
u2, w2 = pair(g2, y2); R["extra"]["gistemp_1880_2020"] = dict(trend_unw=trend_per_decade(u2, y2), trend_w=trend_per_decade(w2, y2))
S["gistemp_full"] = dict(years=y2, unw=u2, w=w2)

# ---------------------------------------------------------------- 3. GPCP v2.3 precip (2.5 deg)
p = op("gpcp.precip.mon.mean.nc").precip * 1.0
pa = annual(p).sel(year=slice(1980, 2020)) * 365.25; y = pa.year.values  # mm/yr
u, w = pair(pa, y); add("GPCP 降水 全球(含海洋)", "P", "mm/yr", y, "2.5°", u, w)
lmg = land_mask_for(pa.lat.values, pa.lon.values)
u, w = pair(pa, y, mask=lmg); add("GPCP 降水 仅陆地(含南极/格陵兰)", "P", "mm/yr", y, "2.5°", u, w, "Natural Earth land mask (cell centre)")
lmg2 = lmg & (pa.lat.values[:, None] > -60)
u, w = pair(pa, y, mask=lmg2); add("GPCP 降水 仅陆地(不含南极)", "P", "mm/yr", y, "2.5°", u, w)

# ---------------------------------------------------------------- 4. GPCC v2020 land precipitation (1 deg, monthly totals)
gp = op("precip.mon.total.1x1.v2020.nc").precip
gpa = gp.sel(time=slice("1980-01-01", "2019-12-31")).groupby("time.year").sum("time", skipna=False).where(
    gp.sel(time=slice("1980-01-01", "2019-12-31")).groupby("time.year").count("time") == 12)
gpa = gpa.where(gpa > -1)  # keep NaN
y = gpa.year.values
u, w = pair(gpa, y); add("GPCC v2020 降水 陆地(观测网格, 1°)", "P", "mm/yr", y, "1°", u, w, "GPCC 仅陆地, 其余为 NaN")
S["gpcc_valid_frac"] = float(gpa.notnull().isel(year=0).mean())
# ---------------------------------------------------------------- 5. NCEP R1 Gaussian grid latent heat (ET proxy) & precip, land only
lh = annual(op("lhtfl.mon.mean.nc").lhtfl).sel(year=slice(1980, 2020)); y = lh.year.values
lmgs = op("land.sfc.gauss.nc").land.isel(time=0).values == 1
et = lh * (365.25 * 86400 / 2.45e6)   # W/m2 -> mm/yr  (lambda = 2.45e6 J/kg)
u, w = pair(et, y, mask=lmgs); add("NCEP R1 潜热通量→蒸散发代用 陆地(T62 高斯网格)", "ET", "mm/yr", y, "T62 Gaussian", u, w, "反演代用量, 非观测")
lmgs2 = lmgs & (lh.lat.values[:, None] > -60)
u, w = pair(et, y, mask=lmgs2); add("NCEP R1 蒸散发代用 陆地(不含南极)", "ET", "mm/yr", y, "T62 Gaussian", u, w)
pr = annual(op("prate.mon.mean.nc").prate).sel(year=slice(1980, 2020)) * 365.25 * 86400
u, w = pair(pr, y, mask=lmgs2); add("NCEP R1 降水 陆地(不含南极)", "P", "mm/yr", y, "T62 Gaussian", u, w)
# Gaussian weights: exact quadrature vs midpoint-bands vs cos
latg = lh.lat.values
x_, wq = np.polynomial.legendre.leggauss(94)
lat_exact = np.rad2deg(np.arcsin(x_))
order = np.argsort(-lat_exact)
wq = wq[order] / wq.sum(); lat_exact = lat_exact[order]
wb = gaw.area_weights(latg, method="band"); wc = gaw.area_weights(latg, method="cos")
R["extra"]["gaussian"] = dict(max_lat_dev_deg=float(np.max(np.abs(lat_exact - latg))),
                              max_rel_err_band=float(np.max(np.abs(wb - wq) / wq)), max_rel_err_cos=float(np.max(np.abs(wc - wq) / wq)))
res = {}
for nm, wt in [("gauss_exact", wq), ("band", wb), ("cos", wc)]:
    W2 = np.repeat(wt[:, None], et.sizes["lon"], 1)
    res[nm] = float(np.mean(wmean(et, weights=W2, mask=lmgs2).values))
res["unweighted"] = float(np.mean(wmean(et, "none", mask=lmgs2).values))
R["extra"]["gaussian"]["land_ET_mean_by_weights"] = res
print("gaussian", R["extra"]["gaussian"])

# ---------------------------------------------------------------- 6. OISST v2 SST (ocean only) & sea ice
sst = annual(op("oisst.sst.mnmean.nc").sst).sel(year=slice(1982, 2022)); y = sst.year.values
u, w = pair(sst, y); add("OISST v2 海温 全海洋(1°)", "SST", "°C", y, "1°", u, w)
ice = op("oisst.icec.mnmean.nc").icec
ice_a = annual(ice).sel(year=slice(1982, 2022)); 
lsm = op("oisst.lsmask.nc").mask.isel(time=0).values == 1   # 1 = ocean
lat1 = ice_a.lat.values
mN = lsm & (lat1[:, None] >= 50)
u, w = pair(ice_a, y, mask=mN); add("OISST 海冰密集度 年均 50-90°N 海洋", "ICE", "%", y, "1°", u, w)
# NH sea-ice extent, sept: area of cells with conc>=15 (km2 weighted) vs 'cell count x mean cell area'
sep = ice.sel(time=ice.time.dt.month == 9); sep = sep.assign_coords(year=("time", sep.time.dt.year.values)).swap_dims(time="year").drop_vars("time").sel(year=slice(1982, 2022))
cond = (sep >= 15).where(sep.notnull())
fr_u = wmean(cond.astype(float), "none", mask=mN).values; fr_w = wmean(cond.astype(float), "band", mask=mN).values
add("OISST 9月北极海冰覆盖面积占比(>=15%, 50-90°N海洋)", "ICEFRAC", "fraction", sep.year.values, "1°", fr_u, fr_w)
A = gaw.cell_area(ice.lat.values, ice.lon.values, units="km2")
ext = (cond.fillna(0).values * (A * mN)[None]).sum((1, 2)) / 1e6
nai = (cond.fillna(0).values * mN[None]).sum((1, 2)) * A.mean() / 1e6   # naive: count cells x global-mean cell area
R["extra"]["arctic_sep_extent"] = dict(years=sep.year.values.tolist(), extent_area_weighted_Mkm2=ext.tolist(),
                                       extent_cellcount_x_meanarea_Mkm2=nai.tolist(),
                                       mean_ratio=float(np.mean(nai / ext)))
print("Arctic Sep extent mean (area-weighted):", ext.mean(), " naive count*mean-area:", nai.mean())

# ---------------------------------------------------------------- 7. CRU TS 4.09 land
import os
Dd = D + "/derived"
if os.path.exists(Dd + "/cru_tmp_annual.nc"):
    ct = xr.open_dataset(Dd + "/cru_tmp_annual.nc").tmp.sel(year=slice(1980, 2020)); y = ct.year.values
    u, w = pair(ct, y); add("CRU TS4.09 气温 陆地(0.5°, 无南极)", "T", "°C", y, "0.5°", u, w)
    cta = ct - ct.mean("year")
    S["cru_tmp_anom"] = dict(years=y, unw=pair(cta, y)[0], w=pair(cta, y)[1])
if os.path.exists(Dd + "/cru_pre_annual.nc"):
    cp = xr.open_dataset(Dd + "/cru_pre_annual.nc").pre.sel(year=slice(1980, 2020)); y = cp.year.values
    u, w = pair(cp, y); add("CRU TS4.09 降水 陆地(0.5°, 无南极)", "P", "mm/yr", y, "0.5°", u, w)
    cp_full = xr.open_dataset(Dd + "/cru_pre_annual.nc").pre
    arid = (cp_full.sel(year=slice(1981, 2010)).mean("year") < 180)   # <180 mm/yr: the paper's arid mask
    R["extra"]["arid_area_fraction_land"] = float(wmean(arid.astype(float).where(cp_full.isel(year=0).notnull()), "band"))
    R["extra"]["arid_cell_fraction_land"] = float(wmean(arid.astype(float).where(cp_full.isel(year=0).notnull()), "none"))

# ---------------------------------------------------------------- 8. SPEI-12 drought-affected area (the retracted drought paper's statistic)
if os.path.exists(Dd + "/spei12_dec.nc"):
    sp = xr.open_dataset(Dd + "/spei12_dec.nc").spei12.sel(year=slice(1981, 2022)); y = sp.year.values
    valid = sp.isel(year=0).notnull()
    u, w = pair(sp, y); add("SPEI-12 (12月) 陆地平均值", "SPEI", "index", y, "0.5°", u, w)
    for thr, nm in [(-1.0, "SPEI<=-1 (中度及以上干旱)"), (-1.5, "SPEI<=-1.5 (严重干旱)")]:
        c = (sp <= thr).where(sp.notnull()).astype(float)
        u, w = pair(c, y); add(f"SPEI-12 {nm} 面积占比", "DAREA", "fraction", y, "0.5°", u, w)
        if thr == -1.0: S["spei_area"] = dict(years=y, unw=u, w=w)
    if "arid_area_fraction_land" in R["extra"]:
        arid_m = xr.open_dataset(Dd + "/cru_pre_annual.nc").pre.sel(year=slice(1981, 2010)).mean("year") < 180
        # align (CRU lat/lon same grid as SPEI? check)
        try:
            arid_m = arid_m.astype(float).interp(lat=sp.lat, lon=sp.lon, method="nearest") > 0.5
            c = (sp <= -1).where(sp.notnull() & ~arid_m).astype(float)
            u, w = pair(c, y); add("SPEI-12<=-1 面积占比 (排除年降水<180 mm 干旱区)", "DAREA", "fraction", y, "0.5°", u, w)
        except Exception as e:
            print("arid mask failed", e)

json.dump(R, open("/workspace/geoareaweight/results/results.json", "w"), ensure_ascii=False, indent=1, default=float)
def _j(o):
    if isinstance(o, dict): return {k: _j(v) for k, v in o.items()}
    if isinstance(o, np.ndarray): return o.tolist()
    return o
json.dump(_j(S), open("/workspace/geoareaweight/results/series.json", "w"), ensure_ascii=False)
print("saved")
