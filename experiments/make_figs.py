import sys, json
sys.path.insert(0, "/workspace/geoareaweight/experiments")
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs, cartopy.feature as cfeature
import geoareaweight as gaw

for f in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"):
    fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "Noto Sans CJK SC", "axes.unicode_minus": False, "font.size": 10.5,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "figure.dpi": 120})
OUT = "/workspace/geoareaweight/figs/"
RES = json.load(open("/workspace/geoareaweight/results/results.json"))
SER = json.load(open("/workspace/geoareaweight/results/series.json"))
EXT = json.load(open("/workspace/geoareaweight/results/extras.json"))
rows = {r["name"]: r for r in RES["rows"]}
C_U, C_W = "#d1495b", "#1f6f8b"   # unweighted red, weighted blue

# ======================================================================= Fig 1: concept
fig = plt.figure(figsize=(13, 4.6))
ax = fig.add_subplot(1, 3, 1)
cm = plt.get_cmap("YlGnBu")
for i, lat0 in enumerate(range(-90, 90, 30)):
    for j, lon0 in enumerate(range(-180, 180, 30)):
        a = np.cos(np.deg2rad(lat0 + 15)) * 1.0
        ax.add_patch(Rectangle((lon0, lat0), 30, 30, fc=cm(0.15 + 0.8 * a), ec="white", lw=1))
ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_aspect("equal")
ax.set_xlabel("经度 (°)"); ax.set_ylabel("纬度 (°)")
ax.set_title("(a) 经纬网格上每个格子“看起来一样大”\n（颜色 = 实际面积，深色更大）", fontsize=10.5)
ax.text(0, 75, "每格 30°×30°，权重 1:1:1…", ha="center", color="k", fontsize=9)
# (b) orthographic globe with one column of cells highlighted
axb = fig.add_subplot(1, 3, 2, projection=ccrs.Orthographic(30, 25))
axb.set_global(); axb.add_feature(cfeature.LAND, fc="#e8e3d3", ec="none"); axb.add_feature(cfeature.OCEAN, fc="#f4f8fb")
axb.gridlines(draw_labels=False, linewidth=0.3, color="gray", xlocs=range(-180, 181, 10), ylocs=range(-90, 91, 10))
for lat0 in range(-90, 90, 10):
    a = np.cos(np.deg2rad(lat0 + 5))
    lons = [20, 30, 30, 20, 20]; lats = [lat0, lat0, lat0 + 10, lat0 + 10, lat0]
    axb.fill(lons, lats, transform=ccrs.Geodetic(), fc=cm(0.15 + 0.8 * a), ec="k", lw=0.5, alpha=0.95)
axb.set_title("(b) 同一经度条带上的 10°×10° 格子：\n越靠近极地，真实面积越小", fontsize=10.5)
# (c) area vs latitude
axc = fig.add_subplot(1, 3, 3)
phi = np.linspace(0, 90, 400)
axc.plot(phi, np.cos(np.deg2rad(phi)), color=C_W, lw=2.2, label="相对面积 ∝ cos φ")
for res, mk in [(10, "o"), (30, "s")]:
    c = np.arange(res / 2, 90, res)
    w = gaw.band_weights(np.concatenate([-c[::-1], c]), normalize=False)
    w = w[len(c):]; w = w / w[0] * np.cos(np.deg2rad(c[0]))
    axc.plot(c, w, mk, color="k", ms=5, mfc="none", label=f"{res}° 格子：精确带面积 sinφ北−sinφ南")
for p in (30, 60, 80):
    axc.plot([p, p], [0, np.cos(np.deg2rad(p))], color="gray", lw=0.8, ls=":")
    axc.annotate(f"{p}°: {np.cos(np.deg2rad(p)):.2f}", (p, np.cos(np.deg2rad(p))), textcoords="offset points", xytext=(5, 6), fontsize=9)
axc.axhline(1, color="gray", lw=0.6)
axc.set_xlabel("纬度 |φ| (°)"); axc.set_ylabel("格子面积 / 赤道格子面积")
axc.set_title("(c) 格子面积随纬度的变化\n普通 mean 却给每个格子同样的权重", fontsize=10.5)
axc.legend(frameon=False, fontsize=8.8, loc="lower left"); axc.set_xlim(0, 90); axc.set_ylim(0, 1.08)
fig.tight_layout(); fig.savefig(OUT + "fig1_concept_cell_area_vs_lat.png", bbox_inches="tight"); plt.close(fig)

# ======================================================================= Fig 2: Cartopy grid maps
fig = plt.figure(figsize=(14, 8.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.12, wspace=0.12)
# (a) PlateCarree 0.25 deg colored by relative area
lat25 = np.arange(-90, 90.01, 0.25); lon25 = np.arange(-180, 180.01, 0.25)
la_c = 0.5 * (lat25[1:] + lat25[:-1])
rel = gaw.band_weights(la_c, lat_bounds=lat25, normalize=False); rel = rel / rel.max()
Z = np.repeat(rel[:, None], lon25.size - 1, 1)
ax1 = fig.add_subplot(gs[0, 0:2], projection=ccrs.PlateCarree())
m1 = ax1.pcolormesh(lon25, lat25, Z, cmap="viridis", vmin=0, vmax=1, transform=ccrs.PlateCarree(), rasterized=True)
ax1.coastlines(lw=0.5, color="w"); ax1.set_global()
ax1.gridlines(draw_labels=True, linewidth=0.3, color="w", alpha=0.6, xlocs=range(-180, 181, 60), ylocs=range(-90, 91, 30))
ax1.set_title("(a) 0.25°×0.25° 规则网格（1440×720）：每个格子的相对面积（等经纬度投影下“格子一样大”，实际不然）", fontsize=11)
cb = fig.colorbar(m1, ax=ax1, orientation="horizontal", pad=0.08, fraction=0.05, shrink=0.8); cb.set_label("相对面积 (赤道 = 1) = (sinφ北 − sinφ南) 归一化")
# (b) polar stereographic zoom with 2.5deg cells
ax2 = fig.add_subplot(gs[0, 2], projection=ccrs.NorthPolarStereo())
ax2.set_extent([-180, 180, 55, 90], ccrs.PlateCarree())
lat5 = np.arange(55, 90.01, 2.5); lon5 = np.arange(-180, 180.01, 2.5)
for la in lat5[:-1]:
    pass
Z5 = np.repeat((np.sin(np.deg2rad(lat5[1:])) - np.sin(np.deg2rad(lat5[:-1])))[:, None], lon5.size - 1, 1)
Z5 = Z5 / (np.sin(np.deg2rad(2.5)))
m2 = ax2.pcolormesh(lon5, lat5, Z5, cmap="viridis", vmin=0, vmax=1, edgecolor="w", linewidth=0.15, transform=ccrs.PlateCarree(), rasterized=True)
ax2.coastlines(lw=0.5, color="w")
ax2.set_title("(b) 北极放大：2.5°×2.5° 格子（55–90°N）\n格子数量相同、面积急剧缩小", fontsize=10.5)
# (c) Robinson with 5deg cells outlines
ax3 = fig.add_subplot(gs[1, 0:2], projection=ccrs.Robinson())
lat1 = np.arange(-90, 90.01, 1.0); lon1 = np.arange(-180, 180.01, 1.0)
la1 = 0.5 * (lat1[1:] + lat1[:-1]); r1 = gaw.band_weights(la1, lat_bounds=lat1, normalize=False); r1 = r1 / r1.max()
ax3.set_global()
ax3.pcolormesh(lon1, lat1, np.repeat(r1[:, None], lon1.size - 1, 1), cmap="viridis", vmin=0, vmax=1, transform=ccrs.PlateCarree(), rasterized=True)
ax3.coastlines(lw=0.5, color="w")
ax3.gridlines(draw_labels=False, linewidth=0.5, color="k", alpha=0.55, xlocs=range(-180, 181, 15), ylocs=range(-90, 91, 15))
ax3.set_title("(c) Robinson 投影 + 15° 经纬网：网格在投影下向两极收拢（1° 格子相对面积着色）", fontsize=11)
# (d) share of cells vs area per 10 deg band (bar)
ax4 = fig.add_subplot(gs[1, 2])
b = EXT["bands"]; ctr = np.array([(x["lo"] + x["hi"]) / 2 for x in b])
# use ideal 1deg grid shares for the concept (not NCEP's pole rows)
lat_g = np.arange(-89.5, 90, 1.0); wa = gaw.area_weights(lat_g); 
edges = np.arange(-90, 91, 10)
cs = [np.mean((lat_g >= lo) & (lat_g < hi)) for lo, hi in zip(edges[:-1], edges[1:])]
as_ = [wa[(lat_g >= lo) & (lat_g < hi)].sum() for lo, hi in zip(edges[:-1], edges[1:])]
ax4.barh(ctr - 1.9, np.array(cs) * 100, height=3.6, color=C_U, label="格子数占比（普通 mean 的权重）")
ax4.barh(ctr + 1.9, np.array(as_) * 100, height=3.6, color=C_W, label="真实面积占比（面积加权）")
ax4.set_ylabel("纬度带 (°)"); ax4.set_xlabel("占比 (%)"); ax4.set_title("(d) 每 10° 纬度带：格子数 vs 面积", fontsize=10.5)
ax4.legend(frameon=False, fontsize=8.5, loc="lower right")
ax4.set_yticks(np.arange(-90, 91, 30))
fig.savefig(OUT + "fig2_cartopy_grid_area.png", bbox_inches="tight"); plt.close(fig)

# ======================================================================= Fig 3: weighted vs unweighted per dataset
def g(n): return rows[n]
tmean = [("NCEP R1 气温 全球(含海洋)", "NCEP R1\n全球 2.5°"), ("NCEP R1 气温 仅陆地(不含南极)", "NCEP R1\n陆地 2.5°"),
         ("CRU TS4.09 气温 陆地(0.5°, 无南极)", "CRU TS\n陆地 0.5°"), ("OISST v2 海温 全海洋(1°)", "OISST\n海温 1°"),
         ("NCEP R1 气温 60-90N", "NCEP R1\n60–90°N"), ("NCEP R1 气温 30S-30N", "NCEP R1\n30°S–30°N")]
pmean = [("GPCP 降水 全球(含海洋)", "GPCP 全球\n降水"), ("GPCC v2020 降水 陆地(观测网格, 1°)", "GPCC 陆地\n降水"),
         ("CRU TS4.09 降水 陆地(0.5°, 无南极)", "CRU TS 陆地\n降水"), ("NCEP R1 蒸散发代用 陆地(不含南极)", "NCEP R1 陆地\n蒸散发代用"),
         ("OISST 海冰密集度 年均 50-90°N 海洋", "OISST 50–90°N\n海冰密集度"),
         ("SPEI-12 SPEI<=-1 (中度及以上干旱) 面积占比", "SPEI-12≤−1\n干旱面积占比")]
trs = [("NCEP R1 气温 全球(含海洋)", "NCEP R1\n全球"), ("GISTEMP v4 温度异常 全球", "GISTEMP\n全球"), ("CRU TS4.09 气温 陆地(0.5°, 无南极)", "CRU TS\n陆地"),
       ("OISST v2 海温 全海洋(1°)", "OISST\n海洋"), ("NCEP R1 气温 60-90N", "NCEP R1\n60–90°N"), ("NCEP R1 气温 30S-30N", "NCEP R1\n30°S–30°N")]
fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
ax = axs[0]
vals = [g(n)["mean_diff"] for n, _ in tmean if n in rows]; labs = [l for n, l in tmean if n in rows]
bars = ax.bar(range(len(vals)), vals, color=[C_U if v < 0 else C_W for v in vals])
for i, (v, n) in enumerate([(g(n)["mean_diff"], n) for n, _ in tmean if n in rows]):
    ax.text(i, v + (0.2 if v >= 0 else -0.2), f"{v:+.1f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=9)
ax.axhline(0, color="k", lw=0.8); ax.set_xticks(range(len(vals))); ax.set_xticklabels(labs, fontsize=8.3)
ax.set_ylabel("未加权 − 面积加权 (°C)"); ax.set_title("(a) 气温/海温均值：绝对差 (°C)", fontsize=11)
ax.set_ylim(min(vals) - 1.3, 1.5)
ax = axs[1]
vals = [g(n)["mean_rel_pct"] for n, _ in pmean if n in rows]; labs = [l for n, l in pmean if n in rows]
ax.bar(range(len(vals)), vals, color=[C_U if v < 0 else C_W for v in vals])
for i, v in enumerate(vals):
    ax.text(i, v + (3 if v >= 0 else -3), f"{v:+.0f}%" if abs(v) >= 10 else f"{v:+.1f}%", ha="center", va="bottom" if v >= 0 else "top", fontsize=9)
ax.axhline(0, color="k", lw=0.8); ax.set_xticks(range(len(vals))); ax.set_xticklabels(labs, fontsize=8.0)
ax.set_ylabel("(未加权 − 加权) / 加权 (%)"); ax.set_title("(b) 降水/蒸散发/海冰/干旱面积：相对差 (%)", fontsize=11)
ax.set_ylim(min(vals) - 22, max(vals) + 25)
ax = axs[2]
x = np.arange(len(trs)); wd = 0.38
tu = [g(n)["trend_unw"] for n, _ in trs]; tw = [g(n)["trend_w"] for n, _ in trs]
ax.bar(x - wd / 2, tu, wd, color=C_U, label="未加权"); ax.bar(x + wd / 2, tw, wd, color=C_W, label="面积加权")
for i in range(len(trs)):
    rel = 100 * (tu[i] - tw[i]) / abs(tw[i]); ax.text(i, max(tu[i], tw[i]) + 0.03, f"{rel:+.0f}%", ha="center", fontsize=9)
ax.set_xticks(x); ax.set_xticklabels([l for _, l in trs], fontsize=8.3); ax.set_ylabel("线性趋势 (°C/10年)")
ax.set_title("(c) 趋势（1980–2020；SST 1982–2022）", fontsize=11); ax.legend(frameon=False, fontsize=9, loc="upper right"); ax.set_ylim(0, 1.15)
fig.tight_layout(); fig.savefig(OUT + "fig3_weighted_vs_unweighted_by_dataset.png", bbox_inches="tight"); plt.close(fig)

# ======================================================================= Fig 4: time series + trend difference
fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.2))
def anom(s, base=(1981, 2010)):
    yrs = np.array(s["years"]); m = (yrs >= base[0]) & (yrs <= base[1]); return yrs, np.array(s["unw"]) - np.mean(np.array(s["unw"])[m]), np.array(s["w"]) - np.mean(np.array(s["w"])[m])
def tline(ax, yrs, y, c, ls="-"):
    p = np.polyfit(yrs, y, 1); ax.plot(yrs, np.polyval(p, yrs), color=c, lw=1, ls="--"); return p[0] * 10
ax = axs[0, 0]; s = SER["ncep_air_global"]; yrs, u, w = anom(s)
ax.plot(yrs, u, color=C_U, lw=1.6, label="未加权"); ax.plot(yrs, w, color=C_W, lw=1.6, label="面积加权")
tu, tw = tline(ax, yrs, u, C_U), tline(ax, yrs, w, C_W)
ax.set_title(f"(a) NCEP/NCAR R1 近地面气温 全球平均距平（相对1981–2010）\n趋势：未加权 {tu:.2f}  vs  加权 {tw:.2f} °C/10年（高估 {100*(tu-tw)/tw:.0f}%）", fontsize=10.5)
ax.set_ylabel("距平 (°C)"); ax.legend(frameon=False)
ax = axs[0, 1]; s = SER["gistemp_full"]; yrs, u, w = anom(s)
ax.plot(yrs, u, color=C_U, lw=1.4, label="未加权"); ax.plot(yrs, w, color=C_W, lw=1.4, label="面积加权")
m = yrs >= 1980; tu = np.polyfit(yrs[m], u[m], 1)[0] * 10; tw = np.polyfit(yrs[m], w[m], 1)[0] * 10
ax.set_title(f"(b) GISTEMP v4 (2°) 全球平均温度距平 1880–2020\n1980–2020 趋势：未加权 {tu:.2f} vs 加权 {tw:.2f} °C/10年（+{100*(tu-tw)/tw:.0f}%）", fontsize=10.5)
ax.set_ylabel("距平 (°C，相对1981–2010)"); ax.legend(frameon=False)
ax = axs[1, 0]
s = SER["spei_area"]; yrs = np.array(s["years"]); u = 100 * np.array(s["unw"]); w = 100 * np.array(s["w"])
ax.plot(yrs, u, color=C_U, lw=1.6, label="未加权"); ax.plot(yrs, w, color=C_W, lw=1.6, label="面积加权")
tu, tw = np.polyfit(yrs, u, 1)[0] * 10, np.polyfit(yrs, w, 1)[0] * 10
ax.set_title(f"(c) SPEI-12 ≤ −1 的陆地面积占比（CSIC SPEIbase v2.10, 0.5°）\n均值 {u.mean():.1f}% vs {w.mean():.1f}%；趋势 {tu:.2f} vs {tw:.2f} 个百分点/10年——差别很小", fontsize=10.5)
ax.set_ylabel("面积占比 (%)"); ax.legend(frameon=False)
ax = axs[1, 1]
s = SER["GPCC v2020 降水 陆地(观测网格, 1°)"] if "GPCC v2020 降水 陆地(观测网格, 1°)" in SER else None
yrs = np.array(s["years"]); u = np.array(s["unw"]); w = np.array(s["w"])
ax.plot(yrs, u, color=C_U, lw=1.6, label="未加权"); ax.plot(yrs, w, color=C_W, lw=1.6, label="面积加权")
tu, tw = np.polyfit(yrs, u, 1)[0] * 10, np.polyfit(yrs, w, 1)[0] * 10
ax.set_title(f"(d) GPCC v2020 陆地年降水量（1°）\n均值 {u.mean():.0f} vs {w.mean():.0f} mm/yr（低 {100*(w.mean()-u.mean())/w.mean():.0f}%）；趋势 {tu:.1f} vs {tw:.1f} mm/yr/10年", fontsize=10.5)
ax.set_ylabel("mm/yr"); ax.legend(frameon=False)
for a in axs.ravel(): a.set_xlabel("年")
fig.tight_layout(); fig.savefig(OUT + "fig4_timeseries_trend_difference.png", bbox_inches="tight"); plt.close(fig)

# ======================================================================= Fig 5: dependence on latitude / resolution / gradient
fig, axs = plt.subplots(2, 3, figsize=(15, 8.4))
ax = axs[0, 0]; ctr = np.array([(x["lo"] + x["hi"]) / 2 for x in EXT["bands"]])
ax.bar(ctr, [x["contrib_diff_K"] for x in EXT["bands"]], width=8.5, color=[C_U if x["contrib_diff_K"] < 0 else C_W for x in EXT["bands"]])
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("纬度带中心 (°)"); ax.set_ylabel("对 (未加权−加权) 的贡献 (K)")
ax.set_title(f"(a) NCEP R1 气温：各 10° 纬度带对偏差的贡献\n总偏差 {EXT['bands_total_diff_K']:.1f} K（高纬冷、被“多数”→ 偏冷）", fontsize=10.5)
ax = axs[0, 1]
L = [r["L"] for r in EXT["lat_window_sym"]]
dT = [r["T_sym_unw"] - r["T_sym_w"] for r in EXT["lat_window_sym"]]
dP = [100 * (r["P_sym_unw"] - r["P_sym_w"]) / r["P_sym_w"] for r in EXT["lat_window_sym"]]
ax.plot(L, dT, "o-", color=C_U, label="气温（左轴，K）"); ax.set_xlabel("纬度窗口 |φ| ≤ L (°)"); ax.set_ylabel("未加权 − 加权 (K)", color=C_U)
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
ax2.plot(L, dP, "s--", color="#e09f3e", label="降水 (右轴，%)"); ax2.set_ylabel("降水相对差 (%)", color="#e09f3e")
ax.set_title("(b) 偏差随纬度窗口的扩大而增大\n低纬窗口（|φ|≤30°）几乎没有差别", fontsize=10.5)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, frameon=False, loc="lower left", fontsize=9)
ax = axs[0, 2]
Lp = [r["L"] for r in EXT["lat_window_pole"]]
ax.plot(Lp, [r["Ttrend_unw"] for r in EXT["lat_window_pole"]], "o-", color=C_U, label="未加权")
ax.plot(Lp, [r["Ttrend_w"] for r in EXT["lat_window_pole"]], "o-", color=C_W, label="面积加权")
ax.set_xlabel("窗口 φ ≥ L（L→90°N）"); ax.set_ylabel("气温趋势 (°C/10年)")
ax.set_title("(c) 区域越偏向高纬，趋势差越小（窗口缩小到北极时各点权重趋同）\n而半球/全球尺度差别大", fontsize=10.5); ax.legend(frameon=False)
ax = axs[1, 0]
rc = EXT["resolution_cru_T"]; rg = EXT["resolution_gpcc_P"]; rn = EXT["resolution_ncep_interp_T"]
ax.plot([r["res"] for r in rc], [100 * (r["T_unw"] - r["T_w"]) / abs(r["T_w"]) for r in rc], "o-", color=C_U, label="CRU TS 陆地气温 (℃为基准)")
ax.plot([r["res"] for r in rg], [100 * (r["P_unw"] - r["P_w"]) / r["P_w"] for r in rg], "s-", color="#e09f3e", label="GPCC 陆地降水")
ax.plot([r["res"] for r in rn], [100 * (r["T_unw"] - r["T_w"]) / abs(r["T_w"]) for r in rn], "^-", color="gray", label="NCEP R1 全球气温（插值到不同分辨率）")
ax.set_xscale("log"); ax.set_xlabel("网格分辨率 (°)"); ax.set_ylabel("(未加权−加权)/加权 (%)")
ax.set_title("(d) 分辨率不是“解药”：网格变细，偏差并不消失", fontsize=10.5); ax.legend(frameon=False, fontsize=8.5)
ax = axs[1, 1]
bs = [float(k) for k in EXT["synthetic"]["极地放大型趋势 s=0.15+b|φ|/90"]]; syn = EXT["synthetic"]["极地放大型趋势 s=0.15+b|φ|/90"]
ratio = [100 * (syn[str(k)]["unw"] - syn[str(k)]["w"]) / syn[str(k)]["w"] for k in bs]
ax.plot(bs, ratio, "o-", color=C_U); ax.set_xlabel("极地放大强度 b（高纬比赤道多增暖 b °C/10年）"); ax.set_ylabel("未加权趋势高估 (%)")
ax.set_title("(e) 理想实验：趋势 s(φ)=0.15+b·|φ|/90\n趋势空间均匀（b=0）时加权与否无差别", fontsize=10.5)
for k, r in zip(bs, ratio): ax.annotate(f"{r:.0f}%", (k, r), textcoords="offset points", xytext=(4, -12), fontsize=8.5)
ax = axs[1, 2]
G = EXT["gap_vs_latvar"]
xs = [100 * x["lat_var_frac"] for x in G]; ys = [abs(x["gap_over_sigma"]) for x in G]
ax.scatter(xs, ys, s=70, color=C_U, zorder=3)
for x, y, d in zip(xs, ys, G):
    ax.annotate(d["name"], (x, y), textcoords="offset points", xytext=(6, -3 if d["name"] != "GPCP 降水(全球)" else 8), fontsize=8.5)
ax.set_xlabel("纬向平均（仅随纬度变化的部分）解释的空间方差 (%)"); ax.set_ylabel("|未加权均值 − 加权均值| / 空间标准差")
ax.set_title("(f) 变量的纬向结构越强，偏差越大\n（SPEI 已标准化，几乎无纬向结构）", fontsize=10.5)
ax.set_xlim(-3, 100); ax.set_ylim(-0.02, max(ys) * 1.15)
fig.tight_layout(); fig.savefig(OUT + "fig5_difference_vs_latitude_resolution_gradient.png", bbox_inches="tight"); plt.close(fig)
print("figs done")
