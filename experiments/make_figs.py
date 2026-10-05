"""Draw figs/fig1-fig5 from docs/results.json, docs/extras.json and docs/series.json.

Figures 1-3 and 5 need only the JSON files committed in docs/. Figure 4 (time
series) needs docs/series.json, which run_experiments.py writes; it is skipped
if that file is missing. Requires matplotlib and cartopy (Natural Earth
coastlines are downloaded on first use).
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import cartopy.crs as ccrs  # noqa: E402
import cartopy.feature as cfeature  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import geoareaweight as gaw  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DOCS, OUT = ROOT / "docs", ROOT / "figs"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "figure.dpi": 120})
RES = json.load(open(DOCS / "results.json"))
EXT = json.load(open(DOCS / "extras.json"))
SER = json.load(open(DOCS / "series.json")) if (DOCS / "series.json").exists() else None
rows = {r["name"]: r for r in RES["rows"]}
C_U, C_W, C_P = "#d1495b", "#1f6f8b", "#e09f3e"   # unweighted red, weighted blue, precipitation orange
L_U, L_W = "Unweighted", "Area-weighted"


def save(fig, name):
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)


# ======================================================================= Fig 1: concept
fig = plt.figure(figsize=(13, 4.6))
ax = fig.add_subplot(1, 3, 1)
cm = plt.get_cmap("YlGnBu")
for lat0 in range(-90, 90, 30):
    for lon0 in range(-180, 180, 30):
        a = np.cos(np.deg2rad(lat0 + 15))
        ax.add_patch(Rectangle((lon0, lat0), 30, 30, fc=cm(0.15 + 0.8 * a), ec="white", lw=1))
ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_aspect("equal")
ax.set_xlabel("Longitude (°)"); ax.set_ylabel("Latitude (°)")
ax.set_title("(a) On a lat-lon grid every cell 'looks' equal\n(colour = true area, darker = larger)", fontsize=10.5)
ax.text(0, 75, "30°×30° cells: plain mean weights 1:1:1...", ha="center", color="k", fontsize=9)
axb = fig.add_subplot(1, 3, 2, projection=ccrs.Orthographic(30, 25))
axb.set_global()
axb.add_feature(cfeature.LAND, fc="#e8e3d3", ec="none"); axb.add_feature(cfeature.OCEAN, fc="#f4f8fb")
axb.gridlines(draw_labels=False, linewidth=0.3, color="gray", xlocs=range(-180, 181, 10), ylocs=range(-90, 91, 10))
for lat0 in range(-90, 90, 10):
    a = np.cos(np.deg2rad(lat0 + 5))
    axb.fill([20, 30, 30, 20, 20], [lat0, lat0, lat0 + 10, lat0 + 10, lat0], transform=ccrs.Geodetic(),
             fc=cm(0.15 + 0.8 * a), ec="k", lw=0.5, alpha=0.95)
axb.set_title("(b) 10°×10° cells along one meridian:\nthe closer to the pole, the smaller the true area",
              fontsize=10.5)
axc = fig.add_subplot(1, 3, 3)
phi = np.linspace(0, 90, 400)
axc.plot(phi, np.cos(np.deg2rad(phi)), color=C_W, lw=2.2, label="relative area ∝ cos φ")
for res, mk in [(10, "o"), (30, "s")]:
    c = np.arange(res / 2, 90, res)
    w = gaw.band_weights(np.concatenate([-c[::-1], c]), normalize=False)
    w = w[len(c):]; w = w / w[0] * np.cos(np.deg2rad(c[0]))
    axc.plot(c, w, mk, color="k", ms=5, mfc="none", label=f"{res}° cells: exact band area sin φN − sin φS")
for p in (30, 60, 80):
    axc.plot([p, p], [0, np.cos(np.deg2rad(p))], color="gray", lw=0.8, ls=":")
    axc.annotate(f"{p}°: {np.cos(np.deg2rad(p)):.2f}", (p, np.cos(np.deg2rad(p))), textcoords="offset points",
                 xytext=(5, 6), fontsize=9)
axc.axhline(1, color="gray", lw=0.6)
axc.set_xlabel("Latitude |φ| (°)"); axc.set_ylabel("Cell area / equatorial cell area")
axc.set_title("(c) Cell area versus latitude;\na plain mean gives every cell the same weight", fontsize=10.5)
axc.legend(frameon=False, fontsize=8.8, loc="lower left"); axc.set_xlim(0, 90); axc.set_ylim(0, 1.08)
fig.tight_layout()
save(fig, "fig1_concept_cell_area_vs_lat.png")

# ======================================================================= Fig 2: Cartopy grid maps
fig = plt.figure(figsize=(14, 9.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.38, wspace=0.12)
lat25 = np.arange(-90, 90.01, 0.25); lon25 = np.arange(-180, 180.01, 0.25)
la_c = 0.5 * (lat25[1:] + lat25[:-1])
rel = gaw.band_weights(la_c, lat_bounds=lat25, normalize=False); rel = rel / rel.max()
ax1 = fig.add_subplot(gs[0, 0:2], projection=ccrs.PlateCarree())
m1 = ax1.pcolormesh(lon25, lat25, np.repeat(rel[:, None], lon25.size - 1, 1), cmap="viridis", vmin=0, vmax=1,
                    transform=ccrs.PlateCarree(), rasterized=True)
ax1.coastlines(lw=0.5, color="w"); ax1.set_global()
ax1.gridlines(draw_labels=True, linewidth=0.3, color="w", alpha=0.6, xlocs=range(-180, 181, 60),
              ylocs=range(-90, 91, 30))
ax1.set_title("(a) Regular 0.25°×0.25° grid (1440×720): relative area of each cell\n"
              "(cells look equal in the equirectangular projection, but they are not)", fontsize=11)
cb = fig.colorbar(m1, ax=ax1, orientation="horizontal", pad=0.08, fraction=0.05, shrink=0.8)
cb.set_label("Relative area (equator = 1) = normalised (sin φN − sin φS)")
ax2 = fig.add_subplot(gs[0, 2], projection=ccrs.NorthPolarStereo())
ax2.set_extent([-180, 180, 55, 90], ccrs.PlateCarree())
lat5 = np.arange(55, 90.01, 2.5); lon5 = np.arange(-180, 180.01, 2.5)
Z5 = np.repeat((np.sin(np.deg2rad(lat5[1:])) - np.sin(np.deg2rad(lat5[:-1])))[:, None], lon5.size - 1, 1)
Z5 = Z5 / np.sin(np.deg2rad(2.5))
ax2.pcolormesh(lon5, lat5, Z5, cmap="viridis", vmin=0, vmax=1, edgecolor="w", linewidth=0.15,
               transform=ccrs.PlateCarree(), rasterized=True)
ax2.coastlines(lw=0.5, color="w")
ax2.set_title("(b) Arctic close-up: 2.5°×2.5° cells (55–90°N)\nsame number of cells per row, rapidly shrinking area",
              fontsize=10.5)
ax3 = fig.add_subplot(gs[1, 0:2], projection=ccrs.Robinson())
lat1 = np.arange(-90, 90.01, 1.0); lon1 = np.arange(-180, 180.01, 1.0)
la1 = 0.5 * (lat1[1:] + lat1[:-1]); r1 = gaw.band_weights(la1, lat_bounds=lat1, normalize=False); r1 = r1 / r1.max()
ax3.set_global()
ax3.pcolormesh(lon1, lat1, np.repeat(r1[:, None], lon1.size - 1, 1), cmap="viridis", vmin=0, vmax=1,
               transform=ccrs.PlateCarree(), rasterized=True)
ax3.coastlines(lw=0.5, color="w")
ax3.gridlines(draw_labels=False, linewidth=0.5, color="k", alpha=0.55, xlocs=range(-180, 181, 15),
              ylocs=range(-90, 91, 15))
ax3.set_title("(c) Robinson projection with a 15° graticule: meridians converge towards the poles\n"
              "(colour = relative area of 1° cells)", fontsize=11)
ax4 = fig.add_subplot(gs[1, 2])
ctr = np.array([(x["lo"] + x["hi"]) / 2 for x in EXT["bands"]])
lat_g = np.arange(-89.5, 90, 1.0); wa = gaw.area_weights(lat_g)   # ideal 1° grid
edges = np.arange(-90, 91, 10)
cs = [np.mean((lat_g >= lo) & (lat_g < hi)) for lo, hi in zip(edges[:-1], edges[1:])]
as_ = [wa[(lat_g >= lo) & (lat_g < hi)].sum() for lo, hi in zip(edges[:-1], edges[1:])]
ax4.barh(ctr - 1.9, np.array(cs) * 100, height=3.6, color=C_U, label="share of cells (plain-mean weight)")
ax4.barh(ctr + 1.9, np.array(as_) * 100, height=3.6, color=C_W, label="share of area (area-weighted)")
ax4.set_ylabel("Latitude band (°)"); ax4.set_xlabel("Share (%)")
ax4.set_title("(d) Per 10° latitude band: cells vs area", fontsize=10.5)
ax4.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=1); ax4.set_yticks(np.arange(-90, 91, 30))
save(fig, "fig2_cartopy_grid_area.png")

# ======================================================================= Fig 3: weighted vs unweighted per dataset
tmean = [("NCEP R1 air temperature, global (land + ocean)", "NCEP R1\nglobal 2.5°"),
         ("NCEP R1 air temperature, land excl. Antarctica", "NCEP R1\nland 2.5°"),
         ("CRU TS4.09 air temperature, land (0.5°, no Antarctica)", "CRU TS\nland 0.5°"),
         ("OISST v2 SST, global ocean (1°)", "OISST\nSST 1°"),
         ("NCEP R1 air temperature, 60-90N", "NCEP R1\n60–90°N"),
         ("NCEP R1 air temperature, 30S-30N", "NCEP R1\n30°S–30°N")]
pmean = [("GPCP precipitation, global (land + ocean)", "GPCP\nglobal P"),
         ("GPCC v2020 precipitation, land (observed cells, 1°)", "GPCC\nland P"),
         ("CRU TS4.09 precipitation, land (0.5°, no Antarctica)", "CRU TS\nland P"),
         ("NCEP R1 ET proxy, land excl. Antarctica", "NCEP R1\nland ET"),
         ("OISST sea-ice concentration, annual mean, 50-90N ocean", "OISST\nsea ice\n50–90°N"),
         ("SPEI-12 <= -1 (moderate or worse drought), area fraction", "SPEI-12 ≤ −1\narea")]
trs = [("NCEP R1 air temperature, global (land + ocean)", "NCEP R1\nglobal"),
       ("GISTEMP v4 temperature anomaly, global", "GISTEMP\nglobal"),
       ("CRU TS4.09 air temperature, land (0.5°, no Antarctica)", "CRU TS\nland"),
       ("OISST v2 SST, global ocean (1°)", "OISST\nocean"),
       ("NCEP R1 air temperature, 60-90N", "NCEP R1\n60–90°N"),
       ("NCEP R1 air temperature, 30S-30N", "NCEP R1\n30°S–30°N")]
fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
ax = axs[0]
sel = [(rows[n]["mean_diff"], lab) for n, lab in tmean if n in rows]
vals = [v for v, _ in sel]
ax.bar(range(len(vals)), vals, color=[C_U if v < 0 else C_W for v in vals])
for i, v in enumerate(vals):
    ax.text(i, v + (0.2 if v >= 0 else -0.2), f"{v:+.1f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=9)
ax.axhline(0, color="k", lw=0.8); ax.set_xticks(range(len(vals))); ax.set_xticklabels([lab for _, lab in sel], fontsize=8.3)
ax.set_ylabel("Unweighted − area-weighted (°C)"); ax.set_title("(a) Mean air / sea-surface temperature:\nabsolute difference (°C)", fontsize=11)
ax.set_ylim(min(vals) - 1.3, 1.5)
ax = axs[1]
sel = [(rows[n]["mean_rel_pct"], lab) for n, lab in pmean if n in rows]
vals = [v for v, _ in sel]
ax.bar(range(len(vals)), vals, color=[C_U if v < 0 else C_W for v in vals])
for i, v in enumerate(vals):
    ax.text(i, v + (3 if v >= 0 else -3), f"{v:+.0f}%" if abs(v) >= 10 else f"{v:+.1f}%", ha="center",
            va="bottom" if v >= 0 else "top", fontsize=9)
ax.axhline(0, color="k", lw=0.8); ax.set_xticks(range(len(vals))); ax.set_xticklabels([lab for _, lab in sel], fontsize=8.0)
ax.set_ylabel("(Unweighted − weighted) / weighted (%)")
ax.set_title("(b) Precipitation / ET / sea ice / drought area:\nrelative difference (%)", fontsize=11)
ax.set_ylim(min(vals) - 22, max(vals) + 25)
ax = axs[2]
x = np.arange(len(trs)); wd = 0.38
tu = [rows[n]["trend_unw"] for n, _ in trs]; tw = [rows[n]["trend_w"] for n, _ in trs]
ax.bar(x - wd / 2, tu, wd, color=C_U, label=L_U); ax.bar(x + wd / 2, tw, wd, color=C_W, label=L_W)
for i in range(len(trs)):
    rel = 100 * (tu[i] - tw[i]) / abs(tw[i])
    ax.text(i, max(tu[i], tw[i]) + 0.03, f"{rel:+.0f}%" if abs(rel) >= 10 else f"{rel:+.1f}%", ha="center", fontsize=9)
ax.set_xticks(x); ax.set_xticklabels([lab for _, lab in trs], fontsize=8.3); ax.set_ylabel("Linear trend (°C per decade)")
ax.set_title("(c) Trends, 1980–2020 (SST 1982–2022)", fontsize=11)
ax.legend(frameon=False, fontsize=9, loc="upper right"); ax.set_ylim(0, 1.15)
fig.tight_layout()
save(fig, "fig3_weighted_vs_unweighted_by_dataset.png")

# ======================================================================= Fig 4: time series + trend difference
if SER is None:
    print("docs/series.json not found: skipping fig4 (run run_experiments.py first)")
else:
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.2))

    def anom(s, base=(1981, 2010)):
        yrs = np.array(s["years"]); m = (yrs >= base[0]) & (yrs <= base[1])
        u, w = np.array(s["unw"]), np.array(s["w"])
        return yrs, u - np.mean(u[m]), w - np.mean(w[m])

    def tline(ax, yrs, y, c):
        p = np.polyfit(yrs, y, 1); ax.plot(yrs, np.polyval(p, yrs), color=c, lw=1, ls="--"); return p[0] * 10

    ax = axs[0, 0]; yrs, u, w = anom(SER["ncep_air_global"])
    ax.plot(yrs, u, color=C_U, lw=1.6, label=L_U); ax.plot(yrs, w, color=C_W, lw=1.6, label=L_W)
    tu, tw = tline(ax, yrs, u, C_U), tline(ax, yrs, w, C_W)
    ax.set_title(f"(a) NCEP/NCAR R1 near-surface air temperature, global anomaly (vs 1981–2010)\n"
                 f"trend: unweighted {tu:.2f} vs weighted {tw:.2f} °C/decade (+{100 * (tu - tw) / tw:.0f}%)",
                 fontsize=10.5)
    ax.set_ylabel("Anomaly (°C)"); ax.legend(frameon=False)
    ax = axs[0, 1]; yrs, u, w = anom(SER["gistemp_full"])
    ax.plot(yrs, u, color=C_U, lw=1.4, label=L_U); ax.plot(yrs, w, color=C_W, lw=1.4, label=L_W)
    m = yrs >= 1980; tu = np.polyfit(yrs[m], u[m], 1)[0] * 10; tw = np.polyfit(yrs[m], w[m], 1)[0] * 10
    ax.set_title(f"(b) GISTEMP v4 (2°) global temperature anomaly, 1880–2020\n"
                 f"1980–2020 trend: unweighted {tu:.2f} vs weighted {tw:.2f} °C/decade (+{100 * (tu - tw) / tw:.0f}%)",
                 fontsize=10.5)
    ax.set_ylabel("Anomaly (°C, vs 1981–2010)"); ax.legend(frameon=False)
    ax = axs[1, 0]; s = SER["spei_area"]
    yrs = np.array(s["years"]); u = 100 * np.array(s["unw"]); w = 100 * np.array(s["w"])
    ax.plot(yrs, u, color=C_U, lw=1.6, label=L_U); ax.plot(yrs, w, color=C_W, lw=1.6, label=L_W)
    tu, tw = np.polyfit(yrs, u, 1)[0] * 10, np.polyfit(yrs, w, 1)[0] * 10
    ax.set_title(f"(c) Land area fraction with SPEI-12 ≤ −1 (CSIC SPEIbase v2.10, 0.5°)\n"
                 f"mean {u.mean():.1f}% vs {w.mean():.1f}%; trend {tu:.2f} vs {tw:.2f} pp/decade: small difference",
                 fontsize=10.5)
    ax.set_ylabel("Area fraction (%)"); ax.legend(frameon=False)
    ax = axs[1, 1]; s = SER["GPCC v2020 precipitation, land (observed cells, 1°)"]
    yrs = np.array(s["years"]); u = np.array(s["unw"]); w = np.array(s["w"])
    ax.plot(yrs, u, color=C_U, lw=1.6, label=L_U); ax.plot(yrs, w, color=C_W, lw=1.6, label=L_W)
    tu, tw = np.polyfit(yrs, u, 1)[0] * 10, np.polyfit(yrs, w, 1)[0] * 10
    ax.set_title(f"(d) GPCC v2020 annual land precipitation (1°)\n"
                 f"mean {u.mean():.0f} vs {w.mean():.0f} mm/yr ({100 * (u.mean() - w.mean()) / w.mean():.0f}%); "
                 f"trend {tu:.1f} vs {tw:.1f} mm/yr per decade", fontsize=10.5)
    ax.set_ylabel("mm/yr"); ax.legend(frameon=False)
    for a in axs.ravel():
        a.set_xlabel("Year")
    fig.tight_layout()
    save(fig, "fig4_timeseries_trend_difference.png")

# ======================================================================= Fig 5: dependence on latitude / resolution / gradient
fig, axs = plt.subplots(2, 3, figsize=(15, 8.4))
ax = axs[0, 0]; ctr = np.array([(x["lo"] + x["hi"]) / 2 for x in EXT["bands"]])
ax.bar(ctr, [x["contrib_diff_K"] for x in EXT["bands"]], width=8.5,
       color=[C_U if x["contrib_diff_K"] < 0 else C_W for x in EXT["bands"]])
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("Latitude band centre (°)"); ax.set_ylabel("Contribution to (unweighted − weighted) (K)")
ax.set_title(f"(a) NCEP R1 temperature: contribution of each 10° band\n"
             f"total bias {EXT['bands_total_diff_K']:.1f} K (cold high latitudes are over-counted)", fontsize=10.5)
ax = axs[0, 1]
L = [r["L"] for r in EXT["lat_window_sym"]]
dT = [r["T_sym_unw"] - r["T_sym_w"] for r in EXT["lat_window_sym"]]
dP = [100 * (r["P_sym_unw"] - r["P_sym_w"]) / r["P_sym_w"] for r in EXT["lat_window_sym"]]
ax.plot(L, dT, "o-", color=C_U, label="temperature (left, K)")
ax.set_xlabel("Latitude window |φ| ≤ L (°)"); ax.set_ylabel("Unweighted − weighted (K)", color=C_U)
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
ax2.plot(L, dP, "s--", color=C_P, label="precipitation (right, %)"); ax2.set_ylabel("Precipitation relative difference (%)", color=C_P)
ax.set_title("(b) The bias grows with the latitude window;\nlow-latitude windows (|φ| ≤ 30°) are barely affected",
             fontsize=10.5)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, frameon=False, loc="lower left", fontsize=9)
ax = axs[0, 2]
Lp = [r["L"] for r in EXT["lat_window_pole"]]
ax.plot(Lp, [r["Ttrend_unw"] for r in EXT["lat_window_pole"]], "o-", color=C_U, label=L_U)
ax.plot(Lp, [r["Ttrend_w"] for r in EXT["lat_window_pole"]], "o-", color=C_W, label=L_W)
ax.set_xlabel("Window φ ≥ L (up to 90°N)"); ax.set_ylabel("Temperature trend (°C per decade)")
ax.set_title("(c) Small polar caps: weights become similar, trend gap shrinks;\n"
             "hemispheric / global windows: large gap", fontsize=10.5)
ax.legend(frameon=False)
ax = axs[1, 0]
rc, rg, rn = EXT["resolution_cru_T"], EXT["resolution_gpcc_P"], EXT["resolution_ncep_interp_T"]
ax.plot([r["res"] for r in rc], [100 * (r["T_unw"] - r["T_w"]) / abs(r["T_w"]) for r in rc], "o-", color=C_U,
        label="CRU TS land temperature (°C based)")
ax.plot([r["res"] for r in rg], [100 * (r["P_unw"] - r["P_w"]) / r["P_w"] for r in rg], "s-", color=C_P,
        label="GPCC land precipitation")
ax.plot([r["res"] for r in rn], [100 * (r["T_unw"] - r["T_w"]) / abs(r["T_w"]) for r in rn], "^-", color="gray",
        label="NCEP R1 global temperature (interpolated)")
ax.set_xscale("log"); ax.set_xlabel("Grid resolution (°)"); ax.set_ylabel("(Unweighted − weighted) / weighted (%)")
ax.set_title("(d) Resolution is no cure:\nthe bias does not vanish on finer grids", fontsize=10.5)
ax.legend(frameon=False, fontsize=8.5, loc="lower right", bbox_to_anchor=(1.0, 0.12))
ax = axs[1, 1]
syn = EXT["synthetic"]["polar-amplified trend s=0.15+b|φ|/90"]
bs = [float(k) for k in syn]
ratio = [100 * (syn[k]["unw"] - syn[k]["w"]) / syn[k]["w"] for k in syn]
ax.plot(bs, ratio, "o-", color=C_U)
ax.set_xlabel("Polar amplification b (extra warming at the pole, °C/decade)")
ax.set_ylabel("Overestimate of the unweighted trend (%)")
ax.set_title("(e) Idealised trend s(φ) = 0.15 + b·|φ|/90:\nno difference when the trend is uniform (b = 0)",
             fontsize=10.5)
for k, r in zip(bs, ratio):
    ax.annotate(f"{r:.0f}%" if abs(r) >= 0.5 else "0%", (k, r), textcoords="offset points", xytext=(4, -12), fontsize=8.5)
ax = axs[1, 2]
G = EXT["gap_vs_latvar"]
xs = [100 * g["lat_var_frac"] for g in G]; ys = [abs(g["gap_over_sigma"]) for g in G]
ax.scatter(xs, ys, s=70, color=C_U, zorder=3)
for xv, yv, g in zip(xs, ys, G):
    ax.annotate(g["name"], (xv, yv), textcoords="offset points",
                xytext={"GPCP": (6, 8), "SPEI-12": (6, 5), "SPEI <": (6, -11)}.get(
                    next((k for k in ("GPCP", "SPEI-12", "SPEI <") if g["name"].startswith(k)), ""), (6, -3)),
                fontsize=8.5,
                ha="right" if xv > 80 else "left")
ax.set_xlabel("Spatial variance explained by the zonal mean (%)")
ax.set_ylabel("|unweighted − weighted mean| / spatial std")
ax.set_title("(f) The stronger the latitudinal structure, the larger the bias\n"
             "(SPEI is standardised: almost no latitudinal structure)", fontsize=10.5)
ax.set_xlim(-3, 100); ax.set_ylim(-0.02, max(ys) * 1.15)
fig.tight_layout()
save(fig, "fig5_difference_vs_latitude_resolution_gradient.png")
print("figs done")
