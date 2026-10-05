# Experiments: unweighted vs area-weighted means on real datasets

How large is the error when a gridded field is averaged with a plain `mean()`
instead of an area-weighted mean? This page compares the two on common public
datasets. All weighted values were computed with GeoAreaWeight
(`method="band"`, exact spherical cell areas) unless stated otherwise.

* Machine-readable results: [`results.json`](results.json) (main table plus
  method comparison, Gaussian-grid check and Arctic sea-ice extent),
  [`extras.json`](extras.json) (latitude windows, resolution, synthetic fields),
  [`summary_table.csv`](summary_table.csv).
* Scripts: [`../experiments/`](../experiments). `prep.py` → `run_experiments.py`
  → `run_extras.py` → `make_tables.py` → `make_figs.py`. Set the environment
  variable `GAW_DATA_DIR` to the folder holding the downloaded files (default:
  `data/` in the repository root). The scripts also need `xarray`, `netCDF4`,
  `scipy`, `regionmask`, `matplotlib` and `cartopy`.

## Datasets

| Dataset | Variable(s) | Grid | Files expected in `GAW_DATA_DIR` |
|---|---|---|---|
| NCEP/NCAR Reanalysis 1 (NOAA PSL) | 2 m air temperature, land mask | 2.5° (with pole rows) | `air.mon.mean.nc`, `land.nc` |
| NCEP/NCAR Reanalysis 1, surface fluxes | latent heat flux, precipitation rate, land mask | T62 Gaussian | `lhtfl.mon.mean.nc`, `prate.mon.mean.nc`, `land.sfc.gauss.nc` |
| GISTEMP v4 (NASA GISS) | temperature anomaly, 1200 km smoothing | 2° | `gistemp1200_GHCNv4_ERSSTv5.nc` |
| GPCP v2.3 (NOAA PSL) | precipitation | 2.5° | `gpcp.precip.mon.mean.nc` |
| GPCC Full Data v2020 (NOAA PSL) | precipitation, land | 1° | `precip.mon.total.1x1.v2020.nc` |
| NOAA OISST v2 (monthly) | SST, sea-ice concentration, land-sea mask | 1° | `oisst.sst.mnmean.nc`, `oisst.icec.mnmean.nc`, `oisst.lsmask.nc` |
| CRU TS 4.09 (UEA CRU) | temperature, precipitation, land | 0.5° | `cru_ts4.09.1901.2024.{tmp,pre}.dat.nc` |
| SPEIbase v2.10 (CSIC) | SPEI-12 | 0.5° | `spei12.nc` |

Processing: annual means are weighted by the number of days per month and use
complete years only; the SPEI-12 value of December represents the calendar year;
GPCC annual totals require all 12 months. Land masks are evaluated at cell
centres (NCEP and OISST masks from the providers, Natural Earth 1:110m via
`regionmask` for GPCP and the land-share statistic). Trends are ordinary
least-squares slopes per decade.

## Main results

| Dataset / region (unit) | Period | Mean, unweighted | Mean, area-weighted | Difference | Rel. diff. | Trend per decade, unweighted | Trend, area-weighted | Rel. diff. |
|---|---|---|---|---|---|---|---|---|
| NCEP R1 air temperature, global (land + ocean) (°C) | 1980-2020 | 4.78 | 14.16 | -9.38 | -66.2% | 0.331 | 0.189 | 75% |
| NCEP R1 air temperature, land incl. Antarctica (°C) | 1980-2020 | -5.00 | 9.14 | -14.14 | -154.7% | 0.419 | 0.274 | 53% |
| NCEP R1 air temperature, land excl. Antarctica (°C) | 1980-2020 | 8.68 | 13.11 | -4.43 | -33.8% | 0.308 | 0.253 | 22% |
| NCEP R1 air temperature, 60-90N (°C) | 1980-2020 | -10.45 | -7.39 | -3.05 | -41.3% | 0.865 | 0.689 | 26% |
| NCEP R1 air temperature, 30S-30N (°C) | 1980-2020 | 24.04 | 24.13 | -0.09 | -0.4% | 0.119 | 0.119 | 0% |
| NCEP R1 air temperature, 0-90N (°C) | 1980-2020 | 7.64 | 15.22 | -7.58 | -49.8% | 0.432 | 0.258 | 67% |
| GISTEMP v4 temperature anomaly, global (K) | 1980-2020 | 0.585 | 0.521 | 0.064 | 12.3% | 0.261 | 0.195 | 34% |
| GPCP precipitation, global (land + ocean) (mm/yr) | 1980-2020 | 818.4 | 982.6 | -164.2 | -16.7% | 2.45 | 1.58 | 55% |
| GPCP precipitation, land incl. Antarctica/Greenland (mm/yr) | 1980-2020 | 592.7 | 794.0 | -201.2 | -25.3% | 0.560 | 0.354 | 58% |
| GPCP precipitation, land excl. Antarctica (mm/yr) | 1980-2020 | 758.3 | 848.2 | -89.9 | -10.6% | 2.133 | 0.472 | 351% |
| GPCC v2020 precipitation, land (observed cells, 1°) (mm/yr) | 1980-2019 | 773.5 | 901.8 | -128.3 | -14.2% | 5.56 | 5.60 | -1% |
| NCEP R1 ET proxy from latent heat flux, land (T62 Gaussian grid) (mm/yr) | 1980-2020 | 486.9 | 652.4 | -165.5 | -25.4% | -1.47 | -4.22 | 65% |
| NCEP R1 ET proxy, land excl. Antarctica (mm/yr) | 1980-2020 | 623.4 | 697.8 | -74.4 | -10.7% | -3.28 | -5.02 | 35% |
| NCEP R1 precipitation, land excl. Antarctica (mm/yr) | 1980-2020 | 803.6 | 907.2 | -103.5 | -11.4% | 12.17 | 15.50 | -21% |
| OISST v2 SST, global ocean (1°) (°C) | 1982-2022 | 12.42 | 17.66 | -5.25 | -29.7% | 0.097 | 0.108 | -10% |
| OISST sea-ice concentration, annual mean, 50-90N ocean (%) | 1982-2022 | 59.96 | 35.46 | 24.50 | 69.1% | -1.58 | -1.15 | -38% |
| OISST September sea-ice area fraction (conc. >= 15%, 50-90N ocean) (fraction) | 1982-2022 | 0.496 | 0.222 | 0.274 | 123.7% | -0.032 | -0.025 | -29% |
| CRU TS4.09 air temperature, land (0.5°, no Antarctica) (°C) | 1980-2020 | 8.77 | 13.75 | -4.98 | -36.2% | 0.303 | 0.267 | 14% |
| CRU TS4.09 precipitation, land (0.5°, no Antarctica) (mm/yr) | 1980-2020 | 728.5 | 836.9 | -108.5 | -13.0% | 5.72 | 6.32 | -10% |
| SPEI-12 (December), land mean (index) | 1981-2022 | -0.086 | -0.107 | 0.021 | 19.6% | -0.097 | -0.103 | 6% |
| SPEI-12 <= -1 (moderate or worse drought), area fraction (fraction) | 1981-2022 | 0.214 | 0.219 | -0.006 | -2.5% | 0.028 | 0.028 | -1% |
| SPEI-12 <= -1.5 (severe drought), area fraction (fraction) | 1981-2022 | 0.098 | 0.100 | -0.003 | -2.8% | 0.018 | 0.018 | 3% |
| SPEI-12 <= -1 area fraction, excl. arid cells (P < 180 mm/yr) (fraction) | 1981-2022 | 0.193 | 0.195 | -0.003 | -1.5% | 0.008 | 0.007 | 21% |

How to read the table:

* **Relative differences are meaningless when the denominator is close to zero.**
  The GPCP land (excl. Antarctica) trend of 2.13 vs 0.47 mm/yr per decade gives
  "351 %" only because the weighted trend is near zero; the land precipitation
  and ET-proxy trends are small and noisy, so look at absolute differences.
* Relative differences of temperatures in °C depend on the arbitrary zero of the
  Celsius scale; use the absolute difference (K).
* Land means are very sensitive to Antarctica (about 9 % of the land area at
  about −50 °C): the NCEP R1 land mean is 9.1 °C including Antarctica and 13.1 °C
  without it. A "global land mean" must state whether Antarctica and Greenland are
  included. This is a separate issue from weighting, but the two are often mixed up.
* The GISTEMP "mean" is a mean of anomalies; only its trend row is informative.

## Key findings

1. **Global mean temperature (NCEP R1, 1980–2020).** Unweighted 4.78 °C vs
   area-weighted 14.16 °C (−9.4 K). Weighting method: cos φ 14.17, exact
   spherical band 14.16, WGS84 ellipsoid 14.10 °C. The ellipsoid differs from the
   sphere by only about 0.05 K; whether to weight matters about 180 times more
   than how to weight (9.4 K vs 0.05 K).
2. **Trends.** NCEP R1 global temperature 0.331 vs 0.189 °C/decade (the plain
   mean overestimates by 75 %). GISTEMP (2°, 1980–2020) 0.261 vs 0.195 (+34 %);
   CRU TS land 0.303 vs 0.267 (+14 %); 30°S–30°N 0.119 vs 0.119 (no difference).
   With the same data and method, the 1979–2025 NCEP R1 trend was 0.340 vs
   0.196 °C/decade, reproducing the numbers of the spatialworkflow.io blog post
   (that run is not stored in `results.json`).
   * Why the plain mean overestimates: the Arctic warms fastest (60–90°N:
     0.87 °C/decade unweighted, 0.69 weighted), and on the NCEP 2.5° grid the
     rows north of 60°N hold 17.8 % of the cells but only 7.3 % of the area.
   * Why some datasets differ less: CRU covers land without Antarctica, GISTEMP has
     gaps near the poles, and land warming is more uniform in latitude.
3. **Precipitation (means are underestimated).** GPCC land 774 vs 902 mm/yr
   (−14 %); CRU land 728 vs 837 (−13 %); GPCP global 818 vs 983 (−17 %). The
   direction agrees with the retraction notice of Pascolini-Campbell et al.
   (global mean precipitation underestimated, hence evapotranspiration
   underestimated). We did **not** reproduce that paper's GRACE water balance or
   precipitation ensemble, so these numbers show that a bias of this size is
   plausible, not what the paper's error was. Trends behave differently from
   means: GPCC 5.56 vs 5.60 mm/yr per decade (−1 %), CRU 5.72 vs 6.32 (−10 %).
4. **Drought area (the statistic of Gebrechorkos et al.).** SPEIbase SPEI-12 ≤ −1
   land area fraction, 1981–2022: 21.4 % vs 21.9 % (−2.5 %); trend 2.80 vs 2.83
   percentage points per decade; 19.3 % vs 19.5 % when arid cells
   (P < 180 mm/yr) are excluded. SPEI is standardised and has little latitudinal
   structure (the zonal mean explains about 17 % of its spatial variance), so
   weighting matters little here. This is consistent with the authors' statement,
   as reported by Retraction Watch, that corrected values are slightly lower with
   unchanged conclusions, but SPEIbase is **not** the paper's own high-resolution
   SPEI/AED ensemble, so it says nothing quantitative about that paper.
5. **Sea ice** is the most weight-sensitive example because it lies entirely at
   high latitudes: OISST annual-mean sea-ice concentration over the 50–90°N
   ocean is 60.0 % unweighted vs 35.5 % weighted (+69 %); the September fraction
   of cells with ≥ 15 % concentration is 0.50 vs 0.22.
6. **Latitude window** (NCEP R1 climatology over |φ| ≤ L): L = 30°: −0.09 K
   (precipitation −1.3 %); L = 60°: −1.7 K (−1.5 %); L = 80°: −5.9 K (−9.5 %);
   L = 90°: −9.4 K (−16.7 %). Tropical studies are hardly affected; the bias grows
   with the polar extent of the domain.
7. **Resolution is no cure.** NCEP temperature interpolated to 10°…0.25°: the
   difference stays at −8.8 to −8.9 K (−62 % to −63 %); CRU land temperature at
   0.5°–20°: −33 % to −36 %; GPCC land precipitation at 1°–20°: −14 % to −15 %.
   The bias comes from the latitudinal distribution of the weights, not from the
   number of cells.
8. **The bias scales with latitudinal structure.** Share of spatial variance
   explained by the zonal mean and |bias|/σ: NCEP temperature 93 %, 0.67; CRU
   temperature 89 %, 0.38; GPCP precipitation 32 %, 0.23; GPCC precipitation 42 %,
   0.15; SPEI 17 %, 0.05–0.07.
9. **Fields without latitudinal structure are unbiased.** For 500 white-noise
   fields (σ = 1) the unweighted − weighted difference has mean 8 × 10⁻⁵ and
   standard deviation 0.002: the error is systematic, not random, and appears only
   when the field varies with latitude.
10. **Idealised trend** s(φ) = 0.15 + b·|φ|/90 °C/decade: b = 0 → no difference;
    b = 0.4 → the unweighted trend is 18 % too high; b = 0.8 → 25 %.
11. **Land share on a 1° grid:** 33.2 % of the cells but 28.9 % of the area are land
    (Natural Earth mask at cell centres, Antarctica included).
12. **Gaussian grid (T62).** Band weights from mid-point bounds differ from the
    Gaussian quadrature weights by at most 0.96 % (cos φ: 1.8 %). Land ET-proxy
    mean: 697.78 mm/yr with Gaussian weights, 697.78 with band weights, 697.78 with
    cos φ, 623.4 unweighted, so the bounds approximation is negligible for regional
    means.
13. **Totals need true cell areas too.** Mean September Arctic sea-ice extent from
    OISST 1982–2022 is 6.40 million km² with true cell areas, but 34.2 million km²
    (5.4 times too large) if computed as "number of ice cells × global mean cell
    area (7871 km² on a 1° grid)".

## Figures

All figures are 300 dpi PNGs in [`../figs`](../figs); figures 1–3 and 5 are drawn
by `experiments/make_figs.py` from the JSON files in this folder.

* `fig1_concept_cell_area_vs_lat.png`: why equal-looking lat-lon cells have
  unequal areas; area ∝ cos φ vs the exact band area.
* `fig2_cartopy_grid_area.png`: relative cell area of a 0.25° grid
  (equirectangular), Arctic 2.5° cells (polar stereographic), Robinson projection
  with a 15° graticule, and share of cells vs share of area per 10° band.
* `fig3_weighted_vs_unweighted_by_dataset.png`: temperature biases (K), relative
  biases of precipitation / ET / sea ice / drought area, and trends.
* `fig4_timeseries_trend_difference.png`: time series and trends for NCEP R1
  temperature, GISTEMP, SPEI drought area and GPCC land precipitation. The
  underlying annual series were not stored in v0.1.0, so the English version of
  this figure was produced by relabelling the original rendering
  (`experiments/relabel_fig4.py`; the plotted data are unchanged).
  `run_experiments.py` now writes `docs/series.json`, after which `make_figs.py`
  redraws the figure from data.
* `fig5_difference_vs_latitude_resolution_gradient.png`: contribution of each
  latitude band, latitude windows, polar-cap trends, resolution, idealised trends
  and latitudinal structure.

## Limitations

* All weighted results come from GeoAreaWeight itself. The library is verified by
  (a) global area sums of 4πR² (sphere) and 510 065 621.7 km² (WGS84),
  (b) an independent comparison with the EPSG:6933 cylindrical equal-area
  projection (relative error < 10⁻⁹), (c) agreement with xarray's
  `DataArray.weighted`, and (d) 312 reference numbers that agree between the
  Python, R and Octave implementations to a relative difference of 2 × 10⁻¹⁴
  (see `../validation`).
* Land masks use the cell centre and ignore partial coastal cells; pass
  `cell_area × land_fraction` as `weights=` to account for them.
* We did not reproduce the original data or numbers of the two retracted papers;
  the experiments only show how large the unweighted-vs-weighted difference is on
  comparable public data.
