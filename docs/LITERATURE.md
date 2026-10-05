# Literature and fact check

Searched on 2026-10-04 (Asia/Shanghai) by the author. Legend: ✔ = DOI checked
against Crossref (title, authors, year); "quantitative result ✔" = read in the
abstract or full text; ⚠ = not verified, do not cite without checking.

## A. The two *Nature* retractions (verified)

| Item | Paper 1 | Paper 2 |
|---|---|---|
| Title | **A 10 per cent increase in global land evapotranspiration from 2003 to 2019** | **Warming accelerates global drought severity** |
| Authors | Pascolini-Campbell, Reager, Chandanpurkar, Rodell (NASA JPL / GSFC) | Gebrechorkos, Sheffield, Vicente-Serrano, Funk, Miralles, Peng, Dyer, Talib, Beck, Singer, Dadson et al. (11 authors) |
| Journal | *Nature* 593(7860): 543–547, published online 2021-05-26 | *Nature* 642(8068): 628–635, published online 2025-06-04 |
| DOI | 10.1038/s41586-021-03503-5 | 10.1038/s41586-025-09047-2 |
| Retraction note | 10.1038/s41586-022-04525-3 (*Nature* 604, 202) | 10.1038/s41586-026-11027-z (*Nature* 657(8132): 844) |
| Retraction date | **2022-02-24** (Crossref `update-to` of type retraction; Retraction Watch database record 35758) | **2026-09-02** (Crossref `update-to` of type retraction); preceded by an Editor's Note on 2026-05-14 saying that criticisms had been raised and were being considered |
| Reason (from the notice) | "We made an error in calculating the global mean precipitation: we used **arithmetic averaging** to calculate the mean, instead of calculating a **spatially weighted mean to account for the changing grid box size with latitude**. As a result, the magnitudes of the global mean precipitation time series were underestimated. This impacted the subsequent calculation of global mean evapotranspiration, resulting in the mean evapotranspiration values being underestimated and altering some results." The notice thanks Ning Ma (Chinese Academy of Sciences) and others for pointing out the error. | Four issues: (1) "**area-weighted grid cells were not used when calculating study-wide averages**", affecting values in the main text, Fig. 1b and Extended Data Figs. 4 and 8; (2) Extended Data Fig. 4 used GLEAM v3 instead of GLEAM4 (with GLEAM4 and weighting, the 2022 drought area is smaller); (3) a time-step parsing error in the code (1981–2017 / 2018–2022 subsets) affecting Extended Data Fig. 8h; (4) an unclear description of the mask (the < 180 mm arid mask was applied to Fig. 1a, not 1b). The authors plan to resubmit a revised version; all authors agreed to the retraction. |
| Original claim | Global land ET increased by 10 ± 2 % in 2003–2019 (GRACE/GRACE-FO water balance) | Atmospheric evaporative demand (AED) increased drought severity by 40 % on average; the 2018–2022 drought area was 74 % larger than in 1981–2017, 58 % of it attributed to AED |
| Sources | nature.com/articles/s41586-021-03503-5; /s41586-022-04525-3; Retraction Watch, 2022-03-03 | nature.com/articles/s41586-025-09047-2; /s41586-026-11027-z; Retraction Watch, 2026-09-02 |

Precise wording:

* Paper 1 is **Pascolini-Campbell et al. 2021** (*Nature* 593), retracted on
  2022-02-24; it is sometimes misattributed (e.g. as "Pan et al. 2020").
* Paper 2 is **Gebrechorkos et al. 2025**, retracted on 2026-09-02, but **missing
  area weighting was only one of four reasons** (also a data-version error, a
  time-step bug and an unclear mask description). Retraction Watch quotes the
  authors as expecting the corrected percentages to be "slightly lower" with the
  overall conclusions unchanged (second-hand; ⚠ full text not read).
* So the accurate statement is: missing area weighting was *the* reason for the
  first retraction and *one of the* reasons for the second.

## B. Studies quantifying unweighted vs area-weighted differences

| # | Reference | Year | DOI | Key result | Status |
|---|---|---|---|---|---|
| 1 | Wei R., Li Y., Yin J., Ma X., *Comparison of Weighted/Unweighted and Interpolated Grid Data at Regional and Global Scales*, Atmosphere 13(12): 2071 | 2022 | 10.3390/atmos13122071 | Compares cos φ weighting, no weighting and interpolation to an equal-area grid for temperature and precipitation. Global temperature differs by **about 5 K** between weighted and unweighted means (CRU land); small differences for China, large for the USA/Canada (Alaska); largest in the 0–60°N window; asks authors to state how gridded data are averaged. | DOI ✔; abstract ✔; the 5 K figure is from a full-text snippet (⚠ full text not read) |
| 2 | Ju S.-D., Choi W.-J., Song H.-J., *Critical Role of Area Weighting on Estimated Long-Term Global Warming and Heat Wave Trends*, AppliedMath 4(4): 1618–1628 | 2024 | 10.3390/appliedmath4040086 | Climate model output 1850–2300: area weighting raises the global mean temperature by **8.2 °C**; warming trend 0.276 (weighted) vs 0.330 °C/decade (unweighted), i.e. **about 20 % overestimation** without weighting; frequency increase of > 35 °C heat waves overestimated by up to 5.4 % without weighting. | DOI ✔; abstract ✔ |
| 3 | Cowtan K., Jacobs P., Thorne P., Wilkinson R., *Statistical analysis of coverage error in simple global temperature estimators*, Dynamics and Statistics of the Climate System | 2018 | 10.1093/climsys/dzy003 | Separates area weighting from coverage bias: a naive area-weighted mean of incomplete observations still has coverage error, which GLS/kriging reduces; this removes about a third of the HadCRUT4 1998–2012 "slowdown". | DOI ✔; abstract excerpt ✔ (Maynooth repository) |
| 4 | Cowtan K., Way R.G., *Coverage bias in the HadCRUT4 temperature series and its impact on recent temperature trends*, QJRMS | 2014 | 10.1002/qj.2297 | Filling the unobserved (polar) regions increases the recent global warming trend: the spatial coverage / weighting choice materially affects trends. | DOI ✔; quantitative result ⚠ |
| 5 | Zeng, Conti, Xu et al., *Correct Temporal and Spatial Averaging of Atmospheric and Surface Variables for Weather and Climate Studies*, BAMS | 2025 | 10.1175/bams-d-25-0164.1 | Rules for correctly averaging temperature, dew point, specific humidity etc. in time and space; ERA5 global land dew point for July 2024 differs by 5.73 K between correct and default averaging; seasonal differences remain large, trend differences are smaller. (About averaging rules for variables, not cos φ weighting itself.) | DOI ✔; abstract ✔ |
| 6 | ROM SAF Report 10, *Latitudinal Binning and Area-Weighted Averaging of Irregularly Distributed Radio Occultation Data* | — | ⚠ no DOI | For irregularly distributed observations, cos φ weighting implicitly assumes uniform sampling within latitude bins and can bias results; on regular grids cos φ is appropriate; without weighting high latitudes are over-represented and gridded means are too cold. | search snippet only ⚠ |
| 7 | Hansen J., Lebedeff S., *Global trends of measured surface air temperature*, JGR 92(D11): 13345 | 1987 | 10.1029/JD092iD11p13345 | The classic equal-area-box method for global means of station data. | DOI ✔; method description from general knowledge ⚠ |
| 8 | Jones P.D., Osborn T.J., Briffa K.R., *Estimating sampling errors in large-scale temperature averages*, J. Climate 10: 2548 | 1997 | 10.1175/1520-0442(1997)010<2548:ESEILS>2.0.CO;2 | Sampling errors of large-scale averages; systematic discussion of grids and weights. | DOI ✔; quantitative result ⚠ |

## C. Cell geometry and tool documentation

| # | Reference | Year | DOI / URL | Notes | Status |
|---|---|---|---|---|---|
| 9 | Santini M., Taramelli A., Sorichetta A., *ASPHAA: A GIS-Based Algorithm to Calculate Cell Area on a Latitude-Longitude (Geographic) Regular Grid*, Transactions in GIS | 2010 | 10.1111/j.1467-9671.2010.01200.x (corrigendum 10.1111/tgis.12785, 2021) | Cell areas on regular lat-lon grids with an authalic sphere, compared with a numerical ellipsoid solution; some formulas of the original paper are corrected in the 2021 corrigendum. | DOI ✔; content from abstract ⚠ |
| 10 | Kmoch A. et al., *Area and shape distortions in open-source discrete global grid systems*, Big Earth Data | 2022 | 10.1080/20964471.2022.2094926 | Area and shape distortion of discrete global grid systems (equal-area alternatives to lat-lon grids). | DOI ✔; details ⚠ |
| 11 | Karney C.F.F., *Algorithms for geodesics*, J. Geodesy | 2012 (online) | 10.1007/s00190-012-0578-z | Standard algorithms for geodesics and areas on the ellipsoid. GeoAreaWeight uses the closed-form authalic-latitude formula instead and checks it against EPSG:6933. | DOI ✔ |
| 12 | Gortan M., Testa L., Fagiolo G., Lamperti F., *A unified dataset for pre-processed climate indicators weighted by gridded economic activity*, Sci. Data 11: 533 | 2024 | 10.1038/s41597-024-03304-1 | Socio-economic weighting; cites Wei et al. 2022 on the impact of the weighting choice. Only loosely related to area weighting. | DOI ✔; weak relevance |
| 13 | xarray example *Compare weighted and unweighted mean temperature* (Mathias Hauser) | — | https://docs.xarray.dev/en/latest/examples/area_weighted_temperature.html | Recommends `np.cos(np.deg2rad(lat))` with `.weighted(weights).mean(("lon", "lat"))`, noting that on a regular grid cos(lat) approximates the cell area. | read ✔ |
| 14 | NCL `wgt_areaave` / `wgt_areaave2`; CDO `fldmean` (uses cell areas from the grid description / spherical geometry) | — | ncl.ucar.edu/Document/Functions/Built-in/wgt_areaave.shtml; code.mpimet.mpg.de/boards/1/topics/484 | Area weighting in common tools. | search snippets ⚠ |
| 15 | spatialworkflow.io, *How to Take an Area-Weighted Mean over a Latitude–Longitude Grid* (blog, not peer reviewed) | — | https://www.spatialworkflow.io/area-weighted-mean-latlon-grid/ | NCEP R1 2.5°: 1991–2020 global mean temperature 4.92 °C unweighted vs 14.24 °C weighted; 1979–2025 trend 0.340 vs 0.196 °C/decade; cos, exact band and WGS84 differ by about 0.06 °C. | **Independently reproduced with NCEP R1: 1979–2025 trend 0.340 vs 0.196 °C/decade** (see EXPERIMENTS.md) |

## D. Material on the impact of the two retractions

* Retraction Watch, 2022-03-03, *NASA researchers retract Nature paper on climate
  change and evapotranspiration*: summarises the notice, names Ning Ma as the
  person who found the error; some reader comments misread "underestimated" as
  "overestimated". ✔ (page read)
* Retraction Watch, 2026-09-02, *Oxford researchers retract Nature paper on
  'atmospheric thirst' for calculation errors*: the first author calls it a
  "coding error" and says a revised manuscript with additional datasets is in
  preparation. ✔ (content from search results)
* Secondary coverage on climatescience.press (low quality, corroboration only).
* **Not found:** a peer-reviewed "Matters Arising" quantifying the impact of either
  retraction; whether Ning Ma's comment on paper 1 was published ⚠.

## E. Not verified: do not cite

* "Santer et al. 2000, JGR 105: 7337 (10.1029/1999JD901105) uses cos φ weighting":
  the DOI exists, but whether the paper discusses or quantifies area weighting was
  **not** checked.
* Specific provisions on area weighting in IPCC AR6 or CMIP guidance were not
  checked item by item. Method notes of diagnostic packages such as CVDP mention
  cos(lat) weighting (search snippets ⚠).
