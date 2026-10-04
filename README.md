# GeoAreaWeight

**经纬度网格面积加权平均工具 (Python / R / MATLAB) · Area-weighted means on latitude–longitude grids**

仓库 / Repo: https://github.com/GISWLH/GeoAreaWeight

> 只需要纬度，就能得到正确的面积权重。  
> Give it the latitudes – get correct area weights and weighted means.

## 为什么需要它 / Why

规则经纬度网格（如 0.25°、1°、2.5°）上，同一纬度带的格子数量相同，但**格子面积随 cos φ 缩小**：
60°N 的格子只有赤道格子的一半，80°N 只有 17%。直接对所有格子做算术平均（`np.mean`、`nanmean`、`mean(x(:))`），
等于让极地小格子拥有和热带大格子**同样的投票权**。

这个细节已经导致两篇 *Nature* 正刊论文撤稿：

| 论文 | 撤稿公告 | 原因（公告原文意思） |
|---|---|---|
| Pascolini-Campbell et al., *A 10 per cent increase in global land evapotranspiration from 2003 to 2019*, Nature 593, 543–547 (2021), doi:10.1038/s41586-021-03503-5 | 2022-02-24，doi:10.1038/s41586-022-04525-3 | 全球平均降水用了**算术平均**而不是考虑网格随纬度变化的空间加权平均，降水被低估，进而 ET 被低估 |
| Gebrechorkos et al., *Warming accelerates global drought severity*, Nature 642, 628–635 (2025), doi:10.1038/s41586-025-09047-2 | 2026-09-02，doi:10.1038/s41586-026-11027-z | 撤稿的**四条原因之一**：计算全域平均时没有使用面积加权的网格（另三条：Extended Data 图误用 GLEAM v3 而非 GLEAM4、代码时间步解析错误、掩膜说明不清） |

On a regular lat–lon grid every latitude row has the same number of cells, but a cell's area scales with cos φ. A plain
`mean()` therefore over-weights the poles. This detail contributed to two retractions in *Nature* (see table above).

## 功能 / Features

三种语言接口一致 (Python `geoareaweight`, R `geoareaweight`, MATLAB/Octave `gaw_*`)：

| 功能 | Python | R | MATLAB / Octave |
|---|---|---|---|
| cos φ 权重 | `cos_lat_weights` | `cos_lat_weights` | `gaw_cos_lat_weights` |
| 精确球面纬度带面积权重（sin φ 边界；任意不等间距、升/降序纬度、极点格点） | `band_weights` | `band_weights` | `gaw_band_weights` |
| WGS84 椭球（authalic 纬度，精确） | `band_weights(..., ellipsoid=True)` | `band_weights(ellipsoid=TRUE)` | `gaw_band_weights(lat,[],true)` |
| 一维纬度权重 (`cos`/`band`/`ellipsoid`/`none`) | `area_weights` | `area_weights` | `gaw_area_weights` |
| 格点面积 (m²/km²)，可给 lon、lat/lon 边界 | `cell_area` | `cell_area` | `gaw_cell_area` |
| 二维权重场 | `weights_2d` | `weights_2d` | `gaw_weights_2d` |
| **面积加权平均**（时间等维度保留；NaN 处权重逐步重新归一化） | `area_mean` | `area_mean` | `gaw_area_mean` |
| 区域 / 仅陆地 / 掩膜平均 | `area_mean(mask=)`, `region_mean` | `area_mean(mask=)` | `'Mask'`, `gaw_region_mean` |
| 面积积分（通量 → 总量） | `area_integral` | `area_integral` | `gaw_area_integral` |
| 面积占比（如干旱面积 %） | `area_fraction` | `area_fraction` | `gaw_area_fraction` |
| 自带权重（如 CMIP `areacella`） | `weights=` | `weights=` | `'Weights'` |
| xarray DataArray（自动识别 lat/lon、维度顺序） | ✔ | – | – |

## 快速开始 / Quick start

### Python
```bash
pip install ./python            # 需要 numpy；xarray 可选
```
```python
import xarray as xr, geoareaweight as gaw
da = xr.open_dataset("air.mon.mean.nc").air          # (time, lat, lon)
gmean = gaw.area_mean(da)                            # 面积加权（精确球面带面积），返回 DataArray(time)
naive = da.mean(("lat", "lon"))                      # ✘ 错误：算术平均
land  = gaw.area_mean(da, mask=land_mask)            # 仅陆地平均（权重在有效格点上重新归一化）
tot   = gaw.area_integral(pr_m_per_yr, lat, lon)     # 全球总量 m3/yr
frac  = gaw.area_fraction(spei <= -1)                # 干旱面积占比（0–1）
w     = gaw.area_weights(lat, method="cos")          # 只要纬度即可得到 1-D 权重（和为 1）
```
numpy 数组：`gaw.area_mean(x, lat, lon, lat_axis=-2, lon_axis=-1)`。

### R
```r
# 在仓库根目录: R CMD INSTALL R
library(geoareaweight)
# x: lon x lat x time（ncdf4 默认顺序）
m <- area_mean(x, lat, lon)                     # 默认 lon_dim=1, lat_dim=2
m_land <- area_mean(x, lat, lon, mask = land)   # mask: nlat x nlon
```

### MATLAB / Octave
```matlab
addpath('matlab');
x = ncread('air.mon.mean.nc','air');            % lon x lat x time
lat = ncread('air.mon.mean.nc','lat'); lon = ncread('air.mon.mean.nc','lon');
m = gaw_area_mean(x, lat, lon);                 % 默认 LonDim=1, LatDim=2
m = gaw_area_mean(x, lat, lon, 'Method','cos', 'Mask', landMask);
```

## 图示与实验 / Figures & experiments

![cell area vs latitude](figs/fig1_concept_cell_area_vs_lat.png)
![grid](figs/fig2_cartopy_grid_area.png)

* `docs/EXPERIMENTS.md`：在 NCEP R1、GISTEMP、GPCP/GPCC、OISST、CRU TS、SPEIbase 上比较未加权与面积加权（均值、趋势、纬度窗口、分辨率）。
* `docs/LITERATURE.md`：撤稿事实核查与相关文献清单（含“未核实”标注）。
* `experiments/`：生成上述结果与图的脚本（需自行下载公开数据，路径在脚本中硬编码为 `/workspace/geoareaweight/data`，使用前请修改）。

## 权重方法说明 / Methods

* `cos` – 权重 = cos φ（φ 为格点中心纬度）。细网格上近似极好；粗网格上有 O(Δφ²) 误差；极点格点（±90°）权重为 0。
* `band`（默认）– 每个格子的精确球面面积 ∝ `sin φ_north − sin φ_south`。若未提供边界，则取相邻中心的中点，首尾外推半个步长并截断到 ±90°，所以 NCEP 风格带 ±90° 点的网格也能正确处理（极点格点变成极冠半格）。支持**不规则间距**与升/降序纬度。
* `ellipsoid` – WGS84 椭球上精确的“纬度–经度格子”面积，使用 authalic 纬度（等积纬度）：`A = R_q² Δλ (sin β₂ − sin β₁)`，`R_q = 6 371 007.181 m`。与球面相比，全球平均通常只相差 ~0.05 °C（气温）量级（见实验），但在做面积积分/总量时更严谨。
* `none` – 全 1，等价于普通平均（用于对照）。
* 二维纬度（曲线/旋转网格）：只能用 cos φ 近似，请改用 `weights=` 传入真实格点面积（如 CMIP 的 `areacella`）。
* 高斯网格（如 NCEP T62、ERA5 N320 的高斯网格）：中点边界给出的权重与高斯积分权重相差最大约 1%（T62）；对区域平均影响可以忽略，需要极致精度时请传入真实面积。
* 强度量（温度、降水率）求**平均**用 `area_mean`；广延量（总降水量、总径流）求**总量**用 `area_integral`。
* 若格子只部分是陆地（海岸），可以把 `面积 × 陆地占比` 作为 `weights=` 传入。

## 测试状态 / Test status

* Python：20 个 pytest 通过；R：testthat 通过且 `R CMD check` 状态 OK；
* **MATLAB 代码只在 GNU Octave 中测试过（未在 MATLAB 中运行）**。/ **The MATLAB code has only been tested with GNU Octave, not with MATLAB itself.**
* 三种语言对同一组测试数据的 312 个结果最大相对差 ≈ 2×10⁻¹⁴（`validation/`）。

## 测试 / Tests
```bash
cd python && pip install -e .[test] && pytest              # 20 个测试
R CMD INSTALL R && Rscript -e 'testthat::test_dir("R/tests/testthat")'
octave-cli --eval "addpath('matlab'); addpath('matlab/tests'); test_gaw"   # MATLAB 同样可运行
cd validation && python ref_python.py ref_python.csv && Rscript ref_r.R ref_r.csv \
   && octave-cli --eval "addpath('../matlab'); ref_octave('ref_octave.csv')" && python compare.py   # 跨语言数值一致性
```
覆盖：全球面积和 = 4πR²（球面）/ 510 065 621.7 km²（WGS84）；与 EPSG:6933（WGS84 柱面等积投影）独立对照；降序/不规则纬度；NaN 重新归一化；掩膜；轴顺序；用户自带权重；xarray。

## 引用与参考 / References

* Pascolini-Campbell et al. 2022, Retraction Note, *Nature* 604, 202. doi:10.1038/s41586-022-04525-3
* Gebrechorkos et al. 2026, Retraction Note, *Nature* 657, 844. doi:10.1038/s41586-026-11027-z
* Wei, Li, Yin, Ma 2022, *Atmosphere* 13(12), 2071. doi:10.3390/atmos13122071
* Ju, Choi, Song 2024, *AppliedMath* 4(4), 1618–1628. doi:10.3390/appliedmath4040086

## License
MIT © 2026 Longhao Wang
