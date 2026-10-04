# 步骤 1：撤稿事实核查 + 文献调研（经 Crossref / Nature 页面核对）

检索日期：2026-10-04（Asia/Shanghai）。✔ = DOI 已用 Crossref 核对标题/作者/年份；“量化结论”✔ = 摘要或原文中读到；⚠ = 未核实。

## A. 两篇 Nature 撤稿（已核实）

| 项 | 论文 1 | 论文 2 |
|---|---|---|
| 原文标题 | **A 10 per cent increase in global land evapotranspiration from 2003 to 2019** | **Warming accelerates global drought severity** |
| 作者 | Pascolini-Campbell, Reager, Chandanpurkar, Rodell（NASA JPL / GSFC） | Gebrechorkos, Sheffield, Vicente-Serrano, Funk, Miralles, Peng, Dyer, Talib, Beck, Singer, Dadson 等 11 人 |
| 期刊信息 | Nature 593(7860): 543–547，2021-05-26 在线 | Nature 642(8068): 628–635，2025-06-04 在线 |
| 原文 DOI | 10.1038/s41586-021-03503-5 | 10.1038/s41586-025-09047-2 |
| 撤稿公告 DOI | 10.1038/s41586-022-04525-3（Nature 604, 202） | 10.1038/s41586-026-11027-z（Nature 657(8132): 844） |
| 撤稿日期 | **2022-02-24**（Crossref `update-to` 类型 retraction，来源 publisher + Retraction Watch 记录 35758） | **2026-09-02**（Crossref `update-to` 类型 retraction）；此前 2026-05-14 有 Editor's Note“结论存在批评与更正，编辑部正在考虑” |
| 撤稿原因（公告原文） | “We made an error in calculating the global mean precipitation: we used **arithmetic averaging** to calculate the mean, instead of calculating a **spatially weighted mean to account for the changing grid box size with latitude**. As a result, the magnitudes of the global mean precipitation time series were underestimated. This impacted the subsequent calculation of global mean evapotranspiration, resulting in the mean evapotranspiration values being underestimated and altering some results.” 并感谢 Ning Ma（中科院）等人指出错误 | 四点：(1) “**area-weighted grid cells were not used when calculating study-wide averages**”，影响正文数值、Fig. 1b、Extended Data Figs. 4 和 8；(2) Extended Data Fig. 4 误用 GLEAM v3 而非 GLEAM4（且 GLEAM4 加权后 2022 年干旱面积减小）；(3) 代码中时间步解析错误（1981–2017 / 2018–2022 子集），影响 ED Fig. 8h；(4) 掩膜说明不清（<180 mm 干旱区掩膜用于 Fig. 1a 而非 1b）。作者将“修订后重投”，所有作者同意撤稿 |
| 论文原结论 | 2003–2019 全球陆地 ET 增加 10±2 %（GRACE/GRACE-FO 水量平衡） | AED 使干旱严重程度平均增加 40 %；2018–2022 干旱面积较 1981–2017 扩大 74 %，AED 贡献 58 % |
| 来源页面 | nature.com/articles/s41586-021-03503-5；/s41586-022-04525-3；Retraction Watch 2022-03-03 | nature.com/articles/s41586-025-09047-2；/s41586-026-11027-z；Retraction Watch 2026-09-02 |

**对你原始表述的校正**：
* 论文 1 不是 “Pan et al. 2020”，而是 **Pascolini-Campbell et al. 2021（Nature 593）**，2022-02-24 撤稿。
* 论文 2 确为 **Gebrechorkos et al.（2025）**，2026-09-02 撤稿；但**撤稿原因不止面积加权一条**（还有数据版本、时间步 bug、掩膜说明）。Retraction Watch 引述作者称修正后的百分比数值预期“略低”、总体结论不变（该点为转述，⚠ 我未读到全文）。
* 因此“两篇都因为没有面积加权”准确说法是：论文 1 的唯一原因是面积加权；论文 2 的原因之一是面积加权。

## B. 直接量化“未加权 vs 面积加权”差异的文献

| # | 文献 | 年 | DOI | 关键结论 / 量化 | 核实状态 |
|---|---|---|---|---|---|
| 1 | Wei R., Li Y., Yin J., Ma X., *Comparison of Weighted/Unweighted and Interpolated Grid Data at Regional and Global Scales*, Atmosphere 13(12):2071 | 2022 | 10.3390/atmos13122071 | 对全球/区域温度、降水比较 cos φ 加权、未加权与插值到等面积网格；全球温度加权与未加权**差可达约 5 K**（CRU 陆地）；中国区差别小，美国/加拿大因高纬（阿拉斯加等）差别大；0–60°N 窗口差别最大；呼吁在论文中明确网格处理方法 | DOI ✔；摘要 ✔；5 K 数字来自检索到的正文片段（⚠ 我没读全文） |
| 2 | Ju S.-D., Choi W.-J., Song H.-J., *Critical Role of Area Weighting on Estimated Long-Term Global Warming and Heat Wave Trends*, AppliedMath 4(4):1618–1628 | 2024 | 10.3390/appliedmath4040086 | 气候模式 1850–2300：加权后全球平均气温比未加权**高 8.2 °C**；增温趋势加权 0.276 vs 未加权 0.330 °C/10年（**未加权高估约 20 %**）；>35 °C 热浪频率增加：未加权高估最多 5.4 % | DOI ✔；摘要 ✔ |
| 3 | Cowtan K., Jacobs P., Thorne P., Wilkinson R., *Statistical analysis of coverage error in simple global temperature estimators*, Dynamics and Statistics of the Climate System | 2018 | 10.1093/climsys/dzy003 | 区分“面积加权”与“覆盖偏差”：朴素面积加权平均在观测不完整时仍有覆盖误差；用 GLS/克里金可减小；用它能把 HadCRUT4 1998–2012 “暂缓”减少约三分之一 | DOI ✔；摘要片段 ✔（来自 Maynooth 档案摘要） |
| 4 | Cowtan K., Way R.G., *Coverage bias in the HadCRUT4 temperature series and its impact on recent temperature trends*, QJRMS | 2014 | 10.1002/qj.2297 | 未覆盖区（极地）填补后全球变暖趋势更大——说明“空间覆盖/加权”选择对趋势有实质影响 | DOI ✔；量化结论 ⚠（未读全文） |
| 5 | Zeng, Conti, Xu 等, *Correct Temporal and Spatial Averaging of Atmospheric and Surface Variables for Weather and Climate Studies*, BAMS | 2025 | 10.1175/bams-d-25-0164.1 | 讨论温度、露点、比湿等变量**正确的时间/空间平均规则**；ERA5 2024-07 全球陆地露点“正确 vs 默认平均”相差 5.73 K；季节差仍大、趋势差较小。（注意：它主要谈**变量平均规则**，不是 cos φ 加权本身） | DOI ✔；摘要 ✔ |
| 6 | ROM SAF Report 10, *Latitudinal Binning and Area-Weighted Averaging of Irregularly Distributed Radio Occultation Data* | — | ⚠ 无 DOI | 对**不规则分布**的观测，cos φ 加权隐含“按纬度均匀分布”假设，可能引入偏差；规则经纬度网格才可直接用 cos φ；无加权时高纬点被过度代表 → 格点平均偏冷 | 来自检索片段（⚠ 未读全文） |
| 7 | Hansen J., Lebedeff S., *Global trends of measured surface air temperature*, JGR 92(D11):13345 | 1987 | 10.1029/JD092iD11p13345 | 经典的**等面积网格盒**全球平均方法（把站点分配到等面积盒后平均） | DOI ✔；方法描述为常识 ⚠（未重读全文） |
| 8 | Jones P.D., Osborn T.J., Briffa K.R., *Estimating sampling errors in large-scale temperature averages*, J. Climate 10:2548 | 1997 | 10.1175/1520-0442(1997)010<2548:ESEILS>2.0.CO;2 | 大尺度平均的采样误差，格点与权重选择的系统讨论 | DOI ✔；量化结论 ⚠ |

## C. 面积/网格几何与工具文档

| # | 文献 / 文档 | 年 | DOI / URL | 要点 | 状态 |
|---|---|---|---|---|---|
| 9 | Santini M., Taramelli A., Sorichetta A., *ASPHAA: A GIS-Based Algorithm to Calculate Cell Area on a Latitude-Longitude (Geographic) Regular Grid*, Transactions in GIS | 2010 | 10.1111/j.1467-9671.2010.01200.x（勘误 10.1111/tgis.12785, 2021） | 经纬度规则网格格点面积用 authalic（等积）球近似，并与椭球数值解对比；原文部分公式有误，2021 勘误给出修正 | DOI ✔；内容来自检索摘要 ⚠ |
| 10 | Kmoch A. 等, *Area and shape distortions in open-source discrete global grid systems*, Big Earth Data | 2022 | 10.1080/20964471.2022.2094926 | 比较各 DGGS 的面积/形状畸变（等面积网格的替代方案） | DOI ✔；细节 ⚠ |
| 11 | Karney C.F.F., *Algorithms for geodesics*, J. Geodesy | 2012 (online) | 10.1007/s00190-012-0578-z | 椭球上测地线/面积计算的标准算法（本库的 WGS84 方案用 authalic 纬度闭式解，并用 EPSG:6933 做了独立对照） | DOI ✔ |
| 12 | Gortan M., Testa L., Fagiolo G., Lamperti F., *A unified dataset for pre-processed climate indicators weighted by gridded economic activity*, Sci. Data 11:533 | 2024 | 10.1038/s41597-024-03304-1 | 谈**社会经济加权**，并引用 Wei et al. 2022 说明加权选择会显著影响结果；与“面积加权”只是旁系相关 | DOI ✔；相关性弱 |
| 13 | xarray 文档示例 *Compare weighted and unweighted mean temperature*（Mathias Hauser） | — | https://docs.xarray.dev/en/latest/examples/area_weighted_temperature.html | 官方推荐 `np.cos(np.deg2rad(lat))` + `.weighted(weights).mean(("lon","lat"))`；说明规则网格“cos 纬度近似格子面积” | 已读 ✔ |
| 14 | NCL `wgt_areaave` / `wgt_areaave2`；CDO `fldmean`（面积权重取 cell_area/球面几何，非 1/cosφ） | — | ncl.ucar.edu/Document/Functions/Built-in/wgt_areaave.shtml；code.mpimet.mpg.de/boards/1/topics/484 | 工具层面的面积加权做法 | 来自检索摘要 ⚠ |
| 15 | spatialworkflow.io, *How to Take an Area-Weighted Mean over a Latitude–Longitude Grid*（博客，非同行评议） | — | https://www.spatialworkflow.io/area-weighted-mean-latlon-grid/ | NCEP R1 2.5°：1991–2020 全球平均气温未加权 4.92 °C vs 加权 14.24 °C；1979–2025 趋势 0.340 vs 0.196 °C/10年；cos、精确带、WGS84 三者差 ~0.06 °C | **我用 NCEP R1 独立复现：1979–2025 趋势 0.340 vs 0.196 °C/10年，完全一致**（见 EXPERIMENTS.md） |

## D. 关于两次撤稿影响的材料

* Retraction Watch 2022-03-03，*NASA researchers retract Nature paper on climate change and evapotranspiration*（https://retractionwatch.com/2022/03/03/...）：转述撤稿公告，Ning Ma 为指出错误者；评论区有误读（把“低估”读成“高估”）。✔（页面已读）
* Retraction Watch 2026-09-02，*Oxford researchers retract Nature paper on ‘atmospheric thirst’ for calculation errors*：作者（Gebrechorkos）称是“coding error”，已准备含额外数据集的修订稿。✔（检索到的页面内容）
* climatescience.press 二次报道（质量较低，仅作旁证）。
* **没有找到**同行评议的“Matters Arising”专门量化这两次撤稿的影响；Ning Ma 对论文 1 的评论是否公开发表 ⚠ 未找到。

## E. 未能核实 / 不要引用的内容

* “Santer et al. 2000, JGR 105:7337（10.1029/1999JD901105）使用 cos φ 加权”——DOI 存在，但我**没有核实**该文是否讨论面积加权及其量化，不建议在文章中引用。
* IPCC AR6 / CMIP 指南中关于面积加权的**具体条款**未逐条核实。CVDP 等诊断包的方法说明提到 cos(lat) 加权（来自检索片段 ⚠）。
