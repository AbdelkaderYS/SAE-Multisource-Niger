# poverty_sae_satellite

Estimating poverty at the regional level in Niger by combining the 2012 DHS survey, satellite variables (NDVI), and a deep learning model (CNN) on Sentinel-2 imagery, within a Small Area Estimation framework (Fay-Herriot model).

---

## Motivation

In Niger, Demographic and Health Surveys (DHS) are conducted by the National Institute of Statistics (INS) every 5 to 10 years. They produce precise estimates at the national level, but precision drops sharply at the regional level and below, because regional samples are too small (between 700 and 1900 households across 8 regions). DHS surveys are typically powered to Admin 1 (region) level, and direct estimates at finer scales suffer from high variance due to data sparsity [Wakefield et al., 2019; Wu et al., 2021]. Increasing sample size is costly and often infeasible [The DHS Program, SAR15].

Satellite imagery (Sentinel-2, MODIS, VIIRS) is free, regularly updated, and covers 100% of the territory. Jean et al. [2016] demonstrated that a CNN trained on publicly available satellite imagery can explain up to 75% of the variation in local-level economic outcomes across five African countries. Yeh et al. [2020] extended this to 20,000 African villages, achieving 70% explained variance in asset wealth. Deep learning can extract signals about wealth, urbanization, or economic activity from these images — signals that regional aggregates lose.

The Fay-Herriot model [Fay & Herriot, 1979] provides the statistical framework to combine a direct estimate (survey) with auxiliary variables (satellite) to "borrow strength" across regions and reduce estimation variance. It is the most widely used area-level small area estimation model [Rao & Molina, 2015] and has been adopted by the U.S. Census Bureau's SAIPE program [Bell et al., 2016] and the World Bank for poverty mapping [Newhouse et al., 2024; Edochie et al., 2024].

This project implements this survey + satellite + deep learning approach on Niger's 2012 DHS.

---

## Research questions

1. Does adding auxiliary satellite variables (mean NDVI, urbanization rate) in a Fay-Herriot model reduce the variance of regional wealth estimates compared to direct survey estimators?
2. Can a CNN trained on Sentinel-2 patches predict the mean wealth of a DHS cluster from satellite imagery alone?

---

## Hypotheses

- **H1**: The EBLUP estimator from the Fay-Herriot model has a strictly lower coefficient of variation (CV) than the direct estimator for small regions.
  - *Basis*: Newhouse et al. [2024] found that SAE with geospatial covariates yields precision gains equivalent to increasing sample size by a factor of 3–5. Edochie et al. [2024] reported 59–68% median CV reduction in West Africa.
- **H2**: The CNN achieves a strictly positive test R², demonstrating that satellite imagery contains a predictive signal of wealth at the cluster level.
  - *Basis*: Jean et al. [2016] achieved R² up to 0.75; Yeh et al. [2020] reported 70% variance explained across Africa; Pettersson et al. [2023] reached 72%.
- **H3**: Combining both approaches (regional FH + aggregated cluster CNN) yields more stable estimates than the survey alone.
  - *Basis*: The DHS Program [SAR15] demonstrated that geostatistical models yield precision equivalent to a survey three times larger when estimating at sub-national levels.

---

## Data

| Source | File | Usage |
|---|---|---|
| DHS Niger 2012 Household Recode | `data/raw/NIHR61FL.DTA` | Wealth index (HV271), region (HV024), urban/rural (HV025), weight (HV005) |
| DHS GPS clusters | `data/raw/NIGE61FL.shp` | Coordinates (LATNUM, LONGNUM) of 476 clusters for Python pipeline |
| Regional NDVI | simulated | Satellite proxy aggregated by region |
| Sentinel-2 | to download via GEE | 64x64 pixel patches, 4 bands (B2, B3, B4, B8) |

---

## Literature review

### Why DHS struggles at sub-national level

DHS surveys use a stratified two-stage cluster design. They are powered to produce reliable estimates for large domains (national, Admin 1). At finer spatial scales (Admin 2, communes), sample sizes per area become too small for reliable direct estimation [Wu et al., 2021; Wakefield et al., 2019]. This is well documented in Niger: Seidler et al. [2025] found that Niger's DHS data quality indicators show high sampling uncertainty at sub-national levels. The World Bank's poverty map for Niger [Edochie et al., 2024] explicitly notes that direct estimates at the commune level have median CVs exceeding 30%.

### Why satellite imagery can help

The link between satellite-observed environmental features and socioeconomic welfare is well established. NDVI, which measures vegetation greenness from red and near-infrared reflectance [Tucker, 1979], is negatively correlated with poverty intensity in West Africa [Sedda et al., 2015]. Engstrom et al. [2022] showed that Landsat-derived NDVI and built-up indices explain 25–50% of variation in economic activity at commune level in Vietnam. Watmough et al. [2019] found that satellite data predicted the poorest households in rural Kenya with 62% accuracy. Nighttime lights have long been used as a proxy for economic activity [Henderson et al., 2012; Chen & Nordhaus, 2011].

### Deep learning on satellite imagery for poverty prediction

The landmark paper by Jean, Burke, Xie, Davis, Lobell, & Ermon [2016] introduced a transfer learning approach: a CNN pre-trained on ImageNet is fine-tuned to predict nighttime lights from daytime imagery, then the learned features are used in ridge regression to predict consumption expenditure and asset wealth. They achieved up to 75% R² in five African countries.

Yeh et al. [2020] extended this with an end-to-end CNN trained on Landsat multispectral and nighttime lights imagery across 20,000 DHS clusters in Africa, explaining 70% of village-level wealth variation and 83% at district level. Chi et al. [2022] produced microestimates of wealth for all low- and middle-income countries using similar methods.

More recently, Pettersson et al. [2023] incorporated temporal information via LSTM layers to achieve 72% variance explained across held-out countries. The DeepWealth framework [Tonneau et al., 2024] achieved R² = 0.69 across 24 African countries and demonstrated generalization to Madagascar, Brazil, and Japan. Vision transformers have also been shown to outperform CNNs for wealth prediction from Landsat imagery [World Bank, 2024].

### Small Area Estimation: the Fay-Herriot model

The Fay-Herriot model [Fay & Herriot, 1979] is the foundational area-level SAE model. It treats the direct estimate as a noisy observation of the true value, modeled as:

```
direct_estimate_i = x_i' β + v_i + e_i
```

where `v_i ~ N(0, σ²_v)` is area-level random effect and `e_i ~ N(0, ψ_i)` is sampling error with known variance. The EBLUP is a weighted average of the direct estimate and the regression prediction, with weight γ_i = σ²_v / (σ²_v + ψ_i) [Rao & Molina, 2015].

The World Bank has extensively used FH models for poverty mapping [Newhouse et al., 2024; Edochie et al., 2024]. Battese, Harter, & Fuller [1988] were the first to combine satellite data with SAE. The DHS Program's SAR15 report found that geostatistical models yield precision at Admin 2 equivalent to a survey three times larger.

### Combining survey, satellite, and ML

The synthesis of these three approaches is the frontier of poverty measurement research. Jean et al. [2016] demonstrated the ML + satellite component; the SAE component is then applied to produce area-level estimates with proper uncertainty quantification. This project applies this integrated framework to Niger, where the need is acute due to infrequent censuses (last census 2012) and limited survey coverage [Edochie et al., 2024].

---

## Project structure

```
poverty_sae_satellite/
|-- README.md                this file
|-- NOTEBOOK.md              FH theory + step-by-step guide
|-- R/
|   |-- charger_dhs.R        load DHS, prepare variables
|   |-- estimateurs_directs.R  svyby by region
|   |-- modele_fay_herriot.R   eblupFH, EBLUP/SE/gamma
|   |-- visualiser_resultats.R 2 publication figures
|-- python/
|   |-- preparer_clusters.py join DHS wealth + GPS
|   |-- entrainer_cnn.py     CNN 3 convolutions, RMSE/R2
|   |-- gradcam.py           activation maps
|-- data/raw/                symlinks to DHS files
|-- outputs/
    |-- data/                intermediate RDS
    |-- tables/              CSV results
    |-- figures/             PNG figures
```

---

## Prerequisites

R packages (tested on version 4.3.3):
```
install.packages(c("haven", "survey", "sae", "dplyr", "ggplot2", "tidyr"))
```

Python packages (optional, for deep learning):
```
pip install torch torchvision pandas numpy geopandas scikit-learn matplotlib
```

Google Earth Engine access is required to download Sentinel-2 patches (free account at code.earthengine.google.com).

---

## How to run

R statistical pipeline (end-to-end, ~30 seconds):
```r
setwd("poverty_sae_satellite")
source("R/charger_dhs.R")
source("R/estimateurs_directs.R")
source("R/modele_fay_herriot.R")
source("R/visualiser_resultats.R")
```

Python deep learning pipeline (after downloading patches via GEE):
```bash
cd poverty_sae_satellite
python python/preparer_clusters.py
python python/entrainer_cnn.py
python python/gradcam.py
```

---

## Expected outputs

| File | Content |
|---|---|
| `outputs/tables/direct_estimates.csv` | Mean wealth and direct CV by region |
| `outputs/tables/fh_results.csv` | EBLUP, SE, gamma by region after Fay-Herriot |
| `outputs/tables/cnn_metrics.csv` | RMSE and R² of CNN on test set |
| `outputs/figures/01_direct_vs_eblup.png` | Direct vs EBLUP comparison with 95% CI |
| `outputs/figures/02_shrinkage.png` | Shrinkage factor gamma by region |
| `outputs/figures/03_cnn_pred_vs_true.png` | CNN predictions vs ground truth (test) |
| `outputs/figures/04_gradcam_*.png` | Grad-CAM activation maps |

---

## References

- Battese, G.E., Harter, R.M., & Fuller, W.A. (1988). An error-components model for prediction of county crop areas using survey and satellite data. *JASA*, 83(401), 28-36.
- Bell, W.R., Basel, W.W., & Maples, J.J. (2016). An overview of the U.S. Census Bureau's Small Area Income and Poverty Estimates program. In *Analysis of Poverty Data by Small Area Estimation*. Wiley.
- Chen, X. & Nordhaus, W.D. (2011). Using luminosity data as a proxy for economic statistics. *PNAS*, 108(21), 8589-8594.
- Chi, G., Fang, H., Chatterjee, S., & Blumenstock, J.E. (2022). Microestimates of wealth for all low- and middle-income countries. *PNAS*, 119(3), e2113658119.
- The DHS Program (2013). *Enquete Demographique et de Sante au Niger 2012*. INS Niger.
- The DHS Program (2019). The DHS Program Modeled Map Surfaces. *Spatial Analysis Reports No. 15*. ICF.
- Edochie, I., Newhouse, D., et al. (2024). Small Area Estimation of Poverty in Four West African Countries by Integrating Survey and Geospatial Data. *World Bank Policy Research Working Paper*.
- Engstrom, R., Hersh, J., & Newhouse, D. (2022). Can Medium-Resolution Satellite Imagery Measure Economic Activity at Small Geographies? *World Bank Policy Research Working Paper*.
- ESA (2015). Sentinel-2 User Handbook. European Space Agency.
- Fay, R.E. & Herriot, R.A. (1979). Estimates of income for small places. *JASA*, 74(366), 269-277.
- Henderson, J.V., Storeygard, A., & Weil, D.N. (2012). Measuring economic growth from outer space. *AER*, 102(2), 994-1028.
- Jean, N., Burke, M., Xie, M., Davis, W.M., Lobell, D.B., & Ermon, S. (2016). Combining satellite imagery and machine learning to predict poverty. *Science*, 353(6301), 790-794.
- Newhouse, D., et al. (2024). Small Area Estimation of Non-Monetary Poverty with Geospatial Data. *World Bank Policy Research Working Paper*.
- OPHI & UNDP (2025). *Global MPI 2025*. OPHI Methodological Note 61.
- Pettersson, M.B., Kakooei, M., Ortheden, J., Johansson, F.D., & Daoud, A. (2023). Time Series of Satellite Imagery Improve Deep Learning Estimates of Neighborhood-Level Poverty in Africa.
- Rao, J.N.K. & Molina, I. (2015). *Small Area Estimation*. 2nd ed. Wiley.
- Sedda, L., et al. (2015). Poverty and malaria: NDVI-based spatial analysis in West Africa.
- Seidler, V., et al. (2025). Subnational variations in the quality of household survey data in sub-Saharan Africa. *Nature Communications*.
- Tonneau, M., et al. (2024). DeepWealth: A generalizable open-source deep learning framework using satellite images for well-being estimation. *SoftwareX*, 26, 101700.
- Tucker, C.J. (1979). Red and photographic infrared linear combinations for monitoring vegetation. *Remote Sensing of Environment*, 8(2), 127-150.
- Wakefield, J., Fuglstad, G.-A., Riebler, A., Godwin, J., Wilson, K., & Clark, S.J. (2019). Estimating under-five mortality in space and time in sub-Saharan Africa. *arXiv:1910.06512*.
- Watmough, G.R., et al. (2019). Socioecologically informed use of remote sensing data to predict rural household poverty. *PNAS*, 116(4), 1213-1218.
- Wu, Y., Li, Z.R., Mayala, B.K., et al. (2021). Spatial Modeling for Subnational Administrative Level 2 Small-Area Estimation. *DHS Spatial Analysis Reports No. 21*. ICF.
- Xie, M., Jean, N., Burke, M., Lobell, D.B., & Ermon, S. (2016). Transfer learning from deep features for remote sensing and poverty mapping. *AAAI*.
- Yeh, C., Perez, A., Driscoll, A., Azzari, G., Tang, Z., Lobell, D., Ermon, S., & Burke, M. (2020). Using publicly available satellite imagery and deep learning to understand economic well-being in Africa. *Nature Communications*, 11, 2583.

---

## Author
