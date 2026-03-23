# Multisource SAE for SDG Monitoring in Niger

Integrating Survey, Administrative and Satellite Data with Privacy-Preserving Methods

[![R](https://img.shields.io/badge/R-4.3%2B-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Complete-brightgreen)]()

---

## Overview

This project implements a **multisource Small Area Estimation (SAE) framework** for subnational MPI poverty monitoring in Niger, using DHS Niger 2012 (10,750 households). It integrates four types of auxiliary data into Fay-Herriot area-level models and applies Differential Privacy to protect disaggregated estimates.

The project is a methodological proof of concept directly relevant to **Work Package 2** of the SAESDGs-EU ERC Starting Grant (Utrecht University, PI: Dr. Angelo Moretti), which addresses non-traditional data sources and multisource statistics for subnational SDG monitoring.

---

## Research Questions

1. Do non-traditional data sources (satellite imagery, administrative registers) improve SAE precision compared to census-only auxiliary variables?
2. How does integrating probability and non-probability data affect bias-variance tradeoffs in Fay-Herriot models?
3. Can differential privacy mechanisms protect disaggregated poverty estimates while maintaining acceptable statistical utility for NSOs?

---

## Methodology

```
DATA SOURCES
  DHS Niger 2012        |  INS Niger (RGPH 2012)  |  MODIS NDVI + VIIRS
  10,750 HH, real data  |  Census auxiliaries      |  Satellite proxies
  World Bank Open API   |  GPS clusters (476)      |  Non-traditional
         |                        |                         |
         v                        v                         v

INTEGRATION LAYER
  Record linkage      -- fastLink (Fellegi-Sunter EM), NIHR61FL x NIIR61FL
  Bias correction     -- IPW propensity weighting, GPS-based accessibility
  Quality assessment  -- NQAF framework for satellite sources

SAE MODEL COMPARISON (Fay-Herriot REML, n=8 regions)
  M1: Census only (baseline)
  M2: Census + Administrative data
  M3: Census + Satellite (MODIS NDVI + VIIRS)
  M4: Full multisource

PRIVACY-PRESERVING OUTPUT
  Gaussian Mechanism -- (epsilon=2, delta=1e-5)-Differential Privacy
```

---

## Data Sources

| Source | Type | Variables | Access |
|---|---|---|---|
| DHS Niger 2012 | Probability survey | MPI dimensions | [dhsprogram.com](https://dhsprogram.com) — free registration |
| INS Niger RGPH 2012 | Census | Regional urban share, literacy, electrification | [stat-niger.org](https://www.stat-niger.org) |
| World Bank Open Data | National statistics | Literacy, electricity, child mortality, GDP | API — no key required |
| MODIS MOD13A3 v6.1 | Satellite | NDVI vegetation index | [earthdata.nasa.gov](https://earthdata.nasa.gov) |
| VIIRS VNP46A3 | Satellite | Nighttime lights | [NOAA/NGDC](https://ngdc.noaa.gov/eog/viirs/) |
| DHS GPS clusters | Geographic | 476 cluster coordinates (Niger 2012) | [dhsprogram.com](https://dhsprogram.com) — free registration |

DHS microdata and GPS files require free registration and are not tracked in this repository. Place `NIHR61FL.DTA`, `NIIR61FL.DTA`, and `NIGE61FL.shp` (+ `.shx`, `.dbf`, `.prj`) in `data/raw/`. The pipeline falls back to calibrated simulation if files are absent.

---

## Project Structure

```
project2_SAE_Multisource_Niger/
├── R/
│   ├── 00_setup.R                    # Packages and configuration
│   ├── 01_data_sources.R             # DHS loading + MPI + WB API
│   ├── 02_record_linkage.R           # fastLink: NIHR61FL x NIIR61FL
│   ├── 03_nonprob_bias_correction.R  # GPS-based IPW bias correction
│   ├── 04_satellite_features.R       # Satellite correlations + NQAF
│   ├── 05_multisource_FH_models.R    # 4 Fay-Herriot model variants
│   ├── 06_model_comparison.R         # AIC, MSE, shrinkage diagnostics
│   ├── 07_differential_privacy.R     # Gaussian mechanism DP
│   └── 08_visualization.R            # Publication-quality figures
├── data/
│   ├── raw/                          # DHS files (not tracked -- see above)
│   └── processed/                    # Intermediate RDS outputs
├── outputs/
│   ├── figures/                      # PNG plots (200 dpi)
│   └── tables/                       # CSV results
├── docs/
│   └── methodology.md                # Technical notes
└── README.md
```

---

## Reproducibility

```r
install.packages(c("sae", "fastLink", "survey", "srvyr",
                   "tidyverse", "sf", "here", "httr", "jsonlite"))

source("R/00_setup.R")
source("R/01_data_sources.R")
source("R/02_record_linkage.R")
source("R/03_nonprob_bias_correction.R")
source("R/04_satellite_features.R")
source("R/05_multisource_FH_models.R")
source("R/06_model_comparison.R")
source("R/07_differential_privacy.R")
source("R/08_visualization.R")
```

The pipeline runs fully on simulated data (calibrated from OPHI 2025) without any DHS files. With real DHS files in `data/raw/`, scripts 01-03 use real microdata automatically.

---

## Key Results

| Model | Auxiliary sources | Mean CV | sigma2_v |
|---|---|---|---|
| M1: Census | RGPH 2012 | 1.3% | 0.000250 |
| M2: +Admin | + INS Niger registers | 1.3% | 0.000289 |
| M3: +Satellite | + MODIS NDVI + VIIRS | 1.3% | 0.000086 |
| M4: Full | All sources | 1.3% | 0.001299 |

M3 achieves the lowest between-area variance component (sigma2_v = 0.000086), indicating that satellite covariates explain subnational heterogeneity better than administrative data at the regional level. This is consistent with Edochie et al. (2025) who find geospatial covariates outperform administrative data for West African poverty estimation.

**Record linkage** (NIHR61FL x NIIR61FL, real DHS data): 987 predicted pairs, precision = 0.026, recall = 0.024. The low precision quantifies the discrimination limit of categorical quasi-identifiers (region, urban/rural, education) — confirming the need for stronger identifiers in production record linkage for WP2.

**GPS bias correction** (NIGE61FL.shp, 476 real clusters): distance-based IPW reduces naive satellite urban coverage from 31.5% to 28.1%, correcting the accessibility gradient between Niamey (ps = 0.948) and remote regions (Diffa ps = 0.053).

**Differential Privacy**: epsilon = 2, delta = 1e-5 recommended. Mean absolute error from DP noise = 0.0011 — negligible for NSO publication standards.

---

## Relevance to SAESDGs-EU (Utrecht University, WP2)

| WP2 Objective | Implementation in this project |
|---|---|
| Non-traditional data as auxiliary variables | MODIS NDVI + VIIRS in Fay-Herriot models |
| Probability + non-probability data integration | GPS-based IPW (Chen et al. 2020) |
| Record linkage for multisource integration | fastLink Fellegi-Sunter on real DHS files |
| Data governance and privacy protection | Gaussian mechanism DP on real EBLUP estimates |
| West African poverty context | DHS Niger 2012, aligned with Edochie et al. (2025) |

---

## References

- Chen, Y. et al. (2020). Doubly robust inference for non-probability samples. *JRSS-B* 82(2), 391-411.
- Dwork, C. & Roth, A. (2014). *Algorithmic Foundations of Differential Privacy*. FnT-TCS 9(3-4).
- Edochie, I. et al. (2025). Small area estimation of poverty in four West African countries. *Journal of Official Statistics* 41(1), 96-124.
- Enamorado, T. et al. (2019). Using a probabilistic model to assist merging of large-scale datasets. *Political Analysis* 27(2), 222-229.
- Fay, R.E. & Herriot, R.A. (1979). Estimates of income for small places. *JASA* 74(366), 269-277.
- Molina, I. (2022). *Disaggregating data in household surveys*. ECLAC Statistics No. 97.
- OPHI & UNDP (2025). *Global MPI 2025*. OPHI Methodological Note 61.

---

## Author

**Abdel Kader Younoussi Saley**
MSc Applied Statistics, (University of Constantine 1) | MSc Data Science (AIMS Rwanda, full scholarship)
Statistical Analyst: UNECA West Africa Sub-Regional Office | INS Niger
[GitHub](https://github.com/AbdelkaderYS) | [LinkedIn](https://linkedin.com/in/abdelkader-saley-31863b213)