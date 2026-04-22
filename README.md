# Multisource SAE for SDG Monitoring in Niger

Integrating Survey, Administrative and Satellite Data with Privacy-Preserving Methods

[![R](https://img.shields.io/badge/R-4.3%2B-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Complete-brightgreen)](https://github.com/AbdelkaderYS/SAE-Multisource-Niger)

---

## Overview

This project implements a **multisource Small Area Estimation (SAE) framework**
for subnational MPI poverty monitoring in Niger, using DHS Niger 2012
(10,750 households). It integrates four types of auxiliary data into
Fay-Herriot area-level models and applies Differential Privacy to protect
disaggregated estimates.

The project is structured as an exploratory study operating under a severe
structural constraint: Niger has **D = 8 administrative regions**, which
limits the statistical power of area-level models and produces well-documented
boundary REML behavior. Results are interpreted accordingly, as findings that
motivate methodological extensions rather than as definitive precision gains.

---

## Research Questions

1. Do non-traditional data sources (satellite imagery, administrative registers)
   improve SAE precision compared to census-only auxiliary variables?
2. How does integrating probability and non-probability data affect
   bias-variance tradeoffs in Fay-Herriot models?
3. Can differential privacy mechanisms protect disaggregated poverty estimates
   while maintaining acceptable statistical utility for NSOs?

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

SAE MODEL COMPARISON (Fay-Herriot REML, D=8 regions)
  M1: Census only (baseline)
  M2: Census + Administrative data
  M3: Census + Satellite (MODIS NDVI + VIIRS)
  M4: Full multisource (4 predictors, D/2 constraint)

PRIVACY-PRESERVING OUTPUT
  Gaussian Mechanism -- (epsilon=2, delta=1e-5)-Differential Privacy
```

---

## Data Sources

| Source | Type | Variables | Access |
|---|---|---|---|
| DHS Niger 2012 | Probability survey | MPI dimensions | dhsprogram.com — free registration |
| INS Niger RGPH 2012 | Census | Regional urban share, literacy, electrification | stat-niger.org |
| World Bank Open Data | National statistics | Literacy, electricity, child mortality, GDP | API — no key required |
| MODIS MOD13A3 v6.1 | Satellite | NDVI vegetation index | earthdata.nasa.gov |
| VIIRS VNP46A3 | Satellite | Nighttime lights | NOAA/NGDC |
| DHS GPS clusters | Geographic | 476 cluster coordinates (Niger 2012) | dhsprogram.com — free registration |

> DHS microdata and GPS files require free registration and are not tracked
> in this repository. Place `NIHR61FL.DTA`, `NIIR61FL.DTA`, and
> `NIGE61FL.shp` (+ `.shx`, `.dbf`, `.prj`) in `data/raw/`.
> The pipeline falls back to calibrated simulation if files are absent.

---

## Project Structure

```
SAE-Multisource-Niger/
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
│   ├── raw/                          # DHS files (not tracked)
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

The pipeline runs fully on simulated data calibrated from OPHI 2025 (Niger
MIS 2021: H = 79.9%, MPI = 0.535) without any DHS files. With real DHS files
in `data/raw/`, scripts 01–03 use real microdata automatically.

---

## Key Results

### 1. Model Comparison — Fay-Herriot REML (D = 8 regions)

| Model | Aux sources | N pred. | sigma2\_v | AIC | BIC | Mean CV |
|---|---|---|---|---|---|---|
| M1: Census | RGPH 2012 | 3 | 0.000000 | −35.38 | −34.99 | 1.262% |
| M2: +Admin | + INS registers | 6 | 0.000000 | −31.46 | −30.83 | 1.277% |
| M3: +Satellite | + NDVI, VIIRS | 6 | 0.000000 | −38.47 | −37.83 | 1.279% |
| M4: Full | Selected 4 pred. | 4 | 0.001299 | −22.88 | −22.41 | 1.262% |

**sigma2\_v boundary solutions (M1–M3).**
Models M1, M2, and M3 produce sigma2\_v = 0 under REML. This is a boundary
solution: the fixed-effect covariates fully absorb the between-area variation
at D = 8, leaving no residual for the random component to estimate. This
behavior is well-documented in the SAE literature for small D and is not a
numerical error (Rao & Molina 2015, §6.1.3; Datta & Lahiri 2000). It means
the EBLUP collapses to the synthetic regression estimator for these models,
and shrinkage factors gamma\_i reflect sampling variance structure only.

**Precision at regional level is already high under direct estimation.**
Direct CVs range from 0.42% (Zinder) to 3.98% (Niamey), with a mean of
1.80%. In this high-precision regime, area-level SAE models offer limited
efficiency gains — a known limitation when D is small and direct estimates
are already reliable (Molina 2022, §3.2). Mean CV differences across models
(1.262% to 1.279%) are statistically negligible.

**M3 does not outperform M1.**
Although M3 has the lowest sigma2\_v (0.000000 as in M1), its Mean CV
(1.279%) is marginally higher than M1 (1.262%). Satellite covariates (NDVI,
VIIRS) do not improve precision at the administrative region level in this
configuration. This may reflect the coarseness of regional-level aggregation,
which smooths out the fine-scale spatial heterogeneity that satellite data
capture well. This is consistent with the finding of Edochie et al. (2025)
that geospatial gains are larger at sub-regional levels than at
administrative macro-regions.

**M4 sigma2\_v is non-zero (0.001299).**
With 4 predictors selected for cross-source coverage (urban share,
infant mortality, NDVI, VIIRS), the M4 design matrix is not fully
rank-saturated at D = 8, allowing REML to estimate a non-zero between-area
variance. However, this comes at a cost: M4 has the worst AIC (−22.88)
and BIC (−22.41) of all models, indicating a poor fit-parsimony tradeoff.
The non-zero sigma2\_v in M4 is better interpreted as residual heterogeneity
unabsorbed by the 4-predictor specification than as a genuine signal of
area-level variation beyond the fixed effects.

**CV gains by region are heterogeneous.**
SAE precision gains are concentrated in high-CV domains (Niamey, Agadez)
and negligible elsewhere:

| Region | CV\_direct | CV\_M4 | CV gain |
|---|---|---|---|
| Niamey | 3.98% | 4.11% | −0.13 pp |
| Agadez | 2.58% | 2.63% | −0.05 pp |
| Diffa | 0.98% | 0.99% | −0.01 pp |
| Maradi | 0.46% | 0.45% | +0.001 pp |
| Zinder | 0.42% | 0.42% | −0.001 pp |

For Niamey and Agadez — the two highest-CV regions — M4 yields slightly
*higher* CVs than direct estimation. This reflects M4's elevated sigma2\_v
inflating MSE beyond the sampling variance for these regions, an effect
documented by Rao & Molina (2015, §6.2) when sigma2\_v is estimated near
the boundary from below.

**Core finding on Research Question 1:**
At D = 8 administrative regions, non-traditional data sources do not produce
measurable precision gains over census-only SAE. The constraint is structural:
Niger's regional-level data are already sufficiently precise for direct
estimation, and the small number of areas prevents robust REML estimation of
the between-area variance. The methodological framework is valid; the
structural context limits what it can demonstrate. Extension to D > 30
subnational units (communes or departments) would be required to test the
framework under conditions where SAE gains are expected.

---

### 2. Record Linkage — Probabilistic Matching (fastLink)

Real DHS data: NIHR61FL × NIIR61FL linked on region, urban/rural, and
education level. Results: 987 predicted pairs, precision = 0.026,
recall = 0.024.

The near-random precision quantifies the **discrimination floor** of
categorical quasi-identifiers at the regional level in a household survey
context: region × urban/rural × education is not sufficient to distinguish
individual records. This is not a failure of the fastLink implementation
but a structural result — it identifies precisely why production record
linkage for official statistics requires stronger quasi-identifiers
(individual-level identifiers, GPS coordinates, or name-based matching).
This finding motivates the use of GPS-based bias correction rather than
record linkage as the integration strategy for non-probability sources
(see §3 below).

---

### 3. GPS Bias Correction — Non-Probability Source Reweighting

Real GPS cluster data: NIGE61FL.shp, 476 clusters. Distance-based IPW
corrects for the accessibility gradient in satellite coverage:

- Naive satellite urban coverage estimate: 31.5%
- IPW-corrected estimate: 28.1%
- Niamey propensity score: ps = 0.948 (high accessibility)
- Diffa propensity score: ps = 0.053 (remote, low accessibility)

The 3.4 percentage point correction (10.8% relative reduction) confirms
that satellite data systematically over-represent accessible urban areas in
Niger's territorial coverage. IPW correction using GPS distance as proxy
for inclusion probability is a valid approach for this type of bias
(Chen et al. 2020), though it relies on the assumption that haversine
distance to Niamey adequately captures accessibility — a simplification
that would require ground-truth validation for production use.

---

### 4. Differential Privacy

Configuration: epsilon = 2, delta = 1e-5 (Gaussian mechanism). Results
by region (full table in `outputs/tables/`):

| Region | EBLUP rate | DP estimate | Absolute error | n\_sample |
|---|---|---|---|---|
| Agadez | 0.728 | 0.727 | 0.0014 | 722 |
| Diffa | 0.904 | 0.905 | 0.0014 | 1116 |
| Dosso | 0.969 | 0.969 | 0.0006 | 1264 |
| Maradi | 0.960 | 0.959 | 0.0010 | 1856 |
| Niamey | 0.404 | 0.407 | 0.0023 | 1669 |
| Tahoua | 0.963 | 0.964 | 0.0012 | 1297 |
| Tillaberi | 0.957 | 0.957 | 0.0001 | 1691 |
| Zinder | 0.955 | 0.955 | 0.0006 | 1135 |

Mean absolute error = 0.0011. Maximum error = 0.0023 (Niamey). All errors
are below 0.003, which is negligible relative to regional poverty rates
ranging from 0.40 to 0.97. The recommended epsilon = 2 provides strong
privacy protection with no meaningful loss of statistical utility for NSO
publication purposes, consistent with the calibration in Dwork & Roth
(2014, §3.5.2) for bounded statistics with sensitivity Δf = 1/n\_i.

**Core finding on Research Question 3:** Differential privacy at
(ε=2, δ=1e-5) is compatible with NSO publication standards for regional
poverty rates in Niger given the sample sizes available (n\_i ≥ 722).

---

## Structural Limitations

**D = 8 is the binding constraint throughout this project.** It determines
the boundary REML behavior, limits the number of predictors in each model,
and restricts the bootstrap MSE estimation. This is not a methodological
flaw — it is the statistical reality of Niger's administrative structure
for area-level SAE. The appropriate response is:

1. Working at commune or department level (D ≈ 50–70) where data allow
2. Unit-level models (e.g., Battese-Harter-Fuller) which do not require
   large D
3. Hierarchical Bayesian approaches that handle boundary behavior more
   gracefully

These extensions define the research agenda that this project motivates.

---

## References

- Chen, Y. et al. (2020). Doubly robust inference for non-probability
  samples. *JRSS-B* 82(2), 391–411.
- Datta, G.S. & Lahiri, P. (2000). A unified measure of uncertainty of
  estimated best linear unbiased predictors in small area estimation
  problems. *Statistica Sinica* 10, 613–627.
- Dwork, C. & Roth, A. (2014). *Algorithmic Foundations of Differential
  Privacy*. FnT-TCS 9(3–4).
- Edochie, I. et al. (2025). Small area estimation of poverty in four West
  African countries. *Journal of Official Statistics* 41(1), 96–124.
- Enamorado, T. et al. (2019). Using a probabilistic model to assist
  merging of large-scale datasets. *Political Analysis* 27(2), 222–229.
- Fay, R.E. & Herriot, R.A. (1979). Estimates of income for small places.
  *JASA* 74(366), 269–277.
- Molina, I. (2022). *Disaggregating data in household surveys*.
  ECLAC Statistics Series No. 97.
- OPHI & UNDP (2025). *Global MPI 2025*. OPHI Methodological Note 61.
- Rao, J.N.K. & Molina, I. (2015). *Small Area Estimation* (2nd ed.).
  Wiley. [Added — was in methodology.md but missing from README]

---

## Author

**Abdel Kader Younoussi Saley**
MSc Mathematical Sciences (Data Science) — AIMS Rwanda (full merit scholarship)  
MSc Applied Statistics — University of Constantine 1
Research Collaborator, AIMS Rwanda \enspace\textbf{|}\enspace
Former Statistician Intern, UNECA West Africa (6 months - to February 2026)[Mail](saley.younoussi@aims.ac.rw) | [LinkedIn](https://linkedin.com/in/abdelkader-saley-31863b213)
