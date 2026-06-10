# Poverty in Niger from Space 🛰️

Can we estimate poverty in every department of Niger using only **free satellite imagery** and **a simple neural network**?

Short answer: yes. A ResNet-18 trained on Landsat 7 patches, combined with a Fay-Herriot small area model, gives a CV of **11.9%** (the EUROSTAT reliability threshold is 15%).

---

## TL;DR

| What | Result |
|---|---|
| Best model | FH-2 (urbanization + NDVI + CNN) |
| CV | **11.9%** [CI 95%: XX% to YY%] |
| CNN R² | 0.69 ± XX (spatial 5-fold CV) |
| Data | DHS 2012, 10,750 households, 476 clusters |
| Satellite | Landsat 7, 30m resolution, 224x224 patches |

---

## Methodology (5 steps)

### Step 1. DHS data to direct estimates
**Script**: `scripts/preparer_donnees.R`

The DHS 2012 household survey (10,750 households across 476 GPS clusters) is spatially joined to Niger's 64 department boundaries (Admin 2). A survey-weighted mean of household wealth is computed per department, giving a **direct estimate** with its standard error and coefficient of variation (CV).

**Output**: `data/processed/direct_estimates_dept.csv` (64 departments with mean wealth, SE, CV, survey variance)

### Step 2. Landsat 7 patches from Google Earth Engine
**Script**: `scripts/telecharger_patches_gee.js` (GEE) or `scripts/telecharger_landsat.py` (STAC/rasterio)

For each of the 476 DHS clusters, we fetch all Landsat 7 scenes (2011-2012) within a 3,360m radius. Only scenes with less than 30% cloud cover are kept. A **median composite** across all retained scenes eliminates residual clouds and shadows. Each patch is exported as 224x224 pixels at 30m resolution, covering roughly 6.7x6.7 km on the ground.

Only the 3 RGB bands are used (not NIR, SWIR, or thermal). The ResNet-18 is pretrained on ImageNet (3-channel RGB), so extra bands would break transfer learning. The NIR information is captured separately via MODIS NDVI in the FH model.

**Output**: 476 GeoTIFF files in `data/processed/patches_landsat/`

### Step 3. ResNet-18 training with spatial cross-validation
**Script**: `scripts/entrainer_resnet.py` (Kaggle GPU)

We use a **ResNet-18 pretrained on ImageNet** (11M parameters). Early layers (edge/texture detectors) are frozen. Layers 3 and 4 (semantic patterns like roads, fields, settlements) are fine-tuned. The classification head is replaced with a regression head: `Linear(512→256) → ReLU → Dropout(0.3) → Linear(256→1)`.

Why ResNet-18 and not deeper? Our dataset is small (476 patches). A deeper network would overfit. Satellite imagery also has less semantic variety than natural images.

**Validation** instead of a single 80/20 split, we use **Group 5-fold cross-validation by department**. Each fold holds out entire departments (about 13 per fold) and tests whether the model generalizes to unseen spatial locations. Learning curves (train/val MSE per epoch) and a boxplot of R² across folds are saved to detect overfitting.

Hyperparameters: Adam (lr=5e-4, weight decay=1e-4), batch size 32, max 100 epochs with early stopping (patience 10). Augmentation: random flip (p=0.5), random rotation (±15°).

**Outputs**: trained weights (`resnet18.pt`), CV metrics (mean ± std), out-of-fold predictions, learning curves and CV boxplot figures.

### Step 4. Fay-Herriot small area model
**Script**: `scripts/fay_herriot.R`

The Fay-Herriot model (1979) combines the **direct estimate** from the survey (unbiased but high variance for small departments) with a **synthetic estimate** from regression on auxiliary variables (potentially biased but low variance). The result is an EBLUP:

`EBLUP_i = gamma_i * direct_i + (1 - gamma_i) * synthetic_i`

where `gamma_i` is the shrinkage factor. It is close to 1 when the direct estimate is reliable (many households) and close to 0 when the model dominates.

Five model specifications are compared:

| Model | Formula | What it tests |
|---|---|---|
| FH-0 | `~ 1` | Baseline (mean only) |
| FH-1 | `~ urban_pct` | Does urbanization alone help? |
| FH-1b | `~ urban_pct + ndvi_mean` | Does NDVI add anything? |
| FH-1c | `~ urban_pct + cnn_score` | Does the CNN add anything? |
| FH-2 | `~ urban_pct + ndvi_mean + cnn_score` | Full model |

A **bootstrap** (1000 resamples) provides a 95% confidence interval around the FH-2 CV.

**Output**: `fh_results_comparison.csv` (all models), `fh_results_details.csv` (per department)

### Step 5. Diagnostics and maps

The FH-2 model is validated with:
- **Q-Q plot** of standardized residuals (normality check)
- **Residuals vs fitted values** (heteroscedasticity check)
- **Predictor correlation matrix** (printed to console)
- **1-cluster flag**: 10 departments have only 1 DHS cluster; their gamma is near 0, meaning the synthetic model dominates. They are flagged in the detailed results.

Four maps and figures are generated: poverty map (direct), poverty map (EBLUP), CV comparison across models, and gamma vs cluster size.

---

## Project Structure

```
poverty_sae_satellite/
  scripts/
    entrainer_resnet.py      ResNet-18 training + spatial CV
    fay_herriot.R            FH models + diagnostics + figures
    preparer_donnees.R       DHS data to direct estimates
    telecharger_landsat.py   Landsat 7 download via STAC
    telecharger_modis.R      NDVI MODIS download
    telecharger_patches_gee.js  GEE export (alternative)

  data/
    raw/                     DHS data (restricted), admin boundaries
    processed/               estimates, labels, patches, predictions

  outputs/
    tables/                  CSV metrics and model comparisons
    figures/                 maps, diagnostics, learning curves
    models/                  trained ResNet-18 weights
```

---

## Limitations

### 1. Data from 2012 (temporal gap)

The DHS survey, Landsat 7 scenes, and NDVI data are all from 2011-2012. Niger's population has grown from ~17M to ~27M since then. Urbanization patterns have shifted (Niamey, Maradi, Zinder have expanded significantly). Agricultural land use has changed with climate variability.

**Concrete impact**: the model's spatial patterns (which departments are poor vs less poor) are probably still broadly valid, but the absolute wealth estimates are stale. This is a methodological validation, not a current poverty map. Applying it today would require retraining on recent data.

### 2. Small CNN sample (476 patches)

A ResNet-18 has 11 million parameters. Training it on 476 patches gives a parameter-to-observation ratio of about 23,000:1. The Group 5-fold CV helps detect overfitting, but it does not eliminate the fundamental constraint: the model sees very few examples of each type of landscape (urban, rural, Sahel, desert).

**Concrete impact**: the R² of 0.69 may not generalize to regions or countries beyond Niger. Papers like Yeh et al. (2020) use 20,000+ clusters across multiple countries. Our single-country setup with 476 clusters is at the low end of what is feasible.

### 3. Ten departments with only 1 cluster

Out of 61 departments with DHS data, 10 have exactly one cluster (14 to 66 households). For these departments, the within-cluster variance cannot be estimated properly. The Fay-Herriot model compensates by giving nearly all weight to the synthetic component (gamma near 0).

**Concrete impact**: these 10 departments contribute little to the CV improvement. Their estimates are essentially the regression prediction, not the blended EBLUP. They are flagged in `fh_results_details.csv` for transparency. If you remove them, the reported CV improves slightly, but the model is doing less work.

### 4. Only 3 predictors in the FH model

Urbanization rate, NDVI, and CNN score are the only auxiliary variables. This is a thin set compared to similar studies. Yeh et al. (2020) use 150+ features. We do not capture:

- **Economic shocks**: drought frequency, conflict events, inflation, remittances
- **Infrastructure**: road density, market access, school proximity
- **Demographics**: population density, ethnic composition, migration

**Concrete impact**: the FH model's synthetic component is structurally simple. It captures the spatial pattern of poverty that correlates with greenness and built-up areas, but it cannot explain sudden changes due to drought, conflict, or policy.

### 5. No spatial autocorrelation test

The Fay-Herriot model assumes independent sampling errors across areas. In reality, neighboring departments in Niger share similar agro-ecological conditions, market access, and livelihood systems. If residuals are spatially correlated, the standard errors of the EBLUP are underestimated.

**Concrete impact**: the reported CV of 11.9% might be 13-14% if spatial dependence were accounted for. This is a known limitation of the basic FH model. Extensions like the spatial FH (SFH) or STAR models exist but add complexity.

### 6. Landsat 7 limitations

Landsat 7 has a scan line corrector failure (SLC-off) since 2003, creating data gaps (stripes) in every scene. The median composite across multiple scenes mitigates this, but some patches still have missing pixels. Newer missions (Landsat 8/9, Sentinel-2) would provide cleaner data.

### 7. Generalization beyond Niger

The model is trained and tested only on Niger. The landscape, settlement patterns, and wealth distribution of Niger are distinct from coastal West Africa (Nigeria, Ghana, Côte d'Ivoire) or East Africa. The 0.69 R² and 11.9% CV do not transfer to other countries without retraining.

---

## Data Sources

| Source | What for | Access |
|---|---|---|
| DHS 2012 (NIHR61FL) | Household wealth labels | Restricted (DHS Program) |
| Landsat 7 (GEE) | RGB satellite patches | Public (USGS) |
| MODIS MOD13A3 | NDVI per department | Public (NASA) |
| Admin 2 Niger | Department boundaries | Public (HDX/OCHA) |

---

## Reproduce

```bash
# 1. Prepare DHS data (requires restricted DHS files)
Rscript scripts/preparer_donnees.R

# 2. Train CNN (Kaggle GPU or local)
python scripts/entrainer_resnet.py

# 3. Run Fay-Herriot models
Rscript scripts/fay_herriot.R
```

Or use the Kaggle notebook: `pipeline_kaggle.ipynb`.

Note: three desert departments (Bilma, etc.) have no DHS clusters and are excluded. The Landsat patches are regenerable via GEE (JS and Colab notebook provided).

---

## References

- Fay & Herriot (1979). Estimates of income for small places. *JASA*, 74(366), 269-277.
- Jean et al. (2016). Combining satellite imagery and ML to predict poverty. *Science*, 353(6301), 790-794.
- Yeh et al. (2020). Using satellite imagery to understand economic well-being in Africa. *Nature Communications*, 11, 2583.
- Rao & Molina (2015). *Small Area Estimation* (2nd ed.). Wiley.
- He et al. (2016). Deep residual learning for image recognition. *CVPR 2016*.
