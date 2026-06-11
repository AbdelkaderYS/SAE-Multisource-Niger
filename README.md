# Poverty in Niger from Space

Can we estimate poverty in every department of Niger using only free satellite imagery and a simple neural network?

Short answer: yes. A ResNet-18 trained on Landsat 7 patches, combined with a Fay-Herriot small area model, gives a CV of **11.1%** (the EUROSTAT reliability threshold is 15%). The full pipeline from raw DHS data to maps is documented and reproducible.

---

## TL;DR

| Metric | Value |
|---|---|
| Best FH model | FH-1c (urban_pct + cnn_score) |
| FH-1c CV | **11.1%** |
| Bootstrap IC 95% | [7.4%, 15.2%] |
| CV without 1-cluster depts | 13.1% |
| CNN R² global OOF | **0.482** |
| CNN R² CV (mean ± std) | 0.320 ± 0.156 |
| If random split (like literature) | R² = 0.69 (bias of 0.21) |
| Data | DHS 2012, 10,750 households, 476 clusters |
| Satellite | Landsat 7, 30m RGB, 224x224 patches |

---

## Methodology (5 steps)

### Step 1: DHS data to direct estimates
**Script**: `scripts/preparer_donnees.R`

The DHS 2012 household survey (10,750 households across 476 GPS clusters) is spatially joined to Niger's 64 department boundaries (Admin 2). A survey-weighted mean of household wealth is computed per department, giving a direct estimate with its standard error and coefficient of variation (CV).

**Output**: `data/processed/direct_estimates_dept.csv` (64 departments)

### Step 2: Landsat 7 patches
**Script**: `scripts/telecharger_patches_gee.js` (GEE) or `scripts/telecharger_landsat.py` (STAC/rasterio)

For each DHS cluster, we fetch all Landsat 7 scenes (2011-2012) within a 3,360m radius, filter by cloud cover (<30%), compute a **median composite** across available scenes, and export as 224x224 pixels at 30m resolution. Three RGB bands only (compatible with ImageNet pretraining).

**Output**: 476 GeoTIFF files in `data/processed/patches_landsat/`

### Step 3: ResNet-18 with spatial cross-validation
**Script**: `scripts/entrainer_resnet.py`

A pretrained ResNet-18 is fine-tuned on the 476 patches. Early layers (edge/texture detectors) are frozen. Layers 3 and 4 (semantic patterns) are fine-tuned, with a regression head (`Linear(512→256) → ReLU → Dropout(0.3) → Linear(256→1)`).

Key features:
- **Group 5-fold CV by department**: entire departments are held out per fold (no spatial leakage)
- **Per-fold z-score normalization** of the target variable
- **T.Resize((224,224))**: handles variable patch sizes from GEE export
- **np.nan_to_num**: handles SLC-off NaN pixels from Landsat 7 (51 patches affected, no wealth bias)
- Learning curves, CV boxplot, and OOF predictions saved

**Outputs**: trained weights, CV metrics, OOF predictions, figures

### Step 4: Fay-Herriot small area model
**Script**: `scripts/fay_herriot.R`

Five model specifications are compared via REML:

| Model | Formula | AIC | CV |
|---|---|---|---|
| FH-0 | direct_mean ~ 1 | 1532.5 | 93.7% |
| FH-1 | ~ urban_pct | 1414.0 | 17.8% |
| FH-1b | ~ urban_pct + ndvi_mean | 1416.1 | 18.1% |
| **FH-1c** | **~ urban_pct + cnn_score** | **1355.0** | **11.1%** |
| FH-2 | ~ urban_pct + ndvi_mean + cnn_score | 1356.8 | 11.4% |

FH-1c is selected as best (minimum CV). A bootstrap (1000 resamples, median) gives the 95% confidence interval [7.4%, 15.2%].

**Outputs**: `fh_results_comparison.csv`, `fh_results_details.csv`

### Step 5: Diagnostics and maps

The best model (FH-1c) is validated with Q-Q plots, residual diagnostics, and predictor correlation matrix. Department maps compare direct estimates with EBLUPs. Nine figures are generated.

---

## Project Structure

```
poverty_sae_satellite/
  scripts/
    entrainer_resnet.py      ResNet-18 training + spatial CV
    fay_herriot.R            FH models + diagnostics + maps
    preparer_donnees.R       DHS data to direct estimates
    telecharger_landsat.py   Landsat 7 download via STAC
    telecharger_modis.R      NDVI MODIS download
    telecharger_patches_gee.js  GEE export (alternative)

  data/
    raw/                     DHS data (restricted), admin boundaries
    processed/               estimates, labels, patches, predictions

  outputs/
    tables/                  CSV metrics and model comparisons
    figures/                 maps, diagnostics, learning curves (9 files)
    models/                  trained ResNet-18 weights
    report/                  LaTeX report (10 pages)

  exploration/               R Markdown notebooks for EDA
  pipeline_colab.ipynb       Colab launcher (T4 GPU, ~2h)
  pipeline_kaggle.ipynb      Kaggle launcher (GPU)
```

---

## Limitations

1. **Data from 2012 (temporal gap)**. DHS, Landsat, and NDVI are all from 2011-2012. Niger's population has grown from ~17M to ~27M. The spatial patterns are probably still valid, but absolute wealth estimates are stale.

2. **Small CNN sample (476 patches)**. A ResNet-18 with 11M parameters trained on 476 patches has a high parameter-to-observation ratio. The Group 5-fold CV detects overfitting but doesn't eliminate this constraint.

3. **Nine departments with 1 cluster**. These have psi=0, gamma=1 (no shrinkage). EBLUP = direct estimate, no SAE gain. Flagged in `fh_results_details.csv`.

4. **Only 3 FH predictors**. Urbanization rate, NDVI, and CNN score are a thin set compared to similar studies. No shocks, infrastructure, or demographic variables.

5. **No spatial autocorrelation test**. The FH model assumes independent sampling errors. Neighboring departments share similar conditions, so standard errors may be underestimated.

6. **Landsat 7 SLC-off**. Scan line gaps affect 51 patches (10.7%). Managed by `np.nan_to_num` with no wealth bias (p=0.997).

7. **Generalization limited to Niger**. The landscape and wealth distribution of Niger are distinct from other African countries. R² and CV do not transfer without retraining.

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
# Full pipeline (requires restricted DHS files)
Rscript scripts/preparer_donnees.R
python scripts/entrainer_resnet.py   # Colab/Kaggle recommended (GPU)
Rscript scripts/fay_herriot.R

# Or use the notebooks:
#   pipeline_colab.ipynb  (T4 GPU, ~2h)
#   pipeline_kaggle.ipynb (Kaggle GPU)
```

---

## References

- Fay & Herriot (1979). Estimates of income for small places. JASA, 74(366), 269-277.
- Jean et al. (2016). Combining satellite imagery and ML to predict poverty. Science, 353(6301), 790-794.
- Yeh et al. (2020). Using satellite imagery to understand economic well-being in Africa. Nature Communications, 11, 2583.
- Perez et al. (2017). Poverty prediction with satellite imagery and ML. NIPS 2017 workshop.
- Rao & Molina (2015). Small Area Estimation (2nd ed.). Wiley.
- He et al. (2016). Deep residual learning for image recognition. CVPR 2016.
