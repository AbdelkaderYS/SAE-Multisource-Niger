# Poverty Estimation in Niger with Satellite Imagery and Small Area Estimation


> ### *Predict poverty across Niger's 64 departments from free satellite imagery, without costly household surveys!*

> ### *Cibler les zones de pauvreté au Niger avec l'IA sur des images satellite, une alternative rapide et gratuite aux enquêtes ménages!*



Combine DHS 2012 household surveys with Landsat 7 satellite imagery and a
ResNet-18 neural network, then apply a Fay-Herriot small area model to estimate
poverty at the department level (Admin 2, 64 departments) in Niger.

**Best model** : FH-2 (urbanization + NDVI + CNN score) achieves a coefficient
of variation (CV) of 11.9 percent.

---

## Key Results

| Model          | Predictors                               | AIC              | CV (%)         | gamma           |
| -------------- | ---------------------------------------- | ---------------- | -------------- | --------------- |
| FH-0           | None (mean only)                         | 1532.5           | 93.7           | 0.916           |
| FH-1           | Urbanization rate                        | 1414.0           | 17.8           | 0.709           |
| **FH-2** | Urbanization + NDVI +**CNN score** | **1325.9** | **11.9** | **0.455** |

- **CNN test R squared**: 0.69 (ResNet-18, 380 train / 96 test clusters)
- **CNN test RMSE**: 49,904 (wealth index ranges from -70,000 to +465,000)
- **CV threshold**: 11.9 percent is below 15 percent, meeting the EUROSTAT
  standard for reliable estimates
- **Gamma**: 0.455 means the model contributes 54.5 percent of the weight,
  while the direct survey contributes 45.5 percent. This is a healthy balance
  for small area estimation.

The CNN alone explains 69 percent of the wealth variance between clusters
using only Landsat 7 pixels at 30 meter resolution (RGB bands). Adding these
predictions to the Fay-Herriot model reduces the CV from 17.8 percent
(urbanization alone) to 11.9 percent.

---

## How the ResNet Works

### Why ResNet-18

We use a ResNet-18 (Residual Network with 18 layers) pretrained on ImageNet
(1.2 million photos, 1000 classes). The key innovation of ResNet is the
**skip connection** (or residual connection), which lets the gradient flow
directly through the network during training. This solves the vanishing
gradient problem that plagued deep networks before 2015.

ResNet-18 has 11 million parameters. We chose 18 layers (not 34, 50, or 101)
because:

- Our dataset is small: 476 patches of 224 by 224 pixels. A deeper network
  would overfit.
- Each department has only 1 to 49 clusters (median 7). The model must
  generalize from limited data.
- Satellite imagery has less semantic variety than natural images. The
  features that distinguish a road from a field are simpler than those that
  distinguish a dog from a cat.

### Architecture Details

```
Input: 224 x 224 x 3 RGB patch (normalized to [0, 1])
  |
  |-- Conv1 (7x7, 64 channels, stride 2)
  |-- MaxPool (3x3, stride 2)
  |
  |-- Layer 1 (2 residual blocks, 64 channels)     -- FROZEN
  |-- Layer 2 (2 residual blocks, 128 channels)    -- FROZEN
  |-- Layer 3 (2 residual blocks, 256 channels)    -- FINE-TUNED
  |-- Layer 4 (2 residual blocks, 512 channels)    -- FINE-TUNED
  |
  |-- Average Pool (7x7)
  |-- Fully Connected (512 -> 256 -> 1)
  |-- Output: predicted wealth score
```

Each residual block contains two convolutional layers (3x3 kernels) with batch
normalization and ReLU activation, plus a skip connection that adds the input
directly to the output.

### Transfer Learning Strategy

We freeze layers 1 and 2 (the early layers that detect edges, colors, and
simple textures) because these features are universal and already well learned
from ImageNet. We fine-tune layers 3 and 4 (which capture higher-level
semantic patterns like roads, fields, and settlements) to adapt the network
to satellite imagery.

The original 1000-class classification head is replaced with a regression head:

- Linear(512 -> 256) + ReLU + Dropout(0.3)
- Linear(256 -> 1) for wealth prediction

### Training Setup

| Hyperparameter | Value             | Reason                                                         |
| -------------- | ----------------- | -------------------------------------------------------------- |
| Optimizer      | Adam              | Adaptive learning rate, standard for fine-tuning               |
| Learning rate  | 0.0005            | 10x lower than typical, to avoid destroying pretrained weights |
| Weight decay   | 0.0001            | Light L2 regularization against overfitting                    |
| Batch size     | 32                | Fits T4 GPU memory (16 GB)                                     |
| Max epochs     | 100               | Safety limit                                                   |
| Early stopping | Patience 10       | Stops when validation loss plateaus                            |
| LR scheduler   | ReduceLROnPlateau | Halves LR every 5 epochs without improvement                   |
| Test split     | 20 percent        | 380 train, 96 test clusters                                    |

### Data Augmentation (Training Only)

- Random horizontal flip (50 percent chance): satellite imagery has no
  inherent orientation.
- Random rotation up to 15 degrees: patches may come from different
  satellite paths.
- ImageNet normalization (mean and standard deviation per RGB channel).

### Training Progress

The loss dropped steadily from 6.9 billion (epoch 1) to 0.66 billion
(epoch 100), a 10x reduction. The final checkpoint was saved at epoch 100
when early stopping triggered (10 epochs without validation improvement).

---

## Pipeline Overview

```
Step 1 (R, local)     Prepare DHS data and compute direct estimates
                            |
Step 2 (GEE, cloud)   Export 476 Landsat 7 patches (224x224 RGB)
                            |
Step 3 (Kaggle, GPU)  Train ResNet-18, predict wealth for all clusters
                            |
Step 4 (R, local)     Fay-Herriot model with CNN predictions
                            |
Step 5 (R, local)     Generate maps, tables, and comparison figures
```

### Step 1: Prepare DHS Data

**Script**: `scripts/preparer_donnees.R`

Read the DHS 2012 household recode (NIHR61FL.DTA, 10,750 households),
spatially join 476 GPS cluster coordinates to department boundaries
(Admin 2, OCHA Niger), and compute survey-weighted direct estimates
for each department.

Produces:

- `data/processed/direct_estimates_dept.csv`: 64 departments with mean
  wealth, standard error, CV, and survey variance (psi).
- `data/processed/clusters_gps.csv`: 476 clusters with coordinates,
  department ID, and wealth mean.

### Step 2: Export Landsat 7 Patches from Google Earth Engine

**Script**: `scripts/telecharger_patches_gee.js`

Run this script in the GEE Code Editor (code.earthengine.google.com).
It iterates over all 476 DHS clusters and for each one:

1. Fetches all Landsat 7 scenes (2011-2012) within a 3,360 meter radius.
2. Filters out scenes with more than 30 percent cloud cover.
3. Computes a median composite across all remaining scenes (eliminates
   residual clouds and shadows).
4. Exports a 224 x 224 pixel RGB GeoTIFF to Google Drive.

The `reproject` call ensures the output is in WGS84 at 30 meter resolution.
Each patch covers roughly 6.7 by 6.7 kilometers on the ground.

### Step 3: Train ResNet-18 on Kaggle

**Script**: `scripts/entrainer_resnet.py`

Run on Kaggle with a GPU accelerator (T4 or P100). The input is the
476 patches converted from GeoTIFF to NumPy arrays (224 x 224 x 3,
float32, normalized to [0, 1]) and the labels from `labels.csv`.

Produces:

- `resnet18.pt`: trained model weights (44 MB).
- `cnn_predictions_cluster.csv`: predicted wealth for all 476 clusters.
- `cnn_metrics.csv`: test R squared and RMSE.
- `cnn_pred_vs_true.png` : scatter plot of predictions vs actual values.

### Step 4: Fay-Herriot Small Area Model

**Script**: `scripts/fay_herriot.R`

The Fay-Herriot model (Fay and Herriot, 1979) is a standard area-level
small area estimation method. It combines:

- **Direct estimate** (from DHS survey), unbiased but high variance for
  small departments.
- **Synthetic estimate** (from linear regression on auxiliary variables):
  potentially biased but low variance.

The model produces an Empirical Best Linear Unbiased Predictor (EBLUP)
as a weighted average:

```
EBLUP_i = gamma_i * direct_i + (1 - gamma_i) * synthetic_i
```

where gamma_i is the shrinkage factor, close to 1 when the direct estimate
is reliable (many households), close to 0 when the model dominates.

Four models are compared:

1. FH-0: intercept only (no auxiliary variables).
2. FH-1: urbanization rate (percent urban households per department).
3. FH-2: urbanization rate + NDVI + CNN score.
4. FH-3: urbanization rate + NDVI + CNN score + nightlights *(planned)*.

### Step 5: Outputs and Figures

| File                                         | Description                                  |
| -------------------------------------------- | -------------------------------------------- |
| `outputs/tables/fh_results_comparison.csv` | AIC, CV, gamma for all fitted models         |
| `outputs/tables/fh_results_details.csv`    | Per-department: direct, EBLUP, SE, gamma     |
| `outputs/figures/carte_wealth_direct.png`  | Map of direct survey estimates by department |
| `outputs/figures/carte_wealth_eblup.png`   | Map of EBLUP estimates (FH-2)                |
| `outputs/figures/comparaison_modeles.png`  | Bar chart comparing CV across models         |
| `outputs/figures/shrinkage_vs_taille.png`  | Gamma vs number of clusters per department   |
| `outputs/figures/cnn_pred_vs_true.png`     | ResNet-18 predictions vs actual wealth       |

---

## Project Structure

```
poverty_sae_satellite/
  scripts/               Production scripts
    preparer_donnees.R     DHS data preparation and direct estimates
    telecharger_modis.R    NDVI MODIS download
    telecharger_patches_gee.js  Landsat 7 export (GEE Code Editor)
    entrainer_resnet.py    ResNet-18 training (Kaggle)
    fay_herriot.R          Fay-Herriot models and figures

  exploration/           R Markdown exploratory analysis
    explorer_dhs.Rmd       Clusters, wealth distribution, sample sizes
    explorer_modis.Rmd     NDVI vs wealth, spatial correlation
    explorer_patches.Rmd   Landsat patch visualization (requires Python)

  data/
    raw/                   DHS 2012 data (restricted), Admin 2 boundaries
    processed/             Direct estimates, cluster GPS, CNN predictions,
                           NDVI, Landsat patches (476 .tif files)

  outputs/
    tables/                CSV files with model comparison and details
    figures/               Maps, scatter plots, and comparison charts
    models/                Trained ResNet-18 weights (regenerable)

  pipeline_kaggle.ipynb   Kaggle notebook for ResNet training
  poverty_sae_kaggle.zip  Kaggle dataset archive (patches + scripts)
  gee_export_patches.ipynb  Colab notebook for batch GEE export
  gee_clusters.csv        476 cluster coordinates for GEE upload
  requirements_kaggle.txt  Python dependencies for Kaggle
```

---

## Data Sources

| Source                  | Description                                              | Access                                     |
| ----------------------- | -------------------------------------------------------- | ------------------------------------------ |
| DHS 2012 (NIHR61FL)     | Household survey, 10,750 households, 476 clusters        | Restricted (DHS Program, registered users) |
| GPS clusters (NIGE61FL) | Cluster coordinates with random displacement             | Restricted (DHS Program)                   |
| Admin 2 Niger           | 67 department boundaries                                 | Public (HDX, OCHA, Creative Commons)       |
| Landsat 7 (GEE)         | Surface reflectance, 30m resolution, 2011-2012 composite | Public (USGS, via Google Earth Engine)     |
| MODIS MOD13A3           | Monthly NDVI at 1km, 2011-2012 average                   | Public (NASA, via MODISTools)              |
| VIIRS/DMSP-OLS          | Nightlights, 500m/1km, 2012 annual composite             | Public (NOAA, via Google Earth Engine)     |

---

## Notes on Reproducibility

- DHS data is confidential and must be requested from the DHS Program
  (dhsprogram.com). The raw .DTA and .SHP files are gitignored.
- The 476 Landsat patches are regenerable via GEE (the JavaScript and
  Colab notebook are provided).
- Three departments have no DHS clusters (desert areas: Bilma, etc.) and
  are excluded from estimation. This is standard for small area methods.

---

## References

- Fay, R. E., & Herriot, R. A. (1979). Estimates of income for small places:
  An application of James-Stein procedures to census data. *Journal of the
  American Statistical Association*, 74(366), 269-277.
- Jean, N., Burke, M., Xie, M., Davis, W. M., Lobell, D. B., & Ermon, S.
  (2016). Combining satellite imagery and machine learning to predict poverty.
  *Science*, 353(6301), 790-794.
- Yeh, C., Perez, A., Driscoll, A., Azzari, G., Tang, Z., Lobell, D., ...
  & Ermon, S. (2020). Using publicly available satellite imagery and deep
  learning to understand economic well-being in Africa. *Nature
  Communications*, 11, 2583.
- Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.).
  Wiley.
- He, K., Zhang, X., Ren, S., & Sun, J. (2016). Deep residual learning for
  image recognition. *CVPR 2016*.
- Edochie, I., Newhouse, D., & Cochinard, F. (2024). Small area estimation
  of poverty under high volatility. *World Bank Policy Research Working Paper*.
