# poverty_sae_satellite

Estimation de la pauvrete au niveau regional au Niger en combinant l'enquete DHS 2012, des variables satellite (NDVI) et un modele de deep learning (CNN) sur Sentinel-2, dans un cadre Small Area Estimation (modele de Fay-Herriot).

---

## Motivation de l'etude

Au Niger, les enquetes menage de type Demographic and Health Survey (DHS) sont conduit par l'Institut National de la Statistique (INS) tous les 5 a 10 ans. Elles produisent des estimations precises au niveau national mais leur precision chute fortement au niveau regional et inferieur, car les echantillons par region sont trop petits (entre 700 et 1900 menages sur 8 regions). Or, pour le suivi de l'ODD 1 (elimination de la pauvrete) et pour la planification territoriale, les decideurs ont besoin d'estimations sous-nationales frequentes et precises.

L'imagerie satellite (Sentinel-2, MODIS, VIIRS) est gratuite, mise a jour regulierement, et couvre 100% du territoire. Le deep learning permet d'extraire de ces images des signaux sur la richesse, l'urbanisation ou l'activite economique que les agregats regionaux perdent. Le modele de Fay-Herriot (1979) fournit le cadre statistique pour combiner une estimation directe (l'enquete) avec des variables auxiliaires (satellite) afin de "preter de la force" entre les regions et reduire la variance des estimations.

C'est ce triptyque enquete + satellite + deep learning, integre dans un cadre statistique rigoureux, que ce projet met en oeuvre sur le Niger avec la DHS 2012.

---

## Questions de recherche

1. L'ajout de variables auxiliaires satextractllites (NDVI moyen, taux d'urbanisation) dans un modele de Fay-Herriot reduit-il la variance des estimations regionales de richesse par rapport aux estimateurs directs de l'enquete ?
2. Un reseau de neurones convolutif (CNN) entraine sur des patches Sentinel-2 est-il capable de predire la richesse moyenne d'un cluster DHS a partir de la seule imagerie satellite ?

---

## Hypotheses

- H1 : L'estimateur EBLUP du modele Fay-Herriot, qui combine estimateur direct et prediction par regression, a un coefficient de variation (CV) strictement inferieur a celui de l'estimateur direct pour les regions de petite taille.
- H2 : Le CNN atteint un R carre de test strictement positif, ce qui demontre que l'imagerie satellite contient un signal predictif de la richesse au niveau cluster.
- H3 : La combinaison des deux approches (FH regional + CNN cluster agrege) donne une estimation plus stable que l'enquete seule.

---

## Donnees

| Source | Fichier | Utilisation |
|---|---|---|
| DHS Niger 2012 Household Recode | `data/raw/NIHR61FL.DTA` | Indicateur de richesse (HV271), region (HV024), urbain/rural (HV025), poids (HV005) |
| DHS GPS clusters | `data/raw/NIGE61FL.shp` | Coordonnees (LATNUM, LONGNUM) des 476 clusters pour le pipeline Python |
| NDVI regional | simule | Proxy satellite agrege par region |
| Sentinel-2 | a telecharger via GEE | Patches 64x64 pixels, 4 bandes (B2, B3, B4, B8) |

Les fichiers DHS sont lies par symbole depuis le projet parent.

---

## Structure du projet

```
poverty_sae_satellite/
|-- README.md                ce fichier
|-- NOTEBOOK.md              theorie FH + guide pas-a-pas
|-- R/
|   |-- charger_dhs.R        charger le DHS, preparer les variables
|   |-- estimateurs_directs.R  svyby par region
|   |-- modele_fay_herriot.R   eblupFH, EBLUP/SE/gamma
|   |-- visualiser_resultats.R 2 figures publication
|-- python/
|   |-- preparer_clusters.py joindre DHS wealth + GPS
|   |-- entrainer_cnn.py     CNN 3 convolutions, RMSE/R2
|   |-- gradcam.py           visualisation des activations
|-- data/raw/                liens symboliques vers les fichiers DHS
|-- outputs/
    |-- data/                RDS intermediaires
    |-- tables/              CSV des resultats
    |-- figures/             PNG des figures
```

---

## Prerequis

Packages R (testes en version 4.3.3) :
```
install.packages(c("haven", "survey", "sae", "dplyr", "ggplot2", "tidyr"))
```

Packages Python (optionnels, pour la partie deep learning) :
```
pip install torch torchvision pandas numpy geopandas scikit-learn matplotlib
```

L'acces a Google Earth Engine est necessaire pour telecharger les patches Sentinel-2 (compte gratuit sur code.earthengine.google.com).

---

## Comment lancer

Pipeline statistique R (bout en bout, ~30 secondes) :
```r
setwd("poverty_sae_satellite")
source("R/charger_dhs.R")
source("R/estimateurs_directs.R")
source("R/modele_fay_herriot.R")
source("R/visualiser_resultats.R")
```

Pipeline deep learning Python (apres telechargement des patches via GEE) :
```bash
cd poverty_sae_satellite
python python/preparer_clusters.py
python python/entrainer_cnn.py
python python/gradcam.py
```

---

## Sorties attendues

| Fichier | Contenu |
|---|---|
| `outputs/tables/direct_estimates.csv` | Richesse moyenne et CV direct par region |
| `outputs/tables/fh_results.csv` | EBLUP, SE, gamma par region apres Fay-Herriot |
| `outputs/tables/cnn_metrics.csv` | RMSE et R2 du CNN sur le test set |
| `outputs/figures/01_direct_vs_eblup.png` | Comparaison directe vs EBLUP avec IC 95% |
| `outputs/figures/02_shrinkage.png` | Facteur de shrinkage gamma par region |
| `outputs/figures/03_cnn_pred_vs_true.png` | Predictions CNN vs verite terrain (test) |
| `outputs/figures/04_gradcam_*.png` | Cartes d'activation Grad-CAM |

---

## References

- Fay, R.E. et Herriot, R.A. (1979). Estimates of income for small places. *JASA* 74(366), 269-277.
- Rao, J.N.K. et Molina, I. (2015). *Small Area Estimation*. Wiley, 2e edition.
- OPHI et UNDP (2025). *Global MPI 2025*. OPHI Methodological Note 61.
- DHS Program (2013). *Enquete Demographique et de Sante au Niger 2012*. INS Niger.
- ESA (2015). Sentinel-2 User Handbook. European Space Agency.

---

## Auteur

Projet derive du cadre multisource complet `project2_SAE_Multisource_Niger/`, simplifie pour un usage pedagogique.
