# NOTEBOOK

Carnet de bord theorie et pratique du projet `poverty_sae_satellite/`.

Le but est de comprendre pourquoi et comment on estime la pauvrete au niveau regional au Niger a partir de l'enquete DHS, d'une variable satellite, et d'un modele statistique de type Fay-Herriot. Chaque section pose une question, donne l'intuition, et renvoie au script R ou Python qui la met en oeuvre.

---

## 1. Pourquoi le Small Area Estimation ?

### Le probleme

L'enquete DHS Niger 2012 a interroge 10 750 menages repartis dans 8 regions. La taille d'echantillon par region varie : 722 menages a Agadez, 1856 a Maradi, 1669 a Niamey. Quand on calcule la richesse moyenne par region, la precision depend directement de cette taille d'echantillon.

Les enquetes DHS sont conçues pour etre representatives au niveau national et au premier niveau administratif (Admin 1, c'est-a-dire les regions) [Wakefield et al., 2019; Wu et al., 2021]. En dessous de ce niveau, les tailles d'echantillon sont trop faibles pour produire des estimations directes fiables. Au Niger, Edochie et al. [2024] montrent que les CV medians des estimations directes au niveau commune dépassent 30%, ce qui les rend inutilisables pour la décision publique. Augmenter la taille d'echantillon couterait trop cher [The DHS Program, SAR15].

Le coefficient de variation (CV = SE / moyenne) est l'indicateur standard en statistique officielle. Un CV inferieur a 15-20% est generalement considere comme publiable. Au-dessus, l'estimation est trop imprecise pour guider une decision publique.

### L'intuition du SAE

Plutot que d'estimer chaque region independamment (estimateur direct), on "emprunte de la force" aux autres regions via un modele statistique [Rao & Molina, 2015]. L'idee : si deux regions se ressemblent (meme climat, meme urbanisation), leurs estimations peuvent s'informer mutuellement. Le modele de Fay-Herriot [Fay & Herriot, 1979] formalise cette intuition. Newhouse et al. [2024] ont montre que cette approche permet des gains de precision equivalents a multiplier la taille d'echantillon par 3 a 5.

### Sortie attendue

Apres les scripts R, on observe que les regions les moins echantillonnees beneficient le plus du modele : leur estimation finale combine leur propre donnee (peu fiable seule) avec la prediction du modele (plus stable car basee sur toutes les regions).

Voir `R/estimateurs_directs.R` pour les CV directs, `R/modele_fay_herriot.R` pour les EBLUP.

---

## 2. Le modele de Fay-Herriot (1979)

### Le modele

Pour une region i (i = 1, ..., 8) :

```
y_i  = x_i' beta + v_i + e_i

  y_i  : estimation directe de la richesse (connue, avec son SE)
  x_i  : vecteur de variables auxiliaires (connu, pas d'erreur)
  beta : coefficients de regression (inconnus, a estimer)
  v_i  : effet aleatoire regional, v_i ~ N(0, sigma2_v)
  e_i  : erreur d'echantillonnage, e_i ~ N(0, psi_i = SE_i^2)
```

Ce modele a ete introduit par Fay & Herriot [1979] pour estimer le revenu par habitant dans de petites zones aux Etats-Unis. Il est depuis devenu le modele de reference pour l'estimation sur petits domaines au niveau agregé [Rao & Molina, 2015]. Le U.S. Census Bureau l'utilise dans le programme SAIPE (Small Area Income and Poverty Estimates) pour allouer des fonds federaux [Bell et al., 2016].

Le terme cle est `v_i`. C'est lui qui modelise la "vraie" variabilite regionale une fois enleve l'effet lineaire des variables auxiliaires. Le modele suppose que les regions sont echangeables conditionnellement aux x_i.

### L'estimation

Les coefficients `beta` et la variance `sigma2_v` sont estimes par maximum de vraisemblance restreint (REML). En pratique, on utilise le package R `sae` avec la fonction `eblupFH()`. Voir `R/modele_fay_herriot.R`.

### Limitation

Avec seulement 8 regions, on ne peut pas inclure plus de 2 variables auxiliaires (regle de Harrell : n / 4 = 2). C'est pour cela que le projet complet utilise 4 modeles (M1, M2, M3, M4) pour explorer differents sous-ensembles. Dans cette version simplifiee, on utilise 1 seul modele a 2 predicteurs : `urban_pct` et `ndvi_mean`.

---

## 3. L'EBLUP et le facteur de shrinkage gamma

### Formule de l'EBLUP

Une fois `beta` et `sigma2_v` estimes, l'EBLUP (Empirical Best Linear Unbiased Predictor) pour la region i est :

```
EBLUP_i = gamma_i * y_i + (1 - gamma_i) * x_i' beta_chapeau

  gamma_i = sigma2_v / (sigma2_v + psi_i)  dans [0, 1]
```

L'EBLUP est une application de l'estimateur de James-Stein [James & Stein, 1961] au contexte des petits domaines [Fay & Herriot, 1979]. C'est le meilleur predicteur lineaire sans biais empirique [Rao & Molina, 2015].

### Interpretation de gamma

- `gamma_i` proche de 1 : la variance d'echantillonnage psi_i est faible devant la variance regionale sigma2_v. On fait confiance a l'estimateur direct `y_i`. L'EBLUP est proche de `y_i`.
- `gamma_i` proche de 0 : `psi_i` est grande (region mal echantillonnee). On se fie davantage a la prediction du modele `x_i' beta`. L'EBLUP est "tire" vers la droite de regression.

### Exemple typique

- Une region avec 1500 menages aura un psi_i faible, donc gamma proche de 1.
- Une region avec 200 menages aura psi_i eleve, donc gamma proche de 0, et son EBLUP sera fortement shrinkage vers la prediction du modele.

Voir `R/visualiser_resultats.R` pour le graphique du gamma par region.

---

## 4. Pourquoi ajouter une variable satellite ?

### Le probleme du recensement

Le modele FH classique utilise des variables auxiliaires issues du recensement (taux d'urbanisation, taux d'alphabétisation, etc.). Mais le recensement est conduit tous les 10 ans et devient rapidement obsolete. Au Niger, le dernier recensement general date de 2012, et le suivant est prevu en 2024-2025. Les images satellite, quant à elles, sont mises à jour tous les 5 jours (Sentinel-2) ou quotidiennement (MODIS) et couvrent 100% du territoire.

### L'apport du satellite : ce que dit la litterature

L'indice NDVI (Normalized Difference Vegetation Index) a ete introduit par Tucker [1979]. Il mesure la vigueur de la vegetation a partir de la reflectance dans le rouge et le proche infrarouge. De nombreuses etudes ont montre sa correlation avec la pauvrete :

- Sedda et al. [2015] ont montre que l'intensite de la pauvrete varie inversement avec le NDVI en Afrique de l'Ouest.
- Engstrom et al. [2022] ont trouve que le NDVI et l'indice de bâti (NDBI) extraits de Landsat expliquent 25-50% de la variation de l'activite economique au niveau communal au Vietnam.
- Watmough et al. [2019] ont predit les menages les plus pauvres au Kenya rural avec 62% de precision à partir d'images satellite.
- Xie et al. [2016] ont montre que le NDVI seul, combine a un CNN, predit la pauvrete aussi bien que des images haute resolution.
- Une etude de la Banque Mondiale [Barbier & Hochard, 2016] a montre qu'une augmentation de 10% de la vigueur de la vegetation est associee à une reduction de 0.7 point du taux de pauvrete en Afrique subsaharienne.

Les lumieres nocturnes (VIIRS) sont egalement un proxy classique de l'activite economique [Henderson et al., 2012; Chen & Nordhaus, 2011], mais elles montrent peu de variation dans les zones les plus pauvres [Jean et al., 2016].

### Dans ce projet

On utilise le NDVI moyen par region. Le modele FH est : `richesse ~ urban_pct + ndvi_mean`. Si le beta associe a `ndvi_mean` est negatif et significatif, cela confirme que les regions avec plus de vegetation (rural agricole) ont une richesse moyenne plus elevee que la moyenne conditionnelle a l'urbanisation.

Voir `R/estimateurs_directs.R` pour la construction de la variable NDVI regionale, et `R/modele_fay_herriot.R` pour son integration dans le modele.

---

## 5. Pourquoi un CNN sur Sentinel-2 ?

### Le probleme des agregats regionaux

Le NDVI moyen par region perd toute l'information spatiale intra-regionale. Deux regions peuvent avoir le meme NDVI moyen mais des structures differentes (une grande ville isolee + desert vs petites villes distribuees + savane). Le CNN peut, en principe, capturer ces patterns.

### L'idee

Pour chaque cluster DHS (476 au total), on dispose de :
- la richesse moyenne du cluster (label),
- les coordonnees GPS (centre du cluster).

On telecharge via Google Earth Engine un patch Sentinel-2 de 64x64 pixels (640 m x 640 m a 10 m de resolution) centre sur chaque cluster, sur 4 bandes (B2 bleu, B3 vert, B4 rouge, B8 proche infrarouge). On entraine un CNN a predire la richesse du cluster a partir de ce patch.

### Ce que dit la litterature sur les CNN + satellite

L'article fondateur de Jean, Burke, Xie, Davis, Lobell & Ermon [2016] dans *Science* a montre qu'un CNN entraine sur des images satellite peut expliquer jusqu'a 75% de la variation de la richesse au niveau cluster dans cinq pays africains. Leur approche utilise le transfert d'apprentissage : le reseau est d'abord pre-entraine sur ImageNet, puis fine-tune pour predire les lumieres nocturnes, et enfin les caracteristiques apprises sont utilisees dans une regression ridge.

Yeh et al. [2020] dans *Nature Communications* ont etendu cette approche à 20 000 villages africains, atteignant 70% de variance expliquee pour la richesse et 83% au niveau district. Leur modele utilise en entree des images Landsat multispectrales et des lumieres nocturnes.

Plus recemment, Pettersson et al. [2023] ont ajoute une dimension temporelle avec des couches LSTM, atteignant 72% de variance expliquee. Le framework DeepWealth [Tonneau et al., 2024] a montre une generalisation à 24 pays africains avec un R² de 0.69. Chi et al. [2022] ont produit des microestimations de richesse pour tous les pays à revenu faible et intermediaire.

### Architecture simplifiee

```
Input (4, 64, 64)
  -> Conv2D(32, 3) + ReLU + MaxPool(2)   -> (32, 31, 31)
  -> Conv2D(64, 3) + ReLU + MaxPool(2)   -> (64, 14, 14)
  -> Conv2D(128, 3) + ReLU + MaxPool(2)  -> (128, 6, 6)
  -> Flatten + Linear(256) + ReLU + Dropout(0.5) + Linear(1)
Loss MSE, optimiseur Adam lr=1e-3, 20 epochs, batch 32.
```

C'est volontairement simple (3 convolutions au lieu d'un ResNet-18) pour rester pedagogique. Voir `python/entrainer_cnn.py`.

### Hypothese testee

**H2** : Si le R² de test est strictement positif, alors l'imagerie satellite contient un signal predictif de la richesse. Jean et al. [2016] ont obtenu des R² de 0.37 à 0.55 selon les pays. Yeh et al. [2020] ont atteint 0.70 median. Avec 476 clusters et 80% d'entrainement (380 echantillons), on s'attend a un signal faible mais present.

### Visualisation Grad-CAM

Apres entrainement, on utilise Grad-CAM (Gradient-weighted Class Activation Mapping) [Selvaraju et al., 2017] pour visualiser les zones du patch qui activent le plus le neurone de sortie. Jean et al. [2016] ont montre que les filtres appris detectent des structures semantiques comme les routes, les zones urbaines et les terres agricoles, sans supervision directe. Voir `python/gradcam.py`.

---

## 6. Marche a suivre pas-a-pas

### Etape 1 : charger les donnees

```r
source("R/charger_dhs.R")
```

Ce script :
- charge `data/raw/NIHR61FL.DTA` via `haven::read_dta()`,
- selectionne 5 colonnes : HV001 (cluster), HV005 (poids), HV024 (region), HV025 (urbain/rural), HV271 (indice de richesse),
- renomme en `cluster`, `weight`, `region`, `urban`, `wealth`,
- cree un binaire `urban` (1 = urbain, 0 = rural) a partir de HV025,
- sauvegarde dans `outputs/data/dhs_subset.rds`.

Si le fichier DHS est introuvable, le script bascule en mode simulation calibree sur les valeurs connues du Niger (mean wealth ~ 1.7, std ~ 1.0, 8 regions, 476 clusters).

### Etape 2 : estimateurs directs

```r
source("R/estimateurs_directs.R")
```

Ce script :
- declare le plan de sondage `svydesign(ids = ~cluster, weights = ~weight, data = dhs)`,
- calcule la richesse moyenne par region via `svyby(~wealth, ~region, design, svymean)`,
- extrait la SE et calcule le CV (en %) et psi = SE^2,
- simule un NDVI regional (calibre par region : plus eleve en zone agricole Sahelienne, plus faible en zone desertique du Nord),
- sauvegarde la table fusionnee dans `outputs/data/sae_poverty_data.rds` et la table des estimateurs directs dans `outputs/tables/direct_estimates.csv`.

### Etape 3 : modele Fay-Herriot

```r
source("R/modele_fay_herriot.R")
```

Ce script :
- ajuste le modele `eblupFH(direct_mean ~ urban_pct + ndvi_mean, vardir = psi, data = sae_data)`,
- extrait les EBLUP, SE, et le facteur gamma de shrinkage,
- calcule le gain de CV (CV direct - CV EBLUP),
- sauvegarde dans `outputs/tables/fh_results.csv`.

### Etape 4 : visualisation

```r
source("R/visualiser_resultats.R")
```

Ce script genere 2 figures en 200 dpi :
- `outputs/figures/01_direct_vs_eblup.png` : graphique en points et barres d'erreur, direct vs EBLUP, avec intervalles de confiance 95%, regions ordonnees par richesse directe.
- `outputs/figures/02_shrinkage.png` : graphique en barres du facteur gamma par region, avec ligne de reference a gamma = 0.5.

### Etape 5 : pipeline deep learning

Avant tout, il faut telecharger les patches Sentinel-2 via Google Earth Engine :
- Creer un compte sur https://code.earthengine.google.com
- Utiliser le snippet dans `python/preparer_clusters.py` pour exporter 476 patches de 64x64 pixels en EPSG:4326
- Les placer dans `python/data/patches/` (un fichier `.npy` par cluster)

Ensuite :
```bash
python python/preparer_clusters.py    # joint DHS wealth + GPS
python python/entrainer_cnn.py        # entraine le CNN, sauvegarde metriques
python python/gradcam.py              # genere les cartes d'activation
```

---

## 7. Comment interpreter les resultats ?

### Si tout fonctionne

- Les 8 regions ont des richesses moyennes directes distinctes.
- Le facteur gamma est compris entre 0.5 et 1.0 pour toutes les regions (les 8 regions DHS sont toutes assez bien echantillonnees au niveau regional).
- Le gain de CV (CV direct - CV EBLUP) peut etre faible, voire negatif, ce qui suggere qu'au niveau regional, le SAE ajoute peu. La vraie valeur est au niveau departement (n ~ 30-50 menages, CV > 25%), comme le montrent Edochie et al. [2024] au Niger.
- Le CNN devrait avoir un R² de test > 0 sur les 96 clusters de test.

### Limites a reconnaitre

- 8 regions est un tres petit echantillon pour un modele a 2 predicteurs. Les resultats sont indicatifs, pas definitifs.
- Le NDVI est simule dans cette version, pas telecharge. Pour un vrai travail, utiliser MODIS MOD13A3 ou Google Earth Engine.
- Le CNN ne pourra pas predire la richesse avec une grande precision. C'est attendu et documente dans la litterature : avec ~380 echantillons d'entrainement, la performance sera inferieure aux 20 000 echantillons de Yeh et al. [2020].

---

## References

- Barbier, E.B. & Hochard, J.P. (2016). Does land degradation increase poverty in developing countries? *PLoS ONE*, 11(5), e0152973.
- Battese, G.E., Harter, R.M., & Fuller, W.A. (1988). An error-components model for prediction of county crop areas using survey and satellite data. *JASA*, 83(401), 28-36.
- Bell, W.R., Basel, W.W., & Maples, J.J. (2016). An overview of the U.S. Census Bureau's SAIPE program. In *Analysis of Poverty Data by Small Area Estimation*. Wiley.
- Chen, X. & Nordhaus, W.D. (2011). Using luminosity data as a proxy for economic statistics. *PNAS*, 108(21), 8589-8594.
- Chi, G., Fang, H., Chatterjee, S., & Blumenstock, J.E. (2022). Microestimates of wealth for all low- and middle-income countries. *PNAS*, 119(3), e2113658119.
- The DHS Program (2019). The DHS Program Modeled Map Surfaces. *Spatial Analysis Reports No. 15*. ICF.
- Edochie, I., Newhouse, D., et al. (2024). Small Area Estimation of Poverty in Four West African Countries by Integrating Survey and Geospatial Data. *World Bank Policy Research Working Paper*.
- Engstrom, R., Hersh, J., & Newhouse, D. (2022). Can Medium-Resolution Satellite Imagery Measure Economic Activity at Small Geographies? *World Bank Policy Research Working Paper*.
- Fay, R.E. & Herriot, R.A. (1979). Estimates of income for small places. *JASA*, 74(366), 269-277.
- Henderson, J.V., Storeygard, A., & Weil, D.N. (2012). Measuring economic growth from outer space. *AER*, 102(2), 994-1028.
- James, W. & Stein, C. (1961). Estimation with quadratic loss. *Proc. 4th Berkeley Symp.*, 1, 361-379.
- Jean, N., Burke, M., Xie, M., Davis, W.M., Lobell, D.B., & Ermon, S. (2016). Combining satellite imagery and machine learning to predict poverty. *Science*, 353(6301), 790-794.
- Newhouse, D., et al. (2024). Small Area Estimation of Non-Monetary Poverty with Geospatial Data. *World Bank Policy Research Working Paper*.
- Pettersson, M.B., Kakooei, M., Ortheden, J., Johansson, F.D., & Daoud, A. (2023). Time Series of Satellite Imagery Improve Deep Learning Estimates of Neighborhood-Level Poverty in Africa.
- Rao, J.N.K. & Molina, I. (2015). *Small Area Estimation*. 2nd ed. Wiley.
- Sedda, L., et al. (2015). Poverty and malaria: NDVI-based spatial analysis in West Africa.
- Selvaraju, R.R., Cogswell, M., Das, A., Vedantam, R., Parikh, D., & Batra, D. (2017). Grad-CAM: Visual explanations from deep networks via gradient-based localization. *ICCV*.
- Tonneau, M., et al. (2024). DeepWealth: A generalizable open-source deep learning framework using satellite images for well-being estimation. *SoftwareX*, 26, 101700.
- Tucker, C.J. (1979). Red and photographic infrared linear combinations for monitoring vegetation. *Remote Sensing of Environment*, 8(2), 127-150.
- Wakefield, J., Fuglstad, G.-A., Riebler, A., Godwin, J., Wilson, K., & Clark, S.J. (2019). Estimating under-five mortality in space and time in sub-Saharan Africa. *arXiv:1910.06512*.
- Watmough, G.R., et al. (2019). Socioecologically informed use of remote sensing data to predict rural household poverty. *PNAS*, 116(4), 1213-1218.
- Wu, Y., Li, Z.R., Mayala, B.K., et al. (2021). Spatial Modeling for Subnational Administrative Level 2 Small-Area Estimation. *DHS Spatial Analysis Reports No. 21*. ICF.
- Xie, M., Jean, N., Burke, M., Lobell, D.B., & Ermon, S. (2016). Transfer learning from deep features for remote sensing and poverty mapping. *AAAI*.
- Yeh, C., Perez, A., Driscoll, A., Azzari, G., Tang, Z., Lobell, D., Ermon, S., & Burke, M. (2020). Using publicly available satellite imagery and deep learning to understand economic well-being in Africa. *Nature Communications*, 11, 2583.
