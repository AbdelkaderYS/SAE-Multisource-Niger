# NOTEBOOK

Carnet de bord theorie et pratique du projet `poverty_sae_satellite/`.

Le but est de comprendre pourquoi et comment on estime la pauvrete au niveau regional au Niger a partir de l'enquete DHS, d'une variable satellite, et d'un modele statistique de type Fay-Herriot. Chaque section pose une question, donne l'intuition, et renvoie au script R ou Python qui la met en oeuvre.

---

## 1. Pourquoi le Small Area Estimation ?

### Le probleme

L'enquete DHS Niger 2012 a interroge 10 750 menages repartis dans 8 regions. La taille d'echantillon par region varie : 722 menages a Agadez, 1856 a Maradi, 1669 a Niamey. Quand on calcule la richesse moyenne par region, la precision depend directement de cette taille d'echantillon.

Le coefficient de variation (CV = SE / moyenne) est l'indicateur standard en statistique officielle. Un CV inferieur a 15-20% est generalement considere comme publiable. Au-dessus, l'estimation est trop imprecise pour guider une decision publique.

### L'intuition du SAE

Plutot que d'estimer chaque region independamment (estimateur direct), on "emprunte de la force" aux autres regions via un modele statistique. L'idee : si deux regions se ressemblent (meme climat, meme urbanisation), leurs estimations peuvent s'informer mutuellement. Le modele de Fay-Herriot formalise cette intuition.

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

### Interpetition de gamma

- `gamma_i` proche de 1 (proche de 1) : la variance d'echantillonnage psi_i est faible devant la variance regionale sigma2_v. On fait confiance a l'estimateur direct `y_i`. L'EBLUP est proche de `y_i`.
- `gamma_i` proche de 0 : `psi_i` est grande (region mal echantillonnee). On se fie davantage a la prediction du modele `x_i' beta`. L'EBLUP est "tire" vers la droite de regression.

### Exemple typique

- Une region avec 1500 menages aura un psi_i faible, donc gamma proche de 1.
- Une region avec 200 menages aura psi_i eleve, donc gamma proche de 0, et son EBLUP sera fortement shrinkage vers la prediction du modele.

Voir `R/visualiser_resultats.R` pour le graphique du gamma par region.

---

## 4. Pourquoi ajouter une variable satellite ?

### Le probleme du recensement

Le modele FH classique utilise des variables auxiliaires issues du recensement (taux d'urbanisation, taux d'alphabétisation, etc.). Mais le recensement est conduit tous les 10 ans et devient rapidement obsolete. Au Niger, le dernier recensement general date de 2012, et le suivant est prevu en 2024-2025.

### L'apport du satellite

Les images satellite (Sentinel-2, MODIS, VIIRS) sont mises a jour tous les 5 jours (Sentinel-2) ou quotidiennement (MODIS). Elles sont gratuites et couvrent 100% du territoire. L'indice NDVI (Normalized Difference Vegetation Index) est un proxy classique de l'activite agricole, qui est elle-meme correlee a la richesse rurale. Les lumieres nocturnes VIIRS sont un proxy d'electrification et d'activite economique.

### Dans ce projet

On utilise le NDVI moyen par region (simule pour des raisons pedagogiques, mais la logique est la meme avec des donnees reelles Google Earth Engine). Le modele FH est : `richesse ~ urban_pct + ndvi_mean`. Si le beta associe a `ndvi_mean` est negatif et significatif, cela confirme que les regions avec plus de vegetation (rural agricole) ont une richesse moyenne plus elevee que la moyenne conditionnelle a l'urbanisation. C'est pedagogique et discutable, mais la methodologie est rigoureuse.

Voir `R/estimateurs_directs.R` pour la construction de la variable NDVI regionale, et `R/modele_fay_herriot.R` pour son integration dans le modele.

---

## 5. Pourquoi un CNN sur Sentinel-2 ?

### Le probleme des agregats regionaux

Le NDVI moyen par region perd toute l'information spatiale intra-regionale. Deux regions peuvent avoir le meme NDVI moyen mais des structures differentes (une grande ville isolee + desert vs petites villes distribuees + savane). Le CNN peut, en principe, capturer ces patterns.

### L'idee

Pour chaque cluster DHS (476 au total), on dispose de :
- la richesse moyenne du cluster (label),
- les coordonnees GPS (centre du cluster).

On telecharge via Google Earth Engine un patch Sentinel-2 de 64x64 pixels (640 m x 640 m a 10 m de resolution) centred sur chaque cluster, sur 4 bandes (B2 bleu, B3 vert, B4 rouge, B8 proche infrarouge). On entraine un CNN a predire la richesse du cluster a partir de ce patch.

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

H2 : Si le R2 de test est strictement positif, alors l'imagerie satellite contient un signal predictif de la richesse. Honnetement, avec 476 clusters, 80% d'entrainement (380), on s'attend a un signal faible mais present.

### Visualisation Grad-CAM

Apres entrainement, on utilise Grad-CAM (Gradient-weighted Class Activation Mapping) pour visualiser les zones du patch qui activent le plus le neurone de sortie. Voir `python/gradcam.py`.

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

- Les 8 regions ont des richesse moyennes directes distinctes.
- Le facteur gamma est compris entre 0.5 et 1.0 pour toutes les regions (les 8 regions DHS sont toutes assez bien echantillonnees au niveau regional).
- Le gain de CV (CV direct - CV EBLUP) peut etre faible, voire negatif, ce qui suggerere qu'au niveau regional, le SAE ajoute peu. La vraie valeur est au niveau departement (n ~ 30-50 menages, CV > 25%).
- Le CNN devrait avoir un R2 de test > 0 sur les 96 clusters de test.

### Limites a reconnaitre

- 8 regions est un tres petit echantillon pour un modele a 2 predicteurs. Les resultats sont indicatifs, pas definitifs.
- Le NDVI est simule dans cette version, pas telecharge. Pour un vrai travail, utiliser MODIS MOD13A3 ou Google Earth Engine.
- Le CNN ne pourra pas predire la richesse avec une grande precision. C'est attendu et documente.

### Pour aller plus loin

Une fois cette version comprise, le projet complet `project2_SAE_Multisource_Niger/` ajoute :
- le calcul du MPI Alkire-Foster (10 indicateurs),
- le record linkage entre bases DHS,
- la correction de biais IPW pour les sources non-probabilistes,
- 4 modeles FH au lieu d'un seul,
- differential privacy,
- un ResNet-18 avec Grad-CAM complet.

Mais la logique de base est la meme. Une fois FH compris, le reste est de l'extension.
