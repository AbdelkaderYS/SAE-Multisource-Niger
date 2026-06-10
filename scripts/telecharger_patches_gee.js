/**
 * telecharger_patches_gee.js
 * ==========================
 * 
 * À EXÉCUTER DANS : https://code.earthengine.google.com
 * 
 * Ce script extrait des patches satellite Landsat 7 (RGB 224×224 pixels)
 * pour les 476 clusters DHS du Niger, et les exporte vers Google Drive.
 * 
 * Pourquoi GEE plutôt que STAC/rasterio en local ?
 *   - GEE est accessible depuis n'importe quelle connexion internet
 *   - GEE contient TOUTE l'archive Landsat 7 (pas besoin de chercher le bon catalogue)
 *   - GEE gère le cloud masking, le compositing, et l'export automatiquement
 *   - Pas de dépendances Python lourdes (torch, rasterio, etc.)
 * 
 * Ce que ça produit :
 *   476 fichiers GeoTIFF (un par cluster) dans Google Drive
 *   Format : 224×224 pixels, 3 bandes (RVB), 30m résolution
 * 
 * Une fois téléchargés, ces patches serviront à :
 *   - Entraîner un ResNet-18 (sur Kaggle) pour prédire la richesse
 *   - Extraire des features CNN pour le modèle Fay-Herriot
 *   - Comparer avec le NDVI MODIS comme variable auxiliaire
 */


// ============================================================
// 0. CONFIGURATION — À MODIFIER
// ============================================================

// Nom complet de la table uploadée dans GEE Assets
// Après avoir uploadé gee_clusters.csv, le nom ressemble à :
//   "projects/ee-votrecompte/assets/niger_dhs_clusters"
var ASSET_ID = 'projects/ee-votrecompte/assets/niger_dhs_clusters';

// Dossier de sortie dans Google Drive
var DRIVE_FOLDER = 'poverty_sae_patches';

// Paramètres des patches
var SCALE     = 30;      // Résolution Landsat (mètres/pixel)
var RADIUS    = 3360;    // Buffer autour du cluster (mètres)
var START     = '2011-01-01';
var END       = '2012-12-31';
var MAX_CLOUD = 30;      // Couverture nuageuse max (%)


// ============================================================
// 1. CHARGER LES CLUSTERS
// ============================================================

var clusters = ee.FeatureCollection(ASSET_ID);
print('Nombre de clusters chargés :', clusters.size());


// ============================================================
// 2. FONCTION : EXTRAIRE UN PATCH LANDSAT 7 POUR UN CLUSTER
// ============================================================

function getPatch(feature) {
  var cid = feature.get('cluster');
  var pt  = feature.geometry();
  var box = pt.buffer(RADIUS).bounds();

  // Collection Landsat 7 Level-2 (Surface Reflectance)
  // Filtres : période, zone, nuages
  var col = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterDate(START, END)
    .filterBounds(box)
    .filter(ee.Filter.lt('CLOUD_COVER', MAX_CLOUD))
    .select(['SR_B3', 'SR_B2', 'SR_B1']);  // Rouge, Vert, Bleu

  var nbScenes = col.size();

  // Composite median : pour chaque pixel, on prend la valeur médiane
  // de toutes les scènes disponibles. Ça élimine les nuages résiduels,
  // les ombres, et les artefacts.
  var composite = col.median()
    .clip(box)
    .reproject('EPSG:4326', null, SCALE);

  return composite.set({
    'cluster':      cid,
    'scene_count':  nbScenes,
    'lon':          feature.get('longitude'),
    'lat':          feature.get('latitude')
  });
}


// ============================================================
// 3. EXPORTER CHAQUE PATCH VERS GOOGLE DRIVE
// ============================================================

// Convertir la FeatureCollection en liste pour pouvoir itérer
var features = clusters.toList(clusters.size());
var n = features.size();
print('Export de', n, 'patches en cours...');

// Boucle sur tous les clusters
for (var i = 0; i < n.toInt(); i++) {
  var f     = ee.Feature(features.get(i));
  var patch = getPatch(f);
  var cid   = ee.Number(f.get('cluster')).format('%03d');

  // Export de l'image vers Google Drive
  Export.image.toDrive({
    image:         patch,
    description:   'cluster_' + cid,
    fileNamePrefix: 'cluster_' + cid,
    folder:        DRIVE_FOLDER,
    scale:         SCALE,
    crs:           'EPSG:4326',
    maxPixels:     1e9,
    fileDimensions: [224, 224],
    fileFormat:    'GEO_TIFF'
  });
}

print('Tâches d\'export créées. Va dans l\'onglet Tasks → Run All.');
