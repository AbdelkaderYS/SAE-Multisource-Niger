# ===========================================================================
# preparer_clusters.py
#
# Joint la richesse moyenne par cluster (depuis le sous-ensemble DHS sauve
# par charger_dhs.R) avec les coordonnees GPS (NIGE61FL.shp) pour preparer
# la table d'entree du CNN.
#
# Sortie : python/data/clusters_wealth_gps.csv
#           colonnes : cluster, n_hh, wealth_mean, LATNUM, LONGNUM, region
#
# Etape ulterieure (non automatisee ici) : utiliser Google Earth Engine pour
# telecharger un patch Sentinel-2 (64x64 pixels, 4 bandes B2/B3/B4/B8) centre
# sur chaque cluster, et stocker dans python/data/patches/<cluster>.npy
# ===========================================================================

import os
import sys

import numpy as np
import pandas as pd

try:
    import geopandas as gpd
except ImportError:
    sys.exit("Le package geopandas est requis : pip install geopandas")


def main():
    os.makedirs("python/data", exist_ok=True)
    os.makedirs("python/data/patches", exist_ok=True)

    # 1) Charger les coordonnees GPS des clusters DHS
    shp_path = "data/raw/NIGE61FL.shp"
    if not os.path.exists(shp_path):
        sys.exit(f"Fichier introuvable : {shp_path}\n"
                 "Creer le lien symbolique depuis le projet parent.")
    gps = gpd.read_file(shp_path)
    print(f"Clusters GPS charges : {len(gps)}")
    print("Colonnes GPS :", list(gps.columns))

    # Garder les colonnes strictement necessaires
    keep_cols = ["DHSCLUST", "LATNUM", "LONGNUM", "DHSREGNA", "URBAN_RURA"]
    gps_keep = gps[keep_cols].rename(columns={
        "DHSCLUST": "cluster",
        "LATNUM":   "lat",
        "LONGNUM":  "lon",
        "DHSREGNA": "region",
        "URBAN_RURA": "urban_rura"
    }).copy()
    gps_keep["cluster"] = gps_keep["cluster"].astype(int)

    # 2) Charger le sous-ensemble DHS sauvegarde par R
    dhs_rds = "outputs/data/dhs_subset.rds"
    if not os.path.exists(dhs_rds):
        sys.exit(f"Fichier introuvable : {dhs_rds}\n"
                 "Lancer d'abord R/charger_dhs.R")
    dhs = pd.read_csv(dhs_rds) if False else None  # placeholder
    # Le RDS est un format R ; on passe par pyreadr si dispo, sinon on
    # recree la table minimale a partir d'un export CSV.
    dhs_csv = "python/data/dhs_subset.csv"
    if not os.path.exists(dhs_csv):
        sys.exit(f"Export CSV introuvable : {dhs_csv}\n"
                 "Lancer : Rscript -e 'readRDS(\"outputs/data/dhs_subset.rds\") |> "
                 "write.csv(\"python/data/dhs_subset.csv\", row.names=FALSE)'")
    dhs = pd.read_csv(dhs_csv)
    print(f"DHS subset charge : {len(dhs)} menages")

    # 3) Richesse moyenne et taille par cluster
    cluster_stats = (dhs.groupby("cluster")
                       .agg(wealth_mean=("wealth", "mean"),
                            n_hh=("wealth", "size"))
                       .reset_index())
    cluster_stats["cluster"] = cluster_stats["cluster"].astype(int)

    # 4) Jointure avec GPS
    merged = gps_keep.merge(cluster_stats, on="cluster", how="inner")
    print(f"Clusters apres jointure : {len(merged)} (sur {len(gps_keep)} GPS)")

    # 5) Sauvegarde
    out_csv = "python/data/clusters_wealth_gps.csv"
    merged.to_csv(out_csv, index=False)
    print(f"Sauvegarde : {out_csv}")
    print(merged.head())

    # 6) Note sur les patches Sentinel-2
    n = len(merged)
    print("\n--- Etape suivante (manuelle) ---")
    print(f"Pour {n} clusters, utiliser Google Earth Engine pour extraire un")
    print("patch Sentinel-2 L2A de 64x64 pixels (4 bandes B2/B3/B4/B8),")
    print("median composite sans nuages, en 2019-2020 (plus proche de la DHS 2012")
    print("parmi les archives Sentinel-2 disponibles).")
    print("Sauvegarder un fichier .npy par cluster dans python/data/patches/")
    print("Snippet GEE a executer dans https://code.earthengine.google.com :")
    print("""
    var table = ee.FeatureCollection('users/your_user/clusters');
    var s2 = ee.ImageCollection('COPERNICUS/S2_SR')
              .filterBounds(table)
              .filterDate('2019-01-01', '2020-12-31')
              .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));
    function exportPatch(feat) {
      var geom = feat.geometry();
      var img  = s2.median().clip(geom.buffer(320));
      var rgb  = img.select(['B2','B3','B4','B8']).toFloat();
      Export.image.toDrive({
        image: rgb, description: 'patch_' + feat.get('cluster'),
        scale: 10, region: geom.buffer(320), maxPixels: 1e7
      });
    }
    table.toList(500).map(exportPatch);
    """)


if __name__ == "__main__":
    main()
