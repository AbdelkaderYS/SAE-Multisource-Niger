"""
Téléchargement des patches Landsat 7 pour chaque cluster DHS.
Composite sans nuage 2011-2012, bandes RGB, 224x224 pixels.

Usage:
    # Mode test (50 premiers clusters)
    python scripts/telecharger_landsat.py --test

    # Mode complet (476 clusters)
    python scripts/telecharger_landsat.py

Ce script est conçu pour tourner sur Kaggle (Python 3.10+).
Paquets requis : pystac-client, rasterio, numpy, pandas
"""

import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_true",
                        help="Télécharger seulement 50 clusters pour test")
    args = parser.parse_args()

    CLUSTERS_CSV = "data/processed/clusters_gps.csv"
    OUT_DIR = "data/processed/patches_landsat"
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Charger les clusters
    clusters = pd.read_csv(CLUSTERS_CSV)
    clusters = clusters.dropna(subset=["LONGNUM", "LATNUM"]).copy()
    clusters["cluster"] = clusters["cluster"].astype(int)
    print(f"Clusters chargés : {len(clusters)}")

    if args.test:
        clusters = clusters.head(50)
        print(f"Mode TEST : {len(clusters)} clusters")

    # 2. Vérifier les clusters déjà téléchargés
    existing = set()
    for f in os.listdir(OUT_DIR):
        if f.endswith(".npy"):
            existing.add(int(f.replace(".npy", "")))

    clusters = clusters[~clusters["cluster"].isin(existing)]
    print(f"Clusters à télécharger : {len(clusters)}")

    if len(clusters) == 0:
        print("Tous les patches sont déjà présents.")
        return

    # 3. Importer STAC (fait ici pour que l'import échoue seulement si on exécute)
    try:
        import pystac_client
        import rasterio
        from rasterio.warp import transform_bounds
        try:
            import planetary_computer
            HAS_PC_SIGN = True
        except ImportError:
            HAS_PC_SIGN = False
    except ImportError:
        sys.exit("Paquets manquants : pip install pystac-client rasterio")

    # Planetary Computer uniquement (USGS LandsatLook ne liste pas les collections Landsat 7)
    CATALOG_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
    COLLECTION = "landsat-c2-l2"

    try:
        catalog = pystac_client.Client.open(CATALOG_URL)
        # Valider que le catalogue répond
        catalog.search(collections=[COLLECTION], max_items=1).items()
    except Exception as e:
        sys.exit(f"Impossible d'ouvrir le catalogue STAC : {e}")

    collection = COLLECTION
    filter_l7 = True
    use_pc_sign = True
    print(f"Catalogue STAC : {CATALOG_URL}")

    # Bandes Landsat 7 : RGB
    BANDS = ["red", "green", "blue"]

    # Taille : 224x224 pixels à 30m = 6.72 km
    # Buffer en degrés (~0.03° à 15°N) pour couvrir ±5 km de déplacement DHS
    BUFFER_DEG = 0.03

    success = 0
    for _, row in clusters.iterrows():
        cid = int(row["cluster"])
        lon, lat = float(row["LONGNUM"]), float(row["LATNUM"])

        out_path = os.path.join(OUT_DIR, f"{cid}.npy")
        if os.path.exists(out_path):
            continue

        bbox = [lon - BUFFER_DEG, lat - BUFFER_DEG,
                lon + BUFFER_DEG, lat + BUFFER_DEG]

        try:
            search = catalog.search(
                collections=[collection],
                bbox=bbox,
                datetime="2011-01-01/2012-12-31",
                query={"eo:cloud_cover": {"lt": 30}}
            )
            all_items = list(search.items())
            if filter_l7:
                # Garder seulement Landsat 7
                items = [i for i in all_items
                         if i.properties.get("platform") == "landsat-7"
                         or i.properties.get("landsat:product_id", "").startswith("LE07")]
            else:
                items = all_items
        except Exception as e:
            print(f"  Cluster {cid}: erreur STAC - {e}")
            items = []

        if len(items) == 0:
            print(f"  Cluster {cid}: aucune image Landsat 7 trouvée")
            continue

        # Signer les URLs si Planetary Computer
        if use_pc_sign and HAS_PC_SIGN:
            items = [planetary_computer.sign(i) for i in items]

        # Composite median
        patches = []
        for item in items[:20]:  # max 20 scenes par cluster
            try:
                with rasterio.open(item.assets["red"].href) as src:
                    # Convertir bbox WGS84 → CRS de l'image (UTM)
                    left, bottom, right, top = transform_bounds(
                        "EPSG:4326", src.crs, *bbox
                    )
                    window = src.window(left, bottom, right, top)
                    # Lire les 3 bandes RGB
                    rgb = []
                    for band in BANDS:
                        with rasterio.open(item.assets[band].href) as src_b:
                            arr = src_b.read(1, window=window,
                                             boundless=True, fill_value=0)
                        rgb.append(arr)
                    rgb = np.stack(rgb, axis=-1).astype(np.float32)
                    patches.append(rgb)
            except Exception:
                continue

        if len(patches) == 0:
            print(f"  Cluster {cid}: échec de lecture")
            continue

        patch = np.median(np.stack(patches, axis=0), axis=0)

        # Ajuster à 224x224 (crop si plus grand, pad si plus petit)
        h, w = patch.shape[:2]
        if h != 224 or w != 224:
            target = 224
            if h > target:
                y0 = (h - target) // 2
                patch = patch[y0:y0 + target, :, :]
            elif h < target:
                pad_h = target - h
                patch = np.pad(patch, ((pad_h // 2, pad_h - pad_h // 2),
                                         (0, 0), (0, 0)), mode="reflect")
            h, w = patch.shape[:2]
            if w > target:
                x0 = (w - target) // 2
                patch = patch[:, x0:x0 + target, :]
            elif w < target:
                pad_w = target - w
                patch = np.pad(patch, ((0, 0),
                                         (pad_w // 2, pad_w - pad_w // 2),
                                         (0, 0)), mode="reflect")

        # Normalisation percentile
        p2, p98 = np.percentile(patch, [2, 98])
        patch = np.clip((patch - p2) / (p98 - p2 + 1e-8), 0, 1)

        np.save(out_path, patch.astype(np.float32))
        success += 1
        print(f"  Cluster {cid}: OK ({success}/{len(clusters)})")

    print(f"\nTéléchargement terminé : {success}/{len(clusters)} patches")
    print(f"Patches dans : {OUT_DIR}/")

    # 4. Sauvegarder les labels
    labels = clusters[["cluster", "wealth_mean", "n_hh", "dept_id", "dept_name"]].copy()
    labels.to_csv(os.path.join(OUT_DIR, "labels.csv"), index=False)
    print(f"Labels sauvegardés : {OUT_DIR}/labels.csv")


if __name__ == "__main__":
    main()
