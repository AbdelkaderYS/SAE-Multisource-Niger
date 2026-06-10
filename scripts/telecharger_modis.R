suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

set.seed(42)

ADM2_PATH <- "data/raw/ner_admin2.shp"
CLUSTERS_PATH <- "data/processed/clusters_gps.csv"
OUT_NDVI <- "data/processed/ndvi_dept.csv"

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# 1. Charger les départements
cat("=== 1. Chargement des départements ===\n")
admin2 <- sf::st_read(ADM2_PATH, quiet = TRUE) %>%
  mutate(dept_id = adm2_pcode)
cat(nrow(admin2), "départements\n")

# 2. Essayer MODISTools
cat("\n=== 2. Téléchargement NDVI MODIS ===\n")

if (!requireNamespace("MODISTools", quietly = TRUE)) {
  cat("MODISTools non disponible → fallback GEE\n")
  cat("Utiliser le script GEE manuel ci-dessous.\n")

  gee_script <- '
// Google Earth Engine JavaScript
// Coller ce code dans https://code.earthengine.google.com

var depts = ee.FeatureCollection("projects/ner_admin2");
var modis = ee.ImageCollection("MODIS/061/MOD13A3")
  .filterDate("2011-01-01", "2012-12-31")
  .select("NDVI_1km");

var ndvi_mean = modis.mean();
var ndvi_by_dept = ndvi_mean.reduceRegions({
  collection: depts,
  reducer: ee.Reducer.mean(),
  scale: 1000
});

Export.table.toDrive({
  collection: ndvi_by_dept,
  description: "ndvi_depts_2011_2012",
  fileFormat: "CSV",
  folder: "poverty_sae"
});
'
  cat(gee_script, "\n")
  cat("\nTélécharger manuellement le CSV et le placer dans data/processed/ndvi_dept.csv\n")
  quit(save = "no")
}

dept_centers <- admin2 %>%
  st_drop_geometry() %>%
  select(dept_id, adm2_name, center_lat, center_lon)

# MODISTools par lots pour éviter de surcharger le serveur
ndvi_list <- list()
batch_size <- 10
n_depts <- nrow(dept_centers)
n_batches <- ceiling(n_depts / batch_size)

for (b in seq_len(n_batches)) {
  idx_start <- (b - 1) * batch_size + 1
  idx_end   <- min(b * batch_size, n_depts)
  batch <- dept_centers[idx_start:idx_end, ]

  cat(sprintf("Batch %d/%d : départements %d-%d\n",
              b, n_batches, idx_start, idx_end))

  ndvi_batch <- MODISTools::mt_subset(
    product = "MOD13A3",
    band = "NDVI_1km",
    lat = batch$center_lat,
    lon = batch$center_lon,
    start = "2011-01-01",
    end = "2012-12-31",
    progress = FALSE
  )

  # Agrégation temporelle : moyenne des 24 mois
  ndvi_agg <- ndvi_batch %>%
    group_by(xllcorner, yllcorner) %>%
    summarise(ndvi_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(value = ndvi_mean / 10000)  # MODIS NDVI scale factor

  # Assigner chaque pixel au département le plus proche
  pixel_sf <- st_as_sf(
    ndvi_agg,
    coords = c("xllcorner", "yllcorner"),
    crs = 4326
  )

  nearest <- st_nearest_feature(pixel_sf, admin2)
  pixel_sf$dept_id <- admin2$dept_id[nearest]

  ndvi_dept <- pixel_sf %>%
    st_drop_geometry() %>%
    group_by(dept_id) %>%
    summarise(ndvi_mean = mean(value, na.rm = TRUE), .groups = "drop")

  ndvi_list[[b]] <- ndvi_dept
  Sys.sleep(15)  # Pause pour éviter rate-limiting
}

ndvi_final <- bind_rows(ndvi_list) %>%
  group_by(dept_id) %>%
  summarise(ndvi_mean = mean(ndvi_mean, na.rm = TRUE), .groups = "drop") %>%
  left_join(dept_centers, by = "dept_id") %>%
  select(dept_id, dept_name = adm2_name, ndvi_mean)

# 3. Sauvegarde
cat("\n=== 3. Sauvegarde ===\n")
write.csv(ndvi_final, OUT_NDVI, row.names = FALSE)
cat("->", OUT_NDVI, "\n")
cat("NDVI téléchargé pour", sum(!is.na(ndvi_final$ndvi_mean)), "/", nrow(ndvi_final), "départements\n")
cat("NDVI moyen :", round(mean(ndvi_final$ndvi_mean, na.rm = TRUE), 3), "\n")
cat("NDVI min/max :",
    round(min(ndvi_final$ndvi_mean, na.rm = TRUE), 3),
    "/",
    round(max(ndvi_final$ndvi_mean, na.rm = TRUE), 3), "\n")
cat("\nTerminé.\n")
