suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(sf)
  library(survey)
})

set.seed(42)

DHS_PATH <- "data/raw/NIHR61FL.DTA"
GPS_PATH <- "data/raw/NIGE61FL.shp"
ADM2_PATH <- "data/raw/ner_admin2.shp"
OUT_HOUSEHOLDS <- "data/processed/households_dept.rds"
OUT_DIRECT <- "data/processed/direct_estimates_dept.csv"
OUT_CLUSTERS_GPS <- "data/processed/clusters_gps.csv"

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# 1. Charger DHS
cat("=== 1. Chargement DHS ===\n")
if (!file.exists(DHS_PATH)) {
  stop("Fichier DHS introuvable : ", DHS_PATH)
}
dhs <- haven::read_dta(DHS_PATH)
cat("Ménages chargés :", nrow(dhs), "\n")

dhs_clean <- dhs %>%
  transmute(
    cluster = as.integer(hv001),
    weight  = as.numeric(hv005) / 1e6,
    region  = as.character(haven::as_factor(hv024)),
    urban   = as.integer(hv025 == 1),
    wealth  = as.numeric(hv271)
  ) %>%
  filter(!is.na(wealth), !is.na(cluster), !is.na(weight))
cat("Ménages après nettoyage :", nrow(dhs_clean), "\n")

# 2. Charger GPS clusters
cat("\n=== 2. Chargement GPS ===\n")
gps <- sf::st_read(GPS_PATH, quiet = TRUE)
cat("Clusters GPS :", nrow(gps), "\n")

gps_clean <- gps %>%
  mutate(cluster = as.integer(DHSCLUST)) %>%
  select(cluster, LATNUM, LONGNUM) %>%
  st_drop_geometry()

# 3. Charger Admin 2
cat("\n=== 3. Chargement limites Admin 2 ===\n")
admin2 <- sf::st_read(ADM2_PATH, quiet = TRUE)
cat("Départements :", nrow(admin2), "\n")

# Utiliser adm2_pcode comme clé stable (ex: NE005001)
admin2 <- admin2 %>%
  mutate(
    dept_id   = adm2_pcode,
    dept_name = adm2_name,
    region_name = adm1_name
  ) %>%
  select(dept_id, dept_name, region_name, adm2_pcode, adm1_name, geometry)

# 4. Spatial join : assigner chaque cluster à un département
cat("\n=== 4. Spatial join clusters -> départements ===\n")
gps_sf <- sf::st_as_sf(
  gps_clean,
  coords = c("LONGNUM", "LATNUM"),
  crs = sf::st_crs(admin2),
  remove = FALSE
)

cluster_dept <- sf::st_join(gps_sf, admin2["dept_id"]) %>%
  st_drop_geometry()

# Vérifier s'il y a des clusters non assignés (tomber dans l'océan ~0.1%)
na_count <- sum(is.na(cluster_dept$dept_id))
if (na_count > 0) {
  cat("Clusters non assignés au 1er joint :", na_count, "\n")
  cat("Assignation au plus proche département...\n")
  nearest <- sf::st_nearest_feature(gps_sf[is.na(cluster_dept$dept_id), ], admin2)
  cluster_dept$dept_id[is.na(cluster_dept$dept_id)] <- admin2$dept_id[nearest]
}

cat("Départements couverts par au moins 1 cluster :",
    n_distinct(cluster_dept$dept_id, na.rm = TRUE), "/", nrow(admin2), "\n")

effectifs <- cluster_dept %>%
  group_by(dept_id) %>%
  summarise(n_clusters = n(), .groups = "drop")
cat("Clusters par département : min =", min(effectifs$n_clusters),
    ", max =", max(effectifs$n_clusters),
    ", médiane =", median(effectifs$n_clusters), "\n")

# 5. Joindre les ménages aux départements
cat("\n=== 5. Jointure ménages -> départements ===\n")
dhs_with_dept <- dhs_clean %>%
  left_join(cluster_dept %>% select(cluster, dept_id), by = "cluster")

dept_na <- sum(is.na(dhs_with_dept$dept_id))
if (dept_na > 0) {
  cat("Attention :", dept_na, "ménages sans département assigné\n")
  dhs_with_dept <- dhs_with_dept %>% filter(!is.na(dept_id))
}

# 6. Estimateurs directs par département
cat("\n=== 6. Estimateurs directs par département ===\n")
design <- svydesign(
  ids = ~cluster,
  weights = ~weight,
  data = dhs_with_dept,
  nest = TRUE
)

direct <- svyby(
  formula = ~wealth,
  by = ~dept_id,
  design = design,
  FUN = svymean
) %>%
  rename(direct_mean = wealth, se_direct = se) %>%
  mutate(
    cv_direct = 100 * se_direct / abs(direct_mean),
    psi       = se_direct^2
  ) %>%
  arrange(dept_id)

# Ajouter urban_pct, n_clusters
aux <- dhs_with_dept %>%
  group_by(dept_id) %>%
  summarise(
    urban_pct = 100 * weighted.mean(urban, weight),
    n_hh      = n(),
    .groups   = "drop"
  )

direct <- direct %>%
  left_join(aux, by = "dept_id") %>%
  left_join(effectifs, by = "dept_id") %>%
  left_join(
    admin2 %>% st_drop_geometry() %>% select(dept_id, dept_name, region_name),
    by = "dept_id"
  ) %>%
  select(dept_id, dept_name, region_name,
         n_hh, n_clusters, urban_pct,
         direct_mean, se_direct, cv_direct, psi)

cat("Estimateurs calculés pour", nrow(direct), "départements\n")
cat("CV direct moyen :", round(mean(direct$cv_direct, na.rm = TRUE), 1), "%\n")
cat("CV direct médian :", round(median(direct$cv_direct, na.rm = TRUE), 1), "%\n")

# 7. Sauvegardes
cat("\n=== 7. Sauvegardes ===\n")
saveRDS(dhs_with_dept, OUT_HOUSEHOLDS)
write.csv(direct, OUT_DIRECT, row.names = FALSE)

# Export clusters GPS + dépt pour Python
clusters_export <- cluster_dept %>%
  left_join(
    dhs_clean %>%
      group_by(cluster) %>%
      summarise(wealth_mean = mean(wealth), n_hh = n(), .groups = "drop"),
    by = "cluster"
  ) %>%
  left_join(admin2 %>% st_drop_geometry() %>% select(dept_id, dept_name), by = "dept_id")

write.csv(clusters_export, OUT_CLUSTERS_GPS, row.names = FALSE)

cat("->", OUT_HOUSEHOLDS, "\n")
cat("->", OUT_DIRECT, "\n")
cat("->", OUT_CLUSTERS_GPS, "\n")
cat("\nTerminé.\n")
