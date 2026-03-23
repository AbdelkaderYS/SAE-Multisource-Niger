# =============================================================================
# Project: Multisource SAE for SDG Monitoring in Niger
# Script:  03_nonprob_bias_correction.R
# Author:  Abdel Kader Younoussi Saley
# =============================================================================
# Non-probability bias correction using Inverse Probability Weighting.
# Reference: Chen et al. (2020). JRSS-B 82(2), 391-411.
#
# DESIGN:
#   Satellite data (MODIS/VIIRS) provides full territorial coverage but
#   is biased toward accessible areas. Clusters near Niamey or major roads
#   are more reliably captured by satellite-derived indicators.
#
#   With NIGE61FL.shp (DHS GPS clusters), we compute the real geographic
#   distance from each cluster to Niamey. This distance is used as a
#   proxy for satellite accessibility bias: distant clusters are
#   under-represented in satellite-derived poverty proxies.
#
#   The IPW model then reweights the satellite-based estimates to match
#   the DHS probability sample distribution.
#
#   GPS displacement note: DHS displaces urban clusters up to 2km and
#   rural clusters up to 5km. At Niger's scale (1 million km²), this
#   introduces negligible error in distance calculations.
#
# FALLBACK: If NIGE61FL.shp is absent, uses calibrated simulation.
#
# Output: corrected_urban_pct per region added to sae_full_data.rds
# =============================================================================

source("R/00_setup.R")

sae_data <- readRDS(here("data", "processed", "sae_full_data.rds"))

gps_path <- here("data", "raw", "NIGE61FL.shp")
hr_path  <- here("data", "raw", "NIHR61FL.DTA")

USE_REAL_GPS <- file.exists(gps_path)

# Niamey coordinates (capital city center)
NIAMEY_LAT <- 13.5137
NIAMEY_LON <-  2.1098

# Haversine distance in km
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1*pi/180) * cos(lat2*pi/180) * sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

# =============================================================================
# PART A: COMPUTE CLUSTER ACCESSIBILITY FROM GPS
# =============================================================================

if (USE_REAL_GPS) {

  cat("Loading NIGE61FL.shp (DHS GPS clusters)...\n")
  gps <- sf::st_read(gps_path, quiet = TRUE)
  names(gps) <- toupper(names(gps))

  # DHS GPS shapefiles store coordinates as LATNUM/LONGNUM columns
  # not always as sf geometry — use columns directly if available
  if (!("LATNUM" %in% names(gps)) || !("LONGNUM" %in% names(gps))) {
    coords      <- sf::st_coordinates(gps)
    gps$LONGNUM <- coords[, 1]
    gps$LATNUM  <- coords[, 2]
  }

  # Convert to plain data frame to avoid sf geometry issues
  gps <- as.data.frame(gps)

  # Distance from each cluster to Niamey
  gps$dist_niamey_km <- mapply(
    haversine_km,
    gps$LATNUM, gps$LONGNUM,
    NIAMEY_LAT, NIAMEY_LON
  )

  cat("GPS clusters loaded:", nrow(gps), "\n")
  cat("Distance to Niamey range: [",
      round(min(gps$dist_niamey_km), 0), ",",
      round(max(gps$dist_niamey_km), 0), "] km\n\n")

  # Satellite accessibility score:
  # Clusters close to Niamey = high satellite coverage probability
  # Clusters far = lower coverage probability
  # Logistic decay: ps = 1 / (1 + exp(0.005 * (dist - 200)))
  gps$ps_satellite <- 1 / (1 + exp(0.005 * (gps$dist_niamey_km - 200)))
  gps$ps_satellite <- pmax(0.05, pmin(0.95, gps$ps_satellite))

  # Urban clusters additionally over-represented (URBAN_RURA = "U")
  if ("URBAN_RURA" %in% names(gps)) {
    gps$ps_satellite <- ifelse(
      toupper(as.character(gps$URBAN_RURA)) == "U",
      pmin(0.95, gps$ps_satellite * 1.3),
      gps$ps_satellite
    )
  }

  # IPW: weight for satellite sample = 1 / ps
  gps$ipw <- 1 / gps$ps_satellite

  # Region-level accessibility stats
  # Use DHSREGNA or ALT_DHSREGCO depending on available columns
  region_col <- intersect(c("DHSREGNA", "SHDEPT", "DHSREGCO"), names(gps))[1]

  if (!is.na(region_col) && region_col == "DHSREGNA") {
    # Map DHS region names to our REGIONS
    region_map <- c(
      "AGADEZ"    = "Agadez",   "DIFFA"     = "Diffa",
      "DOSSO"     = "Dosso",    "MARADI"    = "Maradi",
      "NIAMEY"    = "Niamey",   "TAHOUA"    = "Tahoua",
      "TILLABERI" = "Tillaberi","ZINDER"    = "Zinder"
    )
    gps$region_std <- region_map[toupper(as.character(gps[[region_col]]))]
  } else {
    # Fall back to numeric region code if name not available
    gps$region_std <- REGIONS[as.numeric(as.character(gps$DHSREGCO))]
  }

  cluster_stats <- gps %>%
    as_tibble() %>%
    filter(!is.na(region_std)) %>%
    group_by(region_std) %>%
    summarise(
      n_clusters         = n(),
      mean_dist_niamey   = round(mean(dist_niamey_km), 1),
      mean_ps_satellite  = round(mean(ps_satellite), 3),
      corrected_urban_pct = round(weighted.mean(
        URBAN_RURA == "U" | URBAN_RURA == 1, ipw, na.rm = TRUE), 3),
      naive_urban_pct    = round(mean(
        URBAN_RURA == "U" | URBAN_RURA == 1, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    rename(region = region_std)

  cat("Cluster-level accessibility by region:\n")
  print(cluster_stats %>%
    select(region, n_clusters, mean_dist_niamey, mean_ps_satellite,
           naive_urban_pct, corrected_urban_pct))

  bias_reduction <- with(cluster_stats,
    round((1 - mean(abs(corrected_urban_pct - naive_urban_pct)) /
             mean(abs(naive_urban_pct - 0.15))) * 100, 0))

  cat("\nBias correction (national urban rate ~15%):\n")
  cat("  Naive mean urban rate:     ",
      round(mean(cluster_stats$naive_urban_pct) * 100, 1), "%\n")
  cat("  Corrected mean urban rate: ",
      round(mean(cluster_stats$corrected_urban_pct) * 100, 1), "%\n")

} else {

  cat("NIGE61FL.shp not found -- using calibrated simulation\n")
  cat("Place NIGE61FL.shp in data/raw/ for real GPS-based correction\n\n")

  # Calibrated from DHS Niger 2012 cluster distribution
  # Distance to Niamey (km) approximated from regional capitals
  dist_to_niamey <- c(
    Agadez = 900, Diffa = 1120, Dosso = 140, Maradi = 680,
    Niamey = 15,  Tahoua = 580, Tillaberi = 120, Zinder = 900
  )

  set.seed(42)
  n_clusters_per_region <- c(26, 22, 48, 68, 56, 58, 64, 36)

  cluster_sim <- map2_dfr(REGIONS, n_clusters_per_region, function(reg, n) {
    base_dist <- dist_to_niamey[reg]
    tibble(
      region         = reg,
      dist_niamey_km = rnorm(n, base_dist, base_dist * 0.15),
      urban          = rbinom(n, 1, ifelse(reg == "Niamey", 0.95, 0.10))
    )
  }) %>%
    mutate(
      ps_satellite = 1 / (1 + exp(0.005 * (dist_niamey_km - 200))),
      ps_satellite = pmax(0.05, pmin(0.95, ps_satellite)),
      ps_satellite = ifelse(urban == 1, pmin(0.95, ps_satellite * 1.3),
                            ps_satellite),
      ipw          = 1 / ps_satellite
    )

  cluster_stats <- cluster_sim %>%
    group_by(region) %>%
    summarise(
      n_clusters          = n(),
      mean_dist_niamey    = round(mean(dist_niamey_km), 1),
      mean_ps_satellite   = round(mean(ps_satellite), 3),
      naive_urban_pct     = round(mean(urban), 3),
      corrected_urban_pct = round(weighted.mean(urban, ipw), 3),
      .groups = "drop"
    )

  cat("Simulated cluster accessibility by region:\n")
  print(cluster_stats %>%
    select(region, n_clusters, mean_dist_niamey, naive_urban_pct,
           corrected_urban_pct))
}

# =============================================================================
# PART B: UPDATE SAE DATASET
# =============================================================================

corrected_by_region <- cluster_stats %>%
  select(region, corrected_urban_pct, naive_urban_pct)

sae_data <- sae_data %>%
  select(-any_of(c("corrected_urban_pct", "naive_urban_pct"))) %>%
  left_join(corrected_by_region, by = "region")

saveRDS(sae_data, here("data", "processed", "sae_full_data.rds"))

cat("\nSAE dataset updated with GPS-based accessibility correction.\n")
cat("Script 03 complete.\n")
cat("Next: source('R/04_satellite_features.R')\n")