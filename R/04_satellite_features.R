# =============================================================================
# Project: Multisource SAE for SDG Monitoring in Niger
# Script:  04_satellite_features.R
# Author:  Abdel Kader Younoussi Saley
# =============================================================================
# Validates satellite variables as auxiliary variables for SAE.
# Computes correlations with direct MPI estimates and documents
# data quality following the UN-NQAF framework.
#
# Sources:
#   NDVI:  MODIS MOD13A3 v6.1 (NASA Earthdata -- earthdata.nasa.gov, free)
#   VIIRS: VNP46A3 (NOAA/NASA NGDC -- ngdc.noaa.gov/eog/viirs/, free)
#
# Both sources are non-probability: full territorial coverage but
# are proxies only -- not direct poverty measures.
# =============================================================================

source("R/00_setup.R")

sae_data <- readRDS(here("data", "processed", "sae_full_data.rds"))

# =============================================================================
# CORRELATION ANALYSIS
# =============================================================================

cat("Satellite auxiliary variables -- correlation with direct MPI:\n\n")

cors <- sae_data %>%
  summarise(
    cor_ndvi_mpi   = cor(ndvi_mean,      direct_mpi, use = "complete.obs"),
    cor_viirs_mpi  = cor(viirs_mean,     direct_mpi, use = "complete.obs"),
    cor_market_mpi = cor(dist_market_km, direct_mpi, use = "complete.obs")
  )

cat("  NDVI vs MPI:              r =", round(cors$cor_ndvi_mpi, 3), "\n")
cat("  Nightlights vs MPI:       r =", round(cors$cor_viirs_mpi, 3), "\n")
cat("  Market distance vs MPI:   r =", round(cors$cor_market_mpi, 3), "\n\n")

cat("Interpretation:\n")
cat("  NDVI: negative expected (more vegetation = less poverty in Sahel)\n")
cat("  VIIRS: strong negative expected (more lights = less poverty)\n")
cat("  Market distance: positive expected (farther = poorer)\n\n")

# =============================================================================
# NQAF QUALITY ASSESSMENT -- satellite sources
# =============================================================================
# Following UN-NQAF Manual for Official Statistics (UNSD, 2019)

cat("NQAF Quality Assessment -- Satellite Data:\n\n")

nqaf <- tibble(
  Source      = c("MODIS NDVI", "VIIRS Nightlights"),
  Relevance   = c("HIGH -- agricultural activity proxy for income",
                  "HIGH -- economic activity and electrification proxy"),
  Accuracy    = c("MEDIUM -- atmospheric interference in cloudy periods",
                  "MEDIUM -- light pollution and blooming effects"),
  Timeliness  = c("HIGH -- monthly composites available",
                  "HIGH -- monthly composites available"),
  Coverage    = c("FULL -- all Niger territory",
                  "FULL -- all Niger territory"),
  Bias        = c("Rainfall seasonality confounds poverty signal",
                  "Urban bias -- rural economic activity underestimated"),
  Recommended = c("Auxiliary variable in FH model only",
                  "Auxiliary variable in FH model only")
)

for (i in seq_len(nrow(nqaf))) {
  cat("Source:", nqaf$Source[i], "\n")
  cat("  Relevance:  ", nqaf$Relevance[i], "\n")
  cat("  Accuracy:   ", nqaf$Accuracy[i], "\n")
  cat("  Timeliness: ", nqaf$Timeliness[i], "\n")
  cat("  Bias:       ", nqaf$Bias[i], "\n")
  cat("  Use:        ", nqaf$Recommended[i], "\n\n")
}

saveRDS(nqaf, here("data", "processed", "satellite_nqaf.rds"))

cat("Script 04 complete.\n")
cat("Next: source('R/05_multisource_FH_models.R')\n")
