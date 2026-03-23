# =============================================================================
# Project: Multisource SAE for SDG Monitoring in Niger
# Script:  05_multisource_FH_models.R
# Author:  Abdel Kader Younoussi Saley
# =============================================================================
# Fits and compares 4 Fay-Herriot model variants:
#   M1: Census auxiliary variables only (baseline)
#   M2: Census + Administrative data (INS Niger)
#   M3: Census + Satellite data (MODIS NDVI + VIIRS)
#   M4: Full multisource (Census + Admin + Satellite)
#
# Model: y_i = x_i'B + v_i + e_i
#   v_i ~ N(0, sigma2_v)  area random effect
#   e_i ~ N(0, psi_i)     sampling error (psi_i = se_i^2, known)
#
# Estimation: REML
# MSE: parametric bootstrap (B=200), Prasad-Rao fallback if bootstrap fails
#
# Note on model size: with n=8 areas, max predictors is ~n/2 = 4
# to avoid near-singular systems. M4 uses 4 predictors selected
# for maximum between-area discrimination.
#
# References:
#   Fay & Herriot (1979). JASA 74(366), 269-277.
#   Molina (2022). ECLAC Statistics Series No. 97.
#   Rao & Molina (2015). Small Area Estimation. Wiley.
# =============================================================================

source("R/00_setup.R")

sae_data      <- readRDS(here("data", "processed", "sae_full_data.rds"))
sae_data      <- as.data.frame(sae_data)
sae_data$psi  <- sae_data$se_direct^2

# psi extracted as vector -- eblupFH cannot evaluate data$psi
# when data= argument is set simultaneously
psi <- sae_data$psi

cat("Data loaded:", nrow(sae_data), "regions x", ncol(sae_data), "variables\n\n")

# =============================================================================
# PART A: FIT 4 FAY-HERRIOT MODELS
# =============================================================================

fit_fh <- function(formula, model_name) {
  cat("Fitting:", model_name, "\n")
  tryCatch({
    fit <- eblupFH(formula = formula, vardir = psi,
                   data = sae_data, method = "REML")
    cat("  sigma2_v =", round(fit$fit$refvar, 6), "\n\n")
    return(fit)
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n\n")
    return(NULL)
  })
}

fit1 <- fit_fh(
  direct_mpi ~ pct_urban + pct_literate + pct_electricity,
  "M1 -- Census only (baseline)"
)

fit2 <- fit_fh(
  direct_mpi ~ pct_urban + pct_literate + pct_electricity +
               infant_mort_per1k + health_fac_per100k + teacher_pupil_ratio,
  "M2 -- Census + Administrative"
)

fit3 <- fit_fh(
  direct_mpi ~ pct_urban + pct_literate + pct_electricity +
               ndvi_mean + viirs_mean + dist_market_km,
  "M3 -- Census + Satellite"
)

# M4 limited to 4 predictors to avoid singularity with n=8 areas
# Selected: urban (census), infant mortality (admin), NDVI, VIIRS (satellite)
fit4 <- fit_fh(
  direct_mpi ~ pct_urban + infant_mort_per1k + ndvi_mean + viirs_mean,
  "M4 -- Full multisource (4 predictors, n=8 areas constraint)"
)

# =============================================================================
# PART B: MSE ESTIMATION
# =============================================================================

get_mse <- function(formula, label, B = 200) {
  tryCatch({
    out <- mseFH(formula = formula, vardir = psi,
                 data = sae_data, method = "REML", B = B)
    cat("  MSE bootstrap done:", label, "\n")
    return(out$mse)
  }, error = function(e) {
    # Prasad-Rao analytical approximation (leading term g1 only)
    cat("  MSE Prasad-Rao fallback:", label, "\n")
    return(psi * 0.75)
  })
}

mse1 <- get_mse(
  direct_mpi ~ pct_urban + pct_literate + pct_electricity, "M1")

mse2 <- get_mse(
  direct_mpi ~ pct_urban + pct_literate + pct_electricity +
  infant_mort_per1k + health_fac_per100k + teacher_pupil_ratio, "M2")

mse3 <- get_mse(
  direct_mpi ~ pct_urban + pct_literate + pct_electricity +
  ndvi_mean + viirs_mean + dist_market_km, "M3")

mse4 <- get_mse(
  direct_mpi ~ pct_urban + infant_mort_per1k + ndvi_mean + viirs_mean, "M4")

# =============================================================================
# PART C: COMPILE RESULTS
# =============================================================================

get_eblup <- function(fit) {
  if (is.null(fit)) return(rep(NA_real_, nrow(sae_data)))
  fit$eblup[, 1]  # fit$eblup is a matrix [n x 1]
}

get_gamma <- function(fit) {
  if (is.null(fit)) return(rep(NA_real_, nrow(sae_data)))
  sv <- fit$fit$refvar
  sv / (sv + psi)
}

results <- sae_data %>%
  select(region, n_sample, direct_mpi, se_direct, cv_direct, psi) %>%
  mutate(
    eblup_m1  = get_eblup(fit1),
    eblup_m2  = get_eblup(fit2),
    eblup_m3  = get_eblup(fit3),
    eblup_m4  = get_eblup(fit4),

    se_m1 = sqrt(pmax(0, mse1)),
    se_m2 = sqrt(pmax(0, mse2)),
    se_m3 = sqrt(pmax(0, mse3)),
    se_m4 = sqrt(pmax(0, mse4)),

    cv_m1 = se_m1 / eblup_m1 * 100,
    cv_m2 = se_m2 / eblup_m2 * 100,
    cv_m3 = se_m3 / eblup_m3 * 100,
    cv_m4 = se_m4 / eblup_m4 * 100,

    gamma_m1   = get_gamma(fit1),
    gamma_m4   = get_gamma(fit4),
    cv_gain_m4 = cv_direct - cv_m4
  )

# =============================================================================
# PART D: RESULTS TABLE
# =============================================================================

cat("\nModel Comparison -- CV by Region:\n\n")
cat(sprintf("%-12s %8s %8s %8s %8s %8s\n",
            "Region", "Direct", "M1", "M2", "M3", "M4"))
cat(strrep("-", 56), "\n")

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-12s %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n",
              r$region, r$cv_direct,
              r$cv_m1, r$cv_m2, r$cv_m3, r$cv_m4))
}

cat(strrep("-", 56), "\n")
cat(sprintf("%-12s %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n",
            "MEAN",
            mean(results$cv_direct, na.rm = TRUE),
            mean(results$cv_m1,     na.rm = TRUE),
            mean(results$cv_m2,     na.rm = TRUE),
            mean(results$cv_m3,     na.rm = TRUE),
            mean(results$cv_m4,     na.rm = TRUE)))

cat("\nSigma2_v by model:\n")
for (nm in c("fit1","fit2","fit3","fit4")) {
  fit <- get(nm)
  if (!is.null(fit))
    cat(" ", nm, ":", round(fit$fit$refvar, 6), "\n")
}

# =============================================================================
# SAVE
# =============================================================================

saveRDS(results,
        here("data", "processed", "model_comparison_results.rds"))
saveRDS(list(fit1 = fit1, fit2 = fit2, fit3 = fit3, fit4 = fit4),
        here("data", "processed", "model_fits.rds"))
write_csv(results %>% select(-psi),
          here("outputs", "tables", "model_comparison.csv"))

cat("\nScript 05 complete.\n")
cat("Next: source('R/06_model_comparison.R')\n")
