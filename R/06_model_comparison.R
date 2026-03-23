# =============================================================================
# Project: Multisource SAE for SDG Monitoring in Niger
# Script:  06_model_comparison.R
# Author:  Abdel Kader Younoussi Saley
# =============================================================================

source("R/00_setup.R")

results  <- readRDS(here("data", "processed", "model_comparison_results.rds"))
fits     <- readRDS(here("data", "processed", "model_fits.rds"))
sae_data <- as.data.frame(
  readRDS(here("data", "processed", "sae_full_data.rds"))
)

fit1 <- fits$fit1
fit2 <- fits$fit2
fit3 <- fits$fit3
fit4 <- fits$fit4

# =============================================================================
# DIAGNOSTIC 1: Residuals for Model 4
# =============================================================================

cat("Residual Diagnostics -- Model 4:\n\n")

if (!is.null(fit4)) {
  resid_m4 <- sae_data$direct_mpi - fit4$eblup[, 1]

  cat("  Mean:  ", round(mean(resid_m4), 4), "(should be near 0)\n")
  cat("  SD:    ", round(sd(resid_m4), 4), "\n")
  cat("  Range: [", round(min(resid_m4), 3), ",",
      round(max(resid_m4), 3), "]\n")

  sw <- shapiro.test(resid_m4)
  cat("  Shapiro-Wilk p:", round(sw$p.value, 3),
      ifelse(sw$p.value > 0.05,
             "-- normality not rejected",
             "-- normality rejected, check model"), "\n\n")
} else {
  cat("  fit4 is NULL -- run script 05 first\n\n")
}

# =============================================================================
# DIAGNOSTIC 2: Shrinkage factors
# =============================================================================

cat("Shrinkage Factors (gamma_i = sigma2_v / (sigma2_v + psi_i)):\n")
cat("  gamma_i -> 1: trust direct estimate\n")
cat("  gamma_i -> 0: borrow from model\n\n")

results %>%
  select(region, n_sample, gamma_m1, gamma_m4, cv_direct, cv_m4) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

# =============================================================================
# DIAGNOSTIC 3: Model selection (AIC, BIC from REML goodness vector)
# =============================================================================
# sae::eblupFH returns fit$fit$goodness as named numeric vector
# with elements: loglike, AIC, BIC, KIC

cat("\nModel Selection Criteria:\n\n")

extract_stats <- function(fit, label, n_aux) {
  if (is.null(fit)) return(tibble(
    Model = label, sigma2_v = NA, AIC = NA, BIC = NA, N_aux = n_aux))

  g <- fit$fit$goodness
  tibble(
    Model    = label,
    sigma2_v = round(fit$fit$refvar, 6),
    AIC      = if ("AIC" %in% names(g)) round(g["AIC"], 2) else NA_real_,
    BIC      = if ("BIC" %in% names(g)) round(g["BIC"], 2) else NA_real_,
    N_aux    = n_aux
  )
}

model_stats <- bind_rows(
  extract_stats(fit1, "M1: Census",     3),
  extract_stats(fit2, "M2: +Admin",     6),
  extract_stats(fit3, "M3: +Satellite", 6),
  extract_stats(fit4, "M4: Full",       4)
) %>%
  mutate(
    Mean_CV = c(
      mean(results$cv_m1, na.rm = TRUE),
      mean(results$cv_m2, na.rm = TRUE),
      mean(results$cv_m3, na.rm = TRUE),
      mean(results$cv_m4, na.rm = TRUE)
    )
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

print(model_stats)

cat("\nInterpretation:\n")
cat("  Lower AIC/BIC = better fit-parsimony tradeoff\n")
cat("  Lower sigma2_v = auxiliary variables better explain area heterogeneity\n")
cat("  Lower Mean_CV = more precise subnational estimates\n")

write_csv(model_stats, here("outputs", "tables", "model_selection.csv"))

cat("\nScript 06 complete.\n")
cat("Next: source('R/07_differential_privacy.R')\n")
