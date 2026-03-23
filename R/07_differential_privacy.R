# =============================================================================
# Project: Multisource SAE for SDG Monitoring in Niger
# Script:  07_differential_privacy.R
# Author:  Abdel Kader Younoussi Saley
# =============================================================================
# Applies Differential Privacy (Gaussian mechanism) to the real SAE
# estimates produced by script 05. This connects DP directly to the
# pipeline output rather than using hardcoded demonstration values.
#
# (epsilon, delta)-DP Gaussian mechanism:
#   M(D) = f(D) + N(0, sigma^2)
#   sigma = sensitivity * sqrt(2 * log(1.25/delta)) / epsilon
#
# For poverty rates bounded in [0,1]:
#   Global sensitivity = 1/n_i (adding one household changes rate by 1/n_i)
#
# References:
#   Dwork & Roth (2014). Foundations & Trends TCS 9(3-4).
#   Ferrando et al. (2022). Stats & Public Policy.
# =============================================================================

source("R/00_setup.R")

# Load real SAE results from script 05 -- not hardcoded values
results <- readRDS(here("data", "processed", "model_comparison_results.rds"))

sae_results <- results %>%
  select(region, eblup_m4, se_m4, n_sample) %>%
  rename(eblup_rate = eblup_m4, se_eblup = se_m4) %>%
  filter(!is.na(eblup_rate))

cat("DP applied to real SAE estimates from script 05.\n")
cat("Regions:", nrow(sae_results), "\n\n")

# =============================================================================
# PART A: GAUSSIAN MECHANISM
# =============================================================================

gaussian_sigma <- function(sensitivity, epsilon, delta) {
  sensitivity * sqrt(2 * log(1.25 / delta)) / epsilon
}

# =============================================================================
# PART B: UTILITY-PRIVACY TRADEOFF ACROSS EPSILON VALUES
# =============================================================================

epsilon_grid <- c(0.1, 0.5, 1.0, 2.0, 5.0, 10.0)
delta        <- 1e-5

cat("Utility-Privacy Tradeoff (delta =", delta, "):\n\n")

tradeoff <- map_dfr(epsilon_grid, function(eps) {
  noisy <- sae_results %>%
    mutate(
      sensitivity = 1 / n_sample,
      noise_sd    = gaussian_sigma(sensitivity, eps, delta),
      dp_estimate = eblup_rate + rnorm(n(), 0, noise_sd),
      dp_estimate = pmax(0, pmin(1, dp_estimate)),
      abs_error   = abs(dp_estimate - eblup_rate),
      rel_noise   = noise_sd / eblup_rate * 100
    )
  tibble(
    epsilon        = eps,
    mean_mae       = round(mean(noisy$abs_error), 4),
    mean_rel_noise = round(mean(noisy$rel_noise), 2),
    pct_within_5pp = round(mean(noisy$abs_error < 0.05) * 100, 0)
  )
})

print(tradeoff)

# =============================================================================
# PART C: RECOMMENDED CONFIGURATION AND FINAL DP ESTIMATES
# =============================================================================
# epsilon = 2 recommended:
#   - Mean absolute error < 1% of poverty rate
#   - All regions within 5pp of true EBLUP
#   - Acceptable utility for NSO publication standards

eps_rec   <- 2.0
delta_rec <- 1e-5

cat("\nRecommended configuration: epsilon =", eps_rec,
    ", delta =", delta_rec, "\n\n")

dp_results <- sae_results %>%
  mutate(
    sensitivity   = 1 / n_sample,
    noise_sd      = gaussian_sigma(sensitivity, eps_rec, delta_rec),
    dp_estimate   = eblup_rate + rnorm(n(), 0, noise_sd),
    dp_estimate   = pmax(0, pmin(1, dp_estimate)),
    absolute_error = abs(dp_estimate - eblup_rate),
    privacy_loss  = eps_rec
  )

print(dp_results %>%
  select(region, eblup_rate, dp_estimate, noise_sd, absolute_error) %>%
  mutate(across(where(is.numeric), ~round(., 4))))

cat("\nMean absolute error from DP noise:",
    round(mean(dp_results$absolute_error), 4), "\n")
cat("Max absolute error:               ",
    round(max(dp_results$absolute_error), 4), "\n")

# =============================================================================
# PLOT: utility-privacy tradeoff
# =============================================================================

p_dp <- tradeoff %>%
  ggplot(aes(x = epsilon, y = mean_mae)) +
  geom_line(color = "#2196F3", linewidth = 1) +
  geom_point(color = "#2196F3", size = 3) +
  geom_vline(xintercept = eps_rec, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  annotate("text", x = eps_rec + 0.3, y = max(tradeoff$mean_mae) * 0.9,
           label = paste("Recommended: e =", eps_rec),
           hjust = 0, size = 3.5) +
  scale_x_continuous(breaks = epsilon_grid) +
  labs(
    title    = "Differential Privacy: Utility-Privacy Tradeoff",
    subtitle = paste("Gaussian mechanism applied to EBLUP poverty rates |",
                     "delta =", delta),
    x        = "Privacy budget (epsilon)",
    y        = "Mean absolute error",
    caption  = "Lower epsilon = stronger privacy | Higher MAE = lower utility"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4)
  )

ggsave(here("outputs", "figures", "07_dp_tradeoff.png"),
       p_dp, width = 9, height = 5, dpi = 200, bg = "white")

saveRDS(dp_results, here("data", "processed", "dp_results.rds"))
write_csv(dp_results, here("outputs", "tables", "dp_results.csv"))

cat("\nScript 07 complete.\n")
cat("Next: source('R/08_visualization.R')\n")
