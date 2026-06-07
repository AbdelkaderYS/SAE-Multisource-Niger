# ===========================================================================
# visualiser_resultats.R
#
# Genere 2 figures publication (200 dpi) :
#   1. Comparaison directe vs EBLUP avec intervalles de confiance 95%
#   2. Facteur de shrinkage gamma par region
#
# Sorties :
#   outputs/figures/01_direct_vs_eblup.png
#   outputs/figures/02_shrinkage.png
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# --- 1. Charger les resultats FH ----------------------------------------
if (!file.exists("outputs/tables/fh_results.csv")) {
  stop("Lancer d'abord source('R/modele_fay_herriot.R')")
}
fh <- read.csv("outputs/tables/fh_results.csv", stringsAsFactors = FALSE)

# --- 2. Preparer les donnees pour ggplot --------------------------------
fh_long <- fh %>%
  dplyr::select(region, direct_mean, eblup, se_direct, se_eblup) %>%
  tidyr::pivot_longer(
    cols      = c(direct_mean, eblup),
    names_to  = "estimator",
    values_to = "estimate"
  ) %>%
  mutate(
    se = ifelse(estimator == "direct_mean", se_direct, se_eblup),
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    estimator = recode(estimator,
                       direct_mean = "Direct (svyby)",
                       eblup       = "EBLUP (Fay-Herriot)"),
    region = factor(region, levels = fh$region[order(fh$direct_mean)])
  )

# --- 3. Figure 1 : Direct vs EBLUP avec IC 95% --------------------------
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

p1 <- ggplot(fh_long, aes(x = region, y = estimate, colour = estimator, shape = estimator)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.15, linewidth = 0.6,
                position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("Direct (svyby)" = "#1f77b4",
                                  "EBLUP (Fay-Herriot)" = "#d62728")) +
  scale_shape_manual(values = c("Direct (svyby)" = 16,
                                 "EBLUP (Fay-Herriot)" = 17)) +
  labs(
    title    = "Richesse moyenne par region du Niger",
    subtitle = "Estimateur direct (DHS) vs EBLUP (Fay-Herriot) avec IC 95%",
    x = "Region (ordonnee par richesse directe croissante)",
    y = "Indice de richesse (DHS HV271)",
    colour = "Estimateur", shape = "Estimateur",
    caption = "DHS Niger 2012 / Modele Fay-Herriot REML avec urban_pct + ndvi_mean"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top",
    plot.caption = element_text(hjust = 0, face = "italic", size = 8)
  )

ggsave("outputs/figures/01_direct_vs_eblup.png", p1, width = 9, height = 5, dpi = 200)
cat("==> Figure 1 sauvegardee : outputs/figures/01_direct_vs_eblup.png\n")

# --- 4. Figure 2 : Facteur de shrinkage gamma par region ----------------
fh$region <- factor(fh$region, levels = fh$region[order(fh$gamma)])

p2 <- ggplot(fh, aes(x = region, y = gamma)) +
  geom_col(fill = "#2ca02c", alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey30") +
  geom_text(aes(label = sprintf("%.2f", gamma)),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Facteur de shrinkage gamma par region",
    subtitle = "gamma proche de 1 = on fait confiance a l'estimateur direct ; proche de 0 = on se fie au modele",
    x = "Region",
    y = expression(paste("Facteur de shrinkage ", gamma[i])),
    caption = "gamma_i = sigma2_v / (sigma2_v + psi_i) ; psi_i = SE_i^2 direct"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.caption = element_text(hjust = 0, face = "italic", size = 8)
  )

ggsave("outputs/figures/02_shrinkage.png", p2, width = 9, height = 5, dpi = 200)
cat("==> Figure 2 sauvegardee : outputs/figures/02_shrinkage.png\n")
