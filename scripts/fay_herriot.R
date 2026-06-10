suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(sae)
})

set.seed(42)

DIRECT_PATH <- "data/processed/direct_estimates_dept.csv"
NDVI_PATH <- "data/processed/ndvi_dept.csv"
CNN_PATH <- "data/processed/cnn_predictions_cluster.csv"
ADM2_PATH <- "data/raw/ner_admin2.shp"
CLUSTERS_PATH <- "data/processed/clusters_gps.csv"
OUT_RESULTS <- "outputs/tables/fh_results_comparison.csv"
OUT_FIG_DIR <- "outputs/figures"

dir.create(OUT_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Charger les données
cat("=== 1. Chargement des données ===\n")
direct <- read.csv(DIRECT_PATH, stringsAsFactors = FALSE)

ndvi <- NULL
if (file.exists(NDVI_PATH)) {
  ndvi <- read.csv(NDVI_PATH, stringsAsFactors = FALSE)
  cat("NDVI chargé :", nrow(ndvi), "départements\n")
} else {
  cat("NDVI non trouvé, utilisation d'une colonne vide\n")
}

cnn_cluster <- NULL
if (file.exists(CNN_PATH)) {
  cnn_cluster <- read.csv(CNN_PATH, stringsAsFactors = FALSE)
  cat("Prédictions CNN chargées :", nrow(cnn_cluster), "clusters\n")
} else {
  cat("Prédictions CNN non trouvées (le modèle n'a pas encore tourné)\n")
}

admin2 <- sf::st_read(ADM2_PATH, quiet = TRUE) %>%
  mutate(dept_id = adm2_pcode)

clusters <- read.csv(CLUSTERS_PATH, stringsAsFactors = FALSE)

# 2. Agrégation CNN par département
cat("\n=== 2. Agrégation CNN par département ===\n")
cnn_dept <- NULL
if (!is.null(cnn_cluster)) {
  cnn_dept <- clusters %>%
    dplyr::select(cluster, dept_id) %>%
    distinct() %>%
    inner_join(cnn_cluster, by = "cluster") %>%
    group_by(dept_id) %>%
    summarise(
      cnn_score = mean(cnn_prediction, na.rm = TRUE),
      cnn_n     = n(),
      .groups   = "drop"
    )
  cat("CNN agrégé :", nrow(cnn_dept), "départements\n")
}

# 3. Fusion
cat("\n=== 3. Fusion des tables ===\n")
  sae <- direct %>%
  left_join(ndvi[, c("dept_id", "ndvi_mean")], by = "dept_id")

if (!is.null(cnn_dept)) {
  sae <- sae %>% left_join(cnn_dept, by = "dept_id")
}

# Remplacer les NA par la moyenne
sae$ndvi_mean[is.na(sae$ndvi_mean)] <- mean(sae$ndvi_mean, na.rm = TRUE)
if (!is.null(cnn_dept)) {
  sae$cnn_score[is.na(sae$cnn_score)] <- mean(sae$cnn_score, na.rm = TRUE)
} else {
  sae$cnn_score <- 0
}

cat("Table SAE fusionnée :", nrow(sae), "départements\n")

# 4. Modèles FH
cat("\n=== 4. Ajustement des modèles FH ===\n")

fit_model <- function(formula, label) {
  tryCatch({
    fit <- sae::eblupFH(
      formula = formula,
      vardir = psi,
      data = sae,
      method = "REML"
    )
    mse <- sae::mseFH(
      formula = formula,
      vardir = psi,
      data = sae,
      method = "REML"
    )

    sigma2_v <- as.numeric(fit$fit$refvar)
    eblup <- as.numeric(fit$eblup)
    se_eblup <- sqrt(as.numeric(mse$mse))
    gamma_vals <- sigma2_v / (sigma2_v + sae$psi)
    gamma_vals <- pmin(pmax(gamma_vals, 0), 1)

    loglik <- as.numeric(fit$fit$goodness["loglike"])
    aic <- as.numeric(fit$fit$goodness["AIC"])
    bic <- as.numeric(fit$fit$goodness["BIC"])

    cat(sprintf("  %s : sigma2_v=%.4f, AIC=%.1f, BIC=%.1f, gamma_moyen=%.3f\n",
                label, sigma2_v, aic, bic, mean(gamma_vals)))

    list(
      label = label,
      sigma2_v = sigma2_v,
      aic = aic,
      bic = bic,
      gamma_mean = mean(gamma_vals),
      eblup = eblup,
      se_eblup = se_eblup,
      gamma = gamma_vals
    )
  }, error = function(e) {
    cat(sprintf("  %s : ÉCHEC (%s)\n", label, e$message))
    NULL
  })
}

models <- list(
  fit_model(direct_mean ~ 1, "FH-0 (moyenne seule)"),
  fit_model(direct_mean ~ urban_pct, "FH-1 (urban_pct)")
)

# FH-2 : selon disponibilité des données
if (!is.null(cnn_dept) && any(sae$cnn_score != 0)) {
  models[[3]] <- fit_model(
    direct_mean ~ urban_pct + ndvi_mean + cnn_score,
    "FH-2 (urban_pct + ndvi + cnn)"
  )
} else {
  models[[3]] <- fit_model(
    direct_mean ~ urban_pct + ndvi_mean,
    "FH-2 (urban_pct + ndvi)"
  )
}

models <- Filter(Negate(is.null), models)

# 5. Tableau comparatif
cat("\n=== 5. Résultats ===\n")
comparison <- do.call(rbind, lapply(models, function(m) {
  data.frame(
    modele = m$label,
    sigma2_v = round(m$sigma2_v, 4),
    AIC = round(m$aic, 1),
    BIC = round(m$bic, 1),
    gamma_moyen = round(m$gamma_mean, 3),
    cv_moyen_direct = round(mean(sae$cv_direct, na.rm = TRUE), 1),
    cv_moyen_eblup = round(mean(100 * m$se_eblup / abs(m$eblup), na.rm = TRUE), 1),
    stringsAsFactors = FALSE
  )
}))
print(comparison)
write.csv(comparison, OUT_RESULTS, row.names = FALSE)

# 6. Tableau détaillé du meilleur modèle
cat("\n=== 6. Meilleur modèle : détails par département ===\n")
best <- models[[length(models)]]
  results_dept <- sae %>%
  mutate(
    eblup    = best$eblup,
    se_eblup = best$se_eblup,
    gamma    = best$gamma,
    cv_eblup = 100 * se_eblup / abs(eblup),
    cv_gain  = cv_direct - cv_eblup
  ) %>%
  dplyr::select(dept_id, dept_name, n_clusters, urban_pct,
                direct_mean, eblup, se_direct, se_eblup,
                cv_direct, cv_eblup, gamma, cv_gain, ndvi_mean) %>%
  arrange(desc(cv_gain))

print(head(results_dept, 15))
write.csv(results_dept, "outputs/tables/fh_results_details.csv", row.names = FALSE)

# 7. Figures
cat("\n=== 7. Génération des figures ===\n")

# Figure 1 : Carte richesse directe
cat("  Carte richesse directe...\n")
admin2_direct <- admin2 %>%
  left_join(sae, by = "dept_id")

p1 <- ggplot() +
  geom_sf(data = admin2_direct, aes(fill = direct_mean), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#d62728", mid = "#f0f0f0", high = "#2ca02c",
                       midpoint = 0, name = "Richesse") +
  labs(title = "Richesse moyenne par département (estimateur direct)",
       subtitle = "DHS Niger 2012, enquête ménages",
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG_DIR, "carte_wealth_direct.png"), p1,
       width = 8, height = 7, dpi = 200)

# Figure 2 : Carte EBLUP
cat("  Carte richesse EBLUP...\n")
admin2_eblup <- admin2 %>%
  left_join(
    sae %>% mutate(eblup = best$eblup) %>% dplyr::select(dept_id, eblup),
    by = "dept_id"
  )

p2 <- ggplot() +
  geom_sf(data = admin2_eblup, aes(fill = eblup), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#d62728", mid = "#f0f0f0", high = "#2ca02c",
                       midpoint = 0, name = "Richesse") +
  labs(title = "Richesse moyenne par département (EBLUP)",
       subtitle = paste0("Modèle : ", best$label),
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG_DIR, "carte_wealth_eblup.png"), p2,
       width = 8, height = 7, dpi = 200)

# Figure 3 : Comparaison des CV
cat("  Comparaison des CV...\n")
cv_data <- do.call(rbind, lapply(models, function(m) {
  data.frame(
    modele = m$label,
    cv = 100 * m$se_eblup / abs(m$eblup),
    dept_id = sae$dept_id,
    stringsAsFactors = FALSE
  )
}))
cv_summary <- cv_data %>%
  group_by(modele) %>%
  summarise(cv_moyen = mean(cv, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(cv_summary, aes(x = modele, y = cv_moyen, fill = modele)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", cv_moyen)), vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Comparaison du CV moyen entre modèles",
       subtitle = "CV = SE / |moyenne| × 100",
       x = "", y = "CV moyen (%)") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(OUT_FIG_DIR, "comparaison_modeles.png"), p3,
       width = 7, height = 5, dpi = 200)

# Figure 4 : gamma vs nombre de clusters
cat("  Gamma vs taille...\n")
gamma_data <- data.frame(
  n_clusters = sae$n_clusters,
  gamma = best$gamma,
  dept_name = sae$dept_name
)

p4 <- ggplot(gamma_data, aes(x = n_clusters, y = gamma)) +
  geom_point(size = 2.5, alpha = 0.7, color = "#1f77b4") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#d62728", alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "#2ca02c", linewidth = 0.8) +
  labs(title = "Facteur de shrinkage gamma vs nombre de clusters",
       subtitle = "gamma proche de 1 = confiance dans l'estimateur direct",
       x = "Nombre de clusters par département",
       y = expression(gamma[i])) +
  theme_bw()
ggsave(file.path(OUT_FIG_DIR, "shrinkage_vs_taille.png"), p4,
       width = 7, height = 5, dpi = 200)

cat("\nTerminé.\n")
cat("Résultats :", OUT_RESULTS, "\n")
cat("Détails   : outputs/tables/fh_results_details.csv\n")
cat("Figures   :", OUT_FIG_DIR, "/*.png (4 figures)\n")
