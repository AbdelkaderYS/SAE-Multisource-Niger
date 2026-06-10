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
OUT_DETAILS <- "outputs/tables/fh_results_details.csv"
OUT_FIG_DIR <- "outputs/figures"

dir.create(OUT_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Charger les donnees
cat("=== 1. Chargement des donnees ===\n")
direct <- read.csv(DIRECT_PATH, stringsAsFactors = FALSE)

# NDVI
ndvi <- NULL
if (file.exists(NDVI_PATH)) {
  ndvi <- read.csv(NDVI_PATH, stringsAsFactors = FALSE)
  cat("NDVI charge :", nrow(ndvi), "departements\n")
} else {
  cat("NDVI non trouve\n")
}

# CNN predictions
cnn_cluster <- NULL
if (file.exists(CNN_PATH)) {
  cnn_cluster <- read.csv(CNN_PATH, stringsAsFactors = FALSE)
  cat("Predictions CNN chargees :", nrow(cnn_cluster), "clusters\n")
} else {
  cat("Predictions CNN non trouvees\n")
}

admin2 <- sf::st_read(ADM2_PATH, quiet = TRUE) %>%
  mutate(dept_id = adm2_pcode)

clusters <- read.csv(CLUSTERS_PATH, stringsAsFactors = FALSE)

# 2. Agregation CNN par departement
cat("\n=== 2. Agregation CNN par departement ===\n")
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
  cat("CNN agrege :", nrow(cnn_dept), "departements\n")
}

# 3. Fusion
cat("\n=== 3. Fusion des tables ===\n")
sae <- direct %>%
  left_join(ndvi[, c("dept_id", "ndvi_mean")], by = "dept_id")

if (!is.null(cnn_dept)) {
  sae <- sae %>% left_join(cnn_dept, by = "dept_id")
}

# Remplacer NA par la moyenne
sae$ndvi_mean[is.na(sae$ndvi_mean)] <- mean(sae$ndvi_mean, na.rm = TRUE)
if (!is.null(cnn_dept)) {
  sae$cnn_score[is.na(sae$cnn_score)] <- mean(sae$cnn_score, na.rm = TRUE)
} else {
  sae$cnn_score <- 0
}

cat("Table SAE fusionnee :", nrow(sae), "departements\n")

# Flag departements a 1 seul cluster
sae$flag_1cluster <- ifelse(sae$n_clusters == 1, 1, 0)
cat("Departements avec 1 cluster :", sum(sae$flag_1cluster), "/", nrow(sae), "\n")

# Matrice de correlation entre les predicteurs
cat("\n=== Matrice de correlation (predicteurs) ===\n")
cor_vars <- sae %>% dplyr::select(urban_pct, ndvi_mean, cnn_score) %>% na.omit()
cor_matrix <- cor(cor_vars)
print(round(cor_matrix, 3))

# 4. Modeles FH
cat("\n=== 4. Ajustement des modeles FH ===\n")

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

    aic <- as.numeric(fit$fit$goodness["AIC"])
    bic <- as.numeric(fit$fit$goodness["BIC"])

    cv_eblup <- 100 * se_eblup / abs(eblup)

    cat(sprintf("  %s : sigma2_v=%.4f, AIC=%.1f, BIC=%.1f, gamma_moyen=%.3f, CV=%.1f%%\n",
                label, sigma2_v, aic, bic, mean(gamma_vals), mean(cv_eblup, na.rm = TRUE)))

    list(
      label = label,
      sigma2_v = sigma2_v,
      aic = aic,
      bic = bic,
      gamma_mean = mean(gamma_vals),
      eblup = eblup,
      se_eblup = se_eblup,
      gamma = gamma_vals,
      cv_eblup = cv_eblup,
      resid = direct$direct_mean - eblup
    )
  }, error = function(e) {
    cat(sprintf("  %s : ECHEC (%s)\n", label, e$message))
    NULL
  })
}

# Definir tous les modeles
models <- list()

models[[1]] <- fit_model(direct_mean ~ 1, "FH-0 (moyenne seule)")
models[[2]] <- fit_model(direct_mean ~ urban_pct, "FH-1 (urban_pct)")
models[[3]] <- fit_model(direct_mean ~ urban_pct + ndvi_mean, "FH-1b (urban_pct + ndvi)")

if (!is.null(cnn_dept) && any(sae$cnn_score != 0)) {
  models[[4]] <- fit_model(direct_mean ~ urban_pct + cnn_score, "FH-1c (urban_pct + cnn)")
  models[[5]] <- fit_model(
    direct_mean ~ urban_pct + ndvi_mean + cnn_score,
    "FH-2 (urban_pct + ndvi + cnn)"
  )
} else {
  models[[4]] <- fit_model(direct_mean ~ urban_pct + ndvi_mean, "FH-2 (urban_pct + ndvi)")
}

models <- Filter(Negate(is.null), models)

# 5. Bootstrap IC pour le CV du meilleur modele
cat("\n=== 5. Bootstrap IC pour le CV ===\n")
best_idx <- length(models)
n_boot <- 1000
boot_cv <- numeric(n_boot)

for (b in 1:n_boot) {
  idx <- sample(nrow(sae), replace = TRUE)
  boot_data <- sae[idx, ]
  tryCatch({
    fit <- sae::eblupFH(
      formula = direct_mean ~ urban_pct + ndvi_mean + cnn_score,
      vardir = psi,
      data = boot_data,
      method = "REML"
    )
    mse <- sae::mseFH(
      formula = direct_mean ~ urban_pct + ndvi_mean + cnn_score,
      vardir = psi,
      data = boot_data,
      method = "REML"
    )
    eblup <- as.numeric(fit$eblup)
    se <- sqrt(as.numeric(mse$mse))
    cv <- 100 * se / abs(eblup)
    boot_cv[b] <- mean(cv, na.rm = TRUE)
  }, error = function(e) {
    boot_cv[b] <- NA
  })
}

boot_cv <- na.omit(boot_cv)
cv_ci <- quantile(boot_cv, probs = c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("Bootstrap (%d iterations) :\n", length(boot_cv)))
cat(sprintf("  CV moyen = %.1f%%\n", mean(boot_cv)))
cat(sprintf("  IC 95%%  = [%.1f%%, %.1f%%]\n", cv_ci[1], cv_ci[2]))

# CV moyen sans les departs a 1 cluster
sae_no1 <- sae %>% filter(flag_1cluster == 0)
if (nrow(sae_no1) > 0) {
  tryCatch({
    fit_no1 <- sae::eblupFH(
      direct_mean ~ urban_pct + ndvi_mean + cnn_score,
      vardir = psi,
      data = sae_no1,
      method = "REML"
    )
    mse_no1 <- sae::mseFH(
      direct_mean ~ urban_pct + ndvi_mean + cnn_score,
      vardir = psi,
      data = sae_no1,
      method = "REML"
    )
    eblup_no1 <- as.numeric(fit_no1$eblup)
    se_no1 <- sqrt(as.numeric(mse_no1$mse))
    cv_no1 <- mean(100 * se_no1 / abs(eblup_no1), na.rm = TRUE)
    cat(sprintf("\nCV sans %d departements 1-cluster = %.1f%%\n",
                sum(sae$flag_1cluster), cv_no1))
  }, error = function(e) {
    cat("Impossible de calculer le CV sans les 1-cluster\n")
  })
}

# 6. Tableau comparatif
cat("\n=== 6. Tableau comparatif ===\n")
comparison <- do.call(rbind, lapply(models, function(m) {
  data.frame(
    modele = m$label,
    sigma2_v = round(m$sigma2_v, 4),
    AIC = round(m$aic, 1),
    BIC = round(m$bic, 1),
    gamma_moyen = round(m$gamma_mean, 3),
    cv_moyen_eblup = round(mean(m$cv_eblup, na.rm = TRUE), 1),
    stringsAsFactors = FALSE
  )
}))
# Ajouter l'IC du bootstrap au meilleur modele
comparison$cv_ic_bas <- NA
comparison$cv_ic_haut <- NA
comparison$cv_ic_bas[nrow(comparison)] <- round(cv_ci[1], 1)
comparison$cv_ic_haut[nrow(comparison)] <- round(cv_ci[2], 1)

print(comparison)
write.csv(comparison, OUT_RESULTS, row.names = FALSE)

# 7. Tableau detaille du meilleur modele
cat("\n=== 7. Meilleur modele : details par departement ===\n")
best <- models[[best_idx]]
results_dept <- sae %>%
  mutate(
    eblup    = best$eblup,
    se_eblup = best$se_eblup,
    gamma    = best$gamma,
    cv_eblup = 100 * se_eblup / abs(eblup),
    cv_gain  = cv_direct - cv_eblup,
    resid    = best$resid
  ) %>%
  dplyr::select(dept_id, dept_name, n_clusters, flag_1cluster,
                urban_pct, ndvi_mean,
                direct_mean, eblup, se_direct, se_eblup,
                cv_direct, cv_eblup, gamma, cv_gain, resid) %>%
  arrange(desc(cv_gain))

print(head(results_dept, 15))
write.csv(results_dept, OUT_DETAILS, row.names = FALSE)

# 8. Diagnostics
cat("\n=== 8. Diagnostics ===\n")

# QQ-plot des residus standardises
resid_std <- scale(best$resid)
df_qq <- data.frame(
  theoretical = qnorm(ppoints(length(resid_std))),
  sample = sort(resid_std)
)
p_qq <- ggplot(df_qq, aes(x = theoretical, y = sample)) +
  geom_point(color = "#1f77b4", size = 2.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Q-Q plot of standardized residuals (FH-2)",
       x = "Theoretical quantiles (normal)",
       y = "Standardized residuals") +
  theme_bw()
ggsave(file.path(OUT_FIG_DIR, "diag_qqplot.png"), p_qq, width = 6, height = 5, dpi = 200)
cat("  Q-Q plot : output/figures/diag_qqplot.png\n")

# Residus vs valeurs ajustees
df_resid <- data.frame(
  fitted = best$eblup,
  resid = best$resid
)
p_resid <- ggplot(df_resid, aes(x = fitted, y = resid)) +
  geom_point(size = 2.5, alpha = 0.7, color = "#1f77b4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs fitted values (FH-2)",
       x = "Fitted values (EBLUP)",
       y = "Residuals") +
  theme_bw()
ggsave(file.path(OUT_FIG_DIR, "diag_resid_fitted.png"), p_resid, width = 6, height = 5, dpi = 200)
cat("  Residus vs fitted : output/figures/diag_resid_fitted.png\n")

# 9. Figures existantes (cartes, comparaison)
cat("\n=== 9. Figures ===\n")

# Carte richesse directe
cat("  Carte richesse directe...\n")
admin2_direct <- admin2 %>%
  left_join(sae, by = "dept_id")

p1 <- ggplot() +
  geom_sf(data = admin2_direct, aes(fill = direct_mean), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f0f0f0", high = "#b2182b",
                       midpoint = 0, name = "Wealth") +
  labs(title = "Mean wealth by department (direct estimate)",
       subtitle = "DHS Niger 2012",
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG_DIR, "carte_wealth_direct.png"), p1,
       width = 8, height = 7, dpi = 200)

# Carte EBLUP
cat("  Carte richesse EBLUP...\n")
admin2_eblup <- admin2 %>%
  left_join(
    sae %>% mutate(eblup = best$eblup) %>% dplyr::select(dept_id, eblup),
    by = "dept_id"
  )

p2 <- ggplot() +
  geom_sf(data = admin2_eblup, aes(fill = eblup), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f0f0f0", high = "#b2182b",
                       midpoint = 0, name = "Wealth") +
  labs(title = "Mean wealth by department (EBLUP)",
       subtitle = paste0("Model: ", best$label),
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG_DIR, "carte_wealth_eblup.png"), p2,
       width = 8, height = 7, dpi = 200)

# Comparaison des CV entre modeles
cat("  Comparaison des CV...\n")
cv_data <- do.call(rbind, lapply(models, function(m) {
  data.frame(
    modele = m$label,
    cv = m$cv_eblup,
    dept_id = sae$dept_id,
    stringsAsFactors = FALSE
  )
}))
cv_summary <- cv_data %>%
  group_by(modele) %>%
  summarise(cv_moyen = mean(cv, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(cv_summary, aes(x = modele, y = cv_moyen, fill = modele)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", cv_moyen)), vjust = -0.5, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Mean CV comparison across models",
       subtitle = "CV = SE / |mean| * 100",
       x = "", y = "Mean CV (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_FIG_DIR, "comparaison_modeles.png"), p3,
       width = 8, height = 5, dpi = 200)

# Gamma vs nombre de clusters
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
  labs(title = "Shrinkage factor gamma vs number of clusters",
       subtitle = "gamma close to 1 = trust in direct estimate",
       x = "Number of clusters per department",
       y = expression(gamma[i])) +
  theme_bw()
ggsave(file.path(OUT_FIG_DIR, "shrinkage_vs_taille.png"), p4,
       width = 7, height = 5, dpi = 200)

cat("\nTermine.\n")
cat("Resultats :", OUT_RESULTS, "\n")
cat("Details   :", OUT_DETAILS, "\n")
cat("Figures   :", OUT_FIG_DIR, "/*.png\n")
