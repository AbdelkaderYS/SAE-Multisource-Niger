# ===========================================================================
# modele_fay_herriot.R
#
# Ajuste un modele de Fay-Herriot (1979) aux estimateurs directs de la
# richesse par region, en utilisant deux variables auxiliaires :
#   - urban_pct : taux d'urbanisation regional (proxy du recensement)
#   - ndvi_mean : NDVI moyen regional (proxy satellite)
#
# Le modele : direct_mean ~ urban_pct + ndvi_mean
# La variance d'echantillonnage psi est fournie comme vardir (connue).
#
# Sortie : outputs/tables/fh_results.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(sae)
})

set.seed(42)

# --- 1. Charger la table SAe -------------------------------------------
if (!file.exists("outputs/data/sae_poverty_data.rds")) {
  stop("Lancer d'abord source('R/estimateurs_directs.R')")
}
sae_data <- as.data.frame(readRDS("outputs/data/sae_poverty_data.rds"))

# --- 2. Ajuster le modele de Fay-Herriot (REML) -------------------------
# vardir = psi : la variance d'echantillonnage, supposee connue.
# Note importante : eblupFH() utilise deparse(substitute(vardir)) pour
# retrouver le nom de la colonne dans `data`. Il faut donc passer le nom
# nu (psi), pas une chaine ("psi") ni une expression (sae_data$psi).
fit <- sae::eblupFH(
  formula  = direct_mean ~ urban_pct + ndvi_mean,
  vardir   = psi,
  data     = sae_data,
  method   = "REML"
)

# --- 3. Extraire EBLUP, SE et shrinkage ---------------------------------
# L'objet fit contient $eblup (vecteur d'EBLUP) et $fit$coef (beta_chapeau)
# Pour la SE et le gamma, on utilise mseFH()
mse <- sae::mseFH(
  formula  = direct_mean ~ urban_pct + ndvi_mean,
  vardir   = psi,
  data     = sae_data,
  method   = "REML"
)

eblup    <- as.numeric(fit$eblup)
se_eblup <- sqrt(as.numeric(mse$mse))
psi      <- sae_data$psi

# Estimation du facteur de shrinkage a partir de la composante de variance
# sigma2_v du modele REML
sigma2_v <- as.numeric(fit$fit$refvar)
gamma    <- sigma2_v / (sigma2_v + psi)
gamma    <- pmin(pmax(gamma, 0), 1)  # borner numeriquement

# --- 4. Table de resultats ----------------------------------------------
fh_results <- sae_data %>%
  mutate(
    eblup    = eblup,
    se_eblup = se_eblup,
    cv_eblup = 100 * abs(se_eblup / eblup),
    gamma    = gamma,
    cv_gain  = cv_direct - 100 * abs(se_eblup / eblup)
  ) %>%
  dplyr::select(region, n_sample, direct_mean, eblup,
                se_direct, se_eblup, cv_direct, cv_eblup, gamma, cv_gain) %>%
  arrange(direct_mean)

# --- 5. Sauvegarde ------------------------------------------------------
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(fh_results, "outputs/tables/fh_results.csv", row.names = FALSE)

cat("==> Modele Fay-Herriot ajuste (REML)\n")
cat("==> sigma2_v (variance de l'effet regional) :", round(sigma2_v, 5), "\n")
cat("==> Log-vraisemblance :", round(as.numeric(fit$fit$goodness["loglike"]), 3), "\n")
cat("==> AIC :", round(as.numeric(fit$fit$goodness["AIC"]), 3), "\n")
cat("==> BIC :", round(as.numeric(fit$fit$goodness["BIC"]), 3), "\n")
cat("==> Coefficients beta :\n")
print(fit$fit$estcoef)
cat("==> Gamma moyen :", round(mean(gamma), 3), "\n")
cat("==> Gain CV moyen (direct - EBLUP) :", round(mean(fh_results$cv_gain), 3), "points de %\n")
if (sigma2_v == 0) {
  cat("\nNote pedagogique : sigma2_v = 0 signifie que la variabilite entre regions\n")
  cat("est entierement expliquee par les predicteurs (urban_pct, ndvi_mean).\n")
  cat("L'EBLUP se confond alors avec la prediction de la droite de regression.\n")
}
cat("==> Sauvegarde : outputs/tables/fh_results.csv\n")
print(fh_results)
