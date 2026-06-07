# ===========================================================================
# estimateurs_directs.R
#
# Calcule les estimateurs directs de la richesse moyenne par region
# a partir du DHS (design de sondage + svyby), puis agrege une variable
# satellite (NDVI moyen regional simule) et sauvegarde une table fusionnee
# prete pour le modele de Fay-Herriot.
#
# Sorties :
#   outputs/data/sae_poverty_data.rds     : table pour Fay-Herriot
#   outputs/tables/direct_estimates.csv   : tableau des estimateurs directs
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
})

set.seed(42)

# --- 1. Charger le sous-ensemble DHS -------------------------------------
if (!file.exists("outputs/data/dhs_subset.rds")) {
  stop("Lancer d'abord source('R/charger_dhs.R')")
}
dhs <- readRDS("outputs/data/dhs_subset.rds")

# --- 2. Declarer le plan de sondage --------------------------------------
# ids = ~cluster (unites primaires de sondage = clusters DHS)
# weights = ~weight (poids de sondage normalises)
design <- svydesign(
  ids = ~cluster,
  weights = ~weight,
  data = dhs,
  nest = TRUE
)

# --- 3. Moyenne de richesse par region ----------------------------------
direct <- svyby(
  formula = ~wealth,
  by = ~region,
  design = design,
  FUN = svymean
) %>%
  rename(direct_mean = wealth, se_direct = se) %>%
  mutate(
    cv_direct = 100 * se_direct / abs(direct_mean),  # CV = |SE/mean| x 100
    psi       = se_direct^2,
    n_sample  = as.integer(unlist(tapply(dhs$weight, dhs$region, length)))
  ) %>%
  arrange(direct_mean)

# --- 4. Variable auxiliaire X1 : taux d'urbanisation regional ------------
aux_urban <- dhs %>%
  group_by(region) %>%
  summarise(urban_pct = 100 * weighted.mean(urban, weight), .groups = "drop")

# --- 5. Variable auxiliaire X2 : NDVI moyen regional (simule) ------------
# Le NDVI est plus eleve dans la zone Sahelienne agricole (sud Maradi, sud Zinder)
# et plus faible dans la zone desertique du Nord (Agadez, nord Tahoua).
# On simule des valeurs realistes dans [0.1, 0.6] a partir du seul taux
# d'urbanisation, avec un peu de bruit. (Le NDVI est negatif dans le desert
# et positif dans la savane ; on reste dans une plage realiste.)
set.seed(123)
ndvi_sim <- function(urban_pct) {
  base <- 0.45 - 0.005 * urban_pct + rnorm(1, 0, 0.04)
  pmax(0.10, pmin(0.65, base))
}
aux_ndvi <- direct %>%
  left_join(aux_urban, by = "region") %>%
  rowwise() %>%
  mutate(ndvi_mean = ndvi_sim(urban_pct)) %>%
  ungroup() %>%
  select(region, ndvi_mean)

# --- 6. Table fusionnee prete pour Fay-Herriot ---------------------------
sae_data <- direct %>%
  left_join(aux_urban, by = "region") %>%
  left_join(aux_ndvi, by = "region") %>%
  select(region, n_sample, direct_mean, se_direct, cv_direct, psi,
         urban_pct, ndvi_mean)

# --- 7. Sauvegardes ------------------------------------------------------
dir.create("outputs/data",   showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

saveRDS(sae_data, "outputs/data/sae_poverty_data.rds")
write.csv(direct %>% select(region, n_sample, direct_mean, se_direct, cv_direct, psi),
          "outputs/tables/direct_estimates.csv", row.names = FALSE)

cat("==> Estimateurs directs calcules pour", nrow(sae_data), "regions\n")
cat("==> CV direct moyen :", round(mean(sae_data$cv_direct), 2), "%\n")
cat("==> CV direct min/max :",
    round(min(sae_data$cv_direct), 2), "/",
    round(max(sae_data$cv_direct), 2), "%\n")
cat("==> Sauvegardes : outputs/data/sae_poverty_data.rds, outputs/tables/direct_estimates.csv\n")
print(sae_data)
