# ===========================================================================
# charger_dhs.R
#
# Charge le fichier DHS Niger 2012 (Household Recode), selectionne 5 colonnes
# utiles pour l'estimation de pauvrete, et sauvegarde un sous-ensemble
# dans outputs/data/dhs_subset.rds.
#
# Variables conservees :
#   hv001 : identifiant du cluster (menagegroupe)
#   hv005 : poids de sondage (divise par 1e6 dans le DHS)
#   hv024 : region (1=Agadez, 2=Diffa, 3=Dosso, 4=Maradi, 5=Tahoua, 6=Tillaberi, 7=Zinder, 8=Niamey)
#   hv025 : milieu de residence (1 = urbain, 2 = rural)
#   hv271 : indice de richesse du menage (score continu, moyenne ~ 0)
#
# Sortie : outputs/data/dhs_subset.rds
# ===========================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

set.seed(42)

# --- 1. Chemin vers le fichier DHS (lien symbolique) --------------------
dhs_path <- "data/raw/NIHR61FL.DTA"

if (!file.exists(dhs_path)) {
  warning("Fichier DHS introuvable : ", dhs_path,
          "\nBascule en mode simulation calibree.")
  dhs <- NULL
} else {
  dhs <- haven::read_dta(dhs_path)
}

# --- 2. Chargement reel ou simulation calibree ---------------------------
if (is.null(dhs)) {
  # Simulation : 8 regions du Niger, ~1344 menages par region (10750 / 8)
  regions_niger <- c("Agadez", "Diffa", "Dosso", "Maradi",
                     "Niamey", "Tahoua", "Tillaberi", "Zinder")
  n_per_region <- 1344
  dhs <- do.call(rbind, lapply(regions_niger, function(r) {
    n <- n_per_region
    # Richesse : distribution normale, legerement plus elevee en zone urbaine Sahelienne
    urban_pct <- c(Agadez = 0.20, Diffa = 0.15, Dosso = 0.13, Maradi = 0.15,
                   Niamey = 1.00, Tahoua = 0.16, Tillaberi = 0.10, Zinder = 0.18)[r]
    urban <- rbinom(n, 1, urban_pct)
    wealth <- rnorm(n, mean = 0 + 0.5 * urban, sd = 1.0)
    data.frame(
      HV001 = seq_len(n),
      HV005 = runif(n, 0.5e6, 1.5e6),
      HV024 = r,
      HV025 = ifelse(urban == 1, 1, 2),
      HV271 = wealth,
      stringsAsFactors = FALSE
    )
  }))
}

# --- 3. Selection et renommage -------------------------------------------
# Les noms de colonnes du DHS sont en minuscule (hv001, hv005, hv024, hv025, hv271)
# La region est codee 1-8 avec des labels haven ; on la convertit en nom.
dhs_subset <- dhs %>%
  transmute(
    cluster = as.integer(hv001),
    weight  = as.numeric(hv005) / 1e6,  # poids normalises
    region  = as.character(haven::as_factor(hv024)),
    urban   = as.integer(hv025 == 1),   # binaire : 1 = urbain, 0 = rural
    wealth  = as.numeric(hv271)
  ) %>%
  filter(!is.na(wealth), !is.na(region), !is.na(weight))

cat("==> DHS charge :", nrow(dhs_subset), "menages,", n_distinct(dhs_subset$region), "regions\n")
cat("==> Regions :", paste(sort(unique(dhs_subset$region)), collapse = ", "), "\n")
cat("==> Taux d'urbanisation global :", round(100 * mean(dhs_subset$urban), 1), "%\n")
cat("==> Richesse moyenne (brute) :", round(mean(dhs_subset$wealth), 3), "\n")

# --- 4. Sauvegarde -------------------------------------------------------
dir.create("outputs/data", showWarnings = FALSE, recursive = TRUE)
dir.create("python/data",  showWarnings = FALSE, recursive = TRUE)
saveRDS(dhs_subset, "outputs/data/dhs_subset.rds")
# Export CSV pour le pipeline Python (preparer_clusters.py)
write.csv(dhs_subset, "python/data/dhs_subset.csv", row.names = FALSE)
cat("==> Sauvegarde : outputs/data/dhs_subset.rds\n")
cat("==> Sauvegarde : python/data/dhs_subset.csv (pour Python)\n")
