# =============================================================================
# Project: Multisource SAE for SDG Monitoring in Niger
# Script:  00_setup.R
# Author:  Abdel Kader Younoussi Saley
# =============================================================================

packages <- c(
  "sae",       # eblupFH, mseFH, eblupBHF
  "fastLink",  # probabilistic record linkage
  "survey",    # complex survey estimation
  "srvyr",     # tidyverse survey interface
  "tidyverse", # data manipulation and plots
  "sf",        # spatial data
  "here",      # relative paths
  "httr",      # HTTP requests for WB API
  "jsonlite"   # JSON parsing
)

library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

new_pkg <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_pkg) > 0) renv::install(new_pkg)

lapply(packages, library, character.only = TRUE)

set.seed(42)
options(scipen = 999)
theme_set(theme_classic(base_size = 12))

REGIONS <- c("Agadez", "Diffa", "Dosso", "Maradi",
             "Tahoua", "Tillaberi", "Zinder", "Niamey")

dirs <- c(
  here("data", "processed"),
  here("data", "raw"),
  here("outputs", "figures"),
  here("outputs", "tables")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

cat("Setup complete.\n")
