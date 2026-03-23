# =============================================================================
# Project: Multisource SAE for SDG Monitoring in Niger
# Script:  01_data_sources.R
# Author:  Abdel Kader Younoussi Saley
# =============================================================================
# Data sources:
#   DHS Niger 2012 (NIHR61FL.DTA) -- dhsprogram.com, free registration
#   World Bank Open Data API       -- no authentication required
#   INS Niger RGPH 2012            -- calibrated from published reports
#   MODIS NDVI + VIIRS             -- calibrated from NASA Earthdata
#
# MPI: Alkire-Foster (2011)
#   Three dimensions: Health (1/3), Education (1/3), Living Standards (1/3)
#   Poverty cutoff k = 1/3
# =============================================================================

source("R/00_setup.R")

USE_REAL_DATA <- file.exists(here("data", "raw", "NIHR61FL.DTA"))

if (USE_REAL_DATA) {
  cat("DHS Niger 2012 found -- loading real data\n\n")
} else {
  cat("DHS file not found in data/raw/\n")
  cat("Place NIHR61FL.DTA in data/raw/ to use real data\n")
  cat("Using simulated data calibrated from OPHI 2025\n\n")
}

# =============================================================================
# PART A: WORLD BANK API -- national indicators for Niger
# =============================================================================
# These replace hardcoded national values with real API data.
# All indicators are national level only (WB does not provide subnational).
# Fallback values used if API unavailable.

get_wb_value <- function(indicator, country = "NER",
                         start = 2010, end = 2022) {
  url <- paste0(
    "https://api.worldbank.org/v2/country/", country,
    "/indicator/", indicator,
    "?format=json&date=", start, ":", end, "&per_page=100"
  )
  tryCatch({
    resp <- GET(url, timeout(10))
    if (status_code(resp) != 200) return(NA_real_)
    data <- content(resp, "parsed", simplifyVector = FALSE)
    if (length(data) < 2 || is.null(data[[2]])) return(NA_real_)
    vals <- keep(data[[2]], ~!is.null(.x$value))
    if (length(vals) == 0) return(NA_real_)
    vals <- vals[order(sapply(vals, function(x) x$date), decreasing = TRUE)]
    as.numeric(vals[[1]]$value)
  }, error = function(e) NA_real_)
}

cat("Fetching World Bank indicators for Niger...\n")

wb_indicators <- list(
  literacy_pct      = get_wb_value("SE.ADT.LITR.ZS"),
  electricity_pct   = get_wb_value("EG.ELC.ACCS.ZS"),
  safe_water_pct    = get_wb_value("SH.H2O.BASW.ZS"),
  child_mort_per1k  = get_wb_value("SH.DYN.MORT"),
  gdp_per_capita    = get_wb_value("NY.GDP.PCAP.CD"),
  poverty_pct       = get_wb_value("SI.POV.NAHC")
)

# Fallback values calibrated from World Bank Niger data portal
# Used when API is unavailable
fallback <- list(
  literacy_pct      = 37.3,
  electricity_pct   = 19.5,
  safe_water_pct    = 59.8,
  child_mort_per1k  = 85.2,
  gdp_per_capita    = 597.0,
  poverty_pct       = 40.8
)

for (nm in names(fallback)) {
  if (is.na(wb_indicators[[nm]])) {
    wb_indicators[[nm]] <- fallback[[nm]]
    cat("  API unavailable for", nm, "-- using fallback\n")
  } else {
    cat("  WB API:", nm, "=", round(wb_indicators[[nm]], 1), "\n")
  }
}

# =============================================================================
# PART B: DHS DATA LOADING AND MPI CONSTRUCTION
# =============================================================================

if (USE_REAL_DATA) {

  hh_raw <- haven::read_dta(here("data", "raw", "NIHR61FL.DTA"))
  names(hh_raw) <- toupper(names(hh_raw))
  cat("\nDHS loaded:", nrow(hh_raw), "households\n")

  # Aggregate member-level variables to household level
  educ_cols   <- grep("^HV106_", names(hh_raw), value = TRUE)
  age_cols    <- grep("^HV105_", names(hh_raw), value = TRUE)
  attend_cols <- grep("^HV121_", names(hh_raw), value = TRUE)

  hh_raw$HV106_HEAD <- as.numeric(hh_raw$HV106_01)
  hh_raw$HV110_HEAD <- as.numeric(hh_raw$HV110_01)

  educ_mat <- hh_raw %>%
    select(all_of(educ_cols)) %>%
    mutate(across(everything(), as.numeric))
  hh_raw$any_no_educ <- as.integer(
    apply(educ_mat, 1, function(x) any(x == 0, na.rm = TRUE))
  )

  age_mat    <- hh_raw %>% select(all_of(age_cols))    %>%
    mutate(across(everything(), as.numeric))
  attend_mat <- hh_raw %>% select(all_of(attend_cols)) %>%
    mutate(across(everything(), as.numeric))

  n_members <- min(ncol(age_mat), ncol(attend_mat))
  school_not_attending <- apply(
    cbind(age_mat[, 1:n_members], attend_mat[, 1:n_members]),
    1,
    function(x) {
      n       <- length(x) / 2
      ages    <- x[1:n]
      attends <- x[(n + 1):(2 * n)]
      school_age <- !is.na(ages) & ages >= 6 & ages <= 17
      if (!any(school_age)) return(0L)
      as.integer(any(school_age & !is.na(attends) & attends == 0))
    }
  )
  hh_raw$school_not_attending <- as.integer(school_not_attending)

  region_labels <- c(
    "1" = "Agadez", "2" = "Diffa",    "3" = "Dosso",   "4" = "Maradi",
    "5" = "Tahoua", "6" = "Tillaberi","7" = "Zinder",  "8" = "Niamey"
  )

  dhs_data <- hh_raw %>%
    mutate(
      weight = HV005 / 1000000,
      region = recode(as.character(HV024), !!!region_labels),
      # Health dimension
      h1_child_mortality = as.integer(HV270 == 1),
      h2_nutrition       = as.integer(HV270 <= 2),
      h3_healthcare      = as.integer(HV244 == 1),
      # Education dimension
      e1_years_schooling  = as.integer(HV106_HEAD == 0),
      e2_school_attendance = school_not_attending,
      e3_literacy         = any_no_educ,
      # Living standards dimension
      l1_electricity  = as.integer(HV206 == 0),
      l2_water     = as.integer(HV201 %in% c(30, 32, 40, 42, 43, 96)),
      l3_sanitation = as.integer(HV205 %in% c(14, 23, 31, 42, 43, 96)),
      #l2_water        = as.integer(HV201 %in% c(32, 42, 43, 61, 62, 96)),
      #l3_sanitation   = as.integer(HV205 %in% c(22, 23, 31, 96)),
      l4_cooking_fuel = as.integer(HV226 %in% 6:11 & HV226 != 95),
      l5_assets       = as.integer(HV270 <= 1)
    ) %>%
    filter(!is.na(region)) %>%
    select(region, weight,
           h1_child_mortality, h2_nutrition, h3_healthcare,
           e1_years_schooling, e2_school_attendance, e3_literacy,
           l1_electricity, l2_water, l3_sanitation, l4_cooking_fuel, l5_assets)

  cat("Real DHS processed:", nrow(dhs_data), "households\n\n")

} else {

  # Simulated data calibrated from OPHI 2025 (Niger MIS 2021)
  # H = 79.9%, MPI = 0.535
  n_per_region <- c(850, 780, 1250, 1980, 1650, 1720, 1931, 1000)
  mpi_by_region <- c(
    Agadez = 0.82, Diffa = 0.87, Dosso = 0.84, Maradi = 0.89,
    Tahoua = 0.86, Tillaberi = 0.85, Zinder = 0.88, Niamey = 0.40
  )

  dhs_data <- map2_dfr(REGIONS, n_per_region, function(reg, n) {
    p <- mpi_by_region[reg]
    tibble(
      region               = reg,
      weight               = runif(n, 0.85, 1.15),
      h1_child_mortality   = rbinom(n, 1, p * 0.35),
      h2_nutrition         = rbinom(n, 1, p * 0.55),
      h3_healthcare        = rbinom(n, 1, p * 0.50),
      e1_years_schooling   = rbinom(n, 1, p * 0.78),
      e2_school_attendance = rbinom(n, 1, p * 0.72),
      e3_literacy          = rbinom(n, 1, p * 0.80),
      l1_electricity       = rbinom(n, 1, p * 0.92),
      l2_water             = rbinom(n, 1, p * 0.62),
      l3_sanitation        = rbinom(n, 1, p * 0.83),
      l4_cooking_fuel      = rbinom(n, 1, p * 0.88),
      l5_assets            = rbinom(n, 1, p * 0.75)
    )
  })
}

# =============================================================================
# PART C: MPI CONSTRUCTION (Alkire-Foster 2011)
# Cutoff k = 1/3, equal weights across dimensions
# =============================================================================

dhs_data <- dhs_data %>%
  mutate(
    dep_health = (h1_child_mortality + h2_nutrition + h3_healthcare) / 3,
    dep_educ   = (e1_years_schooling + e2_school_attendance + e3_literacy) / 3,
    dep_living = (l1_electricity + l2_water + l3_sanitation +
                  l4_cooking_fuel + l5_assets) / 5,
    deprivation_score = (dep_health + dep_educ + dep_living) / 3,
    is_MPI_poor       = as.integer(deprivation_score >= 1/3),
    intensity         = if_else(is_MPI_poor == 1, deprivation_score, NA_real_)
  )

nat_H   <- weighted.mean(dhs_data$is_MPI_poor, dhs_data$weight)
nat_A   <- weighted.mean(dhs_data$intensity,   dhs_data$weight, na.rm = TRUE)

cat("National MPI Statistics:\n")
cat("  H (incidence): ", round(nat_H * 100, 1), "%\n")
cat("  A (intensity): ", round(nat_A * 100, 1), "%\n")
cat("  MPI = H x A:   ", round(nat_H * nat_A, 3), "\n")

# =============================================================================
# PART D: DIRECT SURVEY ESTIMATES BY REGION
# =============================================================================

dhs_design <- svydesign(
  ids     = ~1,
  strata  = ~region,
  weights = ~weight,
  data    = dhs_data
)

direct <- svyby(
  formula = ~is_MPI_poor + deprivation_score,
  by      = ~region,
  design  = dhs_design,
  FUN     = svymean,
  na.rm   = TRUE
) %>%
  as_tibble() %>%
  rename(
    direct_mpi   = is_MPI_poor,
    direct_score = deprivation_score,
    se_direct    = se.is_MPI_poor,
    se_score     = se.deprivation_score
  ) %>%
  mutate(
    cv_direct = se_direct / direct_mpi * 100,
    psi       = se_direct^2,
    n_sample  = if (USE_REAL_DATA) {
      dhs_data %>%
        count(region) %>%
        arrange(match(region, REGIONS)) %>%
        pull(n)
    } else {
      c(850, 780, 1250, 1980, 1650, 1720, 1931, 1000)
    }
  )

cat("\nDirect Estimates by Region:\n")
print(direct %>%
  select(region, n_sample, direct_mpi, se_direct, cv_direct) %>%
  mutate(across(where(is.numeric), ~round(., 3))))

# =============================================================================
# PART E: AUXILIARY VARIABLES
# =============================================================================
# Census: calibrated from INS Niger RGPH 2012 + DHS 2012 Final Report
# Satellite: calibrated from MODIS MOD13A3 + VIIRS VNP46A3
# Admin: calibrated from INS Niger ministry publications
# National WB indicators: from API above (real or fallback)

census_aux <- tibble(
  region            = REGIONS,
  pct_urban         = c(0.22, 0.14, 0.12, 0.14, 0.12, 0.11, 0.19, 0.98),
  pct_literate      = c(0.35, 0.19, 0.18, 0.17, 0.16, 0.15, 0.18, 0.65),
  pct_electricity   = c(0.20, 0.10, 0.07, 0.08, 0.07, 0.06, 0.10, 0.78),
  pct_safe_water    = c(0.58, 0.46, 0.44, 0.43, 0.45, 0.42, 0.44, 0.88),
  pct_sanitation    = c(0.28, 0.14, 0.10, 0.12, 0.11, 0.09, 0.13, 0.62),
  pct_skilled_birth = c(0.42, 0.25, 0.30, 0.28, 0.29, 0.27, 0.31, 0.92),
  child_mort_per1k  = c(68,   89,   82,   95,   88,   91,   93,   42)
)

satellite_aux <- tibble(
  region         = REGIONS,
  ndvi_mean      = c(0.25, 0.18, 0.42, 0.38, 0.31, 0.39, 0.35, 0.22),
  viirs_mean     = c(0.15, 0.12, 0.28, 0.35, 0.22, 0.25, 0.31, 2.85),
  dist_market_km = c(48, 62, 28, 22, 35, 30, 25, 8)
)

admin_aux <- tibble(
  region              = REGIONS,
  infant_mort_per1k   = c(68, 89, 82, 95, 88, 91, 93, 42),
  teacher_pupil_ratio = c(42, 58, 51, 62, 55, 59, 57, 28),
  health_fac_per100k  = c(3.2, 2.1, 2.8, 2.3, 2.5, 2.2, 2.4, 8.1),
  admin_completeness  = c(78, 65, 71, 68, 70, 67, 69, 95)
)

# Add national WB indicators (broadcast to all regions as national context)
# Regional variation is captured by the region-level census/satellite variables
wb_aux <- tibble(
  region              = REGIONS,
  wb_literacy_nat     = wb_indicators$literacy_pct,
  wb_electricity_nat  = wb_indicators$electricity_pct,
  wb_child_mort_nat   = wb_indicators$child_mort_per1k,
  wb_gdp_per_cap      = wb_indicators$gdp_per_capita
)

# =============================================================================
# PART F: MERGE ALL SOURCES
# =============================================================================

sae_data <- direct %>%
  left_join(census_aux,    by = "region") %>%
  left_join(satellite_aux, by = "region") %>%
  left_join(admin_aux,     by = "region") %>%
  left_join(wb_aux,        by = "region")

cat("\nSAE dataset:", nrow(sae_data), "regions x", ncol(sae_data), "variables\n")
cat("Missing values:", sum(is.na(sae_data)), "\n")

# =============================================================================
# SAVE
# =============================================================================

saveRDS(dhs_data,     here("data", "processed", "dhs_household_data.rds"))
saveRDS(direct,       here("data", "processed", "direct_estimates.rds"))
saveRDS(census_aux,   here("data", "processed", "census_auxiliaries.rds"))
saveRDS(satellite_aux,here("data", "processed", "satellite_auxiliaries.rds"))
saveRDS(admin_aux,    here("data", "processed", "admin_auxiliaries.rds"))
saveRDS(sae_data,     here("data", "processed", "sae_full_data.rds"))

cat("Script 01 complete.\n")
cat("Next: source('R/02_record_linkage.R')\n")
