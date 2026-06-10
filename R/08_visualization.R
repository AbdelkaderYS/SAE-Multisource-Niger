# =============================================================================
# Project: Multisource SAE for SDG Monitoring in Niger
# Script:  08_visualization.R
# Author:  Abdel Kader Younoussi Saley
# =============================================================================

source("R/00_setup.R")

results <- readRDS(here("data", "processed", "model_comparison_results.rds"))
fits    <- readRDS(here("data", "processed", "model_fits.rds"))

pub_theme <- theme_classic(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(color = "grey40", size = 10),
    plot.caption      = element_text(color = "grey50", size = 8),
    axis.line         = element_line(color = "grey30")
  )

MODEL_COLORS <- c(
  "Direct"         = "grey40",
  "M1: Census"     = "#2196F3",
  "M2: +Admin"     = "#4CAF50",
  "M3: +Satellite" = "#FF9800",
  "M4: Full"       = "#9C27B0"
)

# =============================================================================
# PLOT 1: Direct vs EBLUP estimates with 95% CI
# =============================================================================

p1 <- results %>%
  mutate(region = fct_reorder(region, direct_mpi)) %>%
  ggplot() +
  geom_pointrange(
    aes(x = region,
        y = direct_mpi,
        ymin = direct_mpi - 1.96 * se_direct,
        ymax = direct_mpi + 1.96 * se_direct,
        color = "Direct"),
    size = 0.7, linewidth = 1
  ) +
  geom_pointrange(
    aes(x = region,
        y = eblup_m4,
        ymin = eblup_m4 - 1.96 * se_m4,
        ymax = eblup_m4 + 1.96 * se_m4,
        color = "M4: Full"),
    size = 0.7, linewidth = 1,
    position = position_nudge(x = 0.25)
  ) +
  scale_color_manual(values = MODEL_COLORS, name = "Estimator") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.4, 1.05)) +
  coord_flip() +
  labs(
    title   = "MPI Poverty Estimates: Direct vs Full Multisource SAE",
    subtitle = "Niger -- 8 Administrative Regions | 95% Confidence Intervals",
    x       = NULL,
    y       = "MPI Poverty Rate",
    caption = "Data: DHS Niger 2012 | Auxiliaries: RGPH 2012 + INS Niger + MODIS/VIIRS"
  ) +
  pub_theme

ggsave(here("outputs", "figures", "01_direct_vs_eblup.png"),
       p1, width = 9, height = 6, dpi = 200, bg = "white")
cat("Plot 1 saved\n")

# =============================================================================
# PLOT 2: CV comparison -- zoomed on actual range
# Note: CVs are low at regional level by DHS design. SAE is critical
# at the departement level (n~30-50 HH, CV > 25%).
# =============================================================================

p2 <- results %>%
  select(region, cv_direct, cv_m1, cv_m2, cv_m3, cv_m4) %>%
  pivot_longer(-region, names_to = "model", values_to = "cv") %>%
  mutate(
    model = recode(model,
      "cv_direct" = "Direct",
      "cv_m1"     = "M1: Census",
      "cv_m2"     = "M2: +Admin",
      "cv_m3"     = "M3: +Satellite",
      "cv_m4"     = "M4: Full"
    ),
    model  = factor(model, levels = names(MODEL_COLORS)),
    region = fct_reorder(region, -cv)
  ) %>%
  ggplot(aes(x = region, y = cv, fill = model)) +
  geom_col(position = "dodge", alpha = 0.9, width = 0.7) +
  scale_fill_manual(values = MODEL_COLORS, name = "Model") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Coefficient of Variation Across Multisource SAE Models",
    subtitle = "All regions below 20% publishability threshold | SAE critical at departement level",
    x        = "Region",
    y        = "Coefficient of Variation (%)",
    caption  = "DHS 2012 representative at regional level by design"
  ) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here("outputs", "figures", "02_cv_comparison.png"),
       p2, width = 11, height = 6, dpi = 200, bg = "white")
cat("Plot 2 saved\n")

# =============================================================================
# PLOT 3: Shrinkage factors
# =============================================================================

p3 <- results %>%
  select(region, gamma_m1, gamma_m4) %>%
  pivot_longer(-region, names_to = "model", values_to = "gamma") %>%
  mutate(
    model  = recode(model,
      "gamma_m1" = "M1: Census only",
      "gamma_m4" = "M4: Full multisource"
    ),
    region = fct_reorder(region, gamma)
  ) %>%
  ggplot(aes(x = region, y = gamma, fill = model)) +
  geom_col(position = "dodge", alpha = 0.9, width = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  annotate("text", x = 0.6, y = 0.53,
           label = "gamma = 0.5", hjust = 0, size = 3, color = "grey30") +
  scale_fill_manual(
    values = c("M1: Census only" = "#2196F3",
               "M4: Full multisource" = "#9C27B0"),
    name   = "Model"
  ) +
  scale_y_continuous(limits = c(0, 1.05)) +
  coord_flip() +
  labs(
    title    = "Shrinkage Factor by Region and Model",
    subtitle = "gamma -> 1: trust direct | gamma -> 0: borrow from model",
    x        = NULL,
    y        = "Shrinkage factor (gamma_i)",
    caption  = "gamma_i = sigma2_v / (sigma2_v + psi_i)"
  ) +
  pub_theme

ggsave(here("outputs", "figures", "03_shrinkage.png"),
       p3, width = 9, height = 6, dpi = 200, bg = "white")
cat("Plot 3 saved\n")

# =============================================================================
# PLOT 4: Between-area variance component by model
# sigma2_v measures how much between-area heterogeneity the model captures.
# Lower sigma2_v = auxiliary variables better explain spatial variation.
# =============================================================================

sigma2_df <- tibble(
  Model   = c("M1: Census", "M2: +Admin", "M3: +Satellite", "M4: Full"),
  sigma2v = c(
    ifelse(is.null(fits$fit1), NA, fits$fit1$fit$refvar),
    ifelse(is.null(fits$fit2), NA, fits$fit2$fit$refvar),
    ifelse(is.null(fits$fit3), NA, fits$fit3$fit$refvar),
    ifelse(is.null(fits$fit4), NA, fits$fit4$fit$refvar)
  )
) %>%
  filter(!is.na(sigma2v)) %>%
  mutate(Model = factor(Model, levels = c("M1: Census", "M2: +Admin",
                                          "M3: +Satellite", "M4: Full")))

p4 <- sigma2_df %>%
  ggplot(aes(x = reorder(Model, sigma2v), y = sigma2v * 1000, fill = Model)) +
  geom_col(alpha = 0.9, show.legend = FALSE, width = 0.5) +
  geom_text(aes(label = sprintf("%.3f", sigma2v * 1000)),
            hjust = -0.2, size = 4) +
  scale_fill_manual(values = c(
    "M1: Census"     = "#2196F3",
    "M2: +Admin"     = "#4CAF50",
    "M3: +Satellite" = "#FF9800",
    "M4: Full"       = "#9C27B0"
  )) +
  scale_y_continuous(limits = c(0, max(sigma2_df$sigma2v * 1000) * 1.4)) +
  coord_flip() +
  labs(
    title    = "Between-Area Variance Component by Model",
    subtitle = "Lower sigma2_v = auxiliary variables better explain spatial heterogeneity",
    x        = NULL,
    y        = "sigma2_v (x 1000)",
    caption  = "REML estimation | Fay-Herriot area-level model"
  ) +
  pub_theme

ggsave(here("outputs", "figures", "04_sigma2v_by_model.png"),
       p4, width = 9, height = 5, dpi = 200, bg = "white")
cat("Plot 4 saved\n")

cat("\nAll plots saved to outputs/figures/\n")
cat("Project 2 complete.\n")
