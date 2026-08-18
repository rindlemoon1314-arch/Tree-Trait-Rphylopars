# ============================================================
# 1. Load functions
# ============================================================

source(here::here("code/functions/calculate_phy_metrics.R"))
source(here::here("code/functions/plot_grouped_r2_range.R"))

# ============================================================
# 2. Define the output directory
# ============================================================

r2_output_dir <- paste0(
  here::here("Results/"),
  "R2_plots"
)

dir.create(
  r2_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 3. Calculate R2 and draw the grouped range plot
# ============================================================

grouped_r2_range <- plot_grouped_r2_range(
  single_dir = paste0(
    here::here("Results/kfold_cv/"),
    "single_traits_PHY"
  ),
  all8_single_dir = paste0(
    here::here("Results/kfold_cv/"),
    "single_traits_all8_PHY"
  ),
  all8_dir = paste0(
    here::here("Results/kfold_cv/"),
    "all8_PHY"
  ),
  high_correlation_dir = paste0(
    here::here("Results/kfold_cv/"),
    "high_correlation_traits_PHY"
  ),
  folds = c(
    1, 3, 5, 6, 8, 10
  ),
  negative_compression = 0.05,
  output_file = file.path(
    r2_output_dir,
    "R2_grouped_trait_range.png"
  ),
  data_file = file.path(
    r2_output_dir,
    "R2_grouped_trait_range_data.csv"
  ),
  base_font_size = 9,
  title_font_size = 10.5,
  output_width = 8.25,
  minimum_height = 10.5,
  row_height = 0.27,
  font_family = "sans"
)
