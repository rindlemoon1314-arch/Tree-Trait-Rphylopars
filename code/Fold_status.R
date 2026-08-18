# ============================================================
# 1. Load the plot function
# ============================================================

source(here::here("code/functions/plot_fold_status.R"))

# ============================================================
# 2. Draw and save the plot
# ============================================================

PHY_fold_status_plot <- plot_fold_status(
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
  single_cell_width = 0.78,
  single_cell_height = 0.72,
  multi_cell_width = 0.98,
  multi_cell_height = 0.92,
  output_file = paste0(
    here::here("Results/"),
    "Fold_status_plot.png"
  )
)
