# ============================================================
# 1. Load the plotting function
# ============================================================

source(here::here("code/functions/TTF_FTF_r2.R"))

# ============================================================
# 2. Define input and output paths
# ============================================================

r2_data_file <- paste0(
  here::here("Results/R2_plots/"),
  "R2_grouped_trait_range_data.csv"
)

comparison_output_dir <- paste0(
  here::here("Results/"),
  "R2_plots"
)

dir.create(
  comparison_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Draw and save the comparison plots
# ============================================================

ttf_ftf_r2 <- plot_ttf_ftf_r2_comparison(
  r2_data = r2_data_file,
  negative_compression = 0.05,
  delta_fill_limit = 0.05,
  output_file = file.path(
    comparison_output_dir,
    "R2_FTF_vs_TTF_rectangular_dumbbell.png"
  ),
  heatmap_file = file.path(
    comparison_output_dir,
    "R2_FTF_vs_TTF_delta_heatmap.png"
  ),
  combined_file = file.path(
    comparison_output_dir,
    "R2_FTF_vs_TTF_combined.png"
  ),
  detailed_table_file = file.path(
    comparison_output_dir,
    "R2_FTF_vs_TTF_detailed.csv"
  ),
  summary_table_file = file.path(
    comparison_output_dir,
    "R2_FTF_vs_TTF_summary.csv"
  ),
  width = 6.5,
  height = 9.0,
  combined_width = 8.25,
  combined_height = 11.0,
  dpi = 300,
  base_font_size = 7.2,
  title_font_size = 8.5,
  panel_title_font_size = 13.5,
  trait_spacing_multiplier = 1
)


# ============================================================
# 4. Inspect the result tables
# ============================================================

View(ttf_ftf_r2$detailed_table)
View(ttf_ftf_r2$summary_table)
