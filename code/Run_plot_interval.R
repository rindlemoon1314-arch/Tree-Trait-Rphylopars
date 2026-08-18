source(here::here("code/functions/plot_interval.R"))

interval_results <- plot_single_lambda_ttf_interval_calibration(
  single_dir = file.path(
    here::here("Results/kfold_cv"),
    c(
      "single_traits_PHY",
      "single_traits_all8_PHY"
    )
  ),
  folds = c(1, 3, 5, 6, 8, 10),
  interval = 95,
  use_transformed_truth = TRUE,
  traits = c(
    "Bark_thickness",
    "Crown_diameter",
    "Crown_height",
    "Leaf_area",
    "Leaf_density",
    "Leaf_K_per_mass",
    "Leaf_N_per_mass",
    "Leaf_P_per_mass",
    "Leaf_thickness",
    "Leaf_Vcmax_per_dry_mass",
    "Specific_leaf_area",
    "Stem_conduit_diameter",
    "Stem_diameter",
    "Stomatal_conductance",
    "Tree_height",
    "Wood_density",
    "Root_depth",
    "Seed_dry_mass"
  ),
  traits_per_page = 18,
  max_points_per_trait = Inf,
  output_dir = paste0(
    here::here("Results/"),
    "prediction_interval_plots"
  ),
  file_prefix = "single_lambda_TTF_95_interval_calibration_18_traits",
  width = 8.25,
  height = 10.8,
  dpi = 300,
  base_size = 9.5,
  horizontal_padding = 0.005
)

View(interval_results$coverage_summary)
print(interval_results$saved_files)
