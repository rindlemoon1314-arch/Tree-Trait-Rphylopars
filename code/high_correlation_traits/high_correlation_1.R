summary_high_correlation_1_single <- run_single_traits(
  cv_datasets = kfold_cv_high_correlation_1,
  tree = phylogeny,
  traits = c(
    "Leaf_N_per_mass",
    "Specific_leaf_area",
    "Leaf_thickness"
  ),
  result_dir =
    here::here("Results/kfold_cv/single_traits_PHY"),
  checkpoint_dir =
    here::here("checkpoints/single_traits"),
  n_cores = 4
)
