summary_high_correlation_6_single <- run_single_traits(
  cv_datasets = kfold_cv_high_correlation_6,
  tree = phylogeny,
  traits = c(
    "Crown_height",
    "Crown_diameter",
    "Tree_height"
  ),
  result_dir =
    here::here("Results/kfold_cv/single_traits_PHY"),
  checkpoint_dir =
    here::here("checkpoints/single_traits"),
  n_cores = 4
)
