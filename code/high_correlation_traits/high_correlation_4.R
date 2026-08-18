summary_high_correlation_4_single <- run_single_traits(
  cv_datasets = kfold_cv_high_correlation_4,
  tree = phylogeny,
  traits = c(
    "Leaf_density",
    "Wood_density"
  ),
  result_dir =
    here::here("Results/kfold_cv/single_traits_PHY"),
  checkpoint_dir =
    here::here("checkpoints/single_traits"),
  n_cores = 4
)
