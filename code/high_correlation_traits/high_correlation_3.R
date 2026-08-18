summary_high_correlation_3_single <- run_single_traits(
  cv_datasets = kfold_cv_high_correlation_3,
  tree = phylogeny,
  traits = c(
    "Leaf_K_per_mass",
    "Leaf_P_per_mass"
  ),
  result_dir =
    here::here("Results/kfold_cv/single_traits_PHY"),
  checkpoint_dir =
    here::here("checkpoints/single_traits"),
  n_cores = 4
)
