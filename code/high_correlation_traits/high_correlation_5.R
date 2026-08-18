summary_high_correlation_5_single <- run_single_traits(
  cv_datasets = kfold_cv_high_correlation_5,
  tree = phylogeny,
  traits = c(
    "Bark_thickness",
    "Stem_diameter"
  ),
  result_dir =
    here::here("Results/kfold_cv/single_traits_PHY"),
  checkpoint_dir =
    here::here("checkpoints/single_traits"),
  n_cores = 4
)
