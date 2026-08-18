summary_high_correlation_2_single <- run_single_traits(
  cv_datasets = kfold_cv_high_correlation_2,
  tree = phylogeny,
  traits = c(
    "Stem_conduit_diameter",
    "Leaf_Vcmax_per_dry_mass",
    "Stomatal_conductance",
    "Leaf_area"
  ),
  result_dir =
    here::here("Results/kfold_cv/single_traits_PHY"),
  checkpoint_dir =
    here::here("checkpoints/single_traits"),
  n_cores = 4
)
