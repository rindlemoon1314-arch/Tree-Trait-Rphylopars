summary_32 <- run_single_traits(
  cv_datasets = kfold_cv_all8,
  tree = phylogeny,
  traits = c(
    "Wood_density",
    "Specific_leaf_area",
    "Seed_dry_mass",
    "Leaf_P_per_mass",
    "Stem_conduit_diameter",
    "Tree_height",
    "Root_depth",
    "Bark_thickness"
  ),
  result_dir =
    here::here("Results/kfold_cv/single_traits_all8_PHY"),
  checkpoint_dir =
    here::here("checkpoints/single_traits_all8"),
  n_cores = 4
)
