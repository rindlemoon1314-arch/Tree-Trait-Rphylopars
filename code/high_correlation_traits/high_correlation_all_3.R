# ------------------------------------------------------------
# PhyloPars analysis: Brownian motion
# ------------------------------------------------------------

PHY_high_correlation_traits_3_BM <- phylopars(
  trait_high_correlation_3,
  phylogeny,
  model = "BM",
  pheno_error = TRUE,
  phylo_correlated = TRUE,
  pheno_correlated = TRUE,
  EM_verbose = TRUE,
  optim_verbose = TRUE
)

summary.phylopars(
  PHY_high_correlation_traits_3_BM
)

head(
  PHY_high_correlation_traits_3_BM$anc_recon
)

head(
  PHY_high_correlation_traits_3_BM$anc_var
)


# ------------------------------------------------------------
# Cross-validation: Brownian motion
# ------------------------------------------------------------

PHY_cv_high_correlation_traits_3_BM <- phylopars(
  trait_high_correlation_3_cv,
  phylogeny,
  model = "BM",
  pheno_error = TRUE,
  phylo_correlated = TRUE,
  pheno_correlated = TRUE,
  EM_verbose = TRUE,
  optim_verbose = TRUE
)

result_high_correlation_traits_3_BM <-
  calculate_cv_performance(
    PHY_cv = PHY_cv_high_correlation_traits_3_BM,
    trait_name = c(
      "Leaf_K_per_mass",
      "Leaf_P_per_mass"
    ),
    testset = testset_high_correlation_3,
    model_name = "BM"
  )

model_performance_high_correlation_traits_3 <-
  result_high_correlation_traits_3_BM$performance

model_performance_high_correlation_traits_3


# ------------------------------------------------------------
# K-fold cross-validation: Brownian motion
# Configuration: BM_FTF
# ------------------------------------------------------------

PHY_kfoldcv_high_correlation_traits_3_BM_FTF <-
  run_kfold_phylopars_parallel(
    cv_datasets = kfold_cv_high_correlation_3,
    tree = phylogeny,
    model = "BM",
    pheno_error = FALSE,
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    usezscores = FALSE,
    EM_verbose = FALSE,
    optim_verbose = FALSE,
    n_cores = 2,
    checkpoint_dir = paste0(
      here::here("checkpoints/"),
      "kfold_high_correlation_traits_3_BM_FTF"
    ),
    resume = TRUE,
    retry_failed = FALSE,
    worker_stop_grace = 0.5,
    result_components = c(
      "anc_recon",
      "anc_var"
    )
  )


# ------------------------------------------------------------
# Save all four high-correlation-traits-3 PHY results
# Run this section after all four objects have been calculated
# ------------------------------------------------------------

high_correlation_traits_3_objects <- c(
  "PHY_kfoldcv_high_correlation_traits_3_BM_FTF",
  "PHY_kfoldcv_high_correlation_traits_3_BM_TTF",
  "PHY_kfoldcv_high_correlation_traits_3_lambda_FTF",
  "PHY_kfoldcv_high_correlation_traits_3_lambda_TTF"
)

saved_high_correlation_traits_3 <- save_PHY_results(
  object_names = high_correlation_traits_3_objects,
  output_dir =
    here::here("Results/kfold_cv/high_correlation_traits_PHY")
)

saved_high_correlation_traits_3


# ------------------------------------------------------------
# Calculate K-fold performance
# ------------------------------------------------------------

kfold_model_performance_high_correlation_traits_3_BM_FTF <-
  calculate_kfold_performance(
    PHY_kfoldcv_high_correlation_traits_3_BM_FTF,
    data_scale = "scaled"
  )

kfold_model_performance_high_correlation_traits_3_BM_FTF


# ------------------------------------------------------------
# Plot K-fold performance
# ------------------------------------------------------------

cv_plots_high_correlation_traits_3_BM_FTF <- plot_cv(
  kfold_model_performance_high_correlation_traits_3_BM_FTF
)

cv_plots_high_correlation_traits_3_BM_FTF

