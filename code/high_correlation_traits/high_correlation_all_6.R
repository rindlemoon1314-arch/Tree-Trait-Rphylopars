# ------------------------------------------------------------
# PhyloPars analysis: Brownian motion
# ------------------------------------------------------------

PHY_high_correlation_traits_6_BM <- phylopars(
  trait_high_correlation_6,
  phylogeny,
  model = "BM",
  pheno_error = TRUE,
  phylo_correlated = TRUE,
  pheno_correlated = TRUE,
  EM_verbose = TRUE,
  optim_verbose = TRUE
)

summary.phylopars(
  PHY_high_correlation_traits_6_BM
)

head(
  PHY_high_correlation_traits_6_BM$anc_recon
)

head(
  PHY_high_correlation_traits_6_BM$anc_var
)


# ------------------------------------------------------------
# Cross-validation: Brownian motion
# ------------------------------------------------------------

PHY_cv_high_correlation_traits_6_BM <- phylopars(
  trait_high_correlation_6_cv,
  phylogeny,
  model = "BM",
  pheno_error = TRUE,
  phylo_correlated = TRUE,
  pheno_correlated = TRUE,
  EM_verbose = TRUE,
  optim_verbose = TRUE
)

result_high_correlation_traits_6_BM <-
  calculate_cv_performance(
    PHY_cv = PHY_cv_high_correlation_traits_6_BM,
    trait_name = c(
      "Crown_height",
      "Crown_diameter",
      "Tree_height"
    ),
    testset = testset_high_correlation_6,
    model_name = "BM"
  )

model_performance_high_correlation_traits_6 <-
  result_high_correlation_traits_6_BM$performance

model_performance_high_correlation_traits_6


# ------------------------------------------------------------
# K-fold cross-validation: Brownian motion
# Configuration: BM_FTF
# ------------------------------------------------------------

PHY_kfoldcv_high_correlation_traits_6_BM_FTF <-
  run_kfold_phylopars_parallel(
    cv_datasets = kfold_cv_high_correlation_6,
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
      "kfold_high_correlation_traits_6_BM_FTF"
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
# K-fold cross-validation: Brownian motion
# Configuration: BM_TTF
# ------------------------------------------------------------

PHY_kfoldcv_high_correlation_traits_6_BM_TTF <-
  run_kfold_phylopars_parallel(
    cv_datasets = kfold_cv_high_correlation_6,
    tree = phylogeny,
    model = "BM",
    pheno_error = TRUE,
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    usezscores = FALSE,
    EM_verbose = FALSE,
    optim_verbose = FALSE,
    n_cores = 2,
    checkpoint_dir = paste0(
      here::here("checkpoints/"),
      "kfold_high_correlation_traits_6_BM_TTF"
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
# K-fold cross-validation: Pagel's lambda
# Configuration: lambda_FTF
# ------------------------------------------------------------

PHY_kfoldcv_high_correlation_traits_6_lambda_FTF <-
  run_kfold_phylopars_parallel(
    cv_datasets = kfold_cv_high_correlation_6,
    tree = phylogeny,
    model = "lambda",
    pheno_error = FALSE,
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    usezscores = FALSE,
    EM_verbose = FALSE,
    optim_verbose = FALSE,
    n_cores = 2,
    checkpoint_dir = paste0(
      here::here("checkpoints/"),
      "kfold_high_correlation_traits_6_lambda_FTF"
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
# K-fold cross-validation: Pagel's lambda
# Configuration: lambda_TTF
# ------------------------------------------------------------

PHY_kfoldcv_high_correlation_traits_6_lambda_TTF <-
  run_kfold_phylopars_parallel(
    cv_datasets = kfold_cv_high_correlation_6,
    tree = phylogeny,
    model = "lambda",
    pheno_error = TRUE,
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    usezscores = FALSE,
    EM_verbose = FALSE,
    optim_verbose = FALSE,
    n_cores = 2,
    checkpoint_dir = paste0(
      here::here("checkpoints/"),
      "kfold_high_correlation_traits_6_lambda_TTF"
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
# Save all four high-correlation-traits-6 PHY results
# Run this section after all four objects have been calculated
# ------------------------------------------------------------

high_correlation_traits_6_objects <- c(
  "PHY_kfoldcv_high_correlation_traits_6_BM_FTF",
  "PHY_kfoldcv_high_correlation_traits_6_BM_TTF",
  "PHY_kfoldcv_high_correlation_traits_6_lambda_FTF",
  "PHY_kfoldcv_high_correlation_traits_6_lambda_TTF"
)

saved_high_correlation_traits_6 <- save_PHY_results(
  object_names = high_correlation_traits_6_objects,
  output_dir =
    here::here("Results/kfold_cv/high_correlation_traits_PHY")
)

saved_high_correlation_traits_6


# ------------------------------------------------------------
# Calculate K-fold performance
# ------------------------------------------------------------

kfold_model_performance_high_correlation_traits_6_BM_FTF <-
  calculate_kfold_performance(
    PHY_kfoldcv_high_correlation_traits_6_BM_FTF,
    data_scale = "scaled"
  )

kfold_model_performance_high_correlation_traits_6_BM_FTF


# ------------------------------------------------------------
# Plot K-fold performance
# ------------------------------------------------------------

cv_plots_high_correlation_traits_6_BM_FTF <- plot_cv(
  kfold_model_performance_high_correlation_traits_6_BM_FTF
)

cv_plots_high_correlation_traits_6_BM_FTF
