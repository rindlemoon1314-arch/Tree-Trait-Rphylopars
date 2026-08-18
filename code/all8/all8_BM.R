# ------------------------------------------------------------
# PhyloPars analysis: brownian motion
# ------------------------------------------------------------
PHY_all8_BM <- phylopars(trait_all8,
                         phylogeny,
                         model = "BM",
                         pheno_error = TRUE,
                         phylo_correlated = TRUE,
                         pheno_correlated = TRUE,
                         EM_verbose = TRUE,
                         optim_verbose = TRUE)

summary.phylopars(PHY_all8_BM)
head(PHY_all8_BM$anc_recon)
head(PHY_all8_BM$anc_var)

# ------------------------------------------------------------
# Cross-validation: brownian motion
# ------------------------------------------------------------
PHY_cv_all8_BM <- phylopars(trait_all8_cv,
                            phylogeny,
                            model = "BM",
                            pheno_error = TRUE,
                            phylo_correlated = TRUE,
                            pheno_correlated = TRUE,
                            EM_verbose = TRUE,
                            optim_verbose = TRUE)

result_all8_BM <- calculate_cv_performance(
  PHY_cv = PHY_cv_all8_BM,
  trait_name = c(
    "Wood_density",
    "Specific_leaf_area",
    "Seed_dry_mass",
    "Leaf_P_per_mass",
    "Stem_conduit_diameter",
    "Tree_height",
    "Root_depth",
    "Bark_thickness"
  ),
  testset = testset_all8,
  model_name = "BM")

model_performance_all8 <- result_all8_BM$performance
model_performance_all8

# ------------------------------------------------------------
# K-fold cross-validation: brownian motion
# ------------------------------------------------------------

PHY_kfoldcv_all8_BM_FTF <- run_kfold_phylopars_parallel(
  cv_datasets = kfold_cv_all8,
  tree = phylogeny,
  model = "BM",
  pheno_error = FALSE,
  phylo_correlated = TRUE,
  pheno_correlated = FALSE,
  usezscores = FALSE,
  EM_verbose = FALSE,
  optim_verbose = FALSE,
  n_cores = 4,
  checkpoint_dir = here::here("checkpoints/kfold_all8_BM_FTF"),
  resume = TRUE,
  retry_failed = FALSE,
  worker_stop_grace = 0.5
)

all8_objects <- c(
  "PHY_kfoldcv_all8_BM_FTF",
  "PHY_kfoldcv_all8_BM_TTF",
  "PHY_kfoldcv_all8_lambda_FTF",
  "PHY_kfoldcv_all8_lambda_TTF"
)

saved_all8 <- save_PHY_results(
  object_names = all8_objects,
  output_dir =
    here::here("Results/kfold_cv/all8_PHY")
)

kfold_model_performance_all8_BM <- calculate_kfold_performance(PHY_kfoldcv_all8_BM_FTF,
                                                               data_scale = "scaled")
kfold_model_performance_all8_BM

cv_plots_all8_BM <- plot_cv(kfold_model_performance_all8_BM)
cv_plots_all8_BM

# ------------------------------------------------------------
# K-fold cross-validation: brownian motion(fast-fit)
# ------------------------------------------------------------

fast_all8 <- estimate_kfold_phylocov_start_parallel(
  cv_datasets = kfold_cv_all8,
  tree = phylogeny,
  pheno_error = FALSE,
  phylo_correlated = FALSE,
  pheno_correlated = FALSE,
  usezscores = TRUE,
  n_cores = 4
)

full_all8 <- run_kfold_phylopars_parallel(
  cv_datasets = kfold_cv_all8,
  tree = phylogeny,
  pheno_error = FALSE,
  phylo_correlated = TRUE,
  phylocov_start = fast_all8$phylocov_start,
  pheno_correlated = FALSE,
  usezscores = TRUE,
  n_cores = 4
)

performance_all8 <- calculate_kfold_performance(
  phylopars_cv = full_all8,
  data_scale = "scaled"
)
performance_all8
