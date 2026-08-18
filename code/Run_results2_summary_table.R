# ============================================================
# 1. Load the function
# ============================================================

source(here::here("code/functions/results2_summary_table.R"))


# ============================================================
# 2. Check required objects
# ============================================================

stopifnot(
  exists("trait_data", inherits = TRUE),
  exists("phylogeny", inherits = TRUE)
)


# ============================================================
# 3. Define analysis sets
# ============================================================

results2_trait_groups <- list(
  All8 = c(
    "Wood_density",
    "Specific_leaf_area",
    "Seed_dry_mass",
    "Leaf_P_per_mass",
    "Stem_conduit_diameter",
    "Tree_height",
    "Root_depth",
    "Bark_thickness"
  ),
  `Cluster 1` = c(
    "Leaf_N_per_mass",
    "Specific_leaf_area",
    "Leaf_thickness"
  ),
  `Cluster 2` = c(
    "Stem_conduit_diameter",
    "Leaf_Vcmax_per_dry_mass",
    "Stomatal_conductance",
    "Leaf_area"
  ),
  `Cluster 3` = c(
    "Leaf_K_per_mass",
    "Leaf_P_per_mass"
  ),
  `Cluster 4` = c(
    "Leaf_density",
    "Wood_density"
  ),
  `Cluster 5` = c(
    "Bark_thickness",
    "Stem_diameter"
  ),
  `Cluster 6` = c(
    "Crown_height",
    "Crown_diameter",
    "Tree_height"
  )
)

results2_group_labels <- c(
  All8 = "Focal eight traits",
  `Cluster 1` = "Cluster 1",
  `Cluster 2` = "Cluster 2",
  `Cluster 3` = "Cluster 3",
  `Cluster 4` = "Cluster 4",
  `Cluster 5` = "Cluster 5",
  `Cluster 6` = "Cluster 6"
)


# ============================================================
# 4. Define input and output paths
# ============================================================

results2_r2_file <- paste0(
  here::here("Results/R2_plots/"),
  "R2_grouped_trait_range_data.csv"
)

results2_output_dir <- paste0(
  here::here("Results/R2_plots/"),
  "Results2_single_multi_table"
)


# ============================================================
# 5. Build and save all tables
# ============================================================

results2_tables <- build_results2_summary_table(
  trait_data = trait_data,
  tree = phylogeny,
  r2_data = results2_r2_file,
  trait_groups = results2_trait_groups,
  group_labels = results2_group_labels,
  species_column = "species",
  configurations = c(
    "BM_FTF",
    "BM_TTF",
    "lambda_FTF",
    "lambda_TTF"
  ),
  single_label = "Single trait",
  multi_label = "Multi-trait",
  output_dir = results2_output_dir,
  file_prefix = "Results2_single_multi",
  save_xlsx = TRUE
)


# ============================================================
# 6. Inspect the results
# ============================================================

View(results2_tables$main_table)
View(results2_tables$trait_detail)
View(results2_tables$pairwise_overlap)
View(results2_tables$delta_r2_detail)

print(results2_tables$saved_files)
