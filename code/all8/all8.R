# ------------------------------------------------------------
# 6. Select species and trait for PhyloPars analysis
# 6. 选择物种名和性状用于 PhyloPars 分析
# ------------------------------------------------------------
all8 <- extract_traits(
  trait_data = trait_data,
  trait_names = c(
    "Wood_density",
    "Specific_leaf_area",
    "Seed_dry_mass",
    "Leaf_P_per_mass",
    "Stem_conduit_diameter",
    "Tree_height",
    "Root_depth",
    "Bark_thickness"
  )
)

all8_plot <- plot_trait_histograms(
  trait_data = all8,
  trait_names = c(
    "Wood_density",
    "Specific_leaf_area",
    "Seed_dry_mass",
    "Leaf_P_per_mass",
    "Stem_conduit_diameter",
    "Tree_height",
    "Root_depth",
    "Bark_thickness"
  )
)

all8_plot

trait_all8 <- extract_traits(
  trait_data = trait_data,
  trait_names = c(
    "Wood_density",
    "Specific_leaf_area",
    "Seed_dry_mass",
    "Leaf_P_per_mass",
    "Stem_conduit_diameter",
    "Tree_height",
    "Root_depth",
    "Bark_thickness"
  ),
  log_transform = TRUE,
  scale_transform = TRUE
)

head(trait_all8)

# ------------------------------------------------------------
# 7.Cross-validation
# 7.交叉验证
# ------------------------------------------------------------

# Create the cross-validation dataset
cv_all8 <- create_cross_validation(
  trait_data = trait_all8,
  test_proportion = 0.25,
  seed = 123
)

trait_all8_cv <- cv_all8$trait_cv
testset_all8 <- cv_all8$testset

# ------------------------------------------------------------
# 8.K-fold cross-validation
# 8.K-fold 交叉验证
# ------------------------------------------------------------

kfold_cv_all8 <- create_kfold_cv(
  trait_data = trait_all8,
  k = 10,
  seed = 123
)

kfold_cv_all8
