# ------------------------------------------------------------
# 6. Select species and traits for PhyloPars analysis
# 6. 选择物种名和性状用于 PhyloPars 分析
# ------------------------------------------------------------

high_correlation_traits_1 <- extract_traits(
  trait_data = trait_data,
  trait_names = c(
    "Leaf_N_per_mass",
    "Specific_leaf_area",
    "Leaf_thickness"
  )
)

high_correlation_traits_1_plot <- plot_trait_histograms(
  trait_data = high_correlation_traits_1,
  trait_names = c(
    "Leaf_N_per_mass",
    "Specific_leaf_area",
    "Leaf_thickness"
  )
)

high_correlation_traits_1_plot

trait_high_correlation_1 <- extract_traits(
  trait_data = trait_data,
  trait_names = c(
    "Leaf_N_per_mass",
    "Specific_leaf_area",
    "Leaf_thickness"
  ),
  log_transform = TRUE,
  scale_transform = TRUE
)

head(trait_high_correlation_1)


# ------------------------------------------------------------
# 7. Cross-validation
# 7. 交叉验证
# ------------------------------------------------------------

# Create the cross-validation dataset
cv_high_correlation_1 <- create_cross_validation(
  trait_data = trait_high_correlation_1,
  test_proportion = 0.25,
  seed = 123
)

trait_high_correlation_1_cv <-
  cv_high_correlation_1$trait_cv

testset_high_correlation_1 <-
  cv_high_correlation_1$testset


# ------------------------------------------------------------
# 8. K-fold cross-validation
# 8. K-fold 交叉验证
# ------------------------------------------------------------

kfold_cv_high_correlation_1 <- create_kfold_cv(
  trait_data = trait_high_correlation_1,
  k = 10,
  seed = 123
)

kfold_cv_high_correlation_1
