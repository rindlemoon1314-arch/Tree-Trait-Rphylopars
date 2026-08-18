library(Rphylopars)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(iterators)
library(snow)
library(foreach)
library(doSNOW)
library(parallel)
source(here::here("code/functions/extract_trait.R"))
source(here::here("code/functions/create_cv.R"))
source(here::here("code/functions/create_kfold_cv.R"))
source(here::here("code/functions/run_kfold_phylopars_parallel.R"))
source(here::here("code/functions/run_single_traits.R"))
source(here::here("code/functions/save_PHY_results.R"))
source(here::here("code/functions/plot_fold_status.R"))
source(here::here("code/functions/plot_grouped_r2_range.R"))
source(here::here("code/functions/single_multi_r2.R"))
source(here::here("code/functions/TTF_FTF_r2.R"))
source(here::here("code/functions/lambda_BM_r2.R"))
source(here::here("code/functions/calculate_phy_metrics.R"))
source(here::here("code/functions/plot_trait_histograms.R"))
source(here::here("code/functions/plot_interval.R"))

phylogeny <- read.tree(here::here("Tree Trait data/no_monocots_tree.nwk"))
tree_data <- read.csv(here::here("Tree Trait data/TRY_trait_data_cleaned.csv"))
plot(phylogeny, show.tip.label = FALSE)

# ------------------------------------------------------------
# 1. Standardise species names and trait names
# 1. 标准化物种名和性状名
# ------------------------------------------------------------
tree_data_adjusted <- tree_data %>%
  mutate(accepted_bin = gsub(" ", "_", accepted_bin),
         trait = gsub(" ", "_", trait))
head(tree_data_adjusted)

# ------------------------------------------------------------
# 2. Match species between trait data and phylogenetic tree
# 2. 匹配性状数据和系统发育树中的物种
# ------------------------------------------------------------
same_tree <- intersect(unique(tree_data_adjusted$accepted_bin), phylogeny$tip.label)
only_in_tree_data <- setdiff(unique(tree_data_adjusted$accepted_bin), phylogeny$tip.label)
length(only_in_tree_data)

# ------------------------------------------------------------
# 3. Keep only species that exist in both datasets
# 3. 只保留同时存在于两个数据集中的物种
# ------------------------------------------------------------
tree_data_same <- tree_data_adjusted %>%
  filter(accepted_bin %in% same_tree)

# Check
same_tree_test <- intersect(unique(tree_data_same$accepted_bin), phylogeny$tip.label)
only_in_tree_data_test <- setdiff(unique(tree_data_same$accepted_bin), phylogeny$tip.label)
length(only_in_tree_data_test) # This should ideally be 0

# ------------------------------------------------------------
# 4. Convert the trait data format
# 4. 转换性状数据格式
# ------------------------------------------------------------

# Make sure that there are not duplicated trait values
duplicated_trait_values <- tree_data_same %>%
  group_by(accepted_bin, LAT, LON, trait) %>%
  filter(n() > 1) %>%
  ungroup()

# Convert
trait_data <- tree_data_same %>%
  pivot_wider(id_cols = c(accepted_bin, LAT, LON),
              names_from = trait,
              values_from = value)
head(trait_data)

# Check
unique_species_in_same <- tree_data_same %>%
  distinct(accepted_bin, LAT, LON) %>%
  nrow()
unique_species_in_same == length(trait_data$accepted_bin) # This should ideally be TRUE

same_location_count <- tree_data_same %>%
  count(accepted_bin, LAT, LON, name = "n") %>%
  filter(n > 1)
sum(same_location_count$n)

diff_location_count <- tree_data_same %>%
  count(accepted_bin, LAT, LON, name = "n") %>%
  filter(n == 1)
nrow(diff_location_count)

sum(same_location_count$n) + nrow(diff_location_count) == nrow(tree_data_same) # This should ideally be TRUE

# ------------------------------------------------------------
# 5. Rename species column and remove latitude/longitude columns
# 5. 重命名物种列，并删除经纬度列
# ------------------------------------------------------------
trait_data <- trait_data %>%
  rename(species = accepted_bin) %>%
  select(-LAT, -LON)
glimpse(trait_data)

# ------------------------------------------------------------
# 6. Initial assessment of the trait dataset
# 6. 初步分析
# ------------------------------------------------------------

# Create summary for all trait columns
trait_assessment_table <- trait_data %>%
  summarise(
    across(
      -species,
      list(
        data_type = ~ class(.)[1],
        observed_count = ~ sum(!is.na(.)),
        observed_proportion = ~ paste0(round(mean(!is.na(.))*100, 3), "%"),
        min = ~ if (is.numeric(.)) min(., na.rm = TRUE) else NA_real_,
        max = ~ if (is.numeric(.)) max(., na.rm = TRUE) else NA_real_,
        mean = ~ if (is.numeric(.)) mean(., na.rm = TRUE) else NA_real_,
        median = ~ if (is.numeric(.)) median(., na.rm = TRUE) else NA_real_
      )
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("trait", ".value"),
    names_pattern = "(.+)_(data_type|observed_count|observed_proportion|min|max|mean|median)"
  ) %>%
  arrange(desc(observed_count))

print(trait_assessment_table, n = Inf)





