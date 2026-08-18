# ============================================================
# Load the plotting function
# ============================================================

source(here::here("code/functions/RMSE_plot.R"))


# ============================================================
# Define the output directory
# ============================================================

output_dir <- here::here("Results/RMSE_plots")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Define all study traits
# ============================================================

all_study_traits <- c(
  "Wood_density",
  "Specific_leaf_area",
  "Seed_dry_mass",
  "Leaf_P_per_mass",
  "Stem_conduit_diameter",
  "Tree_height",
  "Root_depth",
  "Bark_thickness",
  "Leaf_N_per_mass",
  "Leaf_thickness",
  "Leaf_Vcmax_per_dry_mass",
  "Stomatal_conductance",
  "Leaf_area",
  "Leaf_K_per_mass",
  "Leaf_density",
  "Stem_diameter",
  "Crown_height",
  "Crown_diameter"
)


# ============================================================
# Define the three six-trait groups
# ============================================================

trait_groups <- list(
  B = all_study_traits[1:6],
  C = all_study_traits[7:12],
  D = all_study_traits[13:18]
)

trait_group_numbers <- c(
  B = "(I)",
  C = "(II)",
  D = "(III)"
)

trait_group_file_ids <- c(
  B = "I",
  C = "II",
  D = "III"
)


# ============================================================
# Figures B-D: six traits per figure
# ============================================================

for (panel_name in names(trait_groups)) {

  message(
    "Creating Figure ",
    panel_name,
    "..."
  )

  output_stub <- file.path(
    output_dir,
    paste0(
      "Family_trait_RMSE_",
      trait_group_file_ids[[panel_name]]
    )
  )

  common_arguments <- list(
    trait_names = trait_groups[[panel_name]],

    # Use all 18 traits to select the common species set
    species_selection_traits = all_study_traits,

    model = "lambda",
    configuration = "TTF",

    folds = c(
      1L,
      3L,
      5L,
      6L,
      8L,
      10L
    ),

    plot_title = paste0(
      "Family-level RMSE across selected traits ",
      trait_group_numbers[[panel_name]]
    ),

    width = 8.25,
    height = 10.5,

    radial_compression_power = 1.90,
    ring_width_fraction = 0.026,
    ring_gap_fraction = 0.0045,
    largest_family_labels = 20L,
    family_key_groups = 8L,
    show_complete_family_key = FALSE
  )

  # Save PNG
  do.call(
    plot_family_trait_rmse_phylogeny,
    c(
      common_arguments,
      list(
        output_file = paste0(
          output_stub,
          ".png"
        )
      )
    )
  )

  # Save PDF
  do.call(
    plot_family_trait_rmse_phylogeny,
    c(
      common_arguments,
      list(
        output_file = paste0(
          output_stub,
          ".pdf"
        )
      )
    )
  )

  message(
    "Figure ",
    panel_name,
    " completed."
  )
}


# ============================================================
# Figure A: complete family-name overview
# ============================================================

message(
  "Creating Figure A..."
)

family_overview_arguments <- list(
  trait_names = all_study_traits,
  species_selection_traits = all_study_traits,

  model = "lambda",
  configuration = "TTF",

  folds = c(
    1L,
    3L,
    5L,
    6L,
    8L,
    10L
  ),

  plot_title = NULL,

  width = 8.25,
  height = 11.25,

  radial_compression_power = 1.90,
  ring_width_fraction = 0.0165,
  ring_gap_fraction = 0.0022,

  largest_family_labels = 20L,
  family_key_groups = 8L,
  show_complete_family_key = TRUE
)


# Save Figure A as PNG
do.call(
  plot_family_trait_rmse_phylogeny,
  c(
    family_overview_arguments,
    list(
      output_file = file.path(
        output_dir,
        "Family_names_complete_overview_no_title.png"
      )
    )
  )
)


# Save Figure A as PDF
do.call(
  plot_family_trait_rmse_phylogeny,
  c(
    family_overview_arguments,
    list(
      output_file = file.path(
        output_dir,
        "Family_names_complete_overview_no_title.pdf"
      )
    )
  )
)

message(
  "Figure A completed."
)


# ============================================================
# Report the generated files
# ============================================================

generated_files <- c(
  file.path(output_dir, "Family_names_complete_overview_no_title.png"),
  file.path(output_dir, "Family_names_complete_overview_no_title.pdf"),
  unlist(lapply(c("I", "II", "III"), function(group_id) {
    file.path(
      output_dir,
      paste0("Family_trait_RMSE_", group_id, c(".png", ".pdf"))
    )
  }), use.names = FALSE)
)

if (!all(file.exists(generated_files))) {
  stop(
    "One or more expected output files were not created."
  )
}

message(
  "All four figures were created successfully."
)

print(
  generated_files
)
