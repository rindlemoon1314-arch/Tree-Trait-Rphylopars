# Tree-Trait-Rphylopars

This repository contains the R code and summary outputs used to evaluate phylogenetic trait imputation for tree functional traits using `Rphylopars`.

## Data source

Trait observations were derived from the TRY Plant Trait Database and the tree trait dataset compiled by Maynard et al. (2022), *Global relationships in tree functional traits*, Nature Communications.

The raw TRY data are subject to TRY data-use restrictions and are therefore not redistributed in this repository. Users who wish to reproduce the analyses should request access to the relevant TRY data through the TRY database and prepare the input files locally.

Main data references:

- TRY Plant Trait Database: https://www.try-db.org/
- Maynard et al. (2022): https://doi.org/10.1038/s41467-022-30888-2
- Associated archive: https://doi.org/10.5281/zenodo.6564051

## Code

All analysis scripts are stored in:

```text
code/
```

The code uses project-relative paths through `here::here()`. Therefore, scripts should be run from the repository root, preferably by opening the RStudio project file:

```text
Tree-Trait-Rphylopars.Rproj
```

Before running the scripts, users should check that local file paths match the expected repository structure. 

Some scripts were used for exploratory or intermediate analyses. The final model comparisons reported in the dissertation use:

```r
phylo_correlated = TRUE
pheno_correlated = FALSE
```

Users should check these settings before re-running model-fitting scripts.

## Analysis workflow

The code/ folder contains both main analysis scripts and supporting plotting or utility functions. The main workflow is organised as follows:

- Preparation.R prepares the trait dataset and phylogenetic tree for model fitting. 

- Scripts in code/all8/ fit models for the eight selected traits.

- Scripts in code/high_correlation_traits fit the corresponding univariate and multivariate models for the six functional trait clusters.

- Plotting scripts in code/functions/ generate the figures and summary tables used in the dissertation. 

Users wishing to reproduce the main dissertation results should begin with Preparation.R, then run the model-fitting scripts for the eight selected traits and the six functional trait clusters before running the summary and plotting scripts.

## V.PhyloMaker2

The phylogenetic tree used in the analysis was prepared outside this repository using `V.PhyloMaker2`. Because `V.PhyloMaker2` is not available from CRAN, it should be installed from GitHub, for example:

```r
install.packages("remotes")
remotes::install_github("jinyizju/V.PhyloMaker2")
```

After installation, users should prepare the required species list and generate or prune the phylogeny locally. 

## Computational notes

Model fitting can be computationally intensive, especially for multivariate models and cross-validation. Parallel settings such as the number of cores should be adjusted according to the user’s machine.

Large intermediate files, checkpoints and raw restricted data are not included in the repository. These files should be regenerated locally when needed.

