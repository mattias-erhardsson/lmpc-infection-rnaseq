################################## Script start
print("Starting script 1")

################################## Set seed for reproducibility
set.seed(1337)

################################## Install packages
# Renv bootstrats, but on UPPMAX it fails to install packages correctly with renv::restore()
# Therefore, packages need to be installed manually
# Firstly, pbdZMQ@0.3-9 needs to be installed with devtools since even manual renv install fails on UPPMAX
renv::install("devtools@2.4.5", prompt = FALSE)
library("devtools")
devtools::install_version("pbdZMQ", version = "0.3-9")

# Then restore which installs most of the needed packages
renv::restore()

# Then manual installation of important packages that renv seems to struggle with in UPPMAX
# CRAN packages first
renv::install("plyr@1.8.8", prompt = FALSE)
renv::install("tidyverse@2.0.0", prompt = FALSE)
renv::install("tidyr@1.3.0", prompt = FALSE)
renv::install("readr@2.1.4", prompt = FALSE)
renv::install("purrr@1.0.2", prompt = FALSE)
renv::install("forcats@1.0.0", prompt = FALSE)
renv::install("ggplot2@3.4.2", prompt = FALSE)
renv::install("dplyr@1.1.2", prompt = FALSE)
renv::install("tibble@3.2.1", prompt = FALSE)
renv::install("stringr@1.5.0", prompt = FALSE)
renv::install("vroom@1.6.3", prompt = FALSE)
renv::install("svglite@2.1.1", prompt = FALSE)
renv::install("writexl@1.4.2", prompt = FALSE)
renv::install("styler@1.10.1", prompt = FALSE)
renv::install("viridis@0.6.4", prompt = FALSE)
renv::install("BiocManager@1.30.22", prompt = FALSE)

# Then Bioconductor packages
renv::install("bioc::DESeq2", prompt = FALSE)
renv::install("bioc::IHW", prompt = FALSE)
renv::install("bioc::tximport", prompt = FALSE)
renv::install("bioc::tximportData", prompt = FALSE)
renv::install("bioc::umap", prompt = FALSE)
renv::install("bioc::EnsDb.Mmusculus.v79", prompt = FALSE)
renv::install("bioc::ensembldb", prompt = FALSE)
renv::install("bioc::SetRank", prompt = FALSE)
renv::install("bioc::biomaRt", prompt = FALSE)
renv::install("bioc::org.Mm.eg.db", prompt = FALSE)
renv::install("bioc::reactome.db", prompt = FALSE)
renv::install("bioc::GO.db", prompt = FALSE)
renv::install("bioc::KEGGREST", prompt = FALSE)
renv::install("bioc::STRINGdb", prompt = FALSE)
renv::install("bioc::igraph", prompt = FALSE)

# RCy3 is a special case
# It is a bioconductor package, but there is a newer version on GitHub that I need
# It implements a new function to change label positions
# Important addition was made in commit 2023-08-08
# Using latest available commit 2023-08-12 below
renv::install("cytoscape/RCy3@3015129d026346c1307d06e1eb9d48be6d675318", prompt = FALSE)

################################## Load packages again
## Appears tidyverse does not play nicely with renv, have to call packages I need manually
lapply(
  c(
    "renv", # For project management
    "BiocManager", # For project management
    "plyr", # Data wrangling, part of tidyverse but not automatically loaded with it. Always load plyr before dply to avoid known issues # nolint: error. # nolint
    "ggplot2", # Tidyverse. Data wrangling, processing and presentation.
    "dplyr", # Tidyverse. Data wrangling, processing and presentation.
    "tidyr", # Tidyverse. Data wrangling, processing and presentation.
    "readr", # Tidyverse. Data wrangling, processing and presentation.
    "purrr", # Tidyverse. Data wrangling, processing and presentation.
    "tibble", # Tidyverse. Data wrangling, processing and presentation.
    "stringr", # Tidyverse. Data wrangling, processing and presentation.
    "forcats", # Tidyverse. Data wrangling, processing and presentation.
    "vroom", # Faster data wrangling
    "lubridate", # Working with dates, part of tidyverse but not automatically loaded with it
    "svglite", # To make svg files with ggsave
    "writexl", # Writing excel files
    "DESeq2", # Differential gene expression
    "IHW", # Better power for adjusting p-values of differential gene expression
    "tximport", # Importing RNAseq pipeline data
    "tximportData", # Importing RNAseq pipeline data
    "umap", # For umap, Uniform Manifold Approximation and Projection
    "EnsDb.Mmusculus.v79", # ensdb package for mouse
    "ensembldb", # For getting gene lengths necessary for TPM calculation
    "SetRank", # For the SetRank GSEA
    "parallel", # For parallelisation
    "biomaRt", # For annotation and GO gene sets
    "org.Mm.eg.db", # For GO term annotation, might be used instead of biomaRt
    "reactome.db", # For annotationdbi of reactome
    "GO.db", # For GO term annotation, might be used instead of biomaRt
    "KEGGREST", # For KEGG
    "styler", # R studio addin for interactively adhere to the tidyverse style guide
    "RCy3", # For cytoscape programmatic analysis
    "STRINGdb", # For STRING database annotation
    "igraph", # For RCy3/cytoscape
    "viridis" # Color palette
  ),
  library,
  character.only = TRUE
)

##################################################### Script finished
sessionInfo()

print("Script 1 finished")
