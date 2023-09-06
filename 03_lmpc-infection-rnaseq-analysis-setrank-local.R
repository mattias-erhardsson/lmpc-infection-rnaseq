################################## Script start
print("Starting script 3")

################################## Set seed for reproducibility
set.seed(1337)

################################## Install packages
renv::restore()

################################## Load packages
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

################################## Check existance of required directories
if (!dir.exists("P26010")) {
  stop("You need the P26010 root folder in the same directory as the script")
} else {
  print("./P26010 exists")
}

if (!dir.exists("R_input_files")) {
  stop("You need the R_input_files root folder in the same directory as the script")
} else {
  print("./R_input_files exists")
}

if (!dir.exists("R_intermediate_files")) {
  stop("You need the R_intermediate_files root folder in the same directory as the script")
} else {
  print("./R_intermediate_files exists")
}

if (!dir.exists("R_output_files")) {
  stop("You need the R_output_files root folder in the same directory as the script")
} else {
  print("./R_output_files exists")
}

################################## Create SetRank results subdirectory if it doesn't exist
if (!dir.exists("./R_output_files/Setrank_results")) {
  dir.create("./R_output_files/Setrank_results")
  print("Creating directory ./R_output_files/Setrank_results")
} else {
  print("./R_output_files/Setrank_results exists")
}

################################## Data import
annotationTable <- read_tsv(
  file = "./R_intermediate_files/annotationTable.tsv",
  col_types = c("cccc")
)

referenceSet <- read_tsv(
  file = "./R_intermediate_files/referenceSet.tsv",
  col_types = c("c")
) %>%
  deframe()

geneIDs <- as.vector(read_tsv(
  file = "./R_intermediate_files/geneIDs.tsv",
  col_types = c("c")
)) %>%
  deframe()

################################## Restrict to GO:BP, KEGG and Reactome
# I cant make SetRank work properly with big data sets on my HPC cluster UPPMAX
# Therefore I need to run the code locally
# But with all 5 gene sets, this is too heavy for my computer
# Therefore, I will restrict the analysis to what I think are the most biologically relevant, namely:
# GO:BP, KEGG and Reactome
# Also, the bug on UPPMAX replaces all _ with . so I wonder if thats the bug.
annotationTable <- annotationTable %>% 
  dplyr::filter(dbName %in% c("KEGG",
                              "Reactome",
                              "GO_biological_process"))

referenceSet <- base::intersect(referenceSet, annotationTable$geneID)

geneIDs <- base::intersect(geneIDs, annotationTable$geneID)

########################################## GSEA with SetRank
## Create set collection object for SetRank
paste(
  "Available cores:",
  parallel::detectCores(all.tests = FALSE, logical = TRUE)
)
options(mc.cores = as.integer(parallel::detectCores(all.tests = FALSE, logical = TRUE)) - 1) # Running on all cores can make it slower according to package documentation
collection <- buildSetCollection(annotationTable,
  referenceSet = referenceSet,
  maxSetSize = 500 # Default is 500
)

# Save the collection object to make it easier to re-run analysis
base::save(collection,
           file = "./R_intermediate_files/collection.RData")

## Use SetRank in ranked mode
## CAUTION! Might take several days to complete.
network <- setRankAnalysis(
  geneIDs = geneIDs, # Gene list ranked by adjusted p-value
  setCollection = collection, # SetRank collection from above
  use.ranks = TRUE, # Ranked mode
  setPCutoff = 0.01, # This is default of 0.01
  fdrCutoff = 0.05 # This is default of 0.05
)

# Save the network object to make it easier to re-run analysis
base::save(network,
           file = "./R_intermediate_files/network.RData")

## Export results
exportSingleResult(
  network = network,
  selectedGenes = geneIDs,
  collection = collection,
  networkName = "SetRank_Network",
  IDConverter = NULL,
  outputPath = "./R_output_files/Setrank_results"
)

########################################## SessionInfo
sessionInfo()

print("Script 3 finished")
