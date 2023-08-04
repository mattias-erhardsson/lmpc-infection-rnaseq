################################## Set seed for reproducibility
set.seed(1337)

################################## Load packages
lapply(
  c(
    "devtools", # If new packages of a specific version needs to be installed
    "renv", # For project management
    "BiocManager", # For project management
    "httpgd", # For better figures in interactive mode
    "plyr", #Data wrangling, part of tidyverse but not automatically loaded with it. Always load plyr before dply to avoid known issues # nolint: error. # nolint
    "tidyverse", # Data wrangling, processing and presentation. Does not seem to work with renv, so individual packages need to be loaded in this environment
    "vroom", # Faster data wrangling
    "readr", # Data wrangling
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
    "RCy3" # For cytoscape programmatic analysis
    ),
  library,
  character.only = TRUE
)

sessionInfo()

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
annotationTable <- read_tsv(file = "./R_intermediate_files/annotationTable.tsv",
col_types = c("cccc"))

referenceSet <- read_tsv(file = "./R_intermediate_files/referenceSet.tsv",
col_types = c("c"))

geneIDs <- read_tsv(file = "./R_intermediate_files/geneIDs.tsv",
col_types = c("c"))

########################################## GSEA with SetRank
## Create set collection object for SetRank
paste("Available cores:",
parallel::detectCores(all.tests = FALSE, logical = TRUE))
options(mc.cores = as.integer(parallel::detectCores(all.tests = FALSE, logical = TRUE))- 2) # Adapt to the number of cores you use. I have had problems running on 20 cores when I used all cores.
collection <- buildSetCollection(annotationTable,
                                 referenceSet = referenceSet,
                                 maxSetSize = 500 # Default is 500
)

## Use SetRank in ranked mode
## CAUTION! Might take several days to complete.
network <- setRankAnalysis(
  geneIDs = geneIDs, # Gene list ranked by adjusted p-value
  setCollection = collection, # SetRank collection from above
  use.ranks = TRUE, # Ranked mode
  setPCutoff = 0.01, # This is default of 0.01
  fdrCutoff = 0.05 # This is default of 0.05
)

## Export results
exportSingleResult(
  network = network,
  selectedGenes = geneIDs,
  collection = collection,
  networkName = "SetRank_Network",
  IDConverter = NULL,
  outputPath = "./R_output_files/Setrank_results"
)

print("Script 2 finished, continue by running script 3. As of writing script 3 would probably not work in a non-interactive mode due to interactions with cytoscape")