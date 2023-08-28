################################## important notice about environment
# This script is made to be run on UPPMAX
# renv has been deactivated
# modules R/4.2.1 and R_packages/4.2.1 have been loaded prior to runnin this script

################################## Script start
print("Starting script 4")

################################## Set seed for reproducibility
set.seed(1337)

################################## Load packages
lapply(
  c("tidyverse",
    "SetRank"),
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

########################################## GSEA with SetRank
## Create set collection object for SetRank
paste(
  "Available cores:",
  parallel::detectCores(all.tests = FALSE, logical = TRUE)
)
options(mc.cores = as.integer(parallel::detectCores(all.tests = FALSE, logical = TRUE)) - 2) # Running on all cores can make it slower according to package documentation
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

########################################## SessionInfo
sessionInfo()

print("Script 4 finished")
