################################## Script start
print("Starting script 6")

################################## Set seed for reproducibility
set.seed(1337)

################################## Install packages
renv::restore()

renv::install("readxl@1.4.3", prompt = FALSE)

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
    "viridis", # Color palette
    "readxl" # Loading excel files
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

################################## Create plots results subdirectory if it doesn't exist
if (!dir.exists("./R_output_files/Figures")) {
  dir.create("./R_output_files/Figures")
  print("Creating directory ./R_output_files/Figures")
} else {
  print("./R_output_files/Figures exists")
}

################################## Create tables results subdirectory if it doesn't exist
if (!dir.exists("./R_output_files/Tables")) {
  dir.create("./R_output_files/Tables")
  print("Creating directory ./R_output_files/Tables")
} else {
  print("./R_output_files/Tables exists")
}

################################## Data import
# Potential glycosyltransferases for mucin O-glycans
Glycosyltransferase_activity_genes_to_check <- readxl::read_xlsx(path = "./R_input_files/Glycosyltransferase_activity_genes_to_check.xlsx")

# TPM dataframe
TPM_df <- read_tsv("./R_intermediate_files/TPM_df.tsv")

################################## Sort potential glycosyltransferases based on mean tpm of non-infected mice
Glycosyltransferase_activity_genes_to_check_top_mean_TPM <- Glycosyltransferase_activity_genes_to_check %>% 
  dplyr::filter(JB_recommendation_to_keep == 1) %>% 
  inner_join(TPM_df, by = c("GeneSymbol", "ensembldb_ID")) %>% 
  dplyr::filter(condition == "Non_Infected") %>% 
  dplyr::select(GeneSymbol, ensembldb_ID, User_ID, TPM) %>% 
  group_by(GeneSymbol) %>% 
  dplyr::mutate(Mean_Gene_TPM_Non_Infected = mean(TPM)) %>% 
  ungroup() %>% 
  arrange(desc(Mean_Gene_TPM_Non_Infected)) %>% 
  pivot_wider(names_from = User_ID, values_from = TPM)

writexl::write_xlsx(x = Glycosyltransferase_activity_genes_to_check_top_mean_TPM,
                    path = "./R_intermediate_files/Glycosyltransferase_activity_genes_to_check_top_mean_TPM.xlsx")

########################################## SessionInfo
sessionInfo()

print("Script 6 finished")
