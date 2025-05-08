################################## Script start
print("Starting script 5")

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
Histology_Raw_Data <- read_xlsx(path = "./R_input_files/Histology_Data.xlsx",
                                col_types = c("text",
                                              "text",
                                             "numeric",
                                             "numeric",
                                             "numeric",
                                             "numeric",
                                             "numeric",
                                             "numeric",
                                             "numeric",
                                             "numeric",
                                             "numeric",
                                             "text",
                                             "text",
                                             "logical",
                                             "text"))

################################## Plot histology data
# In case any immediate processing of the raw data becomes necessary
Histology_Data <- Histology_Raw_Data %>% 
  mutate(Cohort = factor(Cohort, levels = c("M", "G1", "G2")))

#Statistical test of histological assays with fdr adjustment

HAI_pvalue <- wilcox.test(x = Histology_Data %>%
              dplyr::filter(condition == "Infected") %>% 
              dplyr::select(HAI_Total) %>% 
              deframe(),
            y = Histology_Data %>%
              dplyr::filter(condition == "Non_Infected") %>% 
              dplyr::select(HAI_Total) %>% 
              deframe(),
            exact = FALSE)[["p.value"]]
print(HAI_pvalue)

Bacteria_pvalue <- wilcox.test(x = Histology_Data %>%
                                 dplyr::filter(condition == "Infected") %>% 
                                 dplyr::select(Bacteria) %>% 
                                 deframe(),
                               y = Histology_Data %>%
                                 dplyr::filter(condition == "Non_Infected") %>% 
                                 dplyr::select(Bacteria) %>% 
                                 deframe(),
                               exact = FALSE)[["p.value"]]
print(Bacteria_pvalue)

Histological_Assay_Tests <- tibble(assay = c("HAI", 
                                             "Bacteria"),
                                   pvalue = c(HAI_pvalue, 
                                              Bacteria_pvalue)) %>% 
  mutate(padj = p.adjust(pvalue, method = "BH"))

#HAI total
#Manually count the overlaps
counted_data_HAI <- Histology_Data %>%
  group_by(Cohort, HAI_Total, condition) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(n = as.factor(n))  # Convert counts to factor for manual sizing

HAI_Total_Plot <- ggplot(counted_data_HAI,
       aes(x = Cohort,
           y = HAI_Total,
           color = condition,
           shape = Cohort)) +
  scale_colour_manual(values = c("Non_Infected" = "#4DC36B", "Infected" = "#440C55")) +
  geom_point(aes(size = n)) +
  scale_y_continuous(breaks = 0:24) +
  geom_text(data = Histology_Data %>%
      dplyr::filter(User_ID %in% c("H1001", "H1008", "H1002", "H1009")),
    aes(label = Manuscript_ID),
    nudge_x = 0.35) +
  scale_size_manual(values = c("1" = 3, "2" = 5)) +
  theme_bw() +
  labs(title = str_c("adjusted p-value = ", as.character(Histological_Assay_Tests %>% 
                                            dplyr::filter(assay == "HAI") %>% 
                                            dplyr::select(padj) %>% 
                                            deframe())))
print(HAI_Total_Plot)
ggsave(filename = "./R_output_files/Figures/HAI_Total_Plot.eps",
       plot = HAI_Total_Plot,
       width = 750*5, #750*2
       height = 1120, #1120*2
       units = "px")
ggsave(filename = "./R_output_files/Figures/HAI_Total_Plot.pdf",
       plot = HAI_Total_Plot,
       width = 750*3, #750*2
       height = 1120, #1120*2
       units = "px")

#Bacteria
#Manually count the overlaps
counted_data_Bacteria <- Histology_Data %>%
  group_by(Cohort, Bacteria, condition) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(n = as.factor(n))  # Convert counts to factor for manual sizing
Bacteria_Plot <- ggplot(counted_data_Bacteria,
       aes(x = Cohort,
           y = Bacteria,
           color = condition,
           shape = Cohort)) +
  scale_colour_manual(values = c("Non_Infected" = "#4DC36B", "Infected" = "#440C55")) +
  geom_point(aes(size = n)) +
  scale_y_continuous(breaks = seq(from = 0,
                                  to = 100,
                                  by = 10)) +
  geom_text(data = Histology_Data %>%
              dplyr::filter(User_ID %in% c("H1001", "H1008", "H1002", "H1009")),
            aes(label = Manuscript_ID),
            nudge_x = 0.35) +
  scale_size_manual(values = c("1" = 3, "2" = 5, "3" = 7, "4" = 9, "5" = 11)) +
  theme_bw() +
  labs(title = str_c("adjusted p-value = ", as.character(Histological_Assay_Tests %>% 
                                            dplyr::filter(assay == "Bacteria") %>% 
                                            dplyr::select(padj) %>% 
                                            deframe())))
print(Bacteria_Plot)
ggsave(filename = "./R_output_files/Figures/Bacteria_Plot.eps",
       plot = Bacteria_Plot,
       width = 750*5, #750*2
       height = 1120, #1120*2
       units = "px")
ggsave(filename = "./R_output_files/Figures/Bacteria_Plot.pdf",
       plot = Bacteria_Plot,
       width = 750*3, #750*2
       height = 1120, #1120*2
       units = "px")


########################################## SessionInfo
sessionInfo()

print("Script 5 finished")

