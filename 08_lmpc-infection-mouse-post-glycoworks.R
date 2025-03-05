################################## Script start
print("Starting script 8")

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
    "tidyr", # Tidyverse. Data wrangling, processing and presentation.
    "dplyr", # Tidyverse. Data wrangling, processing and presentation.
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
    "readxl", # Loading excel files
    "reticulate", # Using python in R
    "pheatmap" # Heatmaps
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
df_canonicalized_mouse_all <- read_tsv("./Python_input_files/df_canonicalized_mouse_all.tsv")
df_glycan_canonicalized_all_data <- readxl::read_xlsx(path = "./R_output_files/tables/df_glycan_canonicalized_all_data.xlsx")
quantified_terminal_motifs <- readxl::read_xlsx(path = "./Python_output_files/Tables/quantified_terminal_motifs.xlsx")
glycans_terminal_motifs_annotation <- readxl::read_xlsx(path = "./Python_output_files/Tables/glycans_terminal_motifs_annotation.xlsx")

################################## Show what the top 10 most abundant glycans were
df_glycan_canonicalized_all_data %>% 
  dplyr::select(Glycan_ID, Canonicalized_Structure, Glycan_Mean_Relative_Abundance) %>% 
  dplyr::distinct() %>% 
  dplyr::slice_max(Glycan_Mean_Relative_Abundance,
                   n = 10)

################################## Calculate top 10 vs rest abundances
# Top 10 glycans summed abundance
df_glycan_canonicalized_all_data %>% 
  dplyr::select(Glycan_ID, Canonicalized_Structure, Glycan_Mean_Relative_Abundance) %>% 
  dplyr::distinct() %>% 
  dplyr::slice_max(Glycan_Mean_Relative_Abundance,
                   n = 10) %>% 
  dplyr::select(Glycan_ID) %>% 
  dplyr::inner_join(df_glycan_canonicalized_all_data, by = "Glycan_ID") %>% 
  dplyr::select(Glycan_Mean_Relative_Abundance) %>% 
  dplyr::distinct() %>% 
  dplyr::summarise(Summed_Abundance_Top_10_Glycans = sum(Glycan_Mean_Relative_Abundance))

# Rest (=90) bottom glycans summed abundance
df_glycan_canonicalized_all_data %>% 
  dplyr::select(Glycan_ID, Canonicalized_Structure, Glycan_Mean_Relative_Abundance) %>% 
  dplyr::distinct() %>% 
  dplyr::slice_min(Glycan_Mean_Relative_Abundance,
                   n = 90) %>% 
  dplyr::select(Glycan_ID) %>% 
  dplyr::inner_join(df_glycan_canonicalized_all_data, by = "Glycan_ID") %>% 
  dplyr::select(Glycan_Mean_Relative_Abundance) %>% 
  dplyr::distinct() %>% 
  dplyr::summarise(Summed_Abundance_Bottom_90_Glycans = sum(Glycan_Mean_Relative_Abundance))

# Max abundance of individual bottom 90 glycans
df_glycan_canonicalized_all_data %>% 
  dplyr::select(Glycan_ID, Canonicalized_Structure, Glycan_Mean_Relative_Abundance) %>% 
  dplyr::distinct() %>% 
  dplyr::slice_min(Glycan_Mean_Relative_Abundance,
                   n = 90) %>% 
  dplyr::select(Glycan_ID) %>% 
  dplyr::inner_join(df_glycan_canonicalized_all_data, by = "Glycan_ID") %>% 
  dplyr::select(Glycan_Mean_Relative_Abundance) %>% 
  dplyr::distinct() %>% 
  dplyr::summarise(Max_Abundance_Bottom_90_Glycans = max(Glycan_Mean_Relative_Abundance))

################################## Show what the top 10 most abundant terminal glycan motifs with 1, 2 or 3 monosackarides were
quantified_terminal_motifs %>% 
  pivot_longer(cols = starts_with("G"),
               names_to = "Sample_ID",
               values_to = "Motif_Relative_Abundance") %>% 
  group_by(Terminal_Motif) %>% 
  mutate("Mean_Motif_Relative_Abundance" = mean(Motif_Relative_Abundance)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Terminal_Motif, 
                Mean_Motif_Relative_Abundance) %>% 
  dplyr::distinct() %>% 
  dplyr::slice_max(Mean_Motif_Relative_Abundance,
                   n = 10) %>% 
  dplyr::arrange(desc(Mean_Motif_Relative_Abundance))

################################## PCA on motif level
df_r_glycomics <- df_glycan_canonicalized_all_data %>% 
  dplyr::group_by(Glycan_Size, Sample_ID) %>% 
  dplyr::mutate("Glycan_Size_Relative_Abundance_Percentage" = sum(Glycan_Relative_Abundance_Percentage)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Glycan_ID, 
                Structure, 
                Canonicalized_Structure, 
                Glycan_Size, 
                Hex, 
                HexNAc, 
                Fucose, 
                Neu5Ac, 
                Neu5Gc, 
                Sulf, 
                Sample_ID, 
                Treatment_Group, 
                Cohort, 
                Glycan_Size_Relative_Abundance_Percentage, 
                Glycan_Relative_Abundance_Percentage) %>% 
  distinct()

## Glycan PCA
# Annotations
annotations <- df_r_glycomics %>%
  dplyr::select(Sample_ID, Treatment_Group, Cohort) %>%
  dplyr::distinct()

# Run PCA
df_r_terminal_motif_wide <- quantified_terminal_motifs %>%
  tibble::column_to_rownames("Terminal_Motif")

matrix_terminal_motif <- as.matrix(df_r_terminal_motif_wide)

# Combine with annotation
pca_analysis <- prcomp(t(matrix_terminal_motif),
                       center = TRUE,
                       scale = TRUE)

pca_df <- as.data.frame(pca_analysis$x) %>% 
  rownames_to_column("Sample_ID") %>% 
  full_join(annotations, by = "Sample_ID")

# Define custom colors for annotations
ann_colors <- list(
  Treatment_Group = c(ShamInfected = "#4DC36B", HpyloriInfected = "#440C55"),
  Cohort = c(G1 = "#E85311", G2 = "#15B8E9")
)

pca_plot <- ggplot(pca_df, aes(
  x = PC1,
  y = PC2,
  color = Treatment_Group,
  shape = Cohort
)) +
  scale_colour_manual(values = c("ShamInfected" = "#4DC36B", "HpyloriInfected" = "#440C55")) +
  geom_label(aes(label = Sample_ID),
             size = 2,
             position = position_nudge(1.1)) +
  geom_point() +
  labs(
    x = paste0(
      "PC1 (captures ",
      (summary(pca_analysis)$importance[2, 1])*100,
      "% of variance)"
    ),
    y = paste0(
      "PC2 (captures ",
      (summary(pca_analysis)$importance[2, 2])*100,
      "% of variance)"
    ),
    title = paste0(
      "PC1 and PC2 together captures ",
      (summary(pca_analysis)$importance[3, 2])*100,
      "% of variance in terminal glycan motifs"
    )
  ) +
  theme_bw()
print(pca_plot)

ggsave(filename = "./R_output_files/Figures/pca_plot.eps",
  plot = pca_plot,
  width = 1484,
  height = 700,
  units = "px"
)

################################## Calculate percentage of terminal motifs that were H antigens
# Relative abundance of glycans with at least 1 H antigen
Glycans_With_H_Antigen <- df_canonicalized_mouse_all %>% 
  dplyr::select(Canonicalized_Structure, Glycan_ID, Glycan_Mean_Relative_Abundance) %>% 
  dplyr::full_join(glycans_terminal_motifs_annotation,
                   by = "Canonicalized_Structure") %>%
  pivot_longer(cols = starts_with("Terminal"),
               names_to = "Motif",
               values_to = "Count") %>% 
  dplyr::mutate("H_Antigen_Motif" = str_detect(Motif, "Terminal_Fuc\\(a...\\)Gal\\(b...\\)GlcNAc")) %>% 
  dplyr::mutate("H_Antigen_Motif_And_Detected_In_Glycan" = H_Antigen_Motif == TRUE & Count > 0) %>% 
  dplyr::select(Glycan_ID, Canonicalized_Structure, Glycan_Mean_Relative_Abundance, H_Antigen_Motif_And_Detected_In_Glycan) %>% 
  dplyr::filter(H_Antigen_Motif_And_Detected_In_Glycan == TRUE) %>% 
  dplyr::distinct()

Glycans_With_H_Antigen %>% 
  summarise("Relative_Abundance_Of_Glycans_With_H_Antigen" = sum(Glycan_Mean_Relative_Abundance))

# Sanity check by comparing relative abundance of all other glycans
df_canonicalized_mouse_all %>% 
  dplyr::anti_join(Glycans_With_H_Antigen, by = "Glycan_ID") %>% 
  summarise("Relative_Abundance_Of_Glycans_Without_H_Antigen" = sum(Glycan_Mean_Relative_Abundance))

########################################## SessionInfo
sessionInfo()

print("Script 8 finished")

