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

################################## PCA on motif level with k-means clustering
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

## K means clustering
# Define the range of clusters to test
k_values <- 1:10

# Calculate total within-cluster sum of square (WSS) for each k
wss <- map_dbl(k_values, function(k){
  kmeans(pca_analysis$x[, 1:2], centers = k, nstart = 25)$tot.withinss
})

# Create a dataframe for plotting the elbow curve
elbow_df <- data.frame(k = k_values, wss = wss)

# Plot the elbow curve
elbow_plot <- ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_point() +
  geom_line() +
  labs(
    x = "Number of Clusters (k)",
    y = "Total Within-Cluster Sum of Squares (WSS)",
    title = "Elbow Method for Optimal k"
  ) +
  scale_x_continuous(breaks = k_values) +
  theme_minimal()

print(elbow_plot)

# Select the optimal k (you can choose based on the elbow plot)
optimal_k <- 5  # Adjust based on the elbow plot

# Perform k-means clustering with the selected k
kmeans_result <- kmeans(pca_analysis$x[, 1:2], centers = optimal_k, nstart = 25)

# Add cluster assignments to the PCA dataframe
pca_df <- pca_df %>%
  mutate(Cluster = as.factor(kmeans_result$cluster))

# Plot PCA with cluster assignments
pca_cluster_plot <- ggplot(pca_df, aes(
  x = PC1,
  y = PC2,
  color = Treatment_Group,
  shape = Cluster
)) +
  geom_point(size = 3) +
  geom_label(aes(label = Sample_ID#, fill = Cohort
                 ),
             size = 2,
             position = position_nudge(1.3)) +
  scale_colour_manual(values = c("ShamInfected" = "#4DC36B", "HpyloriInfected" = "#440C55")) +
  #scale_fill_manual(values = c(G1 = "#E85311", G2 = "#15B8E9")) +
  labs(
    x = paste0(
      "PC1 (captures ",
      round((summary(pca_analysis)$importance[2, 1]) * 100, 1),
      "% of variance)"
    ),
    y = paste0(
      "PC2 (captures ",
      round((summary(pca_analysis)$importance[2, 2]) * 100, 1),
      "% of variance)"
    ),
    title = paste0(
      "K-means Clustering on PC1 and PC2 with k = ", optimal_k
    )
  ) +
  theme_minimal()

print(pca_cluster_plot)


ggsave(
  filename = "./R_output_files/Figures/PCA.eps",
  plot = PCA,
  width = 1484,
  height = 700,
  units = "px"
)

########################################## SessionInfo
sessionInfo()

print("Script 8 finished")

