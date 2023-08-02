################################## Set seed for reproducibility
set.seed(1337)

################################## Load packages
#lapply(
#  c(
#    "renv", # For project management
#    "plyr", # Data wrangling, part of tidyverse but not automatically loaded with it. Always load plyr before dply to avoid known issues
#    "tidyverse", # Data wrangling, processing and presentation
#    "lubridate", # Working with dates, part of tidyverse but not automatically loaded with it
#    "svglite", # To make svg files with ggsave
#    "writexl", # Writing excel files
#    "DESeq2", # Differential gene expression
#    "IHW", # Better power for adjusting p-values of differential gene expression
#    "tximport", # Importing RNAseq pipeline data
#    "tximportData", # Importing RNAseq pipeline data
#    "umap", # For umap, Uniform Manifold Approximation and Projection
#    "EnsDb.Mmusculus.v79", # ensdb package for mouse
#    "ensembldb", # For getting gene lengths necessary for TPM calculation
#    "SetRank", # For the SetRank GSEA
#    "biomaRt", # For annotation and GO gene sets
#    "org.Mm.eg.db", # For GO term annotation, might be used instead of biomaRt
#    "reactome.db", # For annotationdbi of reactome
#    "GO.db", # For GO term annotation, might be used instead of biomaRt
#    "KEGGREST", # For KEGG
#    "styler", # R studio addin for interactively adhere to the tidyverse style guide
#    "grateful" # For citing packages used
#  ),
#  library,
#  character.only = TRUE
#)
#renv::init()

install.packages("grateful", repos = "https://ftp.acc.umu.se/mirror/CRAN/")

################################## Load packages
lapply(
  c(
    "renv", # For project management
    "plyr", # Data wrangling, part of tidyverse but not automatically loaded with it. Always load plyr before dply to avoid known issues
    "tidyverse", # Data wrangling, processing and presentation
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
    "biomaRt", # For annotation and GO gene sets
    "org.Mm.eg.db", # For GO term annotation, might be used instead of biomaRt
    "reactome.db", # For annotationdbi of reactome
    "GO.db", # For GO term annotation, might be used instead of biomaRt
    "KEGGREST", # For KEGG
    "styler", # R studio addin for interactively adhere to the tidyverse style guide
    "grateful" # For citing packages used
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

################################## Create R output directory if it doesn't exist
if (!dir.exists("R_output_files")) {
  dir.create("R_output_files")
  print("Creating directory ./R_output_files")
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

################################## Create grateful results subdirectory if it doesn't exist
if (!dir.exists("./R_output_files/Grateful")) {
  dir.create("./R_output_files/Grateful")
  print("Creating directory ./R_output_files/Grateful")
} else {
  print("./R_output_files/Grateful exists")
}


################################## Data import
# All file reading assumes you are in a root directory above P26010
# Import sample info file from 00-Reports, annotate sample list and convert potentially problematic characters in column names
# Change path to where in your directory the file S.Linden_22_01_sample_info.txt is
# Excluding samples from treatment group as they are outside the scope of this article
Sample_List <- read_tsv(file = "./P26010/00-Reports/S.Linden_22_01_sample_info.txt",
                        col_types = c("ccdd")) %>%
  dplyr::rename(names = "NGI ID") %>%
  dplyr::rename(User_ID = "User ID") %>%
  dplyr::rename(GreaterThan_Q30 = "≥Q30") %>%
  dplyr::filter(!User_ID %in% c(
    "H9_11",
    "H9_12",
    "H9_13",
    "H9_14",
    "H9_15"
  )) %>%
  mutate(Infected = c(
    "Non_Infected", "Non_Infected", "Non_Infected", "Non_Infected", "Non_Infected",
    "Infected", "Infected", "Infected", "Infected"
  )) %>%
  mutate(condition = Infected)

# Importing data about the LMPC and the RNA purification
# Excluding samples from treatment group since they are outside the scope of this article
# Code assumes this data file is in the same directory as the script
LMPC_RNA_Data <- read_tsv(file = "./R_input_files/LMPC_RNA_Data.txt",
                          col_types = c("cifdddccidddd")) %>%
  mutate("Date_harvested" = lubridate::as_date(Date_harvested)) %>% #I couldn't get the date to parse correctly with readr so I change from character to date after reading
  mutate("Date_sectioned" = lubridate::as_date(Date_sectioned)) %>%
  dplyr::filter(!User_ID %in% c(
    "H9_11",
    "H9_12",
    "H9_13",
    "H9_14",
    "H9_15"
  ))

files <- file.path(
  "./P26010/01-RNA-Results/star_salmon",
  Sample_List$names,
  "quant.sf"
)
print(files)
file.exists(files)

# tx2gene from the rnaseq pipeline
tx2gene <- read_tsv(
  file = "./P26010/01-RNA-Results/genome/salmon_tx2gene.tsv",
  col_names = c(
    "TXNAME",
    "GENEID",
    "GeneSymbol"
  ),
  col_types = (c("ccc"))
)

# import star salmon results with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
glimpse(txi$counts)
colnames(txi$counts) <- Sample_List$names
glimpse(txi$counts)

################################## Investigate statistics about LMPC and RNA purification data
## Table of statistics about sectioning time, LMPC area, RIN and V200 values
Summary_LMPC_RNA_data <- LMPC_RNA_Data %>%
  dplyr::select(
    Days_between_harvest_and_sectioning,
    Extracted_area_with_10_micrometer_thick_sections,
    Tapestation_RIN_Value,
    SciLifeLab_RIN_Value,
    DV200,
    Tapestation_RNA_concentration_ng_ul,
    SciLifeLab_RNA_concentration_ng_ul,
    SciLifeLab_RNA_amount_ng
  ) %>%
  rename_all(~ str_replace_all(., "_", "-")) %>%
  dplyr::summarise(across(everything(), list(mean = mean, median = median, min = min, max = max))) %>%
  pivot_longer(cols = everything(), names_to = "summary_stat", values_to = "value") %>%
  separate(summary_stat, into = c("Variable", "summary_stat"), sep = "_") %>%
  pivot_wider(names_from = summary_stat, values_from = value)

# Write table of summary statistics
write_xlsx(
  x = Summary_LMPC_RNA_data,
  path = "./R_output_files/Tables/Summary_LMPC_RNA_data.xlsx"
)
Summary_LMPC_RNA_data

## Any correlation between extracted area and purified RNA?
# Fit a linear model
lm_model <- lm(SciLifeLab_RNA_amount_ng ~ Extracted_area_with_10_micrometer_thick_sections, data = LMPC_RNA_Data)

# Get the model summary
summary_lm <- summary(lm_model)

# Extract the R-squared value
r_squared <- summary_lm$r.squared

# Print the R-squared value
print(r_squared)

# Visualise correlation
LMPC_RNA_Plot <- ggplot(
  LMPC_RNA_Data,
  aes(
    x = Extracted_area_with_10_micrometer_thick_sections,
    y = SciLifeLab_RNA_amount_ng
  )
) +
  geom_point() +
  geom_hline(yintercept = 10) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  expand_limits(y = 0, x = 0) +
  scale_y_continuous(breaks = seq(0, max(LMPC_RNA_Data$SciLifeLab_RNA_amount_ng), by = 10)) +
  scale_x_continuous(breaks = seq(0, max(LMPC_RNA_Data$Extracted_area_with_10_micrometer_thick_sections), by = 1000000)) +
  labs(
    title = paste("R-squared =", round(r_squared, 3)),
    x = "LMPC catapulted area (square micrometers)",
    y = "Purified RNA (nanograms)"
  ) +
  theme_minimal()
LMPC_RNA_Plot
ggsave(
  filename = "./R_output_files/Figures/LMPC_RNA_Plot.svg",
  plot = LMPC_RNA_Plot,
  width = 2250,
  height = 2625,
  units = "px"
)

################################## Create DESeq object
### There are ENSEMBL IDs which are matched to the same gene symbol
### Our raw data have to be imported with Tximport the best way
### The best way I could think of to handle the duplicated symbols is a multi-step process:
### 1 Create DESeq object from tximport.
### 2 Extract counts from DESeq object
### 3 Merge duplicated GeneSymbols
### 4 Re-create DESeq object with the merged counts

## 1 Create DESeq object from tximport.
# Sample annotation data frame for deseq2
sampleTable <- data.frame(
  condition = factor(
    x = Sample_List$condition,
    levels = c(
      "Non_Infected",
      "Infected"
    )
  ),
  row.names = Sample_List$names
)

# Create table matching ENSEMBL gene id to gene symbol for future annotation
Gene_Symbols <- tx2gene %>%
  dplyr::select(-TXNAME) %>%
  distinct() %>%
  arrange(GENEID)

# Prepare deseq2 object from tximport, this is the "original counts and offset" from the vignette
txi_dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

## Extract counts from DESeq object
Counts_All <- data.frame(BiocGenerics::counts(txi_dds)) %>%
  rownames_to_column("GENEID") %>%
  inner_join(Gene_Symbols, # With full join it will be a complete match with tx2gene, but there will be NA values
             by = "GENEID"
  ) %>%
  replace(is.na(.), 0) # Replace NA values with 0

# Merge counts/gene ids for genes which have been annotated with the same symbol
# Every gene should only be matched to one gene symbol, otherwise might become a problem in downstream analysis
# First check if there are any gene gene symbols mapped to multiple gene IDs
Counts_All %>%
  dplyr::select(GENEID, GeneSymbol) %>%
  distinct() %>%
  plyr::count("GeneSymbol") %>% # this is not a typo, it uses the package plyr and not dplyr
  dplyr::filter(freq > 1)
# Yes. mmu-mir-669a-4 even has 9.
# As sanity check since this shouldn't be possible, are there any gene ids mapped to multiple symbols?
Counts_All %>%
  dplyr::select(GENEID, GeneSymbol) %>%
  distinct() %>%
  plyr::count("GENEID") %>% # sic!
  dplyr::filter(freq > 1)
# No

## Merge duplicated GeneSymbols. Approach is to simply take the first of GENEIDs mapped to the symbol.
Counts_All_NoDuplicatedSymbol <- Counts_All %>%
  pivot_longer(
    cols = starts_with("P26010"),
    names_to = "names",
    values_to = "Counts"
  ) %>%
  group_by(GeneSymbol, names) %>%
  dplyr::summarise("Counts" = sum(Counts)) %>%
  #ddply(.(GeneSymbol, names), summarise, "Counts" = sum(Counts)) %>%
  pivot_wider(
    values_from = Counts,
    names_from = names
  ) %>%
  inner_join(Gene_Symbols,
             by = "GeneSymbol",
             multiple = "first" # this takes the first of the multiple matches. Might not work depnding on package verison
  ) %>%
  distinct(GeneSymbol, .keep_all = TRUE) #If the mutiple = "first" doesn't work with your R version this should work
# Double check duplicates were removed
Counts_All_NoDuplicatedSymbol %>%
  dplyr::select(GENEID, GeneSymbol) %>%
  distinct() %>%
  plyr::count("GeneSymbol") %>% # sic!
  dplyr::filter(freq > 1)
# They were removed

# Convert to tidy format (easier to work with)
Tidy_Counts_All_NoDuplicatedSymbol <- Counts_All_NoDuplicatedSymbol %>%
  pivot_longer(
    cols = starts_with("P26010"),
    names_to = "names",
    values_to = "Counts"
  )

# Add variable showing which genes have a count of at least 1 in at least 1 sample ("TRUE"), and which ones we did not detect ("FALSE")
Tidy_Counts_DetectedGenes_NoDuplicatedSymbol <- Tidy_Counts_All_NoDuplicatedSymbol %>%
  ddply(.(GENEID), summarise, "Summed_Counts" = sum(Counts)) %>%
  mutate("Gene_Detected" = ifelse(Summed_Counts > 0, TRUE, FALSE)) %>%
  dplyr::select(-Summed_Counts) %>%
  inner_join(Tidy_Counts_All_NoDuplicatedSymbol,
             by = "GENEID",
             multiple = "all"
  ) %>%
  inner_join(Sample_List %>% dplyr::select(names, User_ID),
             by = "names"
  ) %>%
  dplyr::select(-names)

## 4 Re-create DESeq object with the merged counts
# Now matrix to put back into deseq2
# Since it is easier to understand gene symbols and we now have unique symbols, we will use that primarily
Counts_Matrix_Unique_Symbols <- as.matrix(Tidy_Counts_DetectedGenes_NoDuplicatedSymbol %>%
                                            pivot_wider(
                                              values_from = Counts,
                                              names_from = User_ID
                                            ) %>%
                                            dplyr::select(-GENEID, -Gene_Detected) %>%
                                            column_to_rownames("GeneSymbol"))

mtx_dds <- DESeqDataSetFromMatrix(
  countData = Counts_Matrix_Unique_Symbols,
  colData = Sample_List %>%
    dplyr::select(User_ID, condition) %>%
    mutate(condition = factor(condition,
                              levels = c(
                                "Non_Infected",
                                "Infected"
                              )
    )) %>%
    column_to_rownames("User_ID"),
  design = ~condition
)

## Process with DESeq2
dds <- DESeq(mtx_dds)

################################## PCA and Uniform Manifold Approximation and Projection (umap) of gene expression
## Using principal components as input for umap

## PCA
## Principal Component Analysis(PCA)
# PCA PC1 and PC2
rld <- DESeq2::rlog(dds, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)
pca <- prcomp(t(rld_mat))
df <- cbind(Sample_List, pca$x)

PCA <- ggplot(df, aes(
  x = PC1,
  y = PC2,
  color = condition
)) +
  geom_label(aes(label = User_ID)) +
  labs(
    x = paste0(
      "PC1 (captures ",
      summary(pca)$importance[2, 1],
      " of variance)"
    ),
    y = paste0(
      "PC2 (captures ",
      summary(pca)$importance[2, 2],
      " of variance)"
    ),
    title = paste0(
      "PC1 and PC2 together captures ",
      summary(pca)$importance[3, 2],
      " of variance in gene expresssion"
    )
  )
print(PCA)

ggsave(
  filename = "./R_output_files/Figures/PCA.svg",
  plot = PCA,
  width = 2250,
  height = 2625,
  units = "px"
)

## Visualise how much each PC explains variance in order to select PC for UMAP
summary_pca <- summary(pca)
PC_Variance <- ggplot(
  enframe(summary_pca$importance[2, ],
          name = "PC",
          value = "Proportion_of_variance_captured"
  ) %>%
    mutate(PC = str_replace(PC, "(^[:alpha:]{2})(\\d)$", "\\10\\2")),
  aes(
    x = PC,
    y = Proportion_of_variance_captured
  )
) +
  geom_point()
print(PC_Variance)

ggsave(
  filename = "./R_output_files/Figures/PC_Variance.svg",
  plot = PC_Variance,
  width = 2250,
  height = 2625,
  units = "px"
)

## First 5 PCs seems appropriate to me.
UMAP_Input <- df[, 7:11]

## Set seed for umap, if I understand correctly this is by default random and not the same seed as set before.
custom.config <- umap.defaults
custom.config$random_state <- 1337

## Create umap object from wide counts matrix

UMAP_Output <- umap(UMAP_Input,
                    n_neighbors = 8
) # Highest value possible with 9 samples. Default is 15.

UMAP_mtx <- UMAP_Output[["layout"]]

colnames(UMAP_mtx) <- c("UMAP_1", "UMAP_2")

UMAP_df <- as.data.frame(UMAP_mtx) %>%
  rownames_to_column("User_ID") %>%
  inner_join(Sample_List,
             by = "User_ID"
  )

## Visualize umap
UMAP_Plot <- ggplot(
  data = UMAP_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = condition
  )
) +
  # geom_point() +
  geom_label(aes(label = User_ID)) +
  labs(
    title = "UMAP dimensionality reduction on PCA with PC1 through PC5"
  )
print(UMAP_Plot)

ggsave(
  filename = "./R_output_files/Figures/UMAP_Plot.svg",
  plot = UMAP_Plot,
  width = 2250,
  height = 2625,
  units = "px"
)

################################## Differential Gene Expression (DGE)
## Extract DESeq2 results
res_df <- DESeq2::results(dds,
                          alpha = 0.05,
                          contrast = c("condition", "Infected", "Non_Infected")
)
summary(res_df)
## Independent Hypothesis Weighting http://bioconductor.org/packages/release/bioc/html/IHW.html
## Data frame from deseq results

deRes_res_df <- as.data.frame(res_df)

# IHW adjustment of pvalues
ihwRes_res_df <- ihw(pvalue ~ baseMean,
                     data = deRes_res_df,
                     alpha = 0.05
)

# Number of significant genes
paste(
  "The number of genes with adjusted p-value < 0.05 in Infected_vs_Non_Infected were",
  rejections(ihwRes_res_df),
  "with IHW and",
  sum(deRes_res_df$padj <= 0.05, na.rm = TRUE),
  "with BH"
)

# Dataframe with IHW adjusted p values, fold changes, counts and rlog
Combined_Results_DF <- as.data.frame(ihwRes_res_df,
                                     row.names = rownames(deRes_res_df)
) %>%
  rownames_to_column("GeneSymbol") %>%
  mutate(Comparison = "Infected_vs_Non_Infected") %>%
  inner_join(
    deRes_res_df %>%
      rownames_to_column("GeneSymbol"),
    by = c("GeneSymbol", "pvalue")
  ) %>%
  dplyr::select(GeneSymbol, adj_pvalue, pvalue, log2FoldChange, Comparison) %>%
  inner_join(Tidy_Counts_DetectedGenes_NoDuplicatedSymbol,
             by = "GeneSymbol",
             relationship = "many-to-many"
  ) %>%
  inner_join(
    as.data.frame(assay(rlog(dds, blind = TRUE))) %>%
      rownames_to_column("GeneSymbol") %>%
      pivot_longer(
        cols = starts_with("H9"),
        values_to = "Rlog_Counts",
        names_to = "User_ID"
      ),
    by = c("GeneSymbol", "User_ID")
  ) %>%
  inner_join(
    Sample_List %>%
      dplyr::select(User_ID, condition),
    by = "User_ID"
  )

# Differential gene expression volcano plot
Volcano_Plot_DF <- Combined_Results_DF %>%
  dplyr::filter(Gene_Detected == TRUE) %>%
  dplyr::select(GeneSymbol, adj_pvalue, log2FoldChange, Comparison) %>%
  distinct() %>%
  mutate(adj_pvalue = ifelse(is.na(adj_pvalue), 1, adj_pvalue)) %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange)) %>%
  mutate(log2FoldChange = ifelse(log2FoldChange < -8, -8, log2FoldChange)) %>% # A few non-significant genes have extreme fold changes, this is for better visualisation
  mutate(log2FoldChange = ifelse(log2FoldChange > 8, 8, log2FoldChange)) %>% # A few non-significant genes have extreme fold changes, this is for better visualisation
  mutate(Regulation = factor(case_when(
    log2FoldChange > 0 & adj_pvalue < 0.05 ~ "up",
    log2FoldChange < -0 & adj_pvalue < 0.05 ~ "down",
    TRUE ~ "neither"
  )))

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "neither" = "grey")
sizes <- c("up" = 3, "down" = 3, "neither" = 1)
alphas <- c("up" = 1, "down" = 1, "neither" = 0.5)

Volcano_Plot <- ggplot(
  data = Volcano_Plot_DF,
  aes(
    x = log2FoldChange,
    y = -log10(adj_pvalue),
    size = Regulation,
    colour = Regulation,
    alpha = Regulation
  )
) +
  geom_point() +
  scale_colour_manual(values = cols) +
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  geom_label(
    data = Volcano_Plot_DF %>%
      dplyr::filter(adj_pvalue < 0.05) %>%
      dplyr::filter(log2FoldChange > 2 | log2FoldChange < -2),
    aes(label = GeneSymbol),
    label.size = 1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = c(-2, 2),
    linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = c(seq(-8, 8, 1)),
    limits = c(-8, 8)
  ) +
  scale_y_continuous(
    breaks = c(seq(0, 12, 1), -log10(0.05)),
    limits = c(0, 12),
    labels = as.character(c(seq(0, 12, 1), "adjusted P=0.05"))
  ) +
  labs(
    title = "Gene expression changes in infected vs non-infected samples",
    x = "log2(fold change)",
    y = "-log10(adjusted P-value)"
  ) +
  theme_bw() + # Select theme with a white background
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
print(Volcano_Plot)

ggsave(
  filename = "./R_output_files/Figures/Volcano_Plot_Differential_Gene_Expression_Infection.svg",
  plot = Volcano_Plot,
  width = 2250,
  height = 2625,
  units = "px"
)

## Saving data for all genes
genes_for_annotation <- Combined_Results_DF %>%
  dplyr::select(GENEID, GeneSymbol, Gene_Detected) %>%
  distinct() %>%
  replace_na(list(log2FoldChange = 0, adj_pvalue = 1, pvalue = 1)) # NA values might disrupt downstream analysis

all_detected_genes <- Combined_Results_DF %>%
  dplyr::filter(Gene_Detected == TRUE) %>%
  dplyr::select(GENEID, GeneSymbol, log2FoldChange, adj_pvalue, pvalue, Comparison) %>%
  distinct() %>%
  replace_na(list(log2FoldChange = 0, adj_pvalue = 1, pvalue = 1)) # NA values might disrupt downstream analysis

sig_genes <- all_detected_genes %>%
  dplyr::filter(adj_pvalue < 0.05)
################################## Export
write_tsv(all_detected_genes, "./R_output_files/Tables/All_Detected_Genes.tsv")

write_tsv(
  sig_genes,
  "./R_output_files/Tables/Significant_Genes.tsv"
)

write_xlsx(
  x = sig_genes,
  path = "./R_output_files/Tables/Significant_Genes.xlsx"
)

# top 10 lowest padj
sig_genes %>%
  slice_min(
    order_by = adj_pvalue,
    n = 10
  )

# Significant genes with log2 fold change > 2
sig_genes %>%
  dplyr::filter(log2FoldChange > 2)

# Significant genes with log2 fold change < -2
sig_genes %>%
  dplyr::filter(log2FoldChange < -2)

# Significant genes with absolute log 2 fold change > 4
top_infection_genes_fc4 <- sig_genes %>%
  dplyr::filter(log2FoldChange > 4 | log2FoldChange < -4) %>%
  arrange(log2FoldChange) %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue)

# Exporting significant genes with absolute log 2 fold change
write_xlsx(top_infection_genes_fc4, "./R_output_files/Tables/top_infection_genes_fc4.xlsx")


################################## Calculate TPM (transcripts per kilobase million) for better visualization of within-sample gene expression
ensdb <- EnsDb.Mmusculus.v79

## List all supported keytypes.
keytypes(ensdb)
## List all supported columns for the select and mapIds methods.
columns(ensdb)
## List /real/ database column names.
listColumns(ensdb)
## Retrieve all keys corresponding to keytype
ENTREZID <- keys(ensdb, keytype = "ENTREZID")
length(ENTREZID)
head(ENTREZID)
GENEID <- keys(ensdb, keytype = "GENEID")
length(GENEID)
head(GENEID)
GENENAME <- keys(ensdb, keytype = "GENENAME")
length(GENENAME)
head(GENENAME)
SYMBOL <- keys(ensdb, keytype = "SYMBOL")
length(SYMBOL)
head(SYMBOL)
TXID <- keys(ensdb, keytype = "TXID")
length(TXID)
head(TXID)

## Genename and symbol look identical, are they?
length(SYMBOL) == length(GENENAME)
SYMBOL[!SYMBOL == GENENAME]
## They are identical

## Lets get a mapping table between ensdb and our data

ensdb_gene_counts <- enframe(
  x = mapIds(ensdb, keys = SYMBOL, column = "GENEID", keytype = "SYMBOL"),
  name = "GeneSymbol",
  value = "ensembldb_ID"
) %>%
  inner_join(Tidy_Counts_DetectedGenes_NoDuplicatedSymbol,
             by = "GeneSymbol"
  ) %>%
  inner_join(
    enframe(
      x = mapIds(ensdb, keys = SYMBOL, column = "TXID", keytype = "SYMBOL"),
      name = "GeneSymbol",
      value = "ensembldb_TXID"
    ),
    by = "GeneSymbol"
  )

## Now get transcript lengths and calculate RPK (reads per kilobase)
ensdb_gene_counts_transcripts <- enframe(
  lengthOf(
    x = ensdb,
    of = "tx"
  ),
  name = "ensembldb_TXID",
  value = "Transcript_Length"
) %>%
  inner_join(ensdb_gene_counts,
             by = "ensembldb_TXID"
  ) %>%
  mutate(Transcript_Kbp = Transcript_Length / 1000) %>%
  mutate(RPK = Counts / Transcript_Kbp)

## Calculate scaling factor for TPM (transcripts per million), calculate TPM, and add more annotation
TPM_df <- ensdb_gene_counts_transcripts %>%
  group_by(User_ID) %>%
  dplyr::summarise(TPM_Scaling_Factor = sum(RPK) / 1e6) %>%
  inner_join(ensdb_gene_counts_transcripts,
             by = "User_ID"
  ) %>%
  mutate(TPM = RPK / TPM_Scaling_Factor) %>%
  inner_join(Sample_List,
             by = "User_ID"
  )

########################################## Gene markers for cell types to visualize what cells might be present
## PanglaoDB marker genes download
PanglaoDB_Marker_Genes <- read_tsv(file = "./R_input_files/PanglaoDB_markers_27_Mar_2020.tsv",
                                   col_types = c("ccccdcciccdddd")) %>%
  dplyr::rename("GeneSymbol" = `official gene symbol`) %>%
  dplyr::rename("Cell_Type" = `cell type`) %>%
  dplyr::mutate(GeneSymbol = str_to_title(GeneSymbol))

## Organ filter for cell types that could feasibly exist in the tissue
Organ_Filter <- c(
  "GI tract",
  "Immune system"
)

## Manually adding genes based on Bockerstett et al 2020 https://gut.bmj.com/content/69/6/1027.long
Bockerstett_Gene_List <- read_tsv(file = "./R_input_files/Bockerstett et all 2020 Gut supplemental material marker gene list.tsv",
                                  col_types = c("cic")) %>%
  dplyr::select(GeneSymbol, Cell_Type) %>%
  mutate("organ" = "GI tract") %>%
  distinct()

## Filter for genes presented as marker genes in several cell types in the cell types selected
Gene_Filter <- PanglaoDB_Marker_Genes %>%
  dplyr::filter(organ %in% Organ_Filter) %>%
  full_join(Bockerstett_Gene_List,
            by = c("GeneSymbol", "Cell_Type", "organ")
  ) %>%
  dplyr::count(GeneSymbol) %>%
  dplyr::filter(n != 1) %>%
  dplyr::select(GeneSymbol) %>%
  deframe()

## Cell type filter for cell types that could feasibly exist in the tissue
Cell_Types_Filter <- c(
  "B cells",
  "Basophils",
  "Dendritic cells",
  "Endothelial cells",
  "Enteroendocrine cells",
  "Eosinophils",
  "Foveolar cells",
  "Gastric chief cells",
  "Macrophages",
  "Mast cells",
  "Metaplastic cells",
  "Monocytes",
  "Mucous neck cells",
  "Natural killer T cells",
  "Enterochromaffin cells",
  # "Neuroendocrine cells", #These marker genes overlap with other genes so this cell type wouldn't show up after downstream filters
  "Neutrophils",
  "Parietal cells",
  "Plasma cells",
  "Plasmacytoid dendritic cells",
  "Proliferating Tff2+ cells", # A version of metaplastic cells
  "T cells"
)

## Filter marker genes which show up in multiple cell types in our analysis in order to increase specificity
Gene_Filter <- PanglaoDB_Marker_Genes %>%
  full_join(Bockerstett_Gene_List,
            by = c("GeneSymbol", "Cell_Type", "organ")
  ) %>%
  dplyr::filter(GeneSymbol != "Muc6") %>% # Manual removal from gene filter. Will instead be removed from metaplastic cells to avoid overlap
  dplyr::filter(Cell_Type %in% Cell_Types_Filter) %>%
  dplyr::count(GeneSymbol) %>%
  dplyr::filter(n != 1) %>%
  dplyr::select(GeneSymbol) %>%
  deframe()

## DF of marker genes
TPM_ggplot_df <- PanglaoDB_Marker_Genes %>%
  dplyr::filter(`canonical marker` == 1) %>% # One way to get as specific gene markers as possible
  dplyr::filter(!specificity_mouse > 0) %>% # Another way to make gene markers more specific
  dplyr::filter(species == "Mm Hs" | species == "Mm") %>% # We only want marker genes relevant for mice
  full_join(Bockerstett_Gene_List,
            by = c("GeneSymbol", "Cell_Type", "organ")
  ) %>%
  dplyr::filter(!(GeneSymbol == "Muc6" & Cell_Type == "Metaplastic cells")) %>% # Filtering out Muc6 for metaplastic cells to avoid overlap with mucous neck cells
  dplyr::filter(!GeneSymbol %in% Gene_Filter) %>%
  dplyr::filter(Cell_Type %in% Cell_Types_Filter) %>%
  inner_join(TPM_df,
             by = "GeneSymbol",
             relationship = "many-to-many"
  ) %>%
  dplyr::select(Cell_Type, GeneSymbol, TPM, names, condition, organ, `canonical marker`, specificity_mouse) %>%
  distinct()

## How many marker genes are there for each cell type?
TPM_ggplot_df %>%
  dplyr::select(Cell_Type, GeneSymbol) %>%
  distinct() %>%
  plyr::count("Cell_Type")

## Save table
write_xlsx(
  x = TPM_ggplot_df %>%
    dplyr::select(Cell_Type, GeneSymbol) %>%
    distinct() %>%
    plyr::count("Cell_Type"),
  path = "./R_output_files/Tables/Cell_Type_Markers_Summary.xlsx"
)

## Plot marker genes
Marker_Genes_Plot <- ggplot(
  TPM_ggplot_df,
  aes(
    x = Cell_Type,
    y = TPM,
    colour = condition
  )
) +
  geom_point() +
  theme_bw()
Marker_Genes_Plot

## Genes that stick out are TFF1, MUC5AC, ATP4A, ATP4B, and Tspan8.
## Add manual annotation for these genes

TPM_ggplot_df_annotated <- TPM_ggplot_df %>%
  mutate(High_TPM_Genes = if_else(GeneSymbol %in% c("Muc5ac", "Tff1", "Atp4a", "Atp4b", "Tspan8"), GeneSymbol, "None")) # These are the top 5 genes, and they have an expression greater than 1000TPM

Marker_Genes_Plot <- ggplot(
  TPM_ggplot_df_annotated,
  aes(
    x = Cell_Type,
    y = TPM,
    shape = condition,
    color = High_TPM_Genes
  )
) +
  geom_jitter(width = 0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Marker_Genes_Plot

## Save the plot
ggsave(
  filename = "./R_output_files/Figures/Marker_Genes_Plot.svg",
  plot = Marker_Genes_Plot,
  width = 2250,
  height = 2625,
  units = "px"
)

## Are any of the marker genes among genes with significant differential gene expression?
Significant_Marker_Genes <- sig_genes %>%
  # mutate(GeneSymbol = str_to_upper(GeneSymbol)) %>%
  dplyr::filter(GeneSymbol %in% TPM_ggplot_df_annotated$GeneSymbol)
Significant_Marker_Genes

## Find which cell types these belong to
Significant_Cell_Types <- TPM_ggplot_df_annotated %>%
  dplyr::filter(GeneSymbol %in% Significant_Marker_Genes$GeneSymbol) %>%
  dplyr::select(Cell_Type, GeneSymbol) %>%
  distinct()
Significant_Cell_Types

## How many marker genes do these cell types have in total?
TPM_ggplot_df_annotated %>%
  dplyr::filter(Cell_Type %in% Significant_Cell_Types$Cell_Type) %>%
  dplyr::select(Cell_Type, GeneSymbol) %>%
  distinct() %>%
  plyr::count("Cell_Type")

##################################################### Next comes several sections for annotating the genes with gene sets. It is very messy, to whoever is re-running this script I'm so sorry about this
##################################################### Annotate GO gene sets with bioMart, might be better than annotationdbi
## Annotate with GO.db and Biomart and compare
## First Biomart
## The rnaseq pipeline used ENSEMBL release 81
listEnsemblArchives()
## BiomaRt doesn't have release 81. I will use the latest realease which of writing is 109. It should have better annotation with for example GO terms.
listEnsembl(version = 109)
ensembl <- useEnsembl(
  biomart = "genes",
  version = 109
)
listDatasets(ensembl)
searchDatasets(mart = ensembl, pattern = "mus")
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  host = "www.ensembl.org",
  mirror = "www"
)
head(listFilters(ensembl))
Attributes <- listAttributes(ensembl)
Legend_bioMart <- getBM(
  mart = ensembl,
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name"
  )
)

# Not all ensembl IDs and gene symbols will match perfectly.
# Creating an annotation dataframe that is as complete as possible by taking into account mis-match
# First, the genes where both identifiers match
Matching_ENSEMBL_Symbol <- genes_for_annotation %>%
  mutate("ensembl_gene_id" = GENEID) %>%
  mutate("external_gene_name" = GeneSymbol) %>%
  inner_join(Legend_bioMart,
             by = c("ensembl_gene_id", "external_gene_name")
  ) %>%
  mutate(biomaRt_compatability = "ENSEMBL_and_symbol_matching") %>%
  relocate(c(
    GENEID,
    ensembl_gene_id,
    GeneSymbol,
    external_gene_name,
    biomaRt_compatability
  )) %>%
  distinct()

# Where only ENSEMBL match
Matching_ENSEMBL_Only <- genes_for_annotation %>%
  dplyr::filter(!GENEID %in% Matching_ENSEMBL_Symbol$GENEID) %>%
  mutate("ensembl_gene_id" = GENEID) %>%
  inner_join(Legend_bioMart,
             by = c("ensembl_gene_id")
  ) %>%
  mutate(biomaRt_compatability = "ENSEMBL_only_matching") %>%
  relocate(c(
    GENEID,
    ensembl_gene_id,
    GeneSymbol,
    external_gene_name,
    biomaRt_compatability
  )) %>%
  distinct()

# Where only symbol match
Matching_Symbol_Only <- genes_for_annotation %>%
  dplyr::filter(!GENEID %in% Matching_ENSEMBL_Symbol$GENEID) %>%
  dplyr::filter(!GENEID %in% Matching_ENSEMBL_Only$GENEID) %>%
  mutate("external_gene_name" = GeneSymbol) %>%
  inner_join(Legend_bioMart,
             by = c("external_gene_name"),
             multiple = "first"
  ) %>% # How I handled multiple matches
  mutate(biomaRt_compatability = "Symbol_only_matching") %>%
  relocate(c(
    GENEID,
    ensembl_gene_id,
    GeneSymbol,
    external_gene_name,
    biomaRt_compatability
  )) %>%
  distinct()

# Where neither ENSEMBL nor symbol matched
No_Matching <- genes_for_annotation %>%
  dplyr::filter(!GENEID %in% Matching_ENSEMBL_Symbol$GENEID) %>%
  dplyr::filter(!GENEID %in% Matching_ENSEMBL_Only$GENEID) %>%
  dplyr::filter(!GENEID %in% Matching_Symbol_Only$GENEID) %>%
  mutate("ensembl_gene_id" = GENEID) %>%
  mutate("external_gene_name" = GeneSymbol) %>%
  mutate(biomaRt_compatability = "No_matching") %>%
  relocate(c(
    GENEID,
    ensembl_gene_id,
    GeneSymbol,
    external_gene_name,
    biomaRt_compatability
  )) %>%
  distinct()

GeneSymbol_biomaRt_Annotated_ENSEMBL_Symbol <- rbind(
  Matching_ENSEMBL_Symbol,
  Matching_ENSEMBL_Only,
  Matching_Symbol_Only,
  No_Matching
) %>%
  arrange(GENEID)

# If need to split up:
Entrez <- getBM(
  mart = ensembl,
  attributes = c(
    "ensembl_gene_id",
    "entrezgene_id",
    "entrezgene_accession",
    "entrezgene_description"
  )
)

GO1 <- getBM(
  mart = ensembl,
  attributes = c(
    "ensembl_gene_id",
    "go_id"
  )
)

GO2 <- getBM(
  mart = ensembl,
  attributes = c(
    "go_id", # GO term accession
    "name_1006", # GO term name
    "definition_1006", # GO term definition
    # "go_linkage_type", #GO term evidence code
    "namespace_1003" # GO domain
  )
)

GO <- GO1 %>%
  inner_join(GO2,
             by = "go_id",
             multiple = "all"
  ) %>%
  dplyr::filter(go_id != "")

GO_Entrez <- Entrez %>%
  mutate(entrezgene_id = as.character(entrezgene_id)) %>%
  full_join(GO,
            by = "ensembl_gene_id",
            relationship = "many-to-many"
  ) %>%
  dplyr::mutate(entrezgene_id = na_if(entrezgene_id, "")) %>%
  dplyr::mutate(entrezgene_accession = na_if(entrezgene_accession, "")) %>%
  dplyr::mutate(entrezgene_description = na_if(entrezgene_description, ""))

# Annotation with ENSEMBL, symbols, Entrez and GO all in one
BiomaRt_ENSEMBL_Symbol_Entrez_GO_Annotation <- GeneSymbol_biomaRt_Annotated_ENSEMBL_Symbol %>%
  full_join(GO_Entrez,
            by = "ensembl_gene_id",
            relationship = "many-to-many"
  )

# nrow(BiomaRt_ENSEMBL_Symbol_Entrez_GO_Annotation %>% dplyr::select(GENEID) %>% distinct()) == nrow(GeneSymbol_biomaRt_Annotated_ENSEMBL_Symbol %>% dplyr::select(GENEID) %>% distinct())
Genes_Without_BiomaRt_GO_Annotation <- BiomaRt_ENSEMBL_Symbol_Entrez_GO_Annotation[which(is.na(BiomaRt_ENSEMBL_Symbol_Entrez_GO_Annotation$go_id)), ]
# 19k genes lack GO-annotation in the ensembl database
sig_genes %>%
  dplyr::filter(GENEID %in% Genes_Without_BiomaRt_GO_Annotation$GENEID)
# significant genes lacking GO annotation.
# Browsing the ENSEMBL genome browser, they actually do lack GO annotation in the ENSEMBL database.
# Meanwhile, some of them are part of reactome pathways

Genes_Without_Entrez <- BiomaRt_ENSEMBL_Symbol_Entrez_GO_Annotation %>%
  dplyr::filter(is.na(entrezgene_id)) %>%
  dplyr::select(GENEID) %>%
  distinct()

sig_genes %>%
  dplyr::filter(GENEID %in% Genes_Without_Entrez$GENEID)
# signifcant genes lacking entrez id.

##################################################### Annotate with annotationdbi, ENSEMBLs database which biomart connects to is missing some GO terms
keytypes(org.Mm.eg.db)
head(keys(x = org.Mm.eg.db, keytype = "ENTREZID"))
columns(org.Mm.eg.db)
Annotationdbi_ENSEMBL <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
# How many ENSEMBL ids match?
nrow(genes_for_annotation %>%
       dplyr::filter(GENEID %in% Annotationdbi_ENSEMBL))
# How many symbols match?
Annotationdbi_SYMBOL <- keys(org.Mm.eg.db, keytype = "SYMBOL")
nrow(genes_for_annotation %>%
       dplyr::filter(GeneSymbol %in% Annotationdbi_SYMBOL))
# Annotationdbi works on an entrez framework.
# Compare entrez conversion in annotationdbi with biomart
Annotationdbi_ENTREZ <- keys(x = org.Mm.eg.db, keytype = "ENTREZID")

Annotationdbi_ENTREZ_ENSEMBL_SYMBOL <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = Annotationdbi_ENTREZ,
  columns = c("ENSEMBL", "SYMBOL")
)
nrow(genes_for_annotation %>%
       dplyr::filter(GeneSymbol %in% Annotationdbi_ENTREZ_ENSEMBL_SYMBOL$SYMBOL))
# matches by symbol in annotationdbi
nrow(genes_for_annotation %>%
       dplyr::filter(GENEID %in% Annotationdbi_ENTREZ_ENSEMBL_SYMBOL$ENSEMBL))
# matches by ensembl in annotationdbi
nrow(
  BiomaRt_ENSEMBL_Symbol_Entrez_GO_Annotation %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    dplyr::select(GENEID, GeneSymbol, entrezgene_id) %>%
    distinct()
)
# matches with biomart

# Match by ensembl IDs and symbols
Annotationdbi_Matching_ENSEMBL_Symbol <- genes_for_annotation %>%
  mutate("ENSEMBL" = GENEID) %>%
  mutate("SYMBOL" = GeneSymbol) %>%
  inner_join(Annotationdbi_ENTREZ_ENSEMBL_SYMBOL,
             by = c("ENSEMBL", "SYMBOL")
  ) %>%
  mutate(Annotationdbi_compatability = "ENSEMBL_and_SYMBOL_matching") %>%
  relocate(c(
    GENEID,
    ENSEMBL,
    GeneSymbol,
    SYMBOL,
    ENTREZID,
    Annotationdbi_compatability
  )) %>%
  distinct()

# Where only ENSEMBL match
Annotationdbi_Matching_ENSEMBL_Only <- genes_for_annotation %>%
  dplyr::filter(!GENEID %in% Annotationdbi_Matching_ENSEMBL_Symbol$GENEID) %>%
  mutate("ENSEMBL" = GENEID) %>%
  inner_join(Annotationdbi_ENTREZ_ENSEMBL_SYMBOL,
             by = c("ENSEMBL"),
             relationship = "many-to-many"
  ) %>%
  mutate(Annotationdbi_compatability = "ENSEMBL_matching") %>%
  relocate(c(
    GENEID,
    ENSEMBL,
    GeneSymbol,
    SYMBOL,
    ENTREZID,
    Annotationdbi_compatability
  )) %>%
  distinct()


# Where only symbol match
Annotationdbi_Matching_SYMBOL_Only <- genes_for_annotation %>%
  dplyr::filter(!GENEID %in% Annotationdbi_Matching_ENSEMBL_Symbol$GENEID) %>%
  dplyr::filter(!GENEID %in% Annotationdbi_Matching_ENSEMBL_Only$GENEID) %>%
  mutate("SYMBOL" = GeneSymbol) %>%
  inner_join(Annotationdbi_ENTREZ_ENSEMBL_SYMBOL,
             by = c("SYMBOL"),
             relationship = "many-to-many"
  ) %>%
  mutate(Annotationdbi_compatability = "SYMBOL_matching") %>%
  relocate(c(
    GENEID,
    ENSEMBL,
    GeneSymbol,
    SYMBOL,
    ENTREZID,
    Annotationdbi_compatability
  )) %>%
  distinct()

# Where neither ENSEMBL nor symbol matched
Annotationdbi_Matching_Nothing_Only <- genes_for_annotation %>%
  dplyr::filter(!GENEID %in% Annotationdbi_Matching_ENSEMBL_Symbol$GENEID) %>%
  dplyr::filter(!GENEID %in% Annotationdbi_Matching_ENSEMBL_Only$GENEID) %>%
  dplyr::filter(!GENEID %in% Annotationdbi_Matching_SYMBOL_Only$GENEID) %>%
  mutate("ENSEMBL" = NA) %>%
  mutate("SYMBOL" = NA) %>%
  mutate("ENTREZID" = NA) %>%
  mutate(Annotationdbi_compatability = "No_matching") %>%
  relocate(c(
    GENEID,
    ENSEMBL,
    GeneSymbol,
    SYMBOL,
    ENTREZID,
    Annotationdbi_compatability
  )) %>%
  distinct()

# Bind dataframes together of all possible matches/non-matches
Annotationdbi_Annotated_ENTREZ_ENSEMBL_SYMBOL <- rbind(
  Annotationdbi_Matching_ENSEMBL_Symbol,
  Annotationdbi_Matching_ENSEMBL_Only,
  Annotationdbi_Matching_SYMBOL_Only,
  Annotationdbi_Matching_Nothing_Only
) %>%
  arrange(GENEID)

## Now GO terms with annotationdbi
keytypes(org.Mm.eg.db)
help("GO")
head(keys(x = org.Mm.eg.db, keytype = "GO"))
head(keys(x = org.Mm.eg.db, keytype = "GOALL"))
columns(org.Mm.eg.db)
# Annotationdbi_GO <- keys(org.Mm.eg.db, keytype = "GO")
Annotationdbi_GOALL <- keys(org.Mm.eg.db, keytype = "GOALL")
# GO_ENTREZID <- AnnotationDbi::select(x = org.Mm.eg.db,
#                      keys = Annotationdbi_GO,
#                      column = c("ENTREZID"),
#                      keytype = "GO") %>%
# dplyr::select(-EVIDENCE) %>%
#  distinct()
GOALL_ENTREZID <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = Annotationdbi_GOALL,
  column = c("ENTREZID"),
  keytype = "GOALL"
) %>%
  dplyr::select(-EVIDENCEALL) %>%
  distinct()

# org.Mm.eg.db does not seem to have GO names, have to use GO.db for this
keytypes(GO.db)
help("TERM")
GO_Annotation <- AnnotationDbi::select(
  x = GO.db,
  keys = Annotationdbi_GOALL,
  column = c("TERM", "DEFINITION"),
  keytype = "GOID"
) %>%
  dplyr::rename("GOALL" = GOID)

Annotationdbi_Annotated_ENTREZ_ENSEMBL_SYMBOL_GOALL <- Annotationdbi_Annotated_ENTREZ_ENSEMBL_SYMBOL %>%
  left_join(GOALL_ENTREZID,
            by = "ENTREZID",
            relationship = "many-to-many"
  ) %>%
  left_join(GO_Annotation,
            by = "GOALL"
  )


## Combine annotationdbi and biomart GO annotations
BiomaRt_Annotationdbi_GO <- Annotationdbi_Annotated_ENTREZ_ENSEMBL_SYMBOL_GOALL %>%
  dplyr::rename(
    "go_id" = GOALL,
    "name_1006" = TERM,
    "namespace_1003" = ONTOLOGYALL
  ) %>%
  dplyr::select(GENEID, go_id, name_1006, namespace_1003) %>%
  mutate(namespace_1003 = str_replace(namespace_1003, "BP", "biological_process")) %>%
  mutate(namespace_1003 = str_replace(namespace_1003, "MF", "molecular_function")) %>%
  mutate(namespace_1003 = str_replace(namespace_1003, "CC", "cellular_component")) %>%
  rbind(BiomaRt_ENSEMBL_Symbol_Entrez_GO_Annotation %>%
          dplyr::select(GENEID, go_id, name_1006, namespace_1003)) %>%
  dplyr::filter(!is.na(GENEID)) %>%
  dplyr::filter(!is.na(go_id)) %>%
  distinct() %>%
  full_join(genes_for_annotation,
            by = "GENEID"
  )

########################################## Annotate reactome gene sets with biomart
Reactome_gene <- getBM(
  mart = ensembl,
  attributes = c(
    "ensembl_gene_id",
    "reactome_gene"
  )
)

Reactome_gene_No_Empty <- Reactome_gene %>%
  dplyr::filter(reactome_gene != "") %>%
  distinct()

BiomaRt_Reactome_ENSEMBL_Symbol_Annotation <- GeneSymbol_biomaRt_Annotated_ENSEMBL_Symbol %>%
  full_join(Reactome_gene_No_Empty,
            by = "ensembl_gene_id",
            relationship = "many-to-many"
  ) %>%
  distinct()

Genes_Without_BiomaRt_Reactome_Annotation <- BiomaRt_Reactome_ENSEMBL_Symbol_Annotation[which(is.na(BiomaRt_Reactome_ENSEMBL_Symbol_Annotation$reactome_gene)), ]
# 37k genes lack GO-annotation in the ensembl database
sig_genes %>%
  dplyr::filter(GENEID %in% Genes_Without_BiomaRt_Reactome_Annotation$GENEID)
# significant genes lacking Reactome annotation

########################################## Annotate reactome with annotationdbi
keytypes(reactome.db)
Annotationdbi_Reactome_ENTREZ_ENSEMBL_SYMBOL <- as.data.frame(reactomeEXTID2PATHID) %>%
  dplyr::rename("ENTREZID" = gene_id) %>%
  dplyr::inner_join(Annotationdbi_ENTREZ_ENSEMBL_SYMBOL,
                    by = "ENTREZID",
                    relationship = "many-to-many"
  ) %>%
  dplyr::inner_join(as.data.frame(reactomePATHID2NAME),
                    by = "DB_ID"
  )

BiomaRt_Reactome_Gene_Pairs <- BiomaRt_Reactome_ENSEMBL_Symbol_Annotation %>%
  dplyr::select(GENEID, reactome_gene) %>%
  drop_na() %>%
  distinct() %>%
  mutate(Reactome_Gene_Pair = paste(GENEID, reactome_gene, sep = "_"))

Annotationdbi_Reactome_Gene_Pairs <- Annotationdbi_Reactome_ENTREZ_ENSEMBL_SYMBOL %>%
  dplyr::select(ENSEMBL, DB_ID, path_name) %>%
  drop_na() %>%
  distinct() %>%
  mutate(Reactome_Gene_Pair = paste(ENSEMBL, DB_ID, sep = "_"))

Both_Reactome_Gene_Pairs <- BiomaRt_Reactome_Gene_Pairs %>%
  full_join(Annotationdbi_Reactome_Gene_Pairs,
            by = "Reactome_Gene_Pair"
  )

# Number of reactome mappings present in annotationdbi but not biomart
nrow(Both_Reactome_Gene_Pairs %>%
       dplyr::filter(is.na(GENEID)))

# Number of reactome mappings present in biomart but not annotationdbi
nrow(Both_Reactome_Gene_Pairs %>%
       dplyr::filter(is.na(ENSEMBL)))

# Reactome pathways in biomart but not annotation dbi
Both_Reactome_Gene_Pairs %>%
  dplyr::filter(!reactome_gene %in% Both_Reactome_Gene_Pairs$DB_ID) %>%
  dplyr::select(reactome_gene) %>%
  distinct()

# Reactome pathways in annotationdbi but not biomart
Both_Reactome_Gene_Pairs %>%
  dplyr::filter(!DB_ID %in% Both_Reactome_Gene_Pairs$reactome_gene) %>%
  dplyr::select(DB_ID) %>%
  distinct()

# Combine annotationdbi with biomart reactome annotation
BiomaRt_Annotationdbi_Reactome_ENTREZ_ENSEMBL_SYMBOL <- Annotationdbi_Reactome_ENTREZ_ENSEMBL_SYMBOL %>%
  dplyr::select(ENSEMBL, DB_ID) %>%
  dplyr::rename(
    "GENEID" = ENSEMBL,
    "reactome_gene" = DB_ID
  ) %>%
  rbind(BiomaRt_Reactome_ENSEMBL_Symbol_Annotation %>%
          dplyr::select(GENEID, reactome_gene)) %>%
  drop_na() %>%
  distinct() %>%
  full_join(
    Annotationdbi_Reactome_ENTREZ_ENSEMBL_SYMBOL %>%
      dplyr::rename("reactome_gene" = DB_ID) %>%
      dplyr::select(reactome_gene, path_name) %>%
      distinct(),
    by = "reactome_gene"
  ) %>%
  full_join(genes_for_annotation,
            by = "GENEID"
  ) %>%
  mutate(path_name = ifelse(is.na(path_name), reactome_gene, path_name))

########################################## KEGG annotation
## Get KEGG pathway ID to pathway name
KEGG_1 <- enframe(
  x = keggList(
    database = "pathway",
    organism = "mmu"
  ),
  name = "KEGG_ID",
  value = "KEGG_Name"
)
## Get KEGG pathway ID to ENTREZID
KEGG_2 <- enframe(
  x = keggLink(
    "mmu",
    "pathway"
  ),
  name = "KEGG_ID",
  value = "ENTREZID"
)
## Map KEGG pathway ID and name to ENTREZ id
KEGG_3 <- KEGG_2 %>%
  mutate(KEGG_ID = str_remove(
    KEGG_ID,
    "path:"
  )) %>%
  inner_join(KEGG_1,
             by = "KEGG_ID"
  ) %>%
  mutate(ENTREZID = str_remove(
    ENTREZID,
    "mmu:"
  ))

# Does our annotationdbi entrez annotation contain all entrez ids for the KEGG pathways?
nrow(KEGG_2) == nrow(KEGG_2 %>%
                       mutate(ENTREZID = str_remove(
                         ENTREZID,
                         "mmu:"
                       )) %>%
                       dplyr::filter(ENTREZID %in% Annotationdbi_ENTREZ_ENSEMBL_SYMBOL$ENTREZID))

## Map KEGG to ENSEMBL and SYMBOL
KEGG_Annotated <- Annotationdbi_Annotated_ENTREZ_ENSEMBL_SYMBOL %>%
  full_join(KEGG_3,
            by = "ENTREZID",
            relationship = "many-to-many"
  )

########################################## Combine GO, Reactome and KEGG gene sets
annotationTable <- rbind(
  KEGG_Annotated %>%
    mutate(dbName = "KEGG") %>%
    dplyr::rename(c(
      "geneID" = GeneSymbol,
      "termID" = KEGG_ID,
      "termName" = KEGG_Name
    )) %>%
    dplyr::select(geneID, termID, termName, dbName) %>%
    distinct() %>%
    dplyr::filter(!is.na(termID)) %>%
    dplyr::filter(!is.na(geneID)),
  BiomaRt_Annotationdbi_Reactome_ENTREZ_ENSEMBL_SYMBOL %>%
    mutate(dbName = "Reactome") %>%
    dplyr::rename(c(
      "geneID" = GeneSymbol,
      "termID" = reactome_gene,
      "termName" = path_name
    )) %>%
    dplyr::select(geneID, termID, termName, dbName) %>%
    distinct() %>%
    dplyr::filter(!is.na(termID)) %>%
    dplyr::filter(!is.na(geneID)),
  BiomaRt_Annotationdbi_GO %>%
    dplyr::rename(c(
      "geneID" = GeneSymbol,
      "termID" = go_id,
      "termName" = name_1006,
      "dbName" = namespace_1003
    )) %>%
    dplyr::select(geneID, termID, termName, dbName) %>%
    distinct() %>%
    dplyr::filter(!is.na(termID)) %>%
    dplyr::filter(!is.na(geneID))
) %>%
  dplyr::filter(dbName != "")


########################################## Small scale test of SetRank
## Create background reference gene set for SetRank
## The background will be all genes with at least 1 count in 1 sample
nrow(annotationTable %>%
       dplyr::filter(dbName == "KEGG"))

nrow(annotationTable %>%
       dplyr::filter(dbName == "Reactome"))

nrow(annotationTable %>%
       dplyr::filter(dbName == "cellular_component"))

nrow(annotationTable %>%
       dplyr::filter(dbName == "biological_process"))

nrow(annotationTable %>%
       dplyr::filter(dbName == "molecular_function"))

annotationTableTest <- annotationTable %>%
  dplyr::filter(dbName == "Reactome")

referenceSetTest <- genes_for_annotation %>%
  dplyr::filter(Gene_Detected == TRUE) %>%
  dplyr::filter(GeneSymbol %in% annotationTableTest$geneID) %>%
  #slice_min(order_by = GeneSymbol, n = 1000) %>%
  dplyr::select(GeneSymbol) %>%
  distinct() %>%
  deframe()

## Create set collection object for SetRank
parallel::detectCores(all.tests = FALSE, logical = TRUE)
options(mc.cores = 10) # Adapt to the number of cores you use
collectionTest <- buildSetCollection(annotationTableTest,
                                     referenceSet = referenceSetTest,
                                     maxSetSize = 100 # Default is 500
)

## Use SetRank in ranked mode
## The genes can be ranked in many different ways
## In this case, I use adjusted p-value since this seems to be what is recommended by the paper

## Gene identifiers ranked by adjusted p-value
geneIDsTest <- Combined_Results_DF %>%
  dplyr::filter(GeneSymbol %in% referenceSetTest) %>%
  dplyr::select(
    GeneSymbol,
    adj_pvalue
  ) %>%
  distinct() %>%
  arrange(adj_pvalue) %>%
  dplyr::select(GeneSymbol) %>%
  deframe()

## And now for the actual SetRank analysis.
## CAUTION! Might take several days to complete.
networkTest <- setRankAnalysis(
  geneIDs = geneIDsTest, # Gene list ranked by adjusted p-value
  setCollection = collectionTest, # SetRank collection from above
  use.ranks = TRUE, # Ranked mode
  setPCutoff = 0.01, # This is default of 0.01
  fdrCutoff = 0.05
) # This is default of 0.05

## Export results
exportSingleResult(
  network = networkTest,
  selectedGenes = geneIDsTest,
  collection = collectionTest,
  networkName = "SetRank_NetworkTest",
  IDConverter = NULL,
  outputPath = "./R_output_files/Setrank_results"
)


########################################## GSEA with SetRank
## Create background reference gene set for SetRank
## The background will be all genes with at least 1 count in 1 sample
referenceSet <- genes_for_annotation %>%
  dplyr::filter(Gene_Detected == TRUE) %>%
  dplyr::filter(GeneSymbol %in% annotationTable$geneID) %>%
  dplyr::select(GeneSymbol) %>%
  distinct() %>%
  deframe()

## Create set collection object for SetRank
parallel::detectCores(all.tests = FALSE, logical = TRUE)
options(mc.cores = 20) # Adapt to the number of cores you use
collection <- buildSetCollection(annotationTable,
                                 referenceSet = referenceSet,
                                 maxSetSize = 500 # Default is 500
)

## Use SetRank in ranked mode
## The genes can be ranked in many different ways
## In this case, I use adjusted p-value since this seems to be what is recommended by the paper

## Gene identifiers ranked by adjusted p-value
geneIDs <- Combined_Results_DF %>%
  dplyr::filter(GeneSymbol %in% referenceSet) %>%
  dplyr::select(
    GeneSymbol,
    adj_pvalue
  ) %>%
  distinct() %>%
  arrange(adj_pvalue) %>%
  dplyr::select(GeneSymbol) %>%
  deframe()

geneIDs[duplicated(geneIDs)]

geneIDs[!geneIDs %in% referenceSet]

referenceSet[!referenceSet %in% geneIDs]

## And now for the actual SetRank analysis.
## CAUTION! Might take several days to complete.
network <- setRankAnalysis(
  geneIDs = geneIDs, # Gene list ranked by adjusted p-value
  setCollection = collection, # SetRank collection from above
  use.ranks = TRUE, # Ranked mode
  setPCutoff = 0.01, # This is default of 0.01
  fdrCutoff = 0.05
) # This is default of 0.05

## Export results
exportSingleResult(
  network = network,
  selectedGenes = geneIDs,
  collection = collection,
  networkName = "SetRank_Network",
  IDConverter = NULL,
  outputPath = "./R_output_files/Setrank_results"
)

########################################## Investigate exported gene sets
Gene_Sets <- read_tsv(file = "./R_output_files/Setrank_results/network_pathways.txt",
                      col_types = c("cccidddd"))

Significant_Gene_Sets <- Gene_Sets %>%
  dplyr::filter(pSetRank < 0.05) # SetRank employs several different p-values, see package documentation for more information

########################################## Exploring which significant gene sets the significant genes belong to

# Which significant genes are part of the significant gene sets?
Sig_Sets_Genes <- sig_genes %>%
  left_join(annotationTable, join_by(GeneSymbol == geneID)) %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue, Comparison, termID) %>%
  inner_join(Significant_Gene_Sets,
             join_by(termID == name),
             suffix = c("_Gene", "_GeneSet")
  )

# What about all detected genes which are present in significant gene sets?
Sets_Genes <- all_detected_genes %>%
  left_join(annotationTable, join_by(GeneSymbol == geneID)) %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue, Comparison, termID) %>%
  inner_join(Significant_Gene_Sets,
             join_by(termID == name),
             suffix = c("_Gene", "_GeneSet")
  )

# Which significant genes are part of the Rho GTPase gene set?
Sig_Sets_Genes %>%
  dplyr::filter(termID == "R-MMU-9012999")

# What about all detected genes which are present in the rho gtpase gene set?
Rho_GTPase <- annotationTable %>%
  dplyr::filter(termID == "R-MMU-9012999") %>%
  dplyr::rename("GeneSymbol" = geneID) %>%
  inner_join(all_detected_genes,
             by = "GeneSymbol"
  )

########################################## Cytoscape/STRING clustering annotation
## Import clustering file generated in Cytoscape with STRING and clusermaker app
Cytoscape_clustering <- read_csv("./R_input_files/Cytoscape_clustering.csv",
                                 col_types = c("cidfdddddddddddccdcdclccccccccccccdddddddddddddddd")) %>%
  dplyr::rename("GeneSymbol" = `query term`) %>%
  dplyr::rename("Cluster" = `__mclCluster`) %>%
  dplyr::select(GeneSymbol, Cluster)

Clusters_All <- annotationTable %>%
  full_join(Significant_Gene_Sets,
            join_by(termID == name),
            suffix = c("_Gene", "_GeneSet")
  ) %>%
  # dplyr::filter(pSetRank < 0.05) %>%
  dplyr::rename("GeneSymbol" = geneID) %>%
  full_join(sig_genes, by = "GeneSymbol") %>%
  full_join(Cytoscape_clustering, by = "GeneSymbol") %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue, termID, termName, dbName, pSetRank, Cluster)

Clusters <- annotationTable %>%
  full_join(Significant_Gene_Sets,
            join_by(termID == name),
            suffix = c("_Gene", "_GeneSet")
  ) %>%
  dplyr::filter(pSetRank < 0.05) %>%
  dplyr::rename("GeneSymbol" = geneID) %>%
  full_join(all_detected_genes, by = "GeneSymbol") %>%
  full_join(Cytoscape_clustering, by = "GeneSymbol") %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue, termID, termName, dbName, Cluster) %>%
  dplyr::filter(adj_pvalue < 0.05)

# Which gene sets are present in clusters 1-5, and are they down or up-regulated?
write_xlsx(
  x = Clusters %>%
    dplyr::filter(Cluster %in% 1:5) %>%
    dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
    distinct() %>%
    group_by(Cluster) %>%
    dplyr::summarise("Mean_log2FoldChange" = mean(log2FoldChange)) %>%
    full_join(Clusters, by = "Cluster") %>%
    dplyr::select(Cluster, termName, Mean_log2FoldChange) %>%
    distinct() %>%
    drop_na() %>%
    group_by(Cluster, Mean_log2FoldChange) %>%
    summarize(Significant_gene_sets = paste(termName, collapse = "\n")),
  path = "./R_output_files/Tables/Clusters.xlsx"
)

# I noticed that cluster 1 contains almost all of the most up-regulated genes
# How many up-regulated genes with log2FoldChange >2 belong to cluster 1 as compared to the other gene sets?
Clusters %>%
  dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
  distinct() %>%
  dplyr::filter(log2FoldChange > 2) %>%
  group_by(Cluster) %>%
  tally()

# Export annotated clusters
write_xlsx(
  x = Clusters %>%
    dplyr::select(Cluster, GeneSymbol, log2FoldChange, adj_pvalue) %>%
    distinct(),
  path = "./R_output_files/Tables/Significant_Genes_Cluster_Annotation.xlsx"
)

# Which clusters do the top differentially expressed genes belong to?
Clusters %>%
  group_by(GeneSymbol) %>%
  summarize(Significant_gene_sets = paste(termName, collapse = "\n")) %>%
  inner_join(sig_genes, by = "GeneSymbol") %>%
  dplyr::filter(log2FoldChange > 4 | log2FoldChange < -4) %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(GeneSymbol, Significant_gene_sets) %>%
  inner_join(Clusters, by = "GeneSymbol") %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue, Cluster, Significant_gene_sets) %>%
  distinct()

################################## Cite packages used with grateful
cite_packages(out.format = "docx", out.dir = "./R_output_files/Grateful")
