################################## Set seed for reproducibility
set.seed(1337)

################################## Load packages, commented out what loaded packages for renv
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
## Used in first script
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

tx2gene <- read_tsv(
  file = "./P26010/01-RNA-Results/genome/salmon_tx2gene.tsv",
  col_names = c(
    "TXNAME",
    "GENEID",
    "GeneSymbol"
  ),
  col_types = (c("ccc"))
)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
glimpse(txi$counts)
colnames(txi$counts) <- Sample_List$names
glimpse(txi$counts)

## Results from first script
sig_genes <- read_tsv(file = "./R_output_files/Tables/Significant_Genes.tsv",
col_types = c("ccdddc"))

all_detected_genes <- read_tsv("./R_output_files/Tables/All_Detected_Genes.tsv",
col_types = c("ccdddc"))

TPM_df <- read_tsv("./R_intermediate_files/TPM_df.tsv")

Mean_condition_TPM_sig_genes <- read_tsv("./R_intermediate_files/Mean_condition_TPM_sig_genes.tsv")

## Used in previous script for SetRank
annotationTable <- read_tsv(file = "./R_intermediate_files/annotationTable.tsv",
col_types = c("cccc"))

referenceSet <- read_tsv(file = "./R_intermediate_files/referenceSet.tsv",
col_types = c("c"))

geneIDs <- read_tsv(file = "./R_intermediate_files/geneIDs.tsv",
col_types = c("c"))

## Results from previous script with SetRank
Gene_Sets <- read_tsv(file = "./R_output_files/Setrank_results/SetRank_Network_pathways.txt",
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

# Join with gene sets
Clusters_All <- annotationTable %>%
  full_join(Gene_Sets,
            join_by(termID == name),
            suffix = c("_Gene", "_GeneSet")
  ) %>%
  # dplyr::filter(pSetRank < 0.05) %>%
  dplyr::rename("GeneSymbol" = geneID) %>%
  full_join(sig_genes, by = "GeneSymbol") %>%
  full_join(Cytoscape_clustering, by = "GeneSymbol") %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue, termID, termName, dbName, pSetRank, Cluster)

Clusters_Sig <- annotationTable %>%
  full_join(Significant_Gene_Sets,
            join_by(termID == name),
            suffix = c("_Gene", "_GeneSet")
  ) %>%
  dplyr::filter(pSetRank < 0.05) %>%
  dplyr::rename("GeneSymbol" = geneID) %>%
  full_join(all_detected_genes, by = "GeneSymbol") %>%
  full_join(Cytoscape_clustering, by = "GeneSymbol") %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue, termID, termName, dbName, pSetRank, Cluster) %>%
  dplyr::filter(adj_pvalue < 0.05)

# Which gene sets are present in clusters 1-5, and are they down or up-regulated?
write_xlsx(
  x = Clusters_Sig %>%
    dplyr::filter(Cluster %in% 1:5) %>%
    dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
    distinct() %>%
    group_by(Cluster) %>%
    dplyr::summarise("Mean_log2FoldChange" = mean(log2FoldChange)) %>%
    full_join(Clusters_Sig, by = "Cluster") %>%
    dplyr::select(Cluster, termName, Mean_log2FoldChange) %>%
    distinct() %>%
    drop_na() %>%
    group_by(Cluster, Mean_log2FoldChange) %>%
    summarize(Significant_gene_sets = paste(termName, collapse = "\n")),
  path = "./R_output_files/Tables/Clusters.xlsx"
)

# I noticed that cluster 1 contains almost all of the most up-regulated genes
# How many up-regulated genes with log2FoldChange >2 belong to cluster 1 as compared to the other gene sets?
Clusters_Sig %>%
  dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
  distinct() %>%
  dplyr::filter(log2FoldChange > 2) %>%
  group_by(Cluster) %>%
  tally()

# Export annotated clusters
write_xlsx(
  x = Clusters_Sig %>%
    dplyr::select(Cluster, GeneSymbol, log2FoldChange, adj_pvalue) %>%
    distinct(),
  path = "./R_output_files/Tables/Significant_Genes_Cluster_Annotation.xlsx"
)

# Which clusters do the top differentially expressed genes belong to?
Clusters_Sig %>%
  group_by(GeneSymbol) %>%
  summarize(Significant_gene_sets = paste(termName, collapse = "\n")) %>%
  inner_join(sig_genes, by = "GeneSymbol") %>%
  dplyr::filter(log2FoldChange > 4 | log2FoldChange < -4) %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(GeneSymbol, Significant_gene_sets) %>%
  inner_join(Clusters_Sig, by = "GeneSymbol") %>%
  dplyr::select(GeneSymbol, log2FoldChange, adj_pvalue, Cluster, Significant_gene_sets) %>%
  distinct()

########################################## Combine gene set and clustering annotations for all significant genes
# A sort of master results file with all significant gene sets, all significant genes, as well as clustering info.
Sig_Sets_Clusters_Genes_TPM <- Clusters_Sig %>%
dplyr::select(Cluster,GeneSymbol,log2FoldChange) %>%
dplyr::distinct()  %>% 
dplyr::group_by(Cluster) %>%
dplyr::mutate("Mean_Cluster_log2FoldChange" = mean(log2FoldChange)) %>%
dplyr::ungroup() %>%
dplyr::inner_join(Clusters_Sig, by = c("Cluster", "GeneSymbol", "log2FoldChange")) %>%
    dplyr::select(Cluster,
    Mean_Cluster_log2FoldChange,
    GeneSymbol, 
    log2FoldChange, 
    adj_pvalue,
    termID,
    termName,
    dbName,
    pSetRank) %>%
    distinct() %>%
dplyr::full_join(Mean_condition_TPM_sig_genes, by = "GeneSymbol")

# Double checking all significant genes are with us
# Same number of entries?
nrow(sig_genes %>% distinct) == nrow(Sig_Sets_Clusters_Genes_TPM %>% dplyr::select(GeneSymbol) %>% distinct)
# Which entries do not match?
sig_genes[!deframe(sig_genes %>% dplyr::select(GeneSymbol) %>% arrange(GeneSymbol)) == deframe(Sig_Sets_Clusters_Genes_TPM %>% dplyr::select(GeneSymbol) %>% arrange(GeneSymbol) %>% distinct()),]

write_xlsx(
  x = Sig_Sets_Clusters_Genes_TPM,
  path = "./R_output_files/Tables/summarised-results-of-significant-genes-sets-clusters.xlsx"
)
