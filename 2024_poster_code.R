################################## Script start
print("Starting script 4")

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

################################## Create cytoscape results subdirectory if it doesn't exist
if (!dir.exists("./R_output_files/Cytoscape_results")) {
  dir.create("./R_output_files/Cytoscape_results")
  print("Creating directory ./R_output_files/Cytoscape_results")
} else {
  print("./R_output_files/Cytoscape_results exists")
}


################################## Data import
## Used in first script
Sample_List <- read_tsv(
  file = "./P26010/00-Reports/S.Linden_22_01_sample_info.txt",
  col_types = c("ccdd")
) %>%
  dplyr::rename(names = "NGI ID") %>%
  dplyr::rename(User_ID = "User ID") %>%
  # dplyr::rename(GreaterThan_Q30 = "≥Q30") %>% Removed for robustness, encountered problem where ≥ was misread as = by read_tsv
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

LMPC_RNA_Data <- read_tsv(
  file = "./R_input_files/LMPC_RNA_Data.txt",
  col_types = c("cifdddccidddd")
) %>%
  mutate("Date_harvested" = lubridate::as_date(Date_harvested)) %>% # I couldn't get the date to parse correctly with readr so I change from character to date after reading
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
sig_genes <- read_tsv(
  file = "./R_output_files/Tables/Significant_Genes.tsv",
  col_types = c("ccdddc")
)

all_detected_genes <- read_tsv("./R_output_files/Tables/All_Detected_Genes.tsv",
  col_types = c("ccdddc")
)

TPM_df <- read_tsv("./R_intermediate_files/TPM_df.tsv")

Mean_condition_TPM_sig_genes <- read_tsv("./R_intermediate_files/Mean_condition_TPM_sig_genes.tsv")

## Used in previous script for SetRank
annotationTable <- read_tsv(
  file = "./R_intermediate_files/annotationTable.tsv",
  col_types = c("cccc")
)

referenceSet <- read_tsv(
  file = "./R_intermediate_files/referenceSet.tsv",
  col_types = c("c")
)

geneIDs <- read_tsv(
  file = "./R_intermediate_files/geneIDs.tsv",
  col_types = c("c")
)

## Results from previous script with SetRank
Gene_Sets <- read_tsv(
  file = "./R_output_files/Setrank_results/SetRank_Network_pathways.txt",
  col_types = c("cccidddd")
)

Significant_Gene_Sets <- Gene_Sets %>%
  dplyr::filter(pSetRank < 0.05) # SetRank employs several different p-values, see package documentation for more information

############################################################## Programmatic analysis with cytoscape
## Setup
# Launch cytoscape before running the following command
# Ensures connection with cytoscape
RCy3::cytoscapePing()

# Log cytoscape version. As this script was written I used cytoscape 3.10.1
RCy3::cytoscapeVersionInfo()

# Install STRINGapp for cytoscape
RCy3::installApp("STRINGapp")

# Install Clustermaker2 app for cytoscape
RCy3::installApp("clusterMaker2")

# Install Legend Creator app for cytoscape
RCy3::installApp("Legend Creator")

#
## SetRank results analysis
#

# Load SetRank network
RCy3::importNetworkFromFile("./R_output_files/Setrank_results/SetRank_Network.net.xml")

# Rename network
RCy3::renameNetwork("SetRank")

## Hierarchical layout analysis
# Select nodes with pSetRank < 0.05
RCy3::createColumnFilter("significant_filter", 
                         "pSetRank", 
                         0.05, 
                         "LESS_THAN")

# Create subnetwork of significant nodes
RCy3::createSubnetwork(subnetwork.name = "significant_gene_sets")

# Hierarchical layout
RCy3::layoutNetwork("hierarchical nodeHorizontalSpacing=20 nodeVerticalSpacing=50")

# Interesting pattern with Mus musculus: RHO GTPase cycle at the top with all its edges directed towards it
# Select command for Rho GTPase if you need it
# RCy3::selectNodes("Mus musculus: RHO GTPase cycle", 
#                  by.col = "description")


# Apply the setrank style
RCy3::importVisualStyles("./R_output_files/Setrank_results/setrank.xml")
RCy3::setVisualStyle("SetRank")

## Modify the SetRank style
# Size
RCy3::deleteStyleMapping(
  style.name = "SetRank",
  visual.prop = "NODE_SIZE"
)
RCy3::setNodeSizeDefault(20, style.name = "SetRank")

# Labels
RCy3::setNodeLabelColorDefault("#000000", 
                               style.name = "SetRank")
RCy3::setNodeLabelPositionDefault(
  new.nodeAnchor = "E",
  new.graphicAnchor = "W",
  new.justification = "c",
  new.xOffset = 0.00,
  new.yOffset = 0.00,
  style.name = "SetRank"
)
setNodeLabelPositionBypass(c("mmu00190",
                             "mmu04932",
                             "R-MMU-9711123"),
                           "W,E,c,0.00,0.00")

setNodeLabelPositionBypass(c("R-MMU-156827",
                             "mmu05171",
                             "R-MMU-9006934"),
                           "E,W,c,0.00,0.00")

# Node fill
RCy3::deleteStyleMapping(
  style.name = "SetRank",
  visual.prop = "NODE_FILL_COLOR"
)
RCy3::setNodeColorDefault("#000000", 
                          style.name = "SetRank")
RCy3::setNodeSelectionColorDefault("#63e5ff", 
                                   style.name = "SetRank")
# Map node color to pSetRank in range 0 to 0.05 with viridis plasma palette
setrank_color_scale_values <- seq(from = 0, to = 0.05, length.out = 64)

setrank_viridis_hex_codes <- viridis::plasma(n = 64,
                                      direction = -1) %>%
  stringr::str_replace("FF$","")
RCy3::setNodeColorMapping("pSetRank",
                          table.column.values = setrank_color_scale_values,
                          mapping.type = "c",
                          colors = setrank_viridis_hex_codes,
                          style.name = "SetRank"
)

# Node border
RCy3::deleteStyleMapping(
  style.name = "SetRank",
  visual.prop = "NODE_BORDER_PAINT"
)
RCy3::setNodeBorderColorDefault("#000000",
                                style.name = "SetRank")
RCy3::setNodeBorderWidthDefault(1,
                                style.name = "SetRank")

# Node shapes
RCy3::setNodeShapeDefault("ELLIPSE",
                          style.name = "SetRank"
)
# Arrows
RCy3::deleteStyleMapping(
  style.name = "SetRank",
  visual.prop = "EDGE_TARGET_ARROW_UNSELECTED_PAINT"
)
RCy3::setEdgeTargetArrowColorDefault("#000026",
                                     style.name = "SetRank"
)
RCy3::deleteStyleMapping(
  style.name = "SetRank",
  visual.prop = "EDGE_SOURCE_ARROW_UNSELECTED_PAINT"
)
RCy3::setEdgeSourceArrowColorDefault("#a9a9a9",
                                     style.name = "SetRank"
)

# Edge lines
RCy3::deleteStyleMapping(
  style.name = "SetRank",
  visual.prop = "EDGE_WIDTH"
)
RCy3::setEdgeLineWidthDefault(0.5,
                              style.name = "SetRank"
)
RCy3::deleteStyleMapping(
  style.name = "SetRank",
  visual.prop = "EDGE_STROKE_UNSELECTED_PAINT"
)

# Set white background
setBackgroundColorDefault(
  "#FFFFFF",
  style.name = "SetRank")

# Set a more universally available font than helvetica which seems to be the default
setNodeFontFaceDefault("Calibri", 
                       style.name = "SetRank")

# Legends need to be annotated manually, but Legend Creator app helps with this

# Export image
# If your operating system is set to a language with characters not present in English, this might cause encoding errors with SVG
RCy3::exportImage(
  "./R_output_files/Cytoscape_results/SetRank_Hierchical_Raw.svg",
  type = "SVG"
)




## Cluster analysis of significant genes based on STRING interactions
#

## First we need to construct a network based on the string interactions
## Using STRINGdb package, see the vignette for more information
string_db <- STRINGdb$new(
  version = "12", # Latest version as of writing
  species = 10090, # Mouse
  score_threshold = 400, # Default
  network_type = "full", # Full network
  input_directory = "./R_input_files",
  protocol = "http"
)

# Map to string identifiers
if (!file.exists("./R_input_files/stringdb-v12-10090-protein-aliases.gz")) {
  download.file(
    url = "https://stringdb-downloads.org/download/protein.aliases.v12.0/10090.protein.aliases.v12.0.txt.gz",
    destfile = "./R_input_files/stringdb-v12-10090-protein-aliases.gz"
  )
} else {
  print("./R_input_files/stringdb-v12-10090-protein-aliases.gz exists")
}

mouse_alias_df <- read_tsv("./R_input_files/stringdb-v12-10090-protein-aliases.gz",
  col_types = c("ccc")
) %>%
  dplyr::rename("string_protein_id" = `#string_protein_id`)

# Mapping by genesymbol
String_GeneSymbol_Mapping <- sig_genes %>%
  dplyr::rename("alias" = GeneSymbol) %>%
  left_join(mouse_alias_df,
    by = "alias"
  ) %>%
  dplyr::select(-c(source)) %>%
  dplyr::distinct(alias, .keep_all = TRUE) # Pragmatic way to handle multiple matching, but another option would be to make use of groups in cytoscape

# Plot the string network for mapped significant genes
sig_genes_mapped_string_id <- String_GeneSymbol_Mapping %>%
  dplyr::filter(!is.na(string_protein_id))

# Get string interactions
if (!file.exists("./R_input_files/stringdb-v12-10090-protein-links.gz")) {
  download.file(
    url = "https://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz",
    destfile = "./R_input_files/stringdb-v12-10090-protein-links.gz"
  )
} else {
  print("./R_input_files/stringdb-v12-10090-protein-links.gz exists")
}

mouse_string_interactions_input_df <- readr::read_delim("./R_input_files/stringdb-v12-10090-protein-links.gz",
  delim = " ",
  col_types = c("ccd")
)
mouse_string_interactions_significant_df <- mouse_string_interactions_input_df %>%
  dplyr::filter(protein1 %in% sig_genes_mapped_string_id$string_protein_id) %>%
  dplyr::filter(protein2 %in% sig_genes_mapped_string_id$string_protein_id) %>%
  dplyr::filter(combined_score >= 400) # Default for string database, "medium confidence"

# Using igraph to create network
significant_genes_igraph <- igraph::graph_from_data_frame(mouse_string_interactions_significant_df,
  directed = FALSE,
  vertices = NULL
)

# Simplify igraph, for example to remove loops
significant_genes_igraph <- igraph::simplify(significant_genes_igraph,
  edge.attr.comb = "mean"
)

# Create cytoscape network from igraph
RCy3::createNetworkFromIgraph(
  significant_genes_igraph,
  title = "Significant genes string interactions network",
  collection = "lmpc_rnaseq"
)

# Add annotations for network
# remove generic gene set descriptions " - Mus musculus (house mouse)" from KEGG,
# "Mus musculus: " from reactome,
# "molecular_function", "cellular_component", and "biological_process" from GO
annotationTable_for_cytoscape <- annotationTable %>%
  dplyr::rename("GeneSymbol" = geneID) %>%
  dplyr::mutate(termName = str_remove_all(termName, "- Mus musculus \\(house mouse\\)")) %>%
  dplyr::mutate(termName = str_remove_all(termName, "Mus musculus:")) %>%
  dplyr::mutate(termName = str_remove_all(termName, "molecular_function")) %>%
  dplyr::mutate(termName = str_remove_all(termName, "cellular_component")) %>%
  dplyr::mutate(termName = str_remove_all(termName, "biological_process")) %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise("Concatenated_Gene_Set_Descriptions" = toString(termName)) %>%
  dplyr::mutate(Concatenated_Gene_Set_Descriptions = str_remove_all(Concatenated_Gene_Set_Descriptions, ",")) %>%
  dplyr::mutate(Concatenated_Gene_Set_Descriptions = str_remove_all(Concatenated_Gene_Set_Descriptions, "^ *")) %>%
  dplyr::mutate(Concatenated_Gene_Set_Descriptions = str_replace_all(Concatenated_Gene_Set_Descriptions, " +", " "))

sig_genes_cytoscape_annotation_df <- sig_genes_mapped_string_id %>%
  dplyr::rename("GeneSymbol" = alias) %>%
  dplyr::left_join(annotationTable_for_cytoscape, by = "GeneSymbol")

# In order to get unconnected but significant nodes in to cytoscape, modification below
addCyNodes(sig_genes_cytoscape_annotation_df$string_protein_id)

RCy3::loadTableData(sig_genes_cytoscape_annotation_df,
  data.key.column = "string_protein_id",
  table = "node",
  table.key.column = "name",
  namespace = "default",
  network = NULL
)

# Custom network style
style.name <- "lmpc-rnaseq"
defaults <- list(
  NODE_SHAPE = "ellipse",
  NODE_SIZE = 18,
  NODE_SHAPE = "ELLIPSE",
  NODE_PAINT = "#000000",
  NODE_LABEL_FONT_SIZE = 18,
  NODE_BORDER_WIDTH = 5,
  EDGE_TRANSPARENCY=120,
  NODE_LABEL_POSITION = "E,W,c,0.00,0.00" # See setNodeCustomPosition
)
RCy3::createVisualStyle(style.name, defaults)
RCy3::setNodeLabelMapping("GeneSymbol", style.name = style.name)
RCy3::setNodeColorMapping("log2FoldChange",
                          table.column.values = c(-2, 0, 2),
                          mapping.type = "c",
                          colors = c("#0000FF","#FFFFFF","#FF0000"),
                          style.name = style.name
)
RCy3::setVisualStyle(style.name)

# Cluster with ClusterMaker2
cluster_command <- paste(
  "cluster mcl",
  "createGroups=false",
  "inflation_parameter=2.5"
) # 2.5 is default
RCy3::commandsGET(cluster_command)

# Extract mcl clustering
Cytoscape_clustering <- RCy3::getTableColumns(table = "node",
                                              columns = c("GeneSymbol","__mclCluster")) %>%
  dplyr::rename("Cluster" = "__mclCluster")

GeneSymbols_Clusters_1_2_3 <- Cytoscape_clustering %>%
  dplyr::inner_join(sig_genes_cytoscape_annotation_df,
                    by = "GeneSymbol") %>%
  dplyr::filter(Cluster %in% c("1","2","3")) %>%
  dplyr::select(string_protein_id) %>%
  tibble::deframe()

GeneSymbols_Clusters_1 <- Cytoscape_clustering %>%
  dplyr::inner_join(sig_genes_cytoscape_annotation_df,
                    by = "GeneSymbol") %>%
  dplyr::filter(Cluster %in% c("1")) %>%
  dplyr::select(string_protein_id) %>%
  tibble::deframe()

GeneSymbols_Clusters_2 <- Cytoscape_clustering %>%
  dplyr::inner_join(sig_genes_cytoscape_annotation_df,
                    by = "GeneSymbol") %>%
  dplyr::filter(Cluster %in% c("2")) %>%
  dplyr::select(string_protein_id) %>%
  tibble::deframe()

GeneSymbols_Clusters_3 <- Cytoscape_clustering %>%
  dplyr::inner_join(sig_genes_cytoscape_annotation_df,
                    by = "GeneSymbol") %>%
  dplyr::filter(Cluster %in% c("3")) %>%
  dplyr::select(string_protein_id) %>%
  tibble::deframe()

LFC_Above_2 <- Cytoscape_clustering %>%
  dplyr::inner_join(sig_genes_cytoscape_annotation_df,
                    by = "GeneSymbol") %>%
  dplyr::filter(log2FoldChange > 2)%>%
  dplyr::select(string_protein_id) %>%
  tibble::deframe()

LFC_Below_2 <- Cytoscape_clustering %>%
  dplyr::inner_join(sig_genes_cytoscape_annotation_df,
                    by = "GeneSymbol") %>%
  dplyr::filter(log2FoldChange < -2)%>%
  dplyr::select(string_protein_id) %>%
  tibble::deframe()

# There are 3 obvious clusters, visualize them
RCy3::setNodeFillOpacityDefault(255,
                                style.name)
RCy3::setNodeBorderOpacityDefault(255,
                                style.name)
RCy3::setNodeBorderOpacityDefault(255,
                                  style.name)
RCy3::setNodeLabelOpacityDefault(255,
                                 style.name)
#RCy3::setNodeFillOpacityBypass(GeneSymbols_Clusters_1_2_3,
#                               255)
#RCy3::setNodeBorderOpacityBypass(GeneSymbols_Clusters_1_2_3,
#                               255)
#RCy3::setNodeBorderOpacityBypass(GeneSymbols_Clusters_1_2_3,
#                               255)
#RCy3::setNodeLabelOpacityBypass(GeneSymbols_Clusters_1_2_3,
#                               255)
#setNodeBorderColorBypass(GeneSymbols_Clusters_1,
#                         "#008C00")
#setNodeBorderColorBypass(GeneSymbols_Clusters_2,
#                         "#FF00FD")
#setNodeBorderColorBypass(GeneSymbols_Clusters_3,
#                         "#DB9D00")
setNodeBorderColorBypass(LFC_Above_2,
                         "#FF7000")
setNodeBorderColorBypass(LFC_Below_2,
                         "#0AB5D2")
setNodeShapeBypass(GeneSymbols_Clusters_1,
                         "DIAMOND")
setNodeShapeBypass(GeneSymbols_Clusters_2,
                         "RECTANGLE")
setNodeShapeBypass(GeneSymbols_Clusters_3,
                         "TRIANGLE")
#getLayoutNames()
# Change layout to kamada-kawai
# The higher the confidence score, the closer the nodes will be
#getLayoutPropertyNames("kamada-kawai")
RCy3::layoutNetwork("kamada-kawai edgeAttribute=combined_score randomize=FALSE")

RCy3::setCurrentNetwork("Significant genes string interactions network")
RCy3::selectNodes(c(""),
                  by.col = "__mclCluster",
                  preserve.current.selection = FALSE)
RCy3::exportImage(
  "./2024_Poster/LMPC_2024_Infection_Poster_Network.svg",
  "svg"
)

# List l2fc of unconnected nodes.
sig_unconnected_nodes <- sig_genes_cytoscape_annotation_df %>% 
  dplyr::filter(!string_protein_id %in% c(mouse_string_interactions_significant_df$protein1, mouse_string_interactions_significant_df$protein2)) %>% 
  arrange(log2FoldChange)

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

########################################## Cytoscape/STRING clustering annotation
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

# Which gene sets are present in clusters 1-3, and are they down or up-regulated?
write_xlsx(
  x = Clusters_Sig %>%
    dplyr::filter(Cluster %in% 1:3) %>%
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

# Vector of genes present in cluster 1
Cluster1_genes <- Clusters_Sig %>%
  dplyr::filter(Cluster %in% c(1)) %>%
  dplyr::select(GeneSymbol) %>%
  dplyr::distinct() %>%
  tibble::deframe()

# I noticed that cluster 1 contains almost all of the most up-regulated genes
# How many up-regulated genes with log2FoldChange >2 belong to cluster 1 as compared to the other gene sets?
Clusters_Sig %>%
  dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
  distinct() %>%
  dplyr::filter(log2FoldChange > 2) %>%
  group_by(Cluster) %>%
  tally()

# How about genes with lfch > 0
Clusters_Sig %>%
  dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
  distinct() %>%
  dplyr::filter(log2FoldChange > 0) %>%
  group_by(Cluster) %>%
  tally()

# And specifically, cluster 1 vs everything else
Clusters_Sig %>%
  dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
  distinct() %>%
  dplyr::filter(log2FoldChange > 2) %>% 
  dplyr::mutate("is_cluster_1" = Cluster == 1) %>% 
  dplyr::mutate("is_cluster_1" = replace_na(is_cluster_1, FALSE)) %>% 
  group_by(is_cluster_1) %>%
  tally() %>% 
  mutate("percentage" = n / sum(n)) %>% 
  mutate("sum_n" = sum(n))

Clusters_Sig %>%
  dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
  distinct() %>%
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::mutate("is_cluster_1" = Cluster == 1) %>% 
  dplyr::mutate("is_cluster_1" = replace_na(is_cluster_1, FALSE)) %>% 
  group_by(is_cluster_1) %>%
  tally() %>% 
  mutate("percentage" = n / sum(n)) %>% 
  mutate("sum_n" = sum(n))

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
  dplyr::select(Cluster, GeneSymbol, log2FoldChange) %>%
  dplyr::distinct() %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate("Mean_Cluster_log2FoldChange" = mean(log2FoldChange)) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(Clusters_Sig, by = c("Cluster", "GeneSymbol", "log2FoldChange")) %>%
  dplyr::select(
    Cluster,
    Mean_Cluster_log2FoldChange,
    GeneSymbol,
    log2FoldChange,
    adj_pvalue,
    termID,
    termName,
    dbName,
    pSetRank
  ) %>%
  distinct() %>%
  dplyr::full_join(Mean_condition_TPM_sig_genes, by = "GeneSymbol")

# Double checking all significant genes are with us
# Same number of entries?
nrow(sig_genes %>% distinct()) == nrow(Sig_Sets_Clusters_Genes_TPM %>% dplyr::select(GeneSymbol) %>% distinct())
# Which entries do not match?
sig_genes[!deframe(sig_genes %>% dplyr::select(GeneSymbol) %>% arrange(GeneSymbol)) == deframe(Sig_Sets_Clusters_Genes_TPM %>% dplyr::select(GeneSymbol) %>% arrange(GeneSymbol) %>% distinct()), ]

write_xlsx(
  x = Sig_Sets_Clusters_Genes_TPM,
  path = "./R_output_files/Tables/summarised-results-of-significant-genes-sets-clusters.xlsx"
)

##########################################  Citing packages
#renv::install("knitr")
#library("knitr")
#renv::install("bookdown")
#library("bookdown")

#knitr::write_bib(c(.packages(), "bookdown"), "packages.bib")

########################################## SessionInfo
sessionInfo()

print("Script 4 finished")
