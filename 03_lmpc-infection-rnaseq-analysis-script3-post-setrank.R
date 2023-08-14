################################## Set seed for reproducibility
set.seed(1337)

################################## Load packages, commented out what loaded packages for renv
renv::restore()

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
    "igraph"# For RCy3/cytoscape
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
  #dplyr::rename(GreaterThan_Q30 = "≥Q30") %>% Removed for robustness, encountered problem where ≥ was misread as = by read_tsv
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
  col_types = c("cidfdddddddddddccdcdclccccccccccccdddddddddddddddd")
) %>%
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

############################################################## Planning steps cytoscape
# Install stringapp
# installApp('STRINGapp') 

# How to map with identifier? Do I want to map with GeneSymbols? EnsemblIDs?

# Load SetRank network

# Annotate with StringApp
# for more information on string commands:
# commandsHelp('string')
# commandsHelp('string disease query'
# Use STRINGdb package to create string network?
# https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html

# Set up  layout

# Analyse similar to above, for example with clustering and hierchical layout

# Use group nodes based on for example cluster, and collapse them for better viz?


############################################################## Programmatic analysis with cytoscape
#Launch cytoscape before running the following command
#Ensures connection with cytoscape
RCy3::cytoscapePing()

#Log cytoscape version. As this script was written I used cytoscape 3.10.0
RCy3::cytoscapeVersionInfo()

#Install STRINGapp for cytoscape
RCy3::installApp("STRINGapp")

#Install CLustermaker2 app for cytoscape
RCy3::installApp("clusterMaker2")

#Install AutoAnnotate app for cytoscape
RCy3::installApp("AutoAnnotate")

## SetRank results analysis
#Load SetRank network
RCy3::importNetworkFromFile("./R_output_files/Setrank_results/SetRank_Network.net.xml")

#Apply the setrank style
RCy3::importVisualStyles("./R_output_files/Setrank_results/setrank.xml")

RCy3::setVisualStyle("SetRank")

#Select nodes with pSetRank < 0.05
RCy3::createColumnFilter("significant_filter", "pSetRank", 0.05, "LESS_THAN")

#Create subnetwork of significant nodes
RCy3::createSubnetwork(subnetwork.name = "significant_gene_sets")

#Apply layout Edge-Weighted Spring Embedded - jaccard, recommended in SetRank vignette
#getLayoutNames()
RCy3::getLayoutNameMapping()
#kamada-kawai = Edge-weighted Spring Embedded Layout

#And to set it up with jaccard distances?
RCy3::getLayoutPropertyNames("kamada-kawai") 

#Apply kamada-kawai/edgeweighted spring embedded layout with jaccard as edge attribute
RCy3::layoutNetwork("kamada-kawai edgeAttribute=jaccard")

#Reduce node size so they don't overlap as much
#First need to delete size mapping
getVisualPropertyNames()
RCy3::deleteStyleMapping(style.name = "SetRank", 
                   visual.prop = "NODE_SIZE")

RCy3::setNodeSizeDefault(50, style.name = "SetRank")

#Change text color so it's readable on white background
RCy3::setNodeLabelColorDefault("#000000", style.name = "SetRank")

#Change node color to something light so it's easier to read the text
RCy3::deleteStyleMapping(style.name = "SetRank", 
                   visual.prop = "NODE_FILL_COLOR")

#Don't see any interesting patterns here
#Lets try hierchical layout
RCy3::getLayoutPropertyNames("hierarchical") 

#Improve spacing
RCy3::layoutNetwork("hierarchical nodeHorizontalSpacing=4 nodeVerticalSpacing=2")

#Interesting pattern with Mus musculus: RHO GTPase cycle at the top with all adges directed towards it
#Not sure why mRNA processing ends up so high though, only one node directed there.
#Select Rho GTPase node
RCy3::selectNodes("Mus musculus: RHO GTPase cycle", by.col = "description")

#Change color of selected nodes to make it stand out
RCy3::setNodeSelectionColorDefault("#63e5ff", style.name = "SetRank")

#To make nodes with incoming edges stand out, cchange target arrow color
RCy3::deleteStyleMapping(style.name = "SetRank", 
                         visual.prop = "EDGE_TARGET_ARROW_UNSELECTED_PAINT")
RCy3::setEdgeTargetArrowColorDefault("#000026",
                                     style.name = "SetRank")

#Likewise, to take attention away from outgoing edges change that color to something more discret
RCy3::deleteStyleMapping(style.name = "SetRank", 
                         visual.prop = "EDGE_SOURCE_ARROW_UNSELECTED_PAINT")
RCy3::setEdgeSourceArrowColorDefault("#a9a9a9",
                                     style.name = "SetRank")

#The line widths are based on significantJaccard in a way that doesn't make sense to me
#Lets make line widths the same for all, they are distracting
RCy3::deleteStyleMapping(style.name = "SetRank", 
                         visual.prop = "EDGE_WIDTH")
RCy3::setEdgeLineWidthDefault(1,
                              style.name = "SetRank")

#Change shapes to fit label text
RCy3::setNodeShapeDefault("ROUND_RECTANGLE", 
                          style.name = "SetRank")
RCy3::setNodeWidthDefault(160,
                          style.name = "SetRank")

#Export image
RCy3::exportImage("./R_output_files/Cytoscape_results/SetRank_Hierchical.svg", 
                  "svg")

#To further improve readability I think manual fixes are need
#For example, better node label placements.
#I can't figure out a way to do this in R for this graph
#Edge legends might also have to be done manually

##Now time for cluster analysis of significant genes based on STRING interactions
##First we need to construct a network based on the string interactions
##Using STRINGdb package, see the vignette for more information
string_db <- STRINGdb$new( version="12", #Latest version as of writing
                           species=10090, #Mouse
                           score_threshold=400, #Default
                           network_type="full", #Full network
                           input_directory="./R_input_files",
                           protocol="http")

##Listing methods
STRINGdb$methods()
##Info about get_graph
STRINGdb$help("get_graph")

#Map to string identifiers
#string_db$map(sig_genes, "GeneSymbol", removeUnmappedRows = TRUE)
#Seems like the stringdb-static.org website which stringdb use is down as of writing
#Will do mapping manually instead
if (!file.exists("./R_input_files/mouse_string_alias.gz")) {
  download.file(url = "https://stringdb-downloads.org/download/protein.aliases.v12.0/10090.protein.aliases.v12.0.txt.gz",
                destfile = "./R_input_files/mouse_string_alias.gz")
  mouse_alias_df <- read_tsv("./R_input_files/mouse_string_alias.gz",
                             col_types = c("ccc")) %>%
    dplyr::rename("string_protein_id" = `#string_protein_id`)
} else {
  print("./R_input_files/mouse_string_alias.gz exists")
}


#Explore how mapping by gene symbol vs by ensembl id works
#genesymbol
String_GeneSymbol_Mapping <- sig_genes %>%
  dplyr::rename("alias" = GeneSymbol) %>%
  left_join(mouse_alias_df,
            by = "alias") %>%
  dplyr::select(-c(source)) %>%
  dplyr::distinct(alias, .keep_all = TRUE) # Pragmatic way to handle multiple matching, but another option would be to make use of groups in cytoscape
#Any unmapped?
String_GeneSymbol_Mapping %>%
  dplyr::filter(is.na(string_protein_id))
#ensembl id
String_EnsemblID_Mapping <- sig_genes %>%
  dplyr::rename("alias" = GENEID) %>%
  left_join(mouse_alias_df,
            by = "alias") %>%
  dplyr::select(-c(source)) %>%
  dplyr::distinct(alias, .keep_all = TRUE) # Pragmatic way to handle multiple matching, but another option would be to make use of groups in cytoscape
#Any unmapped?
String_EnsemblID_Mapping %>%
  dplyr::filter(is.na(string_protein_id))
#Mapping by genesymbol was slightly better, got 11 unmatched instead of 12.

#Plot the string network for mapped significant genes
sig_genes_mapped_string_id <- String_GeneSymbol_Mapping %>%
  dplyr::filter(!is.na(string_protein_id))

#The following command only works when the stringdb package website is up
#string_db$get_interactions(sig_genes_mapped_string_id$string_protein_id)

#Which means I have to do it manually
if (!file.exists("./R_input_files/mouse_string_interactions.gz")) {
  download.file(url = "https://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz",
                destfile = "./R_input_files/mouse_string_interactions.gz")
  mouse_string_interactions_input_df <- readr::read_delim("./R_input_files/mouse_string_interactions.gz",
                                                          delim = " ",
                                                          col_types = c("ccd"))
} else {
  print("./R_input_files/mouse_string_interactions.gz exists")
}


mouse_string_interactions_significant_df <- mouse_string_interactions_input_df %>% 
  dplyr::filter(protein1 %in% sig_genes_mapped_string_id$string_protein_id) %>%
  dplyr::filter(protein2 %in% sig_genes_mapped_string_id$string_protein_id) %>%
  dplyr::filter(combined_score >= 400) #Default for string database, "medium confidence"

#Using igraph to create network
significant_genes_igraph <- igraph::graph_from_data_frame(mouse_string_interactions_significant_df, 
                                                  directed = FALSE, 
                                                  vertices = NULL)

#Simplify igraph, for example to remove loops
significant_genes_igraph <- igraph::simplify(significant_genes_igraph,
                                             edge.attr.comb = "mean")

#Create cytoscape network from igraph
createNetworkFromIgraph(
  significant_genes_igraph,
  title = "Significant genes string interactions network",
  collection = "lmpc_rnaseq"
)

# Change layout
# The higher the confidence score, the closer the nodes will be kamada-kawai below
RCy3::layoutNetwork("kamada-kawai edgeAttribute=combined_score")

#Add annotations for network
#for downstream annotation with autoannotate, process annotationtable
#remove generic gene set descriptions " - Mus musculus (house mouse)" from KEGG, 
#"Mus musculus: " from reactome,
#"molecular_function", "cellular_component", and "biological_process" from GO
annotationTable_for_cytoscape <- annotationTable %>%
  dplyr::rename("GeneSymbol" = geneID) %>%
  dplyr::mutate(termName = str_remove_all(termName, "- Mus musculus \\(house mouse\\)")) %>%
  dplyr::mutate(termName = str_remove_all(termName, "Mus musculus:")) %>%
  dplyr::mutate(termName = str_remove_all(termName, "molecular_function")) %>%
  dplyr::mutate(termName = str_remove_all(termName, "cellular_component")) %>%
  dplyr::mutate(termName = str_remove_all(termName, "biological_process")) %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise("Concatenated_Gene_Set_Descriptions" = toString(termName)) %>%
  dplyr::mutate(Concatenated_Gene_Set_Descriptions = str_remove_all(Concatenated_Gene_Set_Descriptions,",")) %>%
  dplyr::mutate(Concatenated_Gene_Set_Descriptions = str_remove_all(Concatenated_Gene_Set_Descriptions,"^ *")) %>%
  dplyr::mutate(Concatenated_Gene_Set_Descriptions = str_replace_all(Concatenated_Gene_Set_Descriptions," +"," "))

sig_genes_cytoscape_annotation_df <- sig_genes_mapped_string_id %>% 
  dplyr::rename("GeneSymbol" = alias) %>%
  dplyr::left_join(annotationTable_for_cytoscape, by = "GeneSymbol")
loadTableData(sig_genes_cytoscape_annotation_df,
              data.key.column = "string_protein_id",
              table = "node",
              table.key.column = "name",
              namespace = "default",
              network = NULL)

#Improve network style
style.name = "lmpc-rnaseq"
defaults <- list(NODE_SHAPE="ellipse",
                 NODE_SIZE=300,
                 #EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="E,W,c,0.00,0.00" #See setNodeCustomPosition
                 )

RCy3::createVisualStyle(style.name, defaults)
RCy3::setVisualStyle(style.name)
#RCy3::setNodeShapeDefault("ROUND_RECTANGLE", style.name = style.name)
RCy3::setNodeShapeDefault("ELLIPSE", style.name = style.name)
RCy3::setNodeColorDefault("#000000", style.name = style.name)
RCy3::setNodeLabelMapping("GeneSymbol", style.name = style.name)
RCy3::setNodeSizeDefault(30, style.name = style.name)
#RCy3::setNodeLabelPositionDefault("S","C","c",0.00,0.00) #This function does not exist yet in this version :(
RCy3::setNodeFontSizeDefault(30, style.name = style.name)
#RCy3::setNodeWidthDefault(160,
#                          style.name = style.name)
#RCy3::deleteVisualStyle(style.name)

#Cluster with AutoAnnotate
# Run the AutoAnnotate command
aa_command <- paste("autoannotate annotate-clusterBoosted",
                    "clusterAlgorithm=MCL",
                    "labelColumn=Concatenated_Gene_Set_Descriptions",
                    "maxWords=5")
print(aa_command)
commandsGET(aa_command)
# Annotate a subnetwork
createSubnetwork(c(1:4),"__mclCluster")
commandsGET(aa_command)



#Great, then cluster analysis
#cluster the network with clustermaker2
commandsGET("string")
commandsGET("cluster mcl")
commandsPOST("cluster mcl")

#Map color by cluster
RCy3::setNodeColorMapping("__mclCluster", 
                          table.column.values = c(1,2,3),
                          mapping.type = "d",
                          colors = c('#5577FF','#FFFFFF','#FF7755'),
                          style.name = style.name)



############################################################## Cytoscape vignette part 1 overview
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Overview-of-RCy3.nb.html
cytoscapePing ()
cytoscapeVersionInfo ()
nodes <- data.frame(id=c("node 0","node 1","node 2","node 3"),
           group=c("A","A","B","B"), # categorical strings
           score=as.integer(c(20,10,15,5)), # integers
           stringsAsFactors=FALSE)
edges <- data.frame(source=c("node 0","node 0","node 0","node 2"),
           target=c("node 1","node 2","node 3","node 3"),
           interaction=c("inhibits","interacts","activates","interacts"),  # optional
           weight=c(5.1,3.0,5.2,9.9), # numeric
           stringsAsFactors=FALSE)

createNetworkFromDataFrames(nodes,edges, title="my first network", collection="DataFrame Example")
setVisualStyle('Marquee')
style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','group','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)
#Pro-tip: if you want to set NODE_WIDTH and NODE_HEIGHT independently, you also need to unlock the node dimensions with…
#lockNodeDimensions(FALSE, style.name)

g = new ('graphNEL', edgemode='directed')
g = graph::addNode ('A', g)
g = graph::addNode ('D', g)
g = graph::addNode ('C', g, edges = list('D'))
g = graph::addNode ('B', g, edges = list(c('A','D','C')))
createNetworkFromGraph (g, title='simple network', collection='GraphNEL Example')

df <- data.frame (moleculeType=c('kinase','TF','cytokine','cytokine'),
                  log2fc=c(1.8,3.0,-1.2,-2.5),
                  row.names = c('A','B','C','D'), # row.names = node names
                  stringsAsFactors = FALSE)       # important when loading strings!
loadTableData (df)

setNodeShapeDefault ('OCTAGON')
setNodeColorDefault ('#AAFF88')
setNodeSizeDefault  (60)
setNodeFontSizeDefault (30)

getNodeShapes ()   # diamond, ellipse, trapezoid, triangle, etc.
column <- 'moleculeType'
values <- c ('kinase',  'TF','cytokine')
shapes <- c ('DIAMOND', 'TRIANGLE', 'RECTANGLE')
setNodeShapeMapping (column, values, shapes)

column <- 'log2fc'
control.points <- c (-3.0, 0.0, 3.0)
colors <-  c ('#5588DD', '#FFFFFF', '#DD8855')
setNodeColorMapping (column, control.points, colors)

control.points <- c (-2.0, 0.0, 2.0)
colors <-  c ('#2255CC', '#5588DD', '#FFFFFF', '#DD8855','#CC5522')
setNodeColorMapping (column, control.points, colors)

control.points = c (-3.0, 2.0, 3.0)
sizes     = c (20, 80, 90)
setNodeSizeMapping (column, control.points, sizes)

selectNodes ('C','name')

getSelectedNodes ()

selectFirstNeighbors ()

node.names <- getSelectedNodes ()

clearSelection()
?clearSelection

saveSession('vignette_session') #.cys

full.path=paste(getwd(),'vignette_image',sep='/')
exportImage(full.path, 'PNG', zoom=200) #.png scaled by 200%
exportImage(full.path, 'PDF') #.pdf
?exportImage

help(package=RCy3)

cyrestAPI()  # CyREST API
commandsAPI()  # Commands API

commandsHelp("help")  

commandsHelp("help network")  

commandsHelp("help network select") 

browseVignettes("RCy3")
############################################################## Cytoscape vignette part 2 importing data
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Importing-data.nb.html
cytoscapePing()

# Always Start with a Network
sif <- system.file("extdata","galFiltered.sif",package="RCy3")
importNetworkFromFile(sif)

#You can import data into Cytoscape from any data.frame in R as long as it contains row.names (or an arbitrary column) that match a Node Table column in Cytoscape. In this example, we are starting with a network with yeast identifiers in the “name” column. We also have a CSV file with gene expression data values keyed by yeast identifiers here:
csv <- system.file("extdata","galExpData.csv", package="RCy3")
data <- read.csv(csv, stringsAsFactors = FALSE)

?mapTableColumn

loadTableData(data,data.key.column="name")

?loadTableData

############################################################## Cytoscape vignette part 3 Network functions and visualization
# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Network-functions-and-visualization.nb.html

cytoscapePing()

## Read a data set.
lesmis <- system.file("extdata","lesmis.txt", package="RCy3")
dataSet <- read.table(lesmis, header = FALSE, sep = "\t")

# Create a graph. Use simplify to ensure that there are no duplicated edges or self loops
gD <- igraph::simplify(igraph::graph.data.frame(dataSet, directed=FALSE))

#Verify the number of nodes (77) and edges (254):
igraph::vcount(gD)
igraph::ecount(gD)

## Common iGraph functions
# Calculate degree for all nodes
degAll <- igraph::degree(gD, v = igraph::V(gD), mode = "all")

# Calculate betweenness for all nodes
betAll <- igraph::betweenness(gD, v = igraph::V(gD), directed = FALSE) / (((igraph::vcount(gD) - 1) * (igraph::vcount(gD)-2)) / 2)
betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
rm(betAll)

# Calculate Dice similarities between all pairs of nodes
dsAll <- igraph::similarity.dice(gD, vids = igraph::V(gD), mode = "all")

## Add attributes to network
# Add new node attributes based on the calculated node properties/similarities
gD <- igraph::set.vertex.attribute(gD, "degree", index = igraph::V(gD), value = degAll)
gD <- igraph::set.vertex.attribute(gD, "betweenness", index = igraph::V(gD), value = betAll.norm)

#Check the attributes. You should see “degree” and “betweeness” now, in addition to “name”.
summary(gD)

# And now for the edge attributes…
F1 <- function(x) {data.frame(V4 = dsAll[which(igraph::V(gD)$name == as.character(x$V1)), which(igraph::V(gD)$name == as.character(x$V2))])}
dataSet.ext <- plyr::ddply(dataSet, .variables=c("V1", "V2", "V3"), function(x) data.frame(F1(x)))

gD <- igraph::set.edge.attribute(gD, "weight", index = igraph::E(gD), value = 0)
gD <- igraph::set.edge.attribute(gD, "similarity", index = igraph::E(gD), value = 0)
# Note: The order of interactions in dataSet.ext is not the same as it is in dataSet or as it is in the edge list and for that reason these values cannot be assigned directly
for (i in 1:nrow(dataSet.ext))
{
  igraph::E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
  igraph::E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)
}
rm(dataSet,dsAll, i, F1)

# Check the edge attributes. You should see “weight” and “similarity” added.
summary(gD)

## Lets check it out in Cytoscape
createNetworkFromIgraph(gD,new.title='Les Miserables')

##Let’s decide on a layout
#A list of available layouts can be accessed from R as follows:
  
getLayoutNames()
#We’ll select the “fruchterman-rheingold” layout. To see properties for the given layout, use:
  
getLayoutPropertyNames("fruchterman-rheingold") 
#We can choose any property we want and provide them as a space-delimited string:
  
layoutNetwork('fruchterman-rheingold gravity_multiplier=1 nIterations=10')
#But that is a crazy layout, so let’s try “force-directed” instead:
  
layoutNetwork('force-directed defaultSpringLength=70 defaultSpringCoefficient=0.000003')
#Next, we can visualize our data

## Next, we can visualize our data
#On nodes…
setNodeColorMapping('degree', c(min(degAll), mean(degAll), max(degAll)), c('#F5EDDD', '#F59777', '#F55333'))
lockNodeDimensions(TRUE)
setNodeSizeMapping('betweenness', c(min(betAll.norm), mean(betAll.norm), max(betAll.norm)), c(30, 60, 100))
#…and edges

setEdgeLineWidthMapping('weight', c(min(as.numeric(dataSet.ext$V3)), mean(as.numeric(dataSet.ext$V3)), max(as.numeric(dataSet.ext$V3))), c(1,3,5))
setEdgeColorMapping('weight', c(min(as.numeric(dataSet.ext$V3)), mean(as.numeric(dataSet.ext$V3)), max(as.numeric(dataSet.ext$V3))), c('#BBEE00', '#77AA00', '#558800'))
#We will define our own default color/size schema after we defined node and edge rules, due to possible issues when using rules

setBackgroundColorDefault('#D3D3D3')
setNodeBorderColorDefault('#000000')
setNodeBorderWidthDefault(3)
setNodeShapeDefault('ellipse')
setNodeFontSizeDefault(20)
setNodeLabelColorDefault('#000000')

#Voila! All done.
## Track versions for your records
cytoscapeVersionInfo()
sessionInfo()

############################################################## Cytoscape vignette part 4 Identifier mapping
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Identifier-mapping.nb.html
## Example: Species specific considerations
#When planning to import data, you need to consider the key columns you have in your network data and in your table data. It’s always recommended that you use proper identifiers as your keys (e.g., from databases like Ensembl and Uniprot-TrEMBL). Relying on conventional symbols and names is not standard and error prone.
#Let’s start with the sample network provided by Cytoscape.
#Caution: Loading a session file will discard your current session. Save first, if you have networks or data you want to keep. Use saveSession(‘path_to_file’).

openSession()  #Closes current session (without saving) and opens a sample session file
#You should now see a network with just over 300 nodes. If you look at the Node Table, you’ll see that there are proper identifiers in the name columns, like “YDL194W”. These are the Ensembl-supported IDs for Yeast.

##Perform identifier mapping
#You need to know a few things about your network in order to run this function, e.g., the species and starting (or source) identifier type. This isn’t usually a problem, but this example highlights a unique case where the Ensembl ID type for a particular species (i.e., Yeast) has a particular format (e.g., YDL194W), rather than the more typical ENSXXXG00001232 format.
#So, with this knowledge, you can run the following function:
mapped.cols <- mapTableColumn('name','Yeast','Ensembl','Entrez Gene')

#We are asking Cytoscape to look in the name column for Yeast Ensembl IDs and then provide a new columns of corresponding Entrez Gene IDs. And if you look back at the Node Table, you’ll see that new column (all the way to the right). That’s it!
#The return value is a data frame of all the mappings between Ensembl and Entrez Gene that were found for your network in case you want those details:
mapped.cols[1:3,] #first three entries

#Note: the row names of the return data frame are the node SUIDs from Cytoscape. These are handy if you want to load the mappings yourself (see last example).

## Example: From proteins to genes
#For this next example, you’ll need the STRING app to access the STRING database from within Cytoscape: * Install the STRING app from https://apps.cytoscape.org/apps/stringapp
#available in Cytoscape 3.7.0 and above
installApp('STRINGapp')  

#Now we can import protein interaction networks with a ton of annotations from the STRING database with a simple commandsGET function, like this:
string.cmd = 'string disease query disease="breast cancer" cutoff=0.9 species="Homo sapiens" limit=150'
commandsGET(string.cmd)

# for more information on string commands:
# commandsHelp('string')
# commandsHelp('string disease query'

## Perform identifier mapping
# Say we have a dataset keyed by Ensembl gene identifiers. Well, then we would want to perform this mapping:
mapped.cols <- mapTableColumn('stringdb::canonical name','Human','Uniprot-TrEMBL','Ensembl')

#Scroll all the way to the right in the Node Table and you’ll see a new column with Ensembl IDs. This example highlights a useful translation from protein to gene identifiers (or vice versa), but is also a caution to be aware of the assumptions involved when making this translation. For example, a typical gene encodes for many proteins, so you may have many-to-one mappings in your results.

## Example: Mixed identifiers
#From time to time, you’ll come across a case where the identifiers in your network are of mixed types. This is a rare scenario, but here is one approach to solving it.
#First, you’ll need the WikiPathways app to access the WikiPathways database. The pathways in WikiPathways are curated by a community of interested researchers and citizen scientists. As such, there are times where authors might use different sources of identifiers. They are valid IDs, just not all from the same source. Future versions of the WikiPathways app will provide pre-mapped columns to a single ID type. But in the meantime (and relevant to other use cases), this example highlights how to handle a source of mixed identifier types.
#Install the WikiPathways app from https://apps.cytoscape.org/apps/wikipathways
#available in Cytoscape 3.7.0 and above
installApp('WikiPathways')  

#Now we can import an Apoptosis Pathway from WikiPathways. Either from the web site (https://wikipathways.org), or from the Network Search Tool in Cytoscape GUI or from the rWikiPathways package, we could identify the pathway as WP254.
wp.cmd = 'wikipathways import-as-pathway id="WP254"'
commandsGET(wp.cmd)

# for more information on wikipathways commands:
# commandsHelp('wikipathways')
# commandsHelp('wikipathways import-as-pathway')
#Take look in the XrefId column and you’ll see a mix of identifier types. The next column over, XrefDatasource, conveniently names each type’s source. Ignoring the metabolites for this example, we just have a mix of Ensembl and Entrez Gene to deal with.

## Perform identifier mapping
#Say we want a column with only Ensembl IDs. The easiest approach is to simply overwrite all the non-Ensembl IDs, i.e., in this case, Entrez Gene IDs. Let’s collect the mappings first:
mapped.cols <- mapTableColumn('XrefId','Human','Entrez Gene','Ensembl')

#Next, we want to remove the values from the Ensembl column in our resulting mapped.cols data frame. We’ll also remove the original source columns (to avoid confusion) and rename our Ensembl column to XrefId to prepare to overwrite. Then we’ll load that into Cytosacpe:
only.mapped.cols <- mapped.cols[complete.cases(mapped.cols), 'Ensembl', drop=FALSE]
colnames(only.mapped.cols) <- 'XrefId'
loadTableData(only.mapped.cols,table.key.column = 'SUID')
#Done! See the updated XrefId column in Cytoscape with all Ensembl IDs.
#Note: you’ll want to either update the XrefDatasource* column as well or simply make a note to ignore it at this point

## More advanced cases
#This identifier mapping function is intended to handle the majority of common ID mapping problems. It has limitation, however.
?mapTableColumn
#If you need an ID mapping solution for species or ID types not covered by this tool, or if you want to connect to alternative sources of mappings, then check out the BridgeDb app: http://apps.cytoscape.org/apps/bridgedb.
#available in Cytoscape 3.7.0 and above
installApp('BridgeDb')  
#And then browse the available function with commandsHelp(‘bridgedb’)

############################################################## Cytoscape vignette part 5 Cytoscape and NDEx
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Cytoscape-and-NDEx.nb.html
## Vignette uses ndexr

## Finding networks
# The Network Data Exchange (NDEx) is a platform for network storage, retrieval and exchange. Developed in close collaboration with Cytoscape, it is a natural partner for Cytoscape network queries and management.
# You can start with NDEx by first establishing a connection.
#ndexcon <- ndex_connect()

#We will use ndexcon throughout the other ndexr calls. For example, a basic search.
#networks <- ndex_find_networks(ndexcon, "Breast Cancer")
#print(networks[,c("name","externalId","nodeCount","edgeCount")])

#That print statement provides a nifty way to browse the search results. You’ll notice that we got results that hit each of the search terms individually, thus including any pathway with “cancer” in the name. That’s perhaps a bit too broad…
#networks <- ndex_find_networks(ndexcon, "BRCA")
#print(networks[,c("name","externalId","nodeCount","edgeCount")])

# Ok. We can work with this list. Let’s use the first hit. Note: you are going to get different hits as this database changes over time, so proceed with any hit you like.
#networkId = networks$externalId[1]
#network = ndex_get_network(ndexcon, networkId)
#print(network)

## Error in RCX:::print.RCX(rcx) : object 'rcx' not found
## I dont think this vignette works
## Skipping this vignette

############################################################## Cytoscape vignette part 6 Group nodes
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Group-nodes.nb.html
cytoscapePing()
cytoscapeVersionInfo ()
## Background
#The ability to group nodes together into “metanodes” and collapse them to a single node in a graph is useful for simplifying views of a complex network.
#The example in this vignette describes application of node grouping functions to data that includes protein-protein interactions and clustered correlations of protein post-translational modifications (Grimes, et al., 2018). This vignette plots five proteins and their modifications, and uses the node grouping functions to manipulate the graph in Cytoscape.

## Example
# First we set up the node and edge data frames.
net.nodes <- c("ALK", "ALK p Y1078", "ALK p Y1096", "ALK p Y1586", "CTNND1", "CTNND1 p Y193", "CTNND1 p Y217", "CTNND1 p Y228", "CTNND1 p Y241", "CTNND1 p Y248", "CTNND1 p Y302", "CTNND1 p Y904", "CTTN", "CTTN ack K107", "CTTN ack K124", "CTTN ack K147", "CTTN ack K161", "CTTN ack K235", "CTTN ack K390", "CTTN ack K87", "CTTN p S113", "CTTN p S224", "CTTN p Y104", "CTTN p Y154", "CTTN p Y162", "CTTN p Y228", "CTTN p Y334", "CTTN p Y421", "IRS1", "IRS1 p Y632", "IRS1 p Y941", "IRS1 p Y989", "NPM1", "NPM1 ack K154", "NPM1 ack K223", "NPM1 p S214", "NPM1 p S218")
net.genes <- sapply(net.nodes,  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1])
parent <- c("", "ALK", "ALK", "ALK", "", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "", "IRS1", "IRS1", "IRS1", "", "NPM1", "NPM1", "NPM1", "NPM1")
nodeType <- c("protein", "modification", "modification", "modification", "protein", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "protein", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "modification", "protein", "modification", "modification", "modification", "protein", "modification", "modification", "modification", "modification")
netnodes.df <- data.frame(id=net.nodes, Gene.Name=net.genes, parent, nodeType, stringsAsFactors = FALSE)

# Define edge data
source.nodes <- c("ALK", "ALK", "ALK", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTNND1", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "CTTN", "IRS1", "IRS1", "IRS1", "NPM1", "NPM1", "NPM1", "NPM1", "ALK p Y1096", "CTNND1 p Y193", "CTNND1 p Y193", "CTNND1 p Y228", "CTNND1 p Y904", "CTNND1 p Y217", "CTNND1 p Y241", "CTNND1 p Y248", "ALK p Y1078", "ALK p Y1096", "ALK p Y1586", "IRS1 p Y941", "CTTN ack K147", "CTTN ack K107", "CTTN ack K235", "CTTN ack K87", "CTTN ack K147", "CTTN ack K124", "CTTN ack K147", "CTTN ack K235", "CTTN ack K161", "CTTN ack K390", "NPM1 ack K223", "NPM1 ack K154", "NPM1 ack K223", "ALK", "CTNND1", "CTNND1", "CTTN", "IRS1")
target.nodes <- c("ALK p Y1078", "ALK p Y1096", "ALK p Y1586", "CTNND1 p Y193", "CTNND1 p Y217", "CTNND1 p Y228", "CTNND1 p Y241", "CTNND1 p Y248", "CTNND1 p Y302", "CTNND1 p Y904", "CTTN ack K107", "CTTN ack K124", "CTTN ack K147", "CTTN ack K161", "CTTN ack K235", "CTTN ack K390", "CTTN ack K87", "CTTN p S113", "CTTN p S224", "CTTN p Y104", "CTTN p Y154", "CTTN p Y162", "CTTN p Y228", "CTTN p Y334", "CTTN p Y421", "IRS1 p Y632", "IRS1 p Y941", "IRS1 p Y989", "NPM1 ack K154", "NPM1 ack K223", "NPM1 p S214", "NPM1 p S218", "ALK p Y1586", "CTNND1 p Y228", "CTNND1 p Y302", "CTNND1 p Y302", "CTTN p Y154", "CTTN p Y162", "CTTN p Y162", "CTTN p Y334", "IRS1 p Y632", "IRS1 p Y989", "IRS1 p Y989", "IRS1 p Y989", "CTTN p S113", "CTTN p S224", "CTTN p S224", "CTTN p S224", "CTTN p Y104", "CTTN p Y228", "CTTN p Y228", "CTTN p Y228", "CTTN p Y421", "CTTN p Y421", "NPM1 p S214", "NPM1 p S218", "NPM1 p S218", "IRS1", "CTTN", "IRS1", "NPM1", "NPM1")
Weight <- c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0.8060606, 0.7575758, 0.7454545, 0.9393939, 0.8949096, 0.7329699, 0.7553845, 0.7866191, 0.775, 0.6969697, 0.7818182, 0.8424242, -0.7714286, -0.8385965, -0.5017544, -0.7473684, -0.5252838, -0.9428571, -0.8285714, -0.6713287, -0.5508772, -0.9428571, -0.8857143, -0.6310881, -0.8285714, 0.6123365, 2.115272, 0.002461723, 0.3354451, 0.5661711)
netedges.df <- data.frame(source=source.nodes, target=target.nodes, Weight, stringsAsFactors = FALSE)

#create network from data frames
net.suid <- createNetworkFromDataFrames(netnodes.df, netedges.df, title=paste(paste("Group Nodes Test"), 1+length(getNetworkList())), collection = "RCy3 Vignettes")

# Make sure nodes are spread out sufficiently
layoutNetwork('force-directed defaultSpringCoefficient=0.00001 defaultSpringLength=50 defaultNodeMass=5')

#Note that for convenience the data frame has defined whether a node is a protein or a modification, and also defined the parent node for each modification.
#The function selectNodes looks by default for the node SUID, which can be retrieved by getTableColumns. Alternatively, the data frame can be used to distinguish proteins and modifications.

nodedata <- getTableColumns("node")
edgedata <- getTableColumns("edge")
genes <- netnodes.df[grep("protein", netnodes.df$nodeType), "id"]

#select by gene SUIDs
geneSUIDs <- nodedata[grep("protein", nodedata$nodeType), 1]
selectNodes(geneSUIDs, preserve.current.selection = FALSE)

# or by names in the "id" column
selectNodes(c("ALK","IRS1"), by.col="id", preserve.current.selection = FALSE)

# or by names based on dataframe subsetting
modifications <- netnodes.df[grep("modification", netnodes.df$nodeType), "id"]
selectNodes(modifications, by='id', pre=FALSE)

# Now select one protein and all its modifications
deltacatnodes <- netnodes.df[grep("CTNND1", netnodes.df$Gene.Name), "id"]
selectNodes(deltacatnodes, by.col="id", preserve=FALSE)

# Let’s create a new group of the selected nodes and collapse it into one node…
createGroup("delta catenin group")
collapseGroup("delta catenin group")

# …then expand it again.
expandGroup("delta catenin group")

# For these data, we can create groups of all proteins together with their modifications. Here we name the groups by their gene names.
deleteGroup("delta catenin group")
for(i in 1:length(genes)) {
  print(genes[i])
  selectNodes(netnodes.df[grep(genes[i], netnodes.df$Gene.Name), "id"], by.col="id", preserve=FALSE)
  createGroup(genes[i])
  collapseGroup(genes[i])
}
groups.1 <- listGroups()
groups.1
# should see 5 group SUIDs reported

# These can all be expanded at once.
expandGroup(genes)

# An alternative method that might be quicker for large networks is to use the input data frame.
deleteGroup(genes)
for(i in 1:length(genes)) {
  print(genes[i])
  createGroup(genes[i], nodes=netnodes.df[grep(genes[i], netnodes.df$Gene.Name), "id"], nodes.by.col = "id")
}
collapseGroup(genes)
expandGroup(genes)

#This can be done more simply using sapply()
deleteGroup(genes)
sapply(genes, function(x) createGroup(x, nodes=netnodes.df[grep(x, netnodes.df$Gene.Name), "id"], nodes.by.col = "id"))
collapseGroup(genes)

#A groups’ information can be retrieved and independently expanded
getGroupInfo("ALK")
expandGroup("ALK")

# Get all groups' info
group.info <- list()
group.info <- lapply(listGroups()$groups, getGroupInfo)
print(group.info)

############################################################## Cytoscape vignette part 7 Custom Graphics and Labels
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Custom-Graphics.nb.html

## Open Sample
# For this tutorial, we will be using the galFiltered sample session file, which includes a yeast network and associated data.
openSession()

## Set style and node color
# First, lets change the style to a simple default and the color of nodes to grey:
setVisualStyle('default')
setNodeColorDefault('#D8D8D8')

## Custom Graphics
## Bar chart
#In this example, we will create a bar chart with the three expression values, gal1RGexp, gal4RGexp and gal80Rexp, available as attributes in the session file.
#Create the Custom Graphic:
setNodeCustomBarChart(c("gal1RGexp","gal4RGexp","gal80Rexp"))

#There are 4 types of Bar Charts and customizable parameters for colors, size, spacing and orientation.
#Position the Bar Chart just below the node. This is an optional step that we are doing here just to make room for subsequent graphics. By specifying both anchors at opposite ends, we can get a lot of space between the node and the graphic.
setNodeCustomPosition(nodeAnchor = "S", graphicAnchor = "N")

## Stripes
# Next we are going to create stripes of gradient mappings using a horizontal “heatmap”" of the same three data columns and position the heatmap right above the node. For this vignette, we need to also specify the slot number to avoid overwriting the Bar Chart:
setNodeCustomHeatMapChart(c("gal1RGexp","gal4RGexp","gal80Rexp"), slot = 2)
setNodeCustomPosition(nodeAnchor = "N", graphicAnchor = "S", slot = 2)

## Pie chart
# Finally, we will create a pie chart with two columns, Radiality and Degree, and place it to the left of the node. Here we’ll use the xOffset parameter to be even more specific about where we want to place the graphic relative to the node.
setNodeCustomPieChart(c("Radiality", "Degree"), slot = 3)
setNodeCustomPosition(nodeAnchor = "W", xOffset = -20, slot = 3)

## Enhanced Graphics
#The nodes in the network are labeled with the corresponding protein names (yeast), but there is additional text information in the Node Table that could be useful to display as labels on the nodes. We are going to use the enhancedGraphics app to create a second node label for the common yeast gene name.
#This involves a new step: Filling a new column with parameters for the enhancedGraphics App. This column is then mapped to a Custom Graphic slot and (optionally) positioned, like in the examples above.

##Install enhancedGraphics
#The enhancedGraphics app is available from the Cytoscape App Store. In Cytoscape 3.7 and above, you can install apps from R with the following function:
#installApp("enhancedGraphics") Pre-installed in Cytoscape 3.10 or later

## Define new label
#The new column values have to follow a specific syntax to be recognized by the enhancedGraphics app. Here, for example, is how you set a label based on another attribute (e.g., the column called “COMMON”), specifying its size, color, outline and background:
#"label: attribute=COMMON labelsize=10 color=red outline=false background=false""
#For more details on the enhancedGraphics format, see the manual.
#First, we define a dataframe with two columns: node names (“name”) and the new label (“my second label”):
all.nodes<-getAllNodes()
label.df<-data.frame(all.nodes, "label: attribute=COMMON labelsize=10 color=red outline=false background=false")
colnames(label.df)<-c("name","my second label")

#Next, we load this dataframe into the Node Table to create and fill a new column:
loadTableData(label.df, data.key.column = "name", table.key.column = "name")

## Map and position label
#We now have a new column, my second label, that we can use for the mapping. This mapping does not come with a custom helper function, se we are going to use two alternative functions to prepare the passthrough mapping property and then update our visual style with the new mapping:
label.map<-mapVisualProperty('node customgraphics 4','my second label','p')
updateStyleMapping('default', label.map)
## Error Error in mapVisualProperty("node customgraphics 4", "my second label",  : Could not find my second label column in defaultnode table.
## Vignette breaks down here, but its the end

############################################################## Cytoscape vignette part 8 Filtering Networks
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Filtering-Networks.nb.html
## This vignette will introduce you to some techniques for filtering a network based on node properties. You will learn to:

#Select a set of nodes based on node degree and attribute filters
#Create a subnetwork based on selected nodes
#Hide a set of nodes based on filters
#For this tutorial, we will use data from the STRING database (https://string-db.org/)
installApp('STRINGapp')

## Get network from STRING
# We are going to query the STRING Disease database for the term “breast cancer”. By default, the app pulls the top 100 human proteins associated with the disease along with edges having an evidence strength of 0.4 or greater:
string.cmd = 'string disease query disease="breast cancer"'
commandsRun(string.cmd)
string.net<-getNetworkSuid()  #grab handle on network for future reference

## Filtering by degree
## Creating a degree filter
#Every node in a network has a Degree property, which corresponds to the number of edges connecting the node to other nodes, either as a target or source. Filtering based on node degree is a useful way to remove nodes with too few (or too many) connections.
#In this example we want to exclude low degree nodes, e.g., those with only 0, 1 or 2 connections:
createDegreeFilter('degree filter', c(0,2), 'IS_NOT_BETWEEN')
# At the bottom of the Select tab, you can see how many edges/nodes where selected.

## Creating a subnetwork from a selection
# We can now create a new network, or subnetwork, from our selected set of nodes and all relevant edges:
createSubnetwork(subnetwork.name ='Breast cancer: highly connected nodes')

## Filtering by attribute
## Creating a column filter
#We could also filter the network based on high disease score. The disease score comes from STRING and indicates the strength of the association to the disease queried.
#Let’s select nodes from the original network with a disease score of greater than 4 (on a scale of 1-5):
createColumnFilter(filter.name='disease score filter', column='stringdb::disease score', 4, 'GREATER_THAN', network = string.net)
#Again, we can create a subnetwork from the selection:
createSubnetwork(subnetwork.name ='Breast cancer: high disease score')

## Combining filters
#But what if we want to combine these two filters? You could apply them sequentially as individual filters, but then you’d need to be careful about the order in which you apply the filters. Alternatively, you can create a composite filter and apply the logic all at once!
#Let’s combine the two filters “degree filter” and “disease score” to produce one filter, then apply it to the original network and create a final subnetwork:
createCompositeFilter('combined filter', c('degree filter','disease score filter'), network = string.net)
createSubnetwork(subnetwork.name ='final subnetwork')

#We can apply a layout to help with interpreting the network:
layoutNetwork('force-directed defaultSpringCoefficient=5E-6')
#This final network obviously contains fewer nodes than the original, but they are the most connected and most highly associated with the disease. If you examine the network you can see several well-known breast cancer oncogenes, for example BRCA1, TP53 and PTEN, near the center of the action.

## Hiding filtered nodes
#As a final example of the filter functions, let’s return to the orignal network once more and apply our “combined filter”. But this time let’s hide the filtered out nodes, rather than forming a selection. This demonstrates the applyFilter function and the hide parameter that is optional for all createXXXFilter functions as well.
applyFilter('combined filter', hide=TRUE, network = string.net)

############################################################## Cytoscape vignette part 9 Network Layout
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/network-layout.nb.html
cytoscapePing()
cytoscapeVersionInfo()
##Applying a Layout Algorithm
#Load the galFiltered.cys session file.
openSession()
#Use “fitContent” to see the whole network:
fitContent()
#The network layout used in this session file is the Prefuse Force Directed Layout. This layout is based on the algorithm implemented as part of the prefuse toolkit. The algorithm is very fast and with the right parameters can provide a very visually pleasing layout. The Force Directed Layout can also use a numeric edge attribute as a weight for the length of the spring.

## Layout Menu
# Available layouts can be found from:
getLayoutNames()

#Note that some Cytoscape apps may add additional layout algorithms to the Layout menu so the listing of available layouts may be influenced by which apps you have loaded.
#In addition to the built-in layout algorithms available by default, a set of yFiles layouts are available for installation. (But yFiles does not support Cytoscape Automation and cannot be controlled from RCy3.)
installApp('yFiles Layout Algorithms')
#These require a license agreement.

## Examples of Layout Algorithms
#Similarly to Prefuse Force Directed, the Edge-weighted spring-embedded layout is also based on a “force-directed” paradigm as implemented by Kamada and Kawai (1988). Network nodes are treated like physical objects that repel each other, such as electrons. The connections between nodes are treated like metal springs attached to the pair of nodes. These springs repel or attract their end points according to a force function. The layout algorithm sets the positions of the nodes in a way that minimizes the sum of forces in the network.
layoutNetwork('kamada-kawai')

#The circular algorithm produces layouts that emphasize group and tree structures within a network. It partitions the network by analyzing its connectivity structure, and arranges the partitions as separate circles. The circles themselves are arranged in a radial tree layout fashion.
layoutNetwork('circular')

#The Compound Spring Embedder (CoSE) layout is based on the traditional force-directed layout algorithm with extensions to handle multi-level nesting (compound nodes), edges between nodes of arbitrary nesting levels and varying node sizes. It is the suggested Cytoscape layout for compound graphs, although it also works very well with noncompound graphs.
layoutNetwork('cose')

##Layout Settings
#To change the settings for a particular algorithm, first you need to check the relevant set of parameters. The following code chunk will display the parameter names.

getLayoutPropertyNames('cose')

# To change the settings for a particular algorithm:
setLayoutProperties('cose', list(incremental='false', idealEdgeLength=50, springStrength=50, repulsionStrength=50, gravityStrength=50, compoundGravityStrength=50, gravityRange=50, compoundGravityRange=50, smartEdgeLengthCalc='true', smartRepulsionRangeCalc='true'))

############################################################## Cytoscape vignette part 10 Loading Networks
## https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/loading-networks.nb.html
#In Cytoscape, network data can be loaded from a variety of sources, and in several different formats. Where you get your network data depends on your biological question and analysis plan. This tutorial outlines how to load network data from several popular sources and formats.
#Public databases
#NDEx
#PSICQUIC
#STRING/STITCH
#WikiPathways
#Local and remote files
#Cytoscape apps (Biopax, KEGG and other formats)

#cytoscapePing()
#cytoscapeVersionInfo()

#Prerequisites
#installApp('stringApp')
#installApp('WikiPathways')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ndexr")
#BiocManager::install("rWikiPathways")

#install.packages("httr")

## Networks from Public Data
#Cytoscape includes a Network Search tool for easy import of public network data. In addition to core apps that are included with your Cytoscape installation (NDEx and PSICQUIC), the resources listed here will depend on which apps you have installed.
#Find out which apps you have installed:
#getInstalledApps()

## NDEx
#The NDEx Project provides an open-source framework where scientists and organizations can share, store, manipulate, and publish biological network knowledge.
#Search NDEx for “TP53 AND BARD1”:
#library(ndexr)
#ndexcon <- ndex_connect()
#networks <- ndex_find_networks(ndexcon, "TP53 AND BARD1")
#print(networks[,c("name","externalId","nodeCount","edgeCount")])
#We can work with this list. Let’s use the first hit. Note: you are going to get different hits as this database changes over time, so proceed with any hit you like.
#networkId = networks$externalId[1]
#network = ndex_get_network(ndexcon, networkId)
#print(network)
#Import the network into Cytoscape:
  
#  importNetworkFromNDEx(networkId)
#For more detailed information about working with NDEx networks in Cytoscape, see the Cytoscape-and-NDEx protocol.
## Does not work

## STRING/STITCH
#STRING is a database of known and predicted protein-protein interactions, and STITCH stored known and predicted interactions between chemicals and proteins. Data types include:
  
#  Genomic Context Predictions
#High-throughput Lab Experiments
#(Conserved) Co-Expression
#Automated Textmining
#Previous Knowledge in Databases
#Search STRING with the disease keyword “ovarian cancer”. The resulting network will load automatically.
string.cmd = 'string disease query disease="ovarian cancer"'
commandsRun(string.cmd)
#STRING networks load with a STRING-specific style, which includes 3D protein structure diagrams.

#Export the image as a png. This will save the png to your current directory.
exportImage('ovarian_cancer', 'PNG')

#STRING networks also include data as node/interaction attributes, that can be used to create a Style. Let’s save the attributes as a dataframe and take a look at the first few rows:
df <- getTableColumns()
head(df)

#The STRING app includes options to change interaction confidence level, expand the network etc.

#Before changing interaction confidence level, let’s find the number of interactions in the network:
getEdgeCount()

#Let’s increase the confidence level to 0.9, from the default 0.4:
string.cmd = 'string change confidence confidence=0.9 network=CURRENT'
commandsRun(string.cmd)

#Now let’s get an edge count after changing interaction confidence level:
getEdgeCount() 

#Again, we can export a figure:
exportImage('before_expand', 'PNG')

#For more detailed information about working with STRING networks in Cytoscape, see the stringApp protocol.

## WikiPathways
