################################## Script start
print("Starting script 6")

################################## Set seed for reproducibility
set.seed(1337)

################################## Install packages
renv::restore()

renv::install("readxl@1.4.3", prompt = FALSE)
renv::install("reticulate@1.35.0", prompt = FALSE)

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
    "readxl", # Loading excel files
    "reticulate" # Using python in R
  ),
  library,
  character.only = TRUE
)

################################## Anaconda python/glycowork environment
# A anaconda enviroment called AnacondaGlycoworks have been set up
# Python interpreter in R studio set to this environment
# Python 3.11.7 installed in the environment
# glycowork[all] v1.1.0 installed from anaconda terminal with pip inside environment 
# https://github.com/BojarLab/glycowork
#reticulate::use_condaenv("base") # Using base environment for now instead of dedicated AnacondaGlycoworks environment

use_python("C:/Users/xerhma/AppData/Local/Programs/Python/Python312/python.exe")

################################## Glycoworks setup
# Use the reticulate package to import the necessary Python modules
glycowork <- import("glycowork")

# Shorten glycoDraw function from glycowork
GlycoDraw <- glycowork$motif$draw$GlycoDraw

plot_glycan_excel <- glycowork$motif$draw$plot_glycans_excel

# Shorten the canonicalize_composition from glycowork
canonicalize_composition <- glycowork$motif$processing$canonicalize_iupac

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

########################################## Wrangle mouse gastric glycan excel file to tidy format
# The input file has the dreaded multi-row header, where the header row above implies information the the row below.
# After searching tidyxl might be possible to use, but this is unsure and the documentation is limited.
# No vignette on how to solve the specific conundrum.
# After deliberating and scratching my hair, my solution is to create a column name vector first based on the header rows
# The names need to consider downstream wrangling to get it into a tidy format
# Then, the more neatly formatted table below the header rows is imported using the column name vector
# When importing, import columns as characters since I already know the cells contain formatting that would break conversion to numeric without previous processing
# tidyr is then used to wrangle into a tidy format
# Finalized with data sanitation and quality control

## Creating column name vector
# Read header rows
Header_Rows_Raw <- read_xlsx(path = "./R_input_files/Mouse corpus MS data adjusted structures by Mattias v2.xlsx",
          col_names = FALSE,
          n_max = 3)

glimpse(Header_Rows_Raw)

# Data QC
# Identify variations of Sample_ID
Sample_IDs <- Header_Rows_Raw %>% 
  dplyr::select(-c(0:13)) %>% 
  mutate("Variable_Type" = c("MS_File_Name", "Sample_ID", "Type")) %>% 
  pivot_longer(cols = -Variable_Type, names_to = "Column_Name", values_to = "Value") %>% 
  pivot_wider(values_from = Value, names_from = Variable_Type) %>% 
  dplyr::select(Sample_ID) %>% 
  dplyr::distinct()
print(Sample_IDs)
# Sample names contain information about treatment
# Formating of treatment information is insistent

# Identify variations of Type
Variable_Types <- Header_Rows_Raw %>% 
  dplyr::select(-c(0:13)) %>% 
  mutate("Variable_Type" = c("MS_File_Name", "Sample_ID", "Type")) %>% 
  pivot_longer(cols = -Variable_Type, names_to = "Column_Name", values_to = "Value") %>% 
  pivot_wider(values_from = Value, names_from = Variable_Type) %>% 
  dplyr::select(Type) %>% 
  dplyr::distinct()
print(Variable_Types)
# No misspelled variable types


# The topmost row are glycan MS file names. There should be one unique name for each sample, lets check that
Header_Rows_Raw %>% 
  dplyr::select(-c(0:13)) %>% 
  mutate("Variable_Type" = c("MS_File_Name", "Sample_ID", "Type")) %>% 
  pivot_longer(cols = -Variable_Type, names_to = "Column_Name", values_to = "Value") %>% 
  pivot_wider(values_from = Value, names_from = Variable_Type) %>% 
  unite("MS_File_Name_Sample_ID", MS_File_Name:Sample_ID, remove = FALSE) %>% 
  dplyr::count(MS_File_Name_Sample_ID) %>% 
  dplyr::count(n)
# 36 samples and as expected they are in triplicates, the assumption checks out.

# Creating the column name vector for sample data columns
# Fixing formatting to be consistent and easier in downstream wrangling
Sample_Data_Column_Names <- Header_Rows_Raw %>% 
  dplyr::select(-c(0:13)) %>% 
  mutate("Variable_Type" = c("MS_File_Name", "Sample_ID", "Type")) %>% 
  pivot_longer(cols = -Variable_Type, names_to = "Column_Name", values_to = "Value") %>% 
  pivot_wider(values_from = Value, names_from = Variable_Type) %>% 
  dplyr::mutate(MS_File_Name = str_replace_all(MS_File_Name, "_", "-")) %>% #_ reserved as separator between variables encoded in column name
  dplyr::mutate(MS_File_Name = str_replace(MS_File_Name, "^", "MS-File-Name-")) %>% #Regex anchor
  dplyr::mutate(Sample_ID = str_replace_all(Sample_ID, "_", "-")) %>% #_ reserved as separator between variables encoded in column name
  dplyr::mutate(Sample_ID = str_replace(Sample_ID, "^", "Sample-ID-")) %>% #Regex anchor
  dplyr::mutate(Sample_ID = str_replace(Sample_ID, "H7-", "H07-")) %>% #Constant length of sample ID
  dplyr::mutate(Sample_ID = str_replace(Sample_ID, "Control", "ShamInfected")) %>% #Consistent naming
  dplyr::mutate(Sample_ID = str_replace(Sample_ID, "Non inf control", "ShamInfected")) %>% #Consistent naming
  dplyr::mutate(Sample_ID = str_replace(Sample_ID, "Infected vehicle", "HpyloriInfected")) %>% #Consistent naming
  dplyr::mutate(Sample_ID = str_replace(Sample_ID, "Infected control", "HpyloriInfected")) %>% #Consistent naming
  dplyr::mutate(Sample_ID = str_replace_all(Sample_ID, " ", "_Treatment-Group-")) %>% #Separator between sample ID and treatment ID
  dplyr::mutate(Type = str_replace(Type, "RT", "Variable-Type-RetentionTime")) %>% #Understandable variable type ID without special characters
  dplyr::mutate(Type = str_replace(Type, "Intensity", "Variable-Type-Intensity")) %>% #Understandable variable type ID without special characters
  dplyr::mutate(Type = str_replace(Type, "Rel. Abundance", "Variable-Type-RelativeAbundance")) %>% #Understandable variable type ID without special characters
  unite("Variable_Name", MS_File_Name:Type, remove = TRUE)
print("Sanity check, there should be 36x3 = 108 unique sample data names")
nrow(Sample_Data_Column_Names %>% 
  distinct()) == (36*3)

# Creating the column name vector for glycan data columns
Glycan_Column_Names <- Header_Rows_Raw %>% 
  dplyr::select(1:13) %>% 
  mutate("Variable_Type" = c("Empty_1", "Empty_2", "Variable_Name")) %>% 
  pivot_longer(cols = -Variable_Type, names_to = "Column_Name", values_to = "Value") %>% 
  pivot_wider(values_from = Value, names_from = Variable_Type) %>% 
  dplyr::select(-c("Empty_1", "Empty_2"))
print("Sanity check, there should be 13 unique glycan variable names")
nrow(Glycan_Column_Names %>% 
       distinct()) == (13)

# Bind together and create final column name vector
Input_Excel_Column_Names <- Glycan_Column_Names %>% 
  full_join(Sample_Data_Column_Names,
            by = c("Column_Name", "Variable_Name")) %>% 
  dplyr::select(Variable_Name) %>% 
  deframe()

Input_Excel_Column_Types <- Glycan_Column_Names %>% 
  full_join(Sample_Data_Column_Names,
            by = c("Column_Name", "Variable_Name")) %>% 
  mutate(Col_Types = c(c("text", #Glycan_ID
                       "numeric", #mz
                       "numeric", #Charge
                       "numeric", #MH
                       "text", #Glycan_Type
                       "text", #IncludedAdductsAndMore
                       "numeric", #Hex
                       "numeric", #HexNAc
                       "numeric", #Fucose 
                       "numeric", #Neu5Ac
                       "numeric", #Neu5Gc
                       "numeric", #Sulf
                       "text" #Structure
                       ),
                       rep("text", 36*3))) %>% 
  dplyr::select(Col_Types) %>% 
  deframe()

# Now import the excel using the created column names
Glycan_Input_Excel_raw <- read_xlsx(path = "./R_input_files/Mouse corpus MS data adjusted structures by Mattias v2.xlsx",
          skip = 3,
          col_names = Input_Excel_Column_Names,
          col_types = Input_Excel_Column_Types,
          n_max = 103)

# NA values to 0 for glycan monosackaride columns
# Pivot longer and split sample names to variables
# Data sanitation
Glycan_Input_Excel_Processed_Long <- Glycan_Input_Excel_raw %>% 
  dplyr::mutate(Hex = replace_na(Hex, 0)) %>% 
  dplyr::mutate(HexNAc = replace_na(HexNAc, 0)) %>% 
  dplyr::mutate(Fucose = replace_na(Fucose, 0)) %>% 
  dplyr::mutate(Neu5Ac = replace_na(Neu5Ac, 0)) %>% 
  dplyr::mutate(Neu5Gc = replace_na(Neu5Gc, 0)) %>% 
  dplyr::mutate(Sulf = replace_na(Sulf, 0)) %>% 
  pivot_longer(cols = starts_with("MS-File-Name"),
               names_to = "File_Sample_Treatment_Variable",
               values_to = "Data") %>% 
  tidyr::separate_wider_delim(col = "File_Sample_Treatment_Variable",
                              delim = "_",
                              names = c("MS_File_Name",
                                        "Sample_ID",
                                        "Treatment_Group",
                                        "Variable_Type")) %>% 
  dplyr::mutate("Variable_Type" = str_remove(Variable_Type, "Variable-Type-")) %>% # To facilitate pivot wider (better colnames)
  pivot_wider(names_from = Variable_Type,
              values_from = Data) %>% 
  dplyr::mutate("DataComment" = if_else(str_detect(RetentionTime, "\\*"), #RT is sometimes commented with a "*" to denote what is described below
                                        "Not verified with MS2, instead verified with MS1 and same RT as another which is verified with MS2",
                                        NA,
                                        missing = NA)) %>% 
  dplyr::mutate(RetentionTime = str_remove(RetentionTime, "\\*")) %>%  #Removing the * used to comment some RT numbers
  pivot_longer(cols = RetentionTime:RelativeAbundance,
               names_to = "Variable_Type",
               values_to = "Data") %>% 
  dplyr::mutate(Data = str_replace(Data, "ND", "0")) %>% #Non-detected glycans marked with ND in what is supposed to be a numeric variable
  dplyr::mutate(Data = str_replace(Data ,"NA", "0")) %>% #Non-detected glycans marked with NA in what is supposed to be a numeric variable
  dplyr::mutate(Data = str_replace(Data ,",", ".")) %>%  #Inconsistent decimals
  dplyr::mutate(Data = as.numeric(Data)) %>%  #Convert Data variable to numeric, should be no NA values here
  mutate(Glycan_ID = str_pad(Glycan_ID, #Padding out ID for more robust sorting and wrangling
                             width = 3,
                             side = "left",
                             pad = "0")) %>% 
  dplyr::mutate(MS_File_Name = str_remove(MS_File_Name, "MS-File-Name-")) %>%  #Better readability
  dplyr::mutate(Sample_ID = str_remove(Sample_ID, "Sample-ID-")) %>%  #Better readability
  dplyr::mutate(Treatment_Group = str_remove(Treatment_Group, "Treatment-Group-")) %>%   #Better readability
  dplyr::mutate(Sample_ID = str_replace(Sample_ID, "-", "")) #Delete special character, creates issues downstream in glycowork
print("Any NA-values in Data, indicating something in the initial data sanitation was missed?")
nrow(Glycan_Input_Excel_Processed_Long %>% 
  dplyr::filter(is.na(Data))) != 0
# Nope =)

# Now lets make it wider, so that each row is one sample/glycan combination
# Then relative abundances can be validated and if need be re-calculated
# Checking for larger differences, not just rounding errors
Glycan_Input_Excel_Processed_Wide_Validation <- Glycan_Input_Excel_Processed_Long %>% 
  pivot_wider(names_from = Variable_Type,
              values_from = Data) %>% 
  group_by(Sample_ID) %>% 
  dplyr::mutate(Total_Intensity_Per_Sample = sum(Intensity)) %>% 
  ungroup() %>% 
  mutate(Recalculated_Relative_Abundance = (Intensity / Total_Intensity_Per_Sample) * 100) %>% 
  mutate(Matching_Relative_Abundance = round(Recalculated_Relative_Abundance, digits = 9) == round(RelativeAbundance, digits = 9))
print("Any relative abundances not matching (unless through rounding error between excel and R)?")
nrow(Glycan_Input_Excel_Processed_Wide_Validation %>% 
  dplyr::filter(Matching_Relative_Abundance == FALSE)) != 0

## Export processed input data
Glycan_Input_Excel_Processed_Output <- Glycan_Input_Excel_Processed_Wide_Validation %>% 
  dplyr::rename("Glycan_Relative_Abundance_Percentage" = Recalculated_Relative_Abundance) %>% 
  dplyr::select(-c("Total_Intensity_Per_Sample", "Matching_Relative_Abundance", "RelativeAbundance")) %>% 
  mutate(Canonicalized_Structure = sapply(Structure, function (x) canonicalize_composition(x))) #Canonicalized structure according to glycoworks IUPAC condensed dialect
# Are all glycan structures unique?
nrow(Glycan_Input_Excel_Processed_Output %>% 
  dplyr::select(Glycan_ID, Structure, Canonicalized_Structure)%>% 
  distinct() %>% 
  group_by(Canonicalized_Structure) %>% 
  summarise(n = n()) %>% 
  dplyr::filter(n != 1)) == 0

write_xlsx(x = Glycan_Input_Excel_Processed_Output,
           path = "./R_intermediate_files/Mattias_Mouse_Gastric_Glycans_All_Fixed_2024_05_31.xlsx")
# Filtering down to the lmpc infection manuscript samples and exporting
Glycan_Input_Excel_Processed_Output_Infection_Manuscript_Input <- Glycan_Input_Excel_Processed_Output %>% 
  dplyr::filter(Treatment_Group == "ShamInfected" | Treatment_Group == "HpyloriInfected") %>% 
  dplyr::select(-Canonicalized_Structure)

write_xlsx(x = Glycan_Input_Excel_Processed_Output_Infection_Manuscript_Input,
           path = "./R_input_files/Glycan_Input.xlsx")

########################################## Import, process, quality control and export complete glycan data
df_glycan_canonicalized_all_data <- read_xlsx(path = "./R_input_files/Glycan_Input.xlsx") %>% 
  mutate(Canonicalized_Structure = sapply(Structure, function (x) canonicalize_composition(x))) %>%  #Canonicalized structure according to glycoworks IUPAC condensed dialect
  relocate(Glycan_ID, Structure, Canonicalized_Structure) %>%
  group_by(Glycan_ID) %>% 
  dplyr::mutate("Glycan_Mean_Relative_Abundance" = mean(Glycan_Relative_Abundance_Percentage)) %>% 
  dplyr::mutate("Glycan_Standard_Deviation_Relative_Abundance" = sd(Glycan_Relative_Abundance_Percentage)) %>% 
  dplyr::mutate("Glycan_Median_Relative_Abundance" = median(Glycan_Relative_Abundance_Percentage)) %>% 
  dplyr::mutate("Glycan_Median_Absolute_Deviation_Relative_Abundance" = mad(Glycan_Relative_Abundance_Percentage)) %>% 
  ungroup()

# QC
# Any duplicate structures?
df_glycan_canonicalized_all_data %>% 
  dplyr::select(Glycan_ID, Canonicalized_Structure) %>% 
  distinct() %>% 
  dplyr::count(Canonicalized_Structure) %>% 
  arrange(n) %>% 
  dplyr::filter(n != 1)

# Export
write_xlsx(x = df_glycan_canonicalized_all_data,
           path = "./R_output_files/tables/df_glycan_canonicalized_all_data.xlsx")
################################## Processing and exporting glycan data for python glycoworks
df_canonicalized_mouse_all <- df_glycan_canonicalized_all_data %>% 
  dplyr::select(-c("MS_File_Name", 
                   "Intensity", 
                   "RetentionTime", 
                   "DataComment",
                   "Treatment_Group")) %>% 
  pivot_wider(names_from = "Sample_ID",
              values_from = "Glycan_Relative_Abundance_Percentage")
write_tsv(x = df_canonicalized_mouse_all,
          file = "./Python_input_files/df_canonicalized_mouse_all.tsv")

df_canonicalized_mouse_minimal_glycan_col_all <- df_canonicalized_mouse_all %>% 
  dplyr::select(Canonicalized_Structure, starts_with(c("H07", "H10"))) %>% 
  pivot_longer(cols = starts_with(c("H07", "H10")),
               names_to = "id",
               values_to = "Glycan_Relative_Abundance_Percentage") %>% 
  pivot_wider(names_from = "Canonicalized_Structure",
              values_from = "Glycan_Relative_Abundance_Percentage")
write_tsv(x = df_canonicalized_mouse_minimal_glycan_col_all,
          file = "./Python_input_files/df_canonicalized_mouse_minimal_glycan_col_all.tsv")

df_canonicalized_mouse_minimal_sample_col_all <- df_canonicalized_mouse_all %>% 
  dplyr::select(Canonicalized_Structure, starts_with(c("H07", "H10"))) %>% 
  dplyr::rename("glycan" = Canonicalized_Structure)
write_tsv(x = df_canonicalized_mouse_minimal_sample_col_all,
          file = "./Python_input_files/df_canonicalized_mouse_minimal_sample_col_all.tsv")

mouse_sample_metadata_all <- df_glycan_canonicalized_all_data %>% 
  dplyr::select(Sample_ID, Treatment_Group) %>% 
  distinct() %>% 
  dplyr::rename("id" = Sample_ID) %>% 
  dplyr::rename("Group" = Treatment_Group) %>% 
  dplyr::mutate(Cohort = str_extract(id, "^..."))
write_tsv(x = mouse_sample_metadata_all,
          file = "./Python_input_files/mouse_sample_metadata_all.tsv")

mouse_infected_vehicle_sample_names <-mouse_sample_metadata_all %>% 
  dplyr::filter(Group == "HpyloriInfected") %>% 
  dplyr::select("id")
write_tsv(x = mouse_infected_vehicle_sample_names,
          file = "./Python_input_files/mouse_infected_vehicle_sample_names.tsv")

mouse_uninfected_vehicle_sample_names <-mouse_sample_metadata_all %>% 
  dplyr::filter(Group == "ShamInfected") %>% 
  dplyr::select("id")
write_tsv(x = mouse_uninfected_vehicle_sample_names,
          file = "./Python_input_files/mouse_uninfected_vehicle_sample_names.tsv")

################################## Drawing glycans
# Creating df
df_glycan_canonicalized <- df_glycan_canonicalized_all_data %>% 
  dplyr::select(Glycan_ID, 
                Canonicalized_Structure, 
                Structure, 
                Glycan_Type, 
                Glycan_Mean_Relative_Abundance, 
                Glycan_Standard_Deviation_Relative_Abundance,
                Glycan_Median_Relative_Abundance,
                Glycan_Median_Absolute_Deviation_Relative_Abundance) %>% 
  distinct()

# Also works as quality control of structure strings
# Canonicalize glycan structures
df_glycodraw <- df_glycan_canonicalized %>% 
  dplyr::select(Glycan_ID,
                Canonicalized_Structure) %>% 
  distinct() %>% 
  unite("File_Name",
        1:2,
        remove = FALSE) %>% 
  mutate(File_Name = str_replace(File_Name,
                                 "$",
                                 ".pdf")) %>% 
  mutate(File_Path = str_replace(File_Name,
                                 "^",
                                 "./R_output_files/Glycan_SNFG/"))

# Draw glycans and save files
for (row in seq_len(nrow(df_glycodraw))) {
  
  file_path <- file.path(df_glycodraw[[row, "File_Path"]])
  
  glycan_structure <- df_glycodraw[[row, "Canonicalized_Structure"]]
  
  # Call GlycoDraw function to draw the glycan structure and save it to SVG file
  GlycoDraw(glycan_structure, filepath = file_path, show_linkage = TRUE)
}

# Check if any glycans failed to draw
# Get created file names
glycan_structure_file_names <- list.files("./R_output_files/Glycan_SNFG/")

# Compare whats missing from expected file names
# Missing files indicates formating errors that canonicalize failed to fix
setdiff(str_replace_all(df_glycodraw$File_Name, "\\?", "_"), str_replace_all(glycan_structure_file_names, "\\?", "_"))

# Draw pixel glycans in excel table format to get quick overview of glycans
plot_glycan_excel(df_glycan_canonicalized,
                  folder_filepath = file.path("./R_output_files/Glycan_SNFG"),
                  glycan_col_num = as.integer(1))

writexl::write_xlsx(x = df_glycan_canonicalized,
                    path = "./R_output_files/tables/df_glycan_canonicalized.xlsx")

################################## Differential glycomics
# Performed in jupyter notebook with the next script. 
# This is done in python and not in R since I noticed issues running these parts of glycoworks with R/reticulate


########################################## SessionInfo
sessionInfo()

print("Script 6 finished")

