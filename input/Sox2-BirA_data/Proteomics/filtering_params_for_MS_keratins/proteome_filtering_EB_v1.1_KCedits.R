##########
# Proteome filtering script v1.1 - Erik Bergstrom, 02/2020
# 
# This was written to supercede the old summaryppc_filter scripts that Team Namrata uses,
# with the goal of being easier to read (and therefore, edit/fix for different use cases).
#
# UPDATE 1.1 fixes
# - This update fixed directory specific filtering issues (left-side vs right-side of the report)
# - UP filtering now uses "unique_peptides" column instead of "numPepsUnique"
# - FQ filtering now uses the left-most columns in the summary,
#   ***this will throw a warning message now, but it can be ignored***
# - Instructions are now updated to reflect these changes
#
# Update 3/25/2020 KC
# - Edited keratin filtering to include uppercase and lowercase keratins
# - Only need to change species when transitioning between human and mouse samples
##########

# install.packages("tidyverse")
# install.packages("magrittr")
library(tidyverse)
library(magrittr)

# Your directory, with forward slashes! If you copy/paste from file window, you will have to change the slashes.
setwd("X:/Kiki/Luo_PPI/Fractions/20200319_SingleShotandFxnwithFT/R_Processing_1p1_KCedits") 

# Set these parameters before running!
input_filename <- "proteinProteinCentricColumnsExport.SGT.1.ssv"
output_filename <- "proteinProteinCentricColumnsExport.SGT.1_FQ_UP2_RC2.ssv"
species_name <- "MOUSE"

#
#
protein_report <- read_delim(input_filename, delim=";") 
keratins <- read_csv("human_keratins_as_downloaded_from_HGNC_06242016.csv") # Put keratin file in working directory

#get list of all keratins, uppercase and lowercase
keratinsCaps = keratins$Approved #uppercase keratin genes
keratinsLowers = substr(keratins$Approved,2,length(keratins$Approved)) 
keratinsLowers = tolower(keratinsLowers)
keratinsLowers = paste0('K', keratinsLowers) #lowercase keratin genes with capital 'K'

keratinGenes = append(keratinsLowers, keratinsCaps) #one list of all keratin genes, upper and lower case, works for mouse and human
print(keratinGenes)

`%ni%` = negate(`%in%`)

##### INSTRUCTIONS #####

# NOTES ON HOW TO USE AND EDIT:
# - After you have checked the directory, filename, species name, and keratin file are all correct and in place,
#   check the FQ filter to make sure your TMT channels are correct (depends on 10, 11, 16).
#   (Note: in a future version, TMT channels will just be a parameter to enter up top.)
# - Then, select everything and click Run (Ctrl+A, Ctrl+Enter).
# - Once you have saved this script to your own directory, you can knock out any step you want.
#   To do this, comment it out (either put a '#' at the start of the line, or highlight and press Ctrl+Shift+C).
#   For example, if you want to remove the keratin step, comment out that chunk and nothing else should change besides the captions.
# 
# DROPPING A CHANNEL
# - To edit "fully quantified" filter, change the strings in "vars()".
#     vars(ends_with("TMT_126_total"):ends_with("TMT_131_total")) 
#         --> Selects and filters the entire range between the two named columns
#     vars(ends_with("TMT_126_total"):ends_with("TMT_131_total"), -ends_with("TMT_129N_total")) 
#         --> Subtracts selected channel from range
# - To drop channels in "ratio counts" filter, add "-contains()" modifiers to vars().
#     vars(contains("numRatios")) --> Searches for every column that contains string "numRatios" in the name
#     vars(contains("numRatios"), -contains("log2_TMT_numRatios_129N")) --> Excludes columns that have string "log2_TMT_numRatios_129N" in their names
# 
# ADDING YOUR OWN FILTERS
# - In the event you want something unique, just copy and paste a chunk into the "pipeline", and edit with your conditions.
#   NOTE: If you copy/paste a chunk onto the end of the pipeline, make sure you remove the last "%>%" operator.
# 
#     If your filter is for one column, use: "filter(LOGICAL EXPRESSION)".
#     If your filter is for multiple columns, use: "filter_at(vars(YOUR COLUMNS),all_vars(EXPRESSION TO BE APPLIED))".
#     If your filter needs a text search through the columns, use: filter_at(vars(contains("YOUR STRING")),all_vars(EXPRESSION))
# 

########################

if ("geneSymbol" %in% names(protein_report)) {
  gene_symbol_col <- sym("geneSymbol")
  cat("geneSymbol detected as gene symbol column", "\n")
} else {
  gene_symbol_col <- sym("gene_symbol")
  cat("gene_symbol detected as gene symbol column", "\n")
}

filtered_PR <- protein_report %T>%
  {cat("Number of unfiltered proteins: ", nrow(.),"\n",  sep ="")} %>%
  
  # Filters for selected species
  filter(species == species_name) %T>% 
  {cat("Number of ", species_name, " proteins: ", nrow(.),"\n",  sep = "")} %>% 
  
  # Filters out keratins from imported keratin file
  filter(!!gene_symbol_col %ni% keratinGenes) %T>% 
  {cat("Number of ", species_name, " proteins, keratins removed: ", nrow(.),"\n",  sep = "")} %>%
  
  # Filters for fully quantified
  # TMT10 = TMT_126_total : TMT_131_total
  # TMT11 = TMT_126C_total : TMT_131C_total
  # TMT16 = TMT_126C_total : TMT_134N_total
  filter_at(vars(ends_with("TMT_126_total"):ends_with("TMT_131_total")), all_vars(!is.na(.))) %T>% 
  {cat("Number of fully quantified ", species_name, " proteins: ", nrow(.),"\n",  sep = "")} %>%
  
  # Filters for >1 unique peptides
  filter_at(vars(contains("unique_peptides")), all_vars(.>1)) %T>%
  {cat("Number of FQ ", species_name, " proteins with at least 2 unique peptides: ", nrow(.),"\n",  sep = "")} %>%
  
  # Filters for >1 ratio count
  #if a channel was dropped (like in TMT8), you must also filter here! See above instructions
  filter_at(vars(contains("numRatios")), all_vars(.>1)) %T>%   
  {cat("Number of FQ ", species_name, " proteins with at least 2 UP and at least 2 ratio counts: ", nrow(.),"\n",  sep = "")}

# Writes your filtered protein report to an SSV file, ready to go for Spectrum Mill process report.
write_delim(filtered_PR, path = output_filename, delim=";", na="")
cut_X1 <- readLines(output_filename)
cut_X1 <- gsub("X1","",cut_X1) # This step cuts out a column header that was added automatically on import.
writeLines(cut_X1,output_filename)
