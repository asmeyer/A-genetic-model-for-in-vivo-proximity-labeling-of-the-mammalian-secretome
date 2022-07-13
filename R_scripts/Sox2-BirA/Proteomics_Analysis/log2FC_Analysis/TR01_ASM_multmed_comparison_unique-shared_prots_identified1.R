# 2021-06-03 
# ASM 
# TR01: BirA/Sox2Cre paper 

# usage: analyse LC-MS/MS data for BirA/Sox2Cre paper

# last updated: 2021-06-03


################################################################## begin ####################################################################

# read in multimedian data; use not normalized for serum (all or nothing condition in serum causes signal to be masked when normalized)
serum_mlt <- read.csv("input/TR01_Serum_AmandaPlex_Multimedian_notNormalized_Unfiltered_correctOrder.csv") # not normalized, still need to filter 
brain_mlt <- read.csv("input/TR01_Brain_MedianNormalized_Multimedian.csv")
kid_mlt <- read.csv("input/TR01_Kidney_mediannormalized_multimedian.csv")
liv_mlt <- read.csv("input/TR01_Liver_MedianNormalized_Multimedian.csv")




# order based on decreasing logFC
serum_mlt <- serum_mlt[order(-serum_mlt$logFC.xBirA.Cre.over.NoBirA), ]
serum_mlt <- serum_mlt[serum_mlt$numPepsUnique >= 2.0, ]
# filter on logFC and adj. p cutoffs 
serum_mlt.ft <- serum_mlt[serum_mlt$logFC.xBirA.Cre.over.NoBirA >= 1.0 & serum_mlt$adj.P.Val.xBirA.Cre.over.NoBirA <= 0.05, ] # could use mKate2 as cutoff
serum_mlt.ft <- serum_mlt[serum_mlt$logFC.xBirA.Cre.over.NoBirA > 0.90 & serum_mlt$adj.P.Val.xBirA.Cre.over.NoBirA <= 0.05, ] # could use mKate2 as cutoff

# order based on decreasing logFC
brain_mlt <- brain_mlt[order(-brain_mlt$logFC.xBirA.over.NoBirA), ]
# filter on logFC and adj. p cutoffs 
brain_mlt.ft <- brain_mlt[brain_mlt$logFC.xBirA.over.NoBirA >= 1.0 & brain_mlt$adj.P.Val.xBirA.over.NoBirA <= 0.05, ] # could use mKate2 as cutoff
brain_mlt.ft <- brain_mlt[brain_mlt$logFC.xBirA.over.NoBirA > 0.90 & brain_mlt$adj.P.Val.xBirA.over.NoBirA <= 0.05, ] # could use mKate2 as cutoff

# order based on decreasing logFC
kid_mlt <- kid_mlt[order(-kid_mlt$logFC.XBirA.over.No_BirA), ]
# filter on logFC and adj. p cutoffs 
kid_mlt.ft <- kid_mlt[kid_mlt$logFC.XBirA.over.No_BirA >= 1.0 & kid_mlt$adj.P.Val.XBirA.over.No_BirA <= 0.05, ] # could use mKate2 as cutoff
kid_mlt.ft <- kid_mlt[kid_mlt$logFC.XBirA.over.No_BirA > 2.45 & kid_mlt$adj.P.Val.XBirA.over.No_BirA <= 0.05, ] # could use mKate2 as cutoff

# order based on decreasing logFC
liv_mlt <- liv_mlt[order(-liv_mlt$logFC.xBirA.over.NoBirA), ]
# filter on logFC and adj. p cutoffs 
liv_mlt.ft <- liv_mlt[liv_mlt$logFC.xBirA.over.NoBirA >= 1.0 & liv_mlt$adj.P.Val.xBirA.over.NoBirA <= 0.05, ] # could use mKate2 as cutoff
liv_mlt.ft <- liv_mlt[liv_mlt$logFC.xBirA.over.NoBirA > 1.8 & liv_mlt$adj.P.Val.xBirA.over.NoBirA <= 0.05, ] # could use mKate2 as cutoff

all_genes.sig <- list(serum = serum_mlt.ft$id.concat, brain = brain_mlt.ft$id.concat, kidney = kid_mlt.ft$id.concat, liver = liv_mlt.ft$id.concat)
all_genes.sig[2]

allGenes <- c(serum_mlt.ft$id.concat, brain_mlt.ft$id.concat, kid_mlt.ft$id.concat, liv_mlt.ft$id.concat)
str(allGenes)

# find shared proteins between all samples. Expect to see ER resident proteins, BirA, non-cell type specific proteins
#dups <- allGenes[any(duplicated(allGenes))]
dups <- allGenes[duplicated(allGenes)]
length(dups)

allTissues <- c(brain_mlt.ft$id.concat, kid_mlt.ft$id.concat, liv_mlt.ft$id.concat)
dups.tissues <- allTissues[duplicated(allTissues)]
length(dups.tissues)
str(dups.tissues)

# check uniques. dups and uniques should all up to total 
unis <- (unique(allGenes))
length(unis)


# look at shared proteins between each tissue and sample 
ser_br <- c(serum_mlt.ft$id.concat, brain_mlt.ft$id.concat)
ser_br.dups <- ser_br[duplicated(ser_br)]
length(ser_br.dups)


ser_kid <- c(serum_mlt.ft$id.concat, kid_mlt.ft$id.concat)
ser_kid.dups <- ser_kid[duplicated(ser_kid)]
length(ser_kid.dups)


ser_liv <- c(serum_mlt.ft$id.concat, liv_mlt.ft$id.concat)
ser_liv.dups <- ser_liv[duplicated(ser_liv)]
length(ser_liv.dups)



# uniques between samples (no serum)
allTissues <- c(brain_mlt.ft$id.concat, kid_mlt.ft$id.concat, liv_mlt.ft$id.concat)
unis.tissues <- unique(allTissues)
length(unis.tissues)

# remove shared proteins observed in all tissues (dups.tissue, above)
brain.ids <- brain_mlt.ft[!grepl(paste(dups.tissues, collapse = "|"), brain_mlt.ft$id.concat), ]
length(brain.ids$id.concat)
brain.ids <- brain.ids$id.concat

kid.ids <- kid_mlt.ft[!grepl(paste(dups.tissues, collapse = "|"), kid_mlt.ft$id.concat), ]
length(kid.ids$id.concat)
kid.ids <- kid.ids$id.concat

liv.ids <- liv_mlt.ft[!grepl(paste(dups.tissues, collapse = "|"), liv_mlt.ft$id.concat), ]
length(liv.ids$id.concat)
liv.ids <- liv.ids$id.concat




### old (doesnt get rid of ER prots)
liv.ids <- c(liv_mlt.ft$id.concat, unis.tissues)
liv.uni <- liv.ids[duplicated(liv.ids)]
length(liv.uni)





# get tissue specific hits 
brain.hits <- brain_mlt.ft[grep(paste(brain.ids, collapse = "|"), brain_mlt.ft$id.concat), ]
kid.hits <- kid_mlt.ft[grep(paste(kid.ids, collapse = "|"), kid_mlt.ft$id.concat), ]
liv.hits <- liv_mlt.ft[grep(paste(liv.ids, collapse = "|"), liv_mlt.ft$id.concat), ]




# PC & NC lists from Claire
PC_CY <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/Mouse_PC_uniprot_secreted.csv")
NC_CY <- read.csv("Box Sync/McMahon Lab/Dissertation/TR01/BirA-Sox2Cre/LC-MSMS/Serum/Data/Mouse_NC_combined.csv")

# PC & NC lists from Rui 
PC_RY <- read.csv("PC-NC_lists/Mouse_PC_uniprot_secreted_RY.csv")
NC_RY <- read.csv("PC-NC_lists/Mouse_NC_combined_updated_RY.csv")

############################################################################# PC/NC ID ###################################################################################
ifelse(serum_A$accession_number == PC_RY$Entry & serum_A$geneSymbol != NC_RY$Gene.Symbol, "PC", "other")
ifelse(serum_A$accession_number != PC_RY$Entry & serum_A$geneSymbol == NC_RY$Gene.Symbol, "NC", "other")

serum_A$PCorNC <- serum_A

brain.hits2 <- ifelse(brain.hits$accession_number == PC_RY$Entry & brain.hits$geneSymbol != NC_RY$Gene.Symbol, "PC", "other")

brain.hits1 <- brain.hits
brain.hits2 <- brain.hits
brain.hits1$PCorNC <- x

PC <- PC_RY$Entry
NC <- NC_RY$Gene.Symbol

brain.hits1$PCorNC <- ifelse(brain.hits[grepl(paste(PC, collapse = "|"), brain.hits$accession_number), ], "PC", "other")

# grepl not specific enough 
brain.hits1$PC <- ifelse(Reduce(`|`, lapply(PC, function(PC) grepl(PC, brain.hits$accession_number))), "PC", "other")
brain.hits1$NC <- ifelse(Reduce(`|`, lapply(NC, function(NC) grepl(NC, brain.hits$geneSymbol))), "NC", "other")


word<-df2$words_to_use
pattern<-paste0('.*',word,'.*', collapse ='|')

ifelse(Reduce(`|`, lapply(word, function(pat) grepl(pat, df$column1))), "Y", "N")

kid.hits 
liv.hits 



















