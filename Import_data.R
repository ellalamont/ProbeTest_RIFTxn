
# 7/8/25

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
# library(dendextend) # May need this for looking at pheatmap clustering
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(edgeR) # for cpm
library(sva) # For ComBat_seq batch correction

# DuffyTools
library(devtools)
# install_github("robertdouglasmorrison/DuffyTools")
library(DuffyTools)
# install_github("robertdouglasmorrison/DuffyNGS")
# BiocManager::install("robertdouglasmorrison/DuffyTools")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase")


cbPalette_1 <- c("#999999", "#E69F00") # Gold and Grey
cbPalette_1.5 <- c("#E69F00", "#999999") # Gold and Grey
cbPalette_2 <- c( "#0072B2", "#999999") # Blue and Grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 <-  c("#bfbfbf", "#56B4E9")
cbPalette3 <-  c("#bfbfbf", "#E69F00")
cbPalette4 <- c("#56B4E9", "#009E73", "#F0E442")
c25 <- c(
  "dodgerblue2", "#E31A1C", "green4",
  "#6A3D9A","#FF7F00","black", "gold1",
  "skyblue2", "#FB9A99","palegreen2","#CAB2D6",
  "#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown"
)
c12 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "palegreen2", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4")
c16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black","gold1", "#FB9A99", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4")


# Stop scientific notation
options(scipen = 999) 
# options(scipen = 0) # To revert back to default


###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############
# Current RIF Txn run then the others for comparison from ProbeTest5


# This has been edited to include more metadata!
RIF_txn_pipeSummary <- read.csv("Data/TBAITRun1_RIF_THP1Spiked/RIF_txn_Pipeline.Summary.Details.csv") 
ProbeTest5_pipeSummary <- read.csv("Data/ProbeTest5/ProbeTest5_Pipeline.Summary.Details_moreTrim.csv")

# Just get the samples I want 
ProbeTest5_pipeSummary_subset <- ProbeTest5_pipeSummary %>% filter(SampleID %in% c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9", "THP1_1e6_1a_S28", "THP1_1e6_2b_S31", "THP1_1e6_3a_S32"))

# Merge the pipeSummaries
All_pipeSummary <- merge(RIF_txn_pipeSummary, ProbeTest5_pipeSummary_subset, all = T)

All_pipeSummary$X <- NULL

# Change the NA in Drug to Untreated
All_pipeSummary <- All_pipeSummary %>% mutate(Drug = replace(Drug, is.na(Drug), "Untreated"))
# All_pipeSummary$Drug[is.na(All_pipeSummary$Drug)] <- "Untreated"

All_pipeSummary$Run[is.na(All_pipeSummary$Run)] <- "RIF_txn"

All_pipeSummary$Drug <- as.character(All_pipeSummary$Drug)
ordered_Drug <- c("Untreated", "RIF")
All_pipeSummary$Drug <- factor(All_pipeSummary$Drug, levels = ordered_Drug)

All_pipeSummary$SampleID <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# All_pipeSummary <- All_pipeSummary %>% mutate(Sputum_Number = str_extract(SampleID, "S_[0-9]+"))



###########################################################
############ IMPORT AND PROCESS ALL TPM VALUES ############
# NOT scaled


RIF_txn_tpm <- read.csv("Data/TBAITRun1_RIF_THP1Spiked/RIF_txn_Mtb.Expression.Gene.Data.TPM.csv")

ProbeTest5_tpm <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.TPM_moreTrim.csv") 
ProbeTest5_tpm_subset <- ProbeTest5_tpm %>% select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9, THP1_1e6_1a_S28, THP1_1e6_2b_S31, THP1_1e6_3a_S32)

# Merge 
All_tpm <- merge(ProbeTest5_tpm_subset, RIF_txn_tpm, all = T)

# Adjust the names so they are slightly shorter
names(All_tpm) <- gsub(x = names(All_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

rownames(All_tpm) <- All_tpm[,1] # add the rownames

# Need to make sure there is a Gene column (gets lots)
All_tpm <- All_tpm %>% 
  rename(Gene = X) 


###########################################################
################### IMPORT BOB's DE DATA ##################

# Don't know why these are only working with the full pathname....
`Broth_RIF.ComparedTo.THP1_RIF` <- read.delim("Data/JOINED_DEG/MTb.MetaResults.RaBroth_RIF.vs.THP1Spiked_RIF/RIF_THP1_1e6.MTb.Meta.JOINED.txt")
`Broth_Untreat.ComparedTo.Broth_RIF` <- read.delim("Data/JOINED_DEG/MTb.MetaResults.RaBroth_Untreated.vs.RaBroth_RIF/RIF_broth.MTb.Meta.JOINED.txt")
`Broth_Untreat.ComparedTo.THP1_RIF` <- read.delim("Data/JOINED_DEG/MTb.MetaResults.RaBroth_Untreated.vs.THP1Spiked_RIF/RIF_THP1_1e6.MTb.Meta.JOINED.txt")
`THP1_Untreat.ComparedTo.Broth_RIF` <- read.delim("/Users/elamont/Documents/RProjects/Sputum/ProbeTest_RIFTxn/Data/JOINED_DEG/MTb.MetaResults.Spiked_Untreated.vs.Broth_RIF/RIF_broth.MTb.Meta.JOINED.txt")
`THP1_Untreat.ComparedTo.THP1_RIF` <- read.delim("/Users/elamont/Documents/RProjects/Sputum/ProbeTest_RIFTxn/Data/JOINED_DEG/MTb.MetaResults.THP1Spiked_Untreated.vs.THP1Spiked_RIF/RIF_THP1_1e6.MTb.Meta.JOINED.txt")

###########################################################
################ MAKE A LIST OF ALL DFs ###################
list_dfs <- list(`Broth_RIF.ComparedTo.THP1_RIF`,
                 `Broth_Untreat.ComparedTo.Broth_RIF`, 
                 `Broth_Untreat.ComparedTo.THP1_RIF`,
                 `THP1_Untreat.ComparedTo.Broth_RIF`, 
                 `THP1_Untreat.ComparedTo.THP1_RIF`)

# Make a list of the names
df_names <- c("Broth_RIF.ComparedTo.THP1_RIF",
              "Broth_Untreat.ComparedTo.Broth_RIF", 
              "Broth_Untreat.ComparedTo.THP1_RIF",
              "THP1_Untreat.ComparedTo.Broth_RIF", 
              "THP1_Untreat.ComparedTo.THP1_RIF")

# Give the df list the correct df names
names(list_dfs) <- df_names


###########################################################
############### ADD COLUMNS OF DE VALUES ##################

# Make a new list to hold dataframes with extra columns
list_dfs_2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs)) {
  
  current_df <- list_dfs[[i]]
  current_df_name <- df_names[i]
  
  # Make the column pointing out which ones are differentially expressed
  current_df$DE <- ifelse(current_df$LOG2FOLD < -1 & current_df$AVG_PVALUE < 0.05, "significant down",
                          ifelse(current_df$LOG2FOLD > 1 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE <- factor(current_df$DE, levels = ordered_DE)
  
  # Make the column with DE gene names for plotting on graph
  current_df$DE_labels <- ifelse(current_df$DE != "not significant", current_df$GENE_NAME, NA)
  
  list_dfs_2[[current_df_name]] <- current_df
  
}

