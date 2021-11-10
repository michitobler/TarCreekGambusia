rm(list = ls())
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(adegenet)
library(gtools)
library(readxl)
library(edgeR)

####################################################################################################################################################################################
#
#
#   1. Import data and do basic manipulations to get it ready for DAPC
#
#
####################################################################################################################################################################################

genes.raw <- read_excel("~/Downloads/05_SupplementalTableS3_GeneExpResults.xlsx", skip = 3, col_names = T)

# We wanted to restrict our analysis to only include differentially expressed genes, so we need to import a list of GeneIDs that were found to be differentially expressed in the other analyses for this manuscript.
deg.raw <- read.csv("~/Downloads/common_DEG_in_all_tissues_and_directions3.csv", header = T)

# Some genes were differentially expressed in more than one tissue, so we need to remove duplicates
deg.dedup <- deg.raw[!duplicated(deg.raw$XmacGene),]

# Subset the raw dataset to only include genes that are in the deduplicated list of differentially expressed genes
deg.exp <- genes.raw[which(genes.raw$GeneID %in% deg.dedup$XmacGene),]

# For DAPC, we only need the counts for each gene in each sample and the gene identifier, we don't need the rest of the dataframe, so we'll subset it here.
deg.subset <- deg.exp[, c(1, 5:57)]

# The output of Stringtie, which is where the genes.raw object came from, has columns representing invdividual samples, and rows representing genes, but 
# DAPC needs columns to represent the genes and rows to represent the samples, so we need to transpose the dataframe.
deg.subset.transposed <- setNames(as.data.frame(t(deg.subset[,-1])), deg.subset$GeneID)

# Separate all samples based on their site and tissue
grps <- unlist(strsplit(rownames(deg.subset.transposed), split = "_R"))[c(T,F)]
deg.subset.transposed$grp <- grps



####################################################################################################################################################################################
#
#
#   2. Running DAPC to visualize expression variation
#
#
####################################################################################################################################################################################

### DAPC on whole, undivided expression dataset to see how expression varies between tissues and among sites within each tissue
# We chose 10 PCs and 3 LDs to retain based on the plots provided by the dapc function.
dapc.all.tissues <- dapc(deg.subset.transposed[,1:ncol(deg.subset.transposed)-1], deg.subset.transposed$grp, var.contrib = TRUE) 
scatter(dapc.all.tissues, scree.pca = T)

### DAPC on dataset divided into each tissue separately
# Divide into tissues
gills <- deg.subset.transposed[which(deg.subset.transposed$grp %in% c("Tar_Gill", "LE_Gill", "Coal_Gill")),]
livers <- deg.subset.transposed[which(deg.subset.transposed$grp %in% c("Tar_Liver", "LE_Liver", "Coal_Liver")),]
brains <- deg.subset.transposed[which(deg.subset.transposed$grp %in% c("Tar_Brain", "LE_Brain", "Coal_Brain")),]

### Run DAPC on each tissue separately
# For all tissues, we retained 8 PCs and 2 DFs
dapc.gill <- dapc(gills[,1:ncol(gills)-1], gills$grp, var.contrib = TRUE)
scatter(dapc.gill, scree.pca = T)

dapc.liver <- dapc(livers[,1:ncol(livers)-1], livers$grp, var.contrib = TRUE)
scatter(dapc.liver, scree.pca = T)

dapc.brain <- dapc(brains[,1:ncol(brains)-1], brains$grp, var.contrib = TRUE)
scatter(dapc.brain, scree.pca = T)

# And here are summaries of DAPC in all three tissues
dapc.gill
dapc.liver
dapc.brain




