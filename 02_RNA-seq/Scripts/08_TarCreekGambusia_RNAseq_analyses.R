rm(list=ls())

# All analyses performed with R 4.0

# Load libraries 
library(limma)
library(DESeq2)
library(edgeR)
library(splines)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(gridExtra)
library(ggpubr)
library(stats)


################################################################################################
#
#           1. Bring in and filter raw counts data and annotation file
#
################################################################################################

# Load gene matrix file generated from stringTie, Xmac annotation file, and Sample names
setwd("~/Desktop/Projects/2018_Tar_Creek_Gambusia_Ionomics_RNAseq_MolEvol/02_Raw_Data/02_RNA-seq_tar_coal_little_elm/")
xmac.annotations <- read.csv("XmacGeneAnnotations.csv", header=T)

data <- read.csv("bwa_subsamples_and_unmodified_gene_matrix.csv", skip = 1, header=T, row.names = 1)
dim(data)

sample_names <- read.csv("sample_names_filtered_reads.csv", header = F)
colnames(sample_names) <- c("ID", "Site", "Tissue", "SampID")

# Trim data: Each row (i.e. transcript) must have greater than 2 counts per million (cpm) and at least 15 of the 53 samples must have counts data
keep <- rowSums(cpm(data)>2) >= 15



# Keep data that meet the filtering criteria
filtered_reads <- data[keep,]
dim(filtered_reads)

# We need to create a lookup table to convert the XmacGene IDs into SwissProt annotations for each gene in the transcriptome
lookup <- data.frame("Annotation" = as.character(xmac.annotations$SwissprotAnnotationID),
                     "XmacGene" = as.character(xmac.annotations$XmacGene))
filtered.annot <- filtered_reads
filtered.annot$Annotation <- lookup$Annotation[match(rownames(filtered.annot), lookup$XmacGene)]


# Now, subset the filtered data into gills, livers, and brains, so analyses can be run on each tissue separately
gill.samps <- sample_names[which(sample_names$Tissue == "Gill"),]
gills = filtered_reads[,which(colnames(filtered_reads) %in% c("LE_Gill_R4","Tar_Gill_R31","LE_Gill_R13", "Tar_Gill_R37","Coal_Gill_R61","Coal_Gill_R70","Tar_Gill_R40","Coal_Gill_R76","LE_Gill_R1","Coal_Gill_R64","Coal_Gill_R67","Tar_Gill_R43","Tar_Gill_R49","LE_Gill_R7","LE_Gill_R10","Coal_Gill_R73"))]
dim(gills)

liver.samps <-sample_names[which(sample_names$Tissue == "Liver"),]
livers = filtered_reads[,which(colnames(filtered_reads) %in% c("Coal_Liver_R68", "LE_Liver_R5", "Coal_Liver_R71", "Tar_Liver_R44", "Tar_Liver_R38", "LE_Liver_R11", "Coal_Liver_R65", "Coal_Liver_R74", "LE_Liver_R14", "LE_Liver_R2", "LE_Liver_R17", "Tar_Liver_R35", "Tar_Liver_R32", "Coal_Liver_R62", "LE_Liver_R8", "Tar_Liver_R47", "Tar_Liver_R41", "Coal_Liver_R77"))]
dim(livers)

brain.samps <- sample_names[which(sample_names$Tissue == "Brain"),]
brain_ids <- brain.samps$ID
brains = filtered_reads[,which(colnames(filtered_reads) %in% brain_ids)]
dim(brains)



#####################################################################################################
#
#       2. Identify candidate gene lists for each tissue separately
#           
#####################################################################################################

# Group samples by population
sites.gills = factor(c("LE","Tar","LE", "Tar","Coal","Coal","Tar","Coal","LE","Coal","Coal","Tar","Tar","LE","LE","Coal"))
sites.livers = liver.samps$Site
sites.brains = brain.samps$Site


# Generate dataframe for each tissue with sample names and site
gill.samples = data.frame(cbind(colnames(gills)), as.character(sites.gills))
colnames(gill.samples) = c("samples","site")

liver.samples = data.frame(cbind(colnames(livers)), as.character(sites.livers))
colnames(liver.samples) = c("samples", "site")

brain.samples = data.frame(cbind(colnames(brains)), as.character(sites.brains))
colnames(brain.samples) = c("samples", "site")


# Create a DGEList object to hold the dataset
samples.gill = DGEList(counts = gills, group = gill.samples$site)
samples.liver = DGEList(counts = livers, group = liver.samples$site)
samples.brain = DGEList(counts = brains, group = brain.samples$site)


# Calculate normalized factors based on raw library sizes
samples.gill = calcNormFactors(samples.gill)
samples.liver = calcNormFactors(samples.liver)
samples.brain = calcNormFactors(samples.brain)


# Create a design matrix (because the datasets are already split in three by tissue, the only predictor is site)
sample.design.gill <- model.matrix(~0+sites.gills)
sample.design.liver <- model.matrix(~0 + sites.livers)
sample.design.brain <- model.matrix(~0 + sites.brains)


# Add column names to the design matrix based on sample names
colnames(sample.design.gill) <- levels(factor(gill.samples$site))
colnames(sample.design.liver) <- levels(factor(liver.samples$site))
colnames(sample.design.brain) <- levels(factor(brain.samples$site))


# Estimate common dispersion and tagwise dispersion
samples.gill <- estimateDisp(samples.gill, sample.design.gill)
samples.liver <- estimateDisp(samples.liver, sample.design.liver)
samples.brain <- estimateDisp(samples.brain, sample.design.brain)


# Given tagwise dispersion and a design matrix, glmFIT fits the negative binomial GLM for each tag. 
# Will then produce a DGEGLM object with new components. 
sample.fit.gill <- glmFit(samples.gill,sample.design.gill)
sample.fit.liver <- glmFit(samples.liver, sample.design.liver)
sample.fit.brain <- glmFit(samples.brain, sample.design.brain)


# Construct contrast matrix of comparisons (compare populations in each tissue)

gill.contrast.tc.le <- makeContrasts(
  Gaffinis = Tar - LE,
  levels = sample.design.gill
)

gill.contrast.tc.cc <- makeContrasts(
  Gaffinis = Tar - Coal,
  levels = sample.design.gill
)



liver.contrast.tc.le <- makeContrasts(
  Gaffinis = Tar - LE, 
  levels = sample.design.liver
)

liver.contrast.tc.cc <- makeContrasts(
  Gaffinis = Tar - Coal,
  levels = sample.design.liver
)


brain.contrast.tc.le <- makeContrasts(
  Gaffinis = Tar - LE, 
  levels = sample.design.brain
)

brain.contrast.tc.cc <- makeContrasts(
  Gaffinis = Tar - Coal,
  levels = sample.design.brain
)



# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons
gill_lrt_tc_le = glmLRT(sample.fit.gill, contrast = gill.contrast.tc.le[,"Gaffinis"])
gill_lrt_tc_cc = glmLRT(sample.fit.gill, contrast = gill.contrast.tc.cc[,"Gaffinis"])

liver_lrt_tc_le = glmLRT(sample.fit.liver, contrast = liver.contrast.tc.le[,"Gaffinis"])
liver_lrt_tc_cc = glmLRT(sample.fit.liver, contrast = liver.contrast.tc.cc[,"Gaffinis"])

brain_lrt_tc_le = glmLRT(sample.fit.brain, contrast = brain.contrast.tc.le[,"Gaffinis"])
brain_lrt_tc_cc = glmLRT(sample.fit.brain, contrast = brain.contrast.tc.cc[,"Gaffinis"])


# Summary of differential expression that was up, down or not significant in each comparison
summary(decideTestsDGE(gill_lrt_tc_le, adjust.method = "BH", p.value = 0.05))
summary(decideTestsDGE(gill_lrt_tc_cc, adjust.method = "BH", p.value = 0.05))

summary(decideTestsDGE(liver_lrt_tc_le, adjust.method = "BH", p.value = 0.05))
summary(decideTestsDGE(liver_lrt_tc_cc, adjust.method = "BH", p.value = 0.05))

summary(decideTestsDGE(brain_lrt_tc_le, adjust.method = "BH", p.value = 0.05))
summary(decideTestsDGE(brain_lrt_tc_cc, adjust.method = "BH", p.value = 0.05))

# Extract the results of the LRT for each gene to see which ones were differentially expressed. This can be found in the 
# "table" element from the lrt objects. We'll also need to add a column for FDR based on the p-value given
gill.le.table <- gill_lrt_tc_le$table
gill.le.table$FDR <- p.adjust(gill_lrt_tc_le$table$PValue, method = "BH")
gill.le.table$GeneID <- rownames(gill.le.table)
colnames(gill.le.table) <- paste("Gill_LE", colnames(gill.le.table), sep = "_")

gill.cc.table <- gill_lrt_tc_cc$table
gill.cc.table$FDR <- p.adjust(gill_lrt_tc_cc$table$PValue, method = "BH")
gill.cc.table$GeneID <- rownames(gill.cc.table)
colnames(gill.cc.table) <- paste("Gill_CC", colnames(gill.cc.table), sep = "_")

liver.le.table <- liver_lrt_tc_le$table
liver.le.table$FDR <- p.adjust(liver_lrt_tc_le$table$PValue, method = "BH")
liver.le.table$GeneID <- rownames(liver.le.table)
colnames(liver.le.table) <- paste("Liver_LE", colnames(liver.le.table), sep = "_")

liver.cc.table <- liver_lrt_tc_cc$table
liver.cc.table$FDR <- p.adjust(liver_lrt_tc_cc$table$PValue, method = "BH")
liver.cc.table$GeneID <- rownames(liver.cc.table)
colnames(liver.cc.table) <- paste("Liver_CC", colnames(liver.cc.table), sep = "_")

brain.le.table <- brain_lrt_tc_le$table
brain.le.table$FDR <- p.adjust(brain_lrt_tc_le$table$PValue, method = "BH")
brain.le.table$GeneID <- rownames(brain.le.table)
colnames(brain.le.table) <- paste("Brain_LE", colnames(brain.le.table), sep = "_")

brain.cc.table <- brain_lrt_tc_cc$table
brain.cc.table$FDR <- p.adjust(brain_lrt_tc_cc$table$PValue, method = "BH")
brain.cc.table$GeneID <- rownames(brain.cc.table)
colnames(brain.cc.table) <- paste("Brain_CC", colnames(brain.cc.table), sep = "_")

# Combine all the results tables above together into the same df
all.tables <- as.data.frame(cbind(gill.le.table, gill.cc.table, liver.le.table, liver.cc.table, brain.le.table, brain.cc.table))

# Check to make sure all the rows are lined up correctly by checking if the GeneID columns all line up 
sum(which(all.tables$Gill_LE_GeneID != all.tables$Gill_LE_GeneID)) # should be zero if they match up
sum(which(all.tables$Gill_LE_GeneID != all.tables$Gill_CC_GeneID)) 
sum(which(all.tables$Gill_LE_GeneID != all.tables$Liver_LE_GeneID))
sum(which(all.tables$Gill_LE_GeneID != all.tables$Liver_CC_GeneID))
sum(which(all.tables$Gill_LE_GeneID != all.tables$Brain_LE_GeneID))
sum(which(all.tables$Gill_LE_GeneID != all.tables$Brain_CC_GeneID))

###### Now let's add this table to the filtered_reads df ######
composite.filtered <- as.data.frame(cbind(filtered_reads, all.tables))

# This composite df is from filtered data, so there are only 18,693 rows instead of 27,266 rows, so we need to match up 
# based on the GeneID column. This "final.data" df should contain the raw counts for all 18,693 genes that were kept after 
# filtering low expressed genes at the beginning of this script, the LRT results for each of those genes, and the annotation.
final.data <- as.data.frame(composite.filtered[match(rownames(data), composite.filtered$Gill_LE_GeneID),])


# Extracting top tags for little elm and dividing into signficantly up- or downregulated
top_gills_tc_le = topTags(gill_lrt_tc_le, n = (1695 + 1383))
top_gills_tc_le.up = rownames(top_gills_tc_le[top_gills_tc_le$table$logFC > 0,])
top_gills_tc_le.down = rownames(top_gills_tc_le[top_gills_tc_le$table$logFC < 0,])

top_livers_tc_le = topTags(liver_lrt_tc_le, n = (497 + 193))
top_livers_tc_le.up = rownames(top_livers_tc_le[top_livers_tc_le$table$logFC > 0,])
top_livers_tc_le.down = rownames(top_livers_tc_le[top_livers_tc_le$table$logFC < 0,])

top_brains_tc_le = topTags(brain_lrt_tc_le, n = (14 + 22))
top_brains_tc_le.up = rownames(top_brains_tc_le[top_brains_tc_le$table$logFC > 0,])
top_brains_tc_le.down = rownames(top_brains_tc_le[top_brains_tc_le$table$logFC < 0,])


# Extracting top tags for Coal Creek and dividing into signficantly up- or downregulated
top_gills_tc_cc = topTags(gill_lrt_tc_cc, n = (2214 + 2414))
top_gills_tc_cc.up = rownames(top_gills_tc_cc[top_gills_tc_cc$table$logFC > 0,])
top_gills_tc_cc.down = rownames(top_gills_tc_cc[top_gills_tc_cc$table$logFC < 0,])

top_livers_tc_cc = topTags(liver_lrt_tc_cc, n = (885 + 442))
top_livers_tc_cc.up = rownames(top_livers_tc_cc[top_livers_tc_cc$table$logFC > 0,])
top_livers_tc_cc.down = rownames(top_livers_tc_cc[top_livers_tc_cc$table$logFC < 0,])

top_brains_tc_cc = topTags(brain_lrt_tc_cc, n = (250 + 236))
top_brains_tc_cc.up = rownames(top_brains_tc_cc[top_brains_tc_cc$table$logFC > 0,])
top_brains_tc_cc.down = rownames(top_brains_tc_cc[top_brains_tc_cc$table$logFC < 0,])



#####################################################################################################
#
#
#     3. Make Venn diagrams and find intersections between both uncontaminated sites
#
#
#####################################################################################################

# Upregulated candidate genes in the gill
universal.gill.up <- unique(c(top_gills_tc_le.up, top_gills_tc_cc.up))
GillA <- universal.gill.up %in% top_gills_tc_le.up
GillB <- universal.gill.up %in% top_gills_tc_cc.up
universal.gill.up.input.df <- data.frame("LE" = GillA, "Coal" = GillB)
head(universal.gill.up.input.df)
venn.gill.up.universal <- vennCounts(universal.gill.up.input.df)
vennDiagram(venn.gill.up.universal, main = "Upregulated Genes in Gill", 
            names = c("Little Elm", "Coal"), cex = 1.5, cex.main = 1.5, 
            circle.col = c("turquoise","steelblue"), mar = c(0, 0, 1, 0))
common.gill.upreg <- universal.gill.up[which(universal.gill.up.input.df["LE"] == T & universal.gill.up.input.df["Coal"] == T )]
common.gill.upreg <- data.frame(common.gill.upreg)
names(common.gill.upreg) <- "XmacGene"
gill.candidates.up <- merge(xmac.annotations[, c("XmacGene", "XmacGeneID", "Evalue", "HumanGene", "SwissprotAnnotationID")], common.gill.upreg, by = "XmacGene") # these are your candidate upregulated genes in the gill


# Downregulated candidate genes in the gill
universal.gill.down <- unique(c(top_gills_tc_le.down, top_gills_tc_cc.down))
Gill.down.A <- universal.gill.down %in% top_gills_tc_le.down
Gill.down.B <- universal.gill.down %in% top_gills_tc_cc.down
universal.gill.down.input.df <- data.frame("LE" = Gill.down.A, "Coal" = Gill.down.B)
head(universal.gill.down.input.df)
venn.gill.down.universal <- vennCounts(universal.gill.down.input.df)
vennDiagram(venn.gill.down.universal, main = "Downregulated Genes in Gill", 
            names = c("Little Elm", "Coal"), cex = 1.5, cex.main = 1.5, 
            circle.col = c("turquoise","steelblue"), mar = c(0, 0, 1, 0))
common.gill.downreg <- universal.gill.down[which(universal.gill.down.input.df["LE"] == T & universal.gill.down.input.df["Coal"] == T )]
common.gill.downreg <- data.frame(common.gill.downreg)
names(common.gill.downreg) <- "XmacGene"
gill.candidates.down <- merge(xmac.annotations[, c("XmacGene", "XmacGeneID", "Evalue", "HumanGene", "SwissprotAnnotationID")], common.gill.downreg, by = "XmacGene") # these are your candidate downregulated genes in the gill


# Upregulated candidate genes in the liver
universal.liver.up <- unique(c(top_livers_tc_le.up, top_livers_tc_cc.up))
liver.up.A <- universal.liver.up %in% top_livers_tc_le.up
liver.up.B <- universal.liver.up %in% top_livers_tc_cc.up
universal.liver.up.input.df <- data.frame("LE" = liver.up.A, "Coal" = liver.up.B)
head(universal.liver.up.input.df)
venn.liver.up.universal <- vennCounts(universal.liver.up.input.df)
vennDiagram(venn.liver.up.universal, main = "Upregulated Genes in liver", 
            names = c("Little Elm", "Coal"), cex = 1.5, cex.main = 1.5, 
            circle.col = c("turquoise","steelblue"), mar = c(0, 0, 1, 0))
common.liver.upreg <- universal.liver.up[which(universal.liver.up.input.df["LE"] == T & universal.liver.up.input.df["Coal"] == T )]
common.liver.upreg <- data.frame(common.liver.upreg)
names(common.liver.upreg) <- "XmacGene"
liver.candidates.up <- merge(xmac.annotations[, c("XmacGene", "XmacGeneID", "Evalue", "HumanGene", "SwissprotAnnotationID")], common.liver.upreg, by = "XmacGene") # these are your candidate upregulated genes in the liver


# Downregulated candidate genes in the liver
universal.liver.down <- unique(c(top_livers_tc_le.down, top_livers_tc_cc.down))
liver.down.A <- universal.liver.down %in% top_livers_tc_le.down
liver.down.B <- universal.liver.down %in% top_livers_tc_cc.down
universal.liver.down.input.df <- data.frame("LE" = liver.down.A, "Coal" = liver.down.B)
head(universal.liver.down.input.df)
venn.liver.down.universal <- vennCounts(universal.liver.down.input.df)
vennDiagram(venn.liver.down.universal, main = "downregulated Genes in liver", 
            names = c("Little Elm", "Coal"), cex = 1.5, cex.main = 1.5, 
            circle.col = c("turquoise","steelblue"), mar = c(0, 0, 1, 0))
common.liver.downreg <- universal.liver.down[which(universal.liver.down.input.df["LE"] == T & universal.liver.down.input.df["Coal"] == T )]
common.liver.downreg <- data.frame(common.liver.downreg)
names(common.liver.downreg) <- "XmacGene"
liver.candidates.down <- merge(xmac.annotations[, c("XmacGene", "XmacGeneID", "Evalue", "HumanGene", "SwissprotAnnotationID")], common.liver.downreg, by = "XmacGene") # these are your candidate downregulated genes in the liver


# Upregulated candidate genes in the brain
universal.brain.up <- unique(c(top_brains_tc_le.up, top_brains_tc_cc.up))
brain.up.A <- universal.brain.up %in% top_brains_tc_le.up
brain.up.B <- universal.brain.up %in% top_brains_tc_cc.up
universal.brain.up.input.df <- data.frame("LE" = brain.up.A, "Coal" = brain.up.B)
head(universal.brain.up.input.df)
venn.brain.up.universal <- vennCounts(universal.brain.up.input.df)
vennDiagram(venn.brain.up.universal, main = "Upregulated Genes in brain", 
            names = c("Little Elm", "Coal"), cex = 1.5, cex.main = 1.5, 
            circle.col = c("turquoise","steelblue"), mar = c(0, 0, 1, 0))
common.brain.upreg <- universal.brain.up[which(universal.brain.up.input.df["LE"] == T & universal.brain.up.input.df["Coal"] == T )]
common.brain.upreg <- data.frame(common.brain.upreg)
names(common.brain.upreg) <- "XmacGene"
brain.candidates.up <- merge(xmac.annotations[, c("XmacGene", "XmacGeneID", "Evalue", "HumanGene", "SwissprotAnnotationID")], common.brain.upreg, by = "XmacGene") # these are your candidate upregulated genes in the brain


# Downregulated candidate genes in the brain
universal.brain.down <- unique(c(top_brains_tc_le.down, top_brains_tc_cc.down))
brain.down.A <- universal.brain.down %in% top_brains_tc_le.down
brain.down.B <- universal.brain.down %in% top_brains_tc_cc.down
universal.brain.down.input.df <- data.frame("LE" = brain.down.A, "Coal" = brain.down.B)
head(universal.brain.down.input.df)
venn.brain.down.universal <- vennCounts(universal.brain.down.input.df)
vennDiagram(venn.brain.down.universal, main = "downregulated Genes in brain", 
            names = c("Little Elm", "Coal"), cex = 1.5, cex.main = 1.5, 
            circle.col = c("turquoise","steelblue"), mar = c(0, 0, 1, 0))
common.brain.downreg <- universal.brain.down[which(universal.brain.down.input.df["LE"] == T & universal.brain.down.input.df["Coal"] == T )]
common.brain.downreg <- data.frame(common.brain.downreg)
names(common.brain.downreg) <- "XmacGene"
brain.candidates.down <- merge(xmac.annotations[, c("XmacGene", "XmacGeneID", "Evalue", "HumanGene", "SwissprotAnnotationID")], common.brain.downreg, by = "XmacGene") # these are your candidate downregulated genes in the brain




#####################################################################################################
#
#
#     4. Filter GO results and identify GO results related to metals
#
#
#####################################################################################################

# Filter GO results
go <- read.csv("~/Desktop/Projects/2018_Tar_Creek_Gambusia_Ionomics_RNAseq_MolEvol/02_Raw_Data/02_RNA-seq_tar_coal_little_elm/GO_enrichment_results_unfiltered.csv")
filtered.go <- go[which(go$FDR.q.value <= 0.05 & go$b >= 5 & go$Enrichment >= 2),]

# Split filtered GO results into tissue-specific datasets
gill.annotations <- filtered.go[filtered.go$tissue == "gill",]
gill.annot.up <- gill.annotations[gill.annotations$direction_of_diff_expression == "upregulated",]
gill.annot.down <- gill.annotations[gill.annotations$direction_of_diff_expression == "downregulated",]

liver.annotations <- filtered.go[filtered.go$tissue == "liver",]
liver.annot.up <- liver.annotations[liver.annotations$direction_of_diff_expression == "upregulated",]
liver.annot.down <- liver.annotations[liver.annotations$direction_of_diff_expression == "downregulated",]

brain.annotations <- filtered.go[filtered.go$tissue == "brain",]
brain.annot.up <- brain.annotations[brain.annotations$direction_of_diff_expression == "upregulated",]
brain.annot.down <- brain.annotations[brain.annotations$direction_of_diff_expression == "downregulated",]

# Import the list of GO terms that came up in a search of AMIGO using the term "metal"
metal.terms <- read.csv("~/Desktop/Projects/2018_Tar_Creek_Gambusia_Ionomics_RNAseq_MolEvol/02_Raw_Data/02_RNA-seq_tar_coal_little_elm/metal_related_go_terms.csv", header = T)
terms <- metal.terms$Term
ids <- metal.terms$GOID

# Find which of the "metal" terms are in the dataframe of GO terms in the DEG
filtered.go[which(filtered.go$Description %in% terms),]
filtered.go[which(filtered.go$GO_term %in% ids),]


##################################################################################################
#
#
#     5. Visualizing shared and unique expression with stacked barplot
#
#
#################################################################################################

expression_stacked <- read.csv("~/Desktop/Projects/2018_Tar_Creek_Gambusia_Ionomics_RNAseq_MolEvol/02_Raw_Data/02_RNA-seq_tar_coal_little_elm/intermediate_steps/expression_stacked_barplot.csv")
#aspect_ratio = 16/9
#pdf("~/Downloads/stacked_expr_barplot.pdf", family = "Open Sans", height = 7, width = 7* aspect_ratio)
ggplot(expression_stacked) + 
  geom_col(aes(y = num_genes, x = tissue, fill = site), width = .5) + 
  theme_classic() + 
  theme(text = element_text(size = 14, family = "Open Sans")) +
  labs(x = "Tissue") + 
  labs(y = "Number of Genes") + 
  geom_hline(yintercept = 0, color = "black", size = .5) + 
  scale_fill_manual(name = "Site", values = c("steelblue", "turquoise", "#FF9933")) + 
  scale_x_discrete(labels=c("agill" = "Gill", "bliver" = "Liver", "cbrain" = "Brain")) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-1)) + 
  theme(legend.position = "none") + 
  ggtitle("Number of Genes Differentially Expressed By Tissue and Site")
#dev.off()



##################################################################################################
#
#
#     6. Calculate deltaK from the multiple runs of NGSadmix to find most likely number of genetic clusters
#
#
#################################################################################################

# Read in data on likelihood estimates from each run of NGSadmix at each level of K
likes <- read.table("~/Downloads/best_likes.txt")
colnames(likes) <- likes[1,]
likes <- likes[-1,]
options(digits = 16)

# DeltaK is calculated as in Evanno et al, 2005, by calculating the mean likelihood and the SD for each value of K, and then 
# dividing mean/SD
likes$likelihood <- as.numeric(likes$likelihood)
likes$abs.like <- abs(likes$likelihood)

# Calculate the average for each value of K
grouped.likes <- aggregate(likes$abs.like, by = list(likes$K), FUN = "mean")
names(grouped.likes) <- c("K", "Avg.likelihood")

# Calculate the SD for each value of K
sds <- aggregate(likes$abs.like, by = list(likes$K), FUN = "sd")
names(sds) <- c("K", "sd")
grouped.likes$sd <- sds$sd

# Calculate deltaK
grouped.likes$delta.K <- grouped.likes$Avg.likelihood / grouped.likes$sd

# The best value of K (most likely number of genetic clusters) is the largest value of deltaK
best.K <- grouped.likes$K[which(grouped.likes$delta.K == max(grouped.likes$delta.K))]
best.K

# Plot of deltaK values (on a log scale because deltaK for k=1 is HUGE)
#pdf("~/Downloads/deltaK.pdf", height = 3.31, width = 5.89, useDingbats = F)
ggplot(data = grouped.likes, aes(x = K, y = log10(delta.K))) + 
  theme_classic() + 
  geom_point(size = 5, color = "#767171") 
#dev.off()


