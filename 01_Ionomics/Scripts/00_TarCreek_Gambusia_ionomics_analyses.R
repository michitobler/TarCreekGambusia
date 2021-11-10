#R 4.0.0
#John L Coffin
#Created 4/30/2020
#Last edited 8/30/2021
#Community-wide ionomics comparisons in Tar Creek vs. Coal Creek in Ottawa County, OK

###set working directory###
rm(list=ls())
setwd("/Users/johncoffin/Desktop/Projects/2018_Tar_Creek_Gambusia_Ionomics_RNAseq_MolEvol/02_Raw_Data/01_ionomics_tar_coal/")

library(car)
library(tidyr)
library(phytools)
library(ggplot2)
library(dplyr)
library(ape)
library(geiger)
library(maps)
library(devtools)
library(plot3D)
library(plotly)
library(convevol)
library(sjstats)
library(mosaic)
library(geometry)
library(ggpubr)
library(gridExtra)
library(stringr)
library(extrafont)
library(showtext)
library(WGCNA)
library(heplots)
library(MASS)
library(matlib)
library(groupedstats)
library(gridExtra)

#################################################################################################
#
#
#   1. Data input and manipulation
#
#
#################################################################################################

########## 1.1. Import data (data is NOT log-transformed here, but is mass-transformed)
data <- read.csv("ionomics_tc-cc_May2018_AllSpecies.csv")[,1:38]


########## 1.2. Subset data to only whole specimens from each species found in both sites
whole_specs <- data[data$Tissue == 'WS',]
rownames(whole_specs) <- c(1:length(rownames(whole_specs)))


########## 1.3. Remove NAs and other spurious data
whole_specs_na <- na_if(whole_specs, "#VALUE!")
WS_minusNAs <- whole_specs_na[ , colSums(is.na(whole_specs_na)) == 0] #removes columns that have NA's in them
species.list <- c("Gambusia_affinis", "Fundulus_notatus", "Pimephales_notatus", "Lepomis_cyanellus", "Lepomis_gulosus", "Lepomis_macrochirus", "Lepomis_megalotis")
WS_minusNAs <- WS_minusNAs[WS_minusNAs$Species %in% species.list,]

# Double check which species are still in the dataset
unique(WS_minusNAs$Species)


########## 1.4. Subset into each species individually
gaff <- WS_minusNAs[WS_minusNAs$Species == "Gambusia_affinis",]
fnot <- WS_minusNAs[WS_minusNAs$Species == "Fundulus_notatus",]
pnot <- WS_minusNAs[WS_minusNAs$Species == "Pimephales_notatus",]
lcya <- WS_minusNAs[WS_minusNAs$Species == "Lepomis_cyanellus",]
lgul <- WS_minusNAs[WS_minusNAs$Species == "Lepomis_gulosus",]
lmac <- WS_minusNAs[WS_minusNAs$Species == "Lepomis_macrochirus",]
lmeg <- WS_minusNAs[WS_minusNAs$Species == "Lepomis_megalotis",]


########## 1.5. Calculate PCscores for each species subset individually
gaff.pca <- prcomp(gaff[,11:34], scale. = T, center = T) # scale. = T means we're using a correlation matrix for the PCA, not a covariance matrix
summary(gaff.pca)
gaff.scores <- cbind(gaff[,c(5, 7)], gaff.pca$x)

fnot.pca <- prcomp(fnot[,11:34], scale. = T, center = T)
summary(fnot.pca)
fnot.scores <- cbind(fnot[,c(5, 7)], fnot.pca$x)

pnot.pca <- prcomp(pnot[,11:34], scale. = T, center = T)
summary(pnot.pca)
pnot.scores <- cbind(pnot[,c(5, 7)], pnot.pca$x)

lcya.pca <- prcomp(lcya[,11:34], scale. = T, center = T)
summary(lcya.pca)
lcya.scores <- cbind(lcya[,c(5, 7)], lcya.pca$x)

lgul.pca <- prcomp(lgul[,11:34], scale. = T, center = T)
summary(lgul.pca)
lgul.scores <- cbind(lgul[,c(5, 7)], lgul.pca$x)

lmac.pca <- prcomp(lmac[,11:34], scale. = T, center = T)
summary(lmac.pca)
lmac.scores <- cbind(lmac[,c(5, 7)], lmac.pca$x)

lmeg.pca <- prcomp(lmeg[,11:34], scale. = T, center = T)
summary(lmeg.pca)
lmeg.scores <- cbind(lmeg[,c(5, 7)], lmeg.pca$x)




############################################################################################
#
#
#   2. Discriminant Function Analyses for Each Species 
#
#
############################################################################################

### For G. affinis
summary(gaff.pca) # For each species, you'll keep the number of PC axes with an eigenvalue greater than 1, which you can find using the summary command
gaff.lda <- lda(Site ~ PC1 + PC2 + PC3 + PC4 + PC5, data = gaff.scores, prior = c(1,1)/2, CV = T)
gaff.table <- table(gaff.scores$Site, gaff.lda$class)
diag(prop.table(gaff.table, 1)) # in each site, what percentage did the lda correctly assign?
sum(diag(prop.table(gaff.table))) # across both sites, what is the percent correctly assigned?

### For F. notatus
summary(fnot.pca)
fnot.lda <- lda(Site ~ ., data = fnot.scores[,c(1, 3:6)], prior = c(1,1)/2, CV = T)
fnot.table <- table(fnot.scores$Site, fnot.lda$class)
diag(prop.table(fnot.table, 1))
sum(diag(prop.table(fnot.table)))

### For P. notatus
summary(pnot.pca)
pnot.lda <- lda(Site ~ ., data = pnot.scores[,c(1, 3:6)], prior = c(1,1)/2, CV = T)
pnot.table <- table(pnot.scores$Site, pnot.lda$class)
diag(prop.table(pnot.table, 1))
sum(diag(prop.table(pnot.table)))

### For L. cyanellus
summary(lcya.pca)
lcya.lda <- lda(Site ~ ., data = lcya.scores[,c(1, 3:7)], prior = c(1,1)/2, CV = T)
lcya.table <- table(lcya.scores$Site, lcya.lda$class)
diag(prop.table(lcya.table, 1))
sum(diag(prop.table(lcya.table)))

### For L. gulosus
summary(lgul.pca)
lgul.lda <- lda(Site ~ ., data = lgul.scores[,c(1, 3:6)], prior = c(1,1)/2, CV = T)
lgul.table <- table(lgul.scores$Site, lgul.lda$class)
diag(prop.table(lgul.table, 1))
sum(diag(prop.table(lgul.table)))

### For L. macrochirus
summary(lmac.pca)
lmac.lda <- lda(Site ~ ., data = lmac.scores[,c(1, 3:6)], prior = c(1,1)/2, CV = T)
lmac.table <- table(lmac.scores$Site, lmac.lda$class)
diag(prop.table(lmac.table, 1))
sum(diag(prop.table(lmac.table)))

### For L. megalotis
summary(lmeg.pca)
lmeg.lda <- lda(Site ~ ., data = lmeg.scores[,c(1, 3:7)], prior = c(1,1)/2, CV = T)
lmeg.table <- table(lmeg.scores$Site, lmeg.lda$class)
diag(prop.table(lmeg.table, 1))
sum(diag(prop.table(lmeg.table)))



############################################################################################
#
#
#   3. Analyzing ionomic differentiation across all species between sites with PCA
#
#
############################################################################################

########## 3.1. Run PCA on raw input data for whole specimens

# PCA run with correlation matrix (scale. = T), which basically calculates z-scores of the input data and then uses the z-scores to run PCA
# In the dataframe, the element concentrations are in columns 11:34
pca <- prcomp(WS_minusNAs[,11:34], center = T, scale. = T)
print(pca)
summary(pca) # first 6 PC axes have eigenvalues greater than 1, so those axes will be retained
plot(pca, type = "l")


########## 3.2. Plot the PC scores (pca$x) and use them as input for future analyses
pca_scores <- as.data.frame(pca$x)
pca_scores$Habitat <- factor(WS_minusNAs$Site)
pca_scores$Spec <- WS_minusNAs$Species
pca_scores$SL <- WS_minusNAs$SL.mm

# Calculate loadings (correlation of linearly transformed variables with original data) as the eigenvector * sqrt(eigenvalue)
eigenvectors <- as.data.frame(pca$rotation)
eigenvalues <- pca$sdev^2

loads <- data.frame(matrix(nrow = 24, ncol = 24))
for (i in 1:length(colnames(eigenvectors))){
  loads[,i] <- eigenvectors[,i] * sqrt(eigenvalues[i])
}

colnames(loads) <- colnames(pca$rotation)
rownames(loads) <- rownames(pca$rotation)


# See which elements load the highest on each axis
pc1 <- c(loads$PC1)
names(pc1) <- rownames(loads)
sort(pc1, decreasing = T)

pc2 <- c(loads$PC2)
names(pc2) <- rownames(loads)
sort(pc2, decreasing = T)

pc3 <- c(loads$PC3)
names(pc3) <- rownames(loads)
sort(pc3, decreasing = T)

# Based on the output of the sorting above, you can see which elements load the most positively and negatively along each axis
high_load_elems_pc1 <- c("Tl", "Bi", "Be", "Zn", "Cr", "S", "Co", "Ca")
high_load_elems_pc2 <- c("Ba", "Sr", "Mn", "Ni", "B", "Si", "Fe", "Cd")
high_load_elems <- c(high_load_elems_pc1, high_load_elems_pc2)
important_loads <- loads[which(rownames(loads) %in% high_load_elems),]

# Now plot the PC scores, using the important_loads to show how those elements load.
#pdf("~/Desktop/Projects/2018_Tar_Creek_Gambusia_Ionomics_RNAseq_MolEvol/04_Manuscript/tar_coal_ionomics_pca_with_Gaff_highlighted.pdf", height = 7, width = 7* 16/9, useDingbats = F)
ggplot(data = pca_scores, aes(x = PC1, y = PC2, color = Habitat)) + 
  geom_point(aes(size = 3)) + 
  geom_point(data = pca_scores[which(pca_scores$Spec == "Gambusia_affinis" & pca_scores$Habitat == "TC"),], fill = "#CB6E37", color = "black", size = 7, pch = 22, stroke = 2) + 
  geom_point(data = pca_scores[which(pca_scores$Spec == "Gambusia_affinis" & pca_scores$Habitat == "CC"),], fill = "steelblue", color = "black", size = 7, pch = 22, stroke = 2) + 
  geom_segment(data = important_loads, 
               aes(x = 0, y = 0, xend = 5*(PC1), yend = 5*(PC2)), 
               arrow = arrow(length = unit(1/2, "picas")), 
               color = "black", 
               size = 1) + 
  theme_classic() + 
  stat_ellipse(size = 2) + 
  annotate("text", 
           x = 5*important_loads$PC1, 
           y = 5*important_loads$PC2, 
           label = rownames(important_loads), 
           size = 6) + 
  scale_color_manual(values = c("steelblue", "#CB6E37"), name = "Site") + 
  ggtitle("PCA of Element Concentration in Coal and Tar Creek Communities") + 
  theme(axis.text.x = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(text = element_text(size = 14)) + 
  theme(axis.line.x = element_line(color="black", size = 0.5)) + 
  theme(axis.text.y = element_text(hjust = 1))
#dev.off()


########## 3.3. Calculate each population pair's means, and just plot those in the first 3 PC axes
poppairs <- aggregate(pca_scores[,1:5], by = list(pca_scores$Spec, pca_scores$Habitat), FUN = "mean")
poppairs$name <- c("F. notatus", "G. affinis", "L. cyanellus", "L. gulosus", "L. macrochirus", "L. megalotis", "P. notatus", "F. notatus", "G. affinis", "L. cyanellus", "L. gulosus", "L. macrochirus", "L. megalotis", "P. notatus")
colnames(poppairs)[2] <- "habitat"

ggplot(data = poppairs, aes(x = PC1, y = PC2, color = habitat, group = name)) + 
  geom_point(aes(size = 3)) + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37")) + 
  theme(legend.position = "none") + 
  geom_text(aes(label = name), hjust = 0.5, vjust = 1.6) + 
  geom_line(color = "black")
ggplot(data = poppairs, aes(x = PC1, y = PC3, color = habitat, group = name)) + 
  geom_point(aes(size = 3)) + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37")) + 
  theme(legend.position = "none") + 
  geom_text(aes(label = name), hjust = 0.5, vjust = 1.6) + 
  geom_line(color = "black")
ggplot(data = poppairs, aes(x = PC2, y = PC3, color = habitat, group = name)) + 
  geom_point(aes(size = 3)) + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37")) + 
  theme(legend.position = "none") + 
  geom_text(aes(label = name), hjust = 0.5, vjust = 1.6) + 
  geom_line(color = "black")



########## 3.4. Take a look at the distribution of each of the first 6 PCs separately
# Plot PC1
ggplot(data = pca_scores, aes(x = PC1, color = Habitat)) + 
  geom_density() + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37"), name = "Site")
# t-test to determine if PC scores from Tar Creek (TC) are significantly different than those from Coal Creek (CC)
t.test(pca_scores$PC1[which(pca_scores$Habitat == "TC")], 
       pca_scores$PC1[which(pca_scores$Habitat == "CC")], 
       var.equal = T)

# PC2
ggplot(data = pca_scores, aes(x = PC2, color = Habitat)) + 
  geom_density() + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37"), name = "Habitat")
t.test(pca_scores$PC2[which(pca_scores$Habitat == "TC")], 
       pca_scores$PC2[which(pca_scores$Habitat == "CC")], 
       var.equal = T)

# PC3
ggplot(data = pca_scores, aes(x = PC3, color = Habitat)) + 
  geom_density() + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37"), name = "Habitat")
t.test(pca_scores$PC3[which(pca_scores$Habitat == "TC")], 
       pca_scores$PC3[which(pca_scores$Habitat == "CC")], 
       var.equal = T)

# PC4
ggplot(data = pca_scores, aes(x = PC4, color = Habitat)) + 
  geom_density() + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37"), name = "Habitat")
t.test(pca_scores$PC4[which(pca_scores$Habitat == "TC")], 
       pca_scores$PC4[which(pca_scores$Habitat == "CC")], 
       var.equal = T)

# PC5
ggplot(data = pca_scores, aes(x = PC5, color = Habitat)) + 
  geom_density() + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37"), name = "Habitat")
t.test(pca_scores$PC5[which(pca_scores$Habitat == "TC")], 
       pca_scores$PC5[which(pca_scores$Habitat == "CC")], 
       var.equal = T)

# PC6
ggplot(data = pca_scores, aes(x = PC6, color = Habitat)) + 
  geom_density() + 
  theme_classic() + 
  scale_color_manual(values = c("steelblue", "#CB6E37"), name = "Habitat")
t.test(pca_scores$PC6[which(pca_scores$Habitat == "TC")], 
       pca_scores$PC6[which(pca_scores$Habitat == "CC")], 
       var.equal = T)

# The first three axes were significantly different between sites, based on the t-tests, so plot those three axes together
plot_ly(pca_scores, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Habitat, colors = c("steelblue", "#CB6E37"))


########## 3.5. Run MANOVA to test for difference in PC scores across species between sites
class(pca_scores$Habitat)
class(pca_scores$Spec)
class(pca_scores$SL)
pca_scores$Spec <- as.factor(pca_scores$Spec)
pca_scores$SL <- as.numeric(pca_scores$SL)

# There are significant differences in average standard length between species, so we need to calculate a measure of 
# relative standard length, by scaling each specimen's SL to its own species' mean and SD. So we will calculate Z scores 
# of each species SL
pca_scores$RelSL <- ave(pca_scores$SL, pca_scores$Spec, FUN = scale)

# The full model has the first 6 PCs as response vars and site, species, Relative SL, and all interactions as predictors
full_mod <- lm(cbind(PC1, PC2, PC3, PC4, PC5, PC6) ~ 
                 Habitat + 
                 Spec + 
                 RelSL + 
                 RelSL*Habitat + 
                 RelSL*Spec + 
                 Habitat*Spec, 
               data = pca_scores)

full_fit <- Manova(full_mod, test.statistic = "Wilks", type = "III")
full_sum <- summary(full_fit, test = "Wilks")
full_fit
etasq(full_fit, partial = T, anova = T, test = "Wilks")

# Updated model removes RelSL, which was not significant in the full model, as well as all interactions containing RelSL
mod1 <- lm(cbind(PC1, PC2, PC3, PC4, PC5, PC6) ~
             Habitat + 
             Spec + 
             Habitat*Spec, 
           data = pca_scores)

fit1 <- Manova(mod1, test.statistic = "Wilks", type = "III")
fit1
fit1_sum <- summary(fit1, test = "Wilks")
etasq(fit1, partial = T, anova = T, test = "Wilks")

########## 4.1. Calculate divergence vector scores by extracting the SSCP matrix for the habitat term in the MANOVA. 
# This tests for convergent ionomic differentiation between habitats across all species

sscp <- as.matrix(as.data.frame(fit1_sum$multivariate.tests$Habitat[1]))

########## 4.2. Run PCA on SSCP matrix
sscp.pca <- prcomp(sscp, center = T, scale. = F) # scale. = F uses a covariance matrix
print(sscp.pca) 
plot(sscp.pca, type = "l")
summary(sscp.pca) # Should have almost all the variance explained on the first PC axis


########## 4.3. Multiply the PC scores by the divergence vector (1st PC)
#divergence vector is just the first PC after PCA on the SSCP matrix
divergence.vector <- sscp.pca$rotation[,1]
dvs <- as.matrix(pca_scores[,1:6]) %*% divergence.vector
pca_scores$DVS <- dvs


########## 4.4. Plot the DVS in each species between the two sites
# add a new column to pca_scores that removes the underscore from the Species column. Also sort it (in the levels = c() part) to whatever order you want the boxplots to be arranged on the y-axis.
pca_scores$name <- factor(gsub("_", " ", pca_scores$Spec), levels = c("Lepomis megalotis", "Lepomis macrochirus", "Lepomis gulosus", "Lepomis cyanellus", "Pimephales notatus", "Fundulus notatus", "Gambusia affinis"), ordered = T)
asp <- 16/9
#pdf("~/Desktop/Projects/2018_Tar_Creek_Gambusia_Ionomics_RNAseq_MolEvol/04_Manuscript/figures/Fig_3_tc-cc-community-dvs.pdf", height = 7, width = 7* asp, useDingbats = F)
#par(mar = c(10,4,4,4))
ggplot(pca_scores, aes(x = name, y = DVS, fill = Habitat)) + 
  geom_boxplot(color = "black", width = .5) + 
  theme_classic() + 
  theme(axis.title.y = element_blank()) + 
  ggtitle("Divergence Vector Scores in Each Population Pair") + 
  scale_fill_manual(values = c("steelblue", "#CB6E37")) + 
  ylab("DVS") + 
  theme(axis.text.x = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(text = element_text(size = 14)) + 
  theme(axis.line.x = element_line(color="black", size = 0.5)) + 
  theme(axis.text.y = element_text(hjust = 1, face = "italic")) + 
  theme(axis.title.x = element_text(vjust = -2)) + 
  theme(legend.position = "bottom") + 
  theme(plot.margin = margin(t = 1, r = 5, b = 1, l = 1, unit = "cm")) + 
  coord_flip()
#dev.off()

# Now use t-tests to compare mean DVS between sites in each species
gaff.dvs <- data.frame(TC = pca_scores[pca_scores$Habitat == "TC" & pca_scores$Spec == "Gambusia_affinis",]$DVS, 
                       CC = pca_scores[pca_scores$Habitat == "CC" & pca_scores$Spec == "Gambusia_affinis",]$DVS)
t.test(gaff.dvs$TC, gaff.dvs$CC, var.equal = T)

fnot.dvs.tc <- pca_scores[pca_scores$Habitat == "TC" & pca_scores$Spec == "Fundulus_notatus",]$DVS
fnot.dvs.cc <- pca_scores[pca_scores$Habitat == "CC" & pca_scores$Spec == "Fundulus_notatus",]$DVS
t.test(fnot.dvs.tc, fnot.dvs.cc, var.equal = T)

lgul.dvs.tc <- pca_scores[pca_scores$Habitat == "TC" & pca_scores$Spec == "Lepomis_gulosus",]$DVS
lgul.dvs.cc <- pca_scores[pca_scores$Habitat == "CC" & pca_scores$Spec == "Lepomis_gulosus",]$DVS
t.test(lgul.dvs.tc, lgul.dvs.cc, var.equal = T)

pnot.dvs.tc <- pca_scores[pca_scores$Habitat == "TC" & pca_scores$Spec == "Pimephales_notatus",]$DVS
pnot.dvs.cc <- pca_scores[pca_scores$Habitat == "CC" & pca_scores$Spec == "Pimephales_notatus",]$DVS
t.test(pnot.dvs.tc, pnot.dvs.cc, var.equal = T)

lcya.dvs.tc <- pca_scores[pca_scores$Habitat == "TC" & pca_scores$Spec == "Lepomis_cyanellus",]$DVS
lcya.dvs.cc <- pca_scores[pca_scores$Habitat == "CC" & pca_scores$Spec == "Lepomis_cyanellus",]$DVS
t.test(lcya.dvs.tc, lcya.dvs.cc, var.equal = T)

lmac.dvs.tc <- pca_scores[pca_scores$Habitat == "TC" & pca_scores$Spec == "Lepomis_macrochirus",]$DVS
lmac.dvs.cc <- pca_scores[pca_scores$Habitat == "CC" & pca_scores$Spec == "Lepomis_macrochirus",]$DVS
t.test(lmac.dvs.tc, lmac.dvs.cc, var.equal = T)

lmeg.dvs.tc <- pca_scores[pca_scores$Habitat == "TC" & pca_scores$Spec == "Lepomis_megalotis",]$DVS
lmeg.dvs.cc <- pca_scores[pca_scores$Habitat == "CC" & pca_scores$Spec == "Lepomis_megalotis",]$DVS
t.test(lmeg.dvs.tc, lmeg.dvs.cc, var.equal = T)




