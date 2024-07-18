#launch needed libraries
library(pacman)
library(gdsfmt)
library(SNPRelate)
library(VariantAnnotation)
library(RColorBrewer)
library(tidyverse)
library(tidyr)
library(dplyr)
library(vcfR)
library(dartR)
library(snpR)
library(poppr)
library(LEA)
library(poppr)
library(ape)
library(ggplot2)

###########################################################################################################################
## noreponly_v2: All 3RADmerged samples with no repeats - 273 individuals ##

data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly_v2"                      # Mac
#data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly"                        # PC

# For filtered through snpFiltR data 
filteredVCF <- paste0(data_dir, "/filteredVCF_ab")
filteredVCF
data_dir <- filteredVCF
data_dir

## Set working directory ##
setwd(data_dir)

################################################################################
## Get in MetaData ## 

## add population data ##
# Population data
md <- read_csv("admixture/popdata/popmap.csv")
md$Ad_cluster

## Adjust Population names ##
pops<-md
pops$Ad_cluster
pops$Monitoring.Unit
pops <- pops %>% select(Sample, Ad_cluster, Morph_Species, sNMF_cluster, Monitoring.Unit)

pops$Ad_cluster <- gsub("CM_90", "moreletti", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("CM_99", "moreletti", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("MCA_90", "acutus_A", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("MCA_99", "acutus_A", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("CA_90", "acutus_B", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("CA_99", "acutus_B", pops$Ad_cluster) 
#pops$Ad_cluster <- gsub(".*Hybrid.*", "Hybrid", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("Acutus_backcross", "Hybrid_acutus", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("F1_Hybrid", "Hybrid", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("Hybrid_MCA_backcross", "Hybrid_acutus_A", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("Hybrid_CA_backcross", "Hybrid_acutus_B", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("Hybrid_CM_backcross", "Hybrid_moreletii", pops$Ad_cluster) 
as.factor(pops$Ad_cluster)

pops$sNMF_cluster
pops$sNMF_cluster <- gsub("CM_90", "moreletti", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("CM_99", "moreletti", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("MCA_90", "acutus_A", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("MCA_99", "acutus_A", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("CA_90", "acutus_B", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("CA_99", "acutus_B", pops$sNMF_cluster) 
#pops$sNMF_cluster <- gsub(".*Hybrid.*", "Hybrid", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("Acutus_backcross", "Hybrid_acutus", pops$sNMF_cluster)
pops$sNMF_cluster <- gsub("F1_Hybrid", "Hybrid", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("Hybrid_MCA_backcross", "Hybrid_acutus_A", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("Hybrid_CA_backcross", "Hybrid_acutus_B", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("Hybrid_CM_backcross", "Hybrid_moreletii", pops$sNMF_cluster) 
as.factor(pops$sNMF_cluster)

pops$Morph_Species <- gsub("MCA", "acutus_Mainland", pops$Morph_Species) 
pops$Morph_Species <- gsub("CA", "acutus_Cayes", pops$Morph_Species) 
pops$Morph_Species <- gsub("CM", "moreletii", pops$Morph_Species) 
pops$Morph_Species <- gsub("HY", "Hybrid", pops$Morph_Species) 
as.factor(pops$Morph_Species)

# Replace NA values in a specific column
#pops <- pops %>%
  #mutate(Monitoring.Unit = ifelse(is.na(Monitoring.Unit), "unknown", Monitoring.Unit))
#pops$Monitoring.Unit

################################################################################
################################################################################
out_dir <- paste0(data_dir, "/PCA")
out_dir

basefile <- "noreponly_v2.75.renamed.LDpruned"

vcf_gds <-paste0(out_dir,"/", basefile,".gds")
vcf_gds

## noreponly_v2 PCA
# genofile used for admixture 
path_vcf <- "/Users/hwsung/Library/CloudStorage/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly_v2/filteredVCF_ab/data/noreponly_v2.75.renamed.LDpruned.vcf.gz"

## convert to GDS format for SNPRelate analysis
#snpgdsVCF2GDS(path_vcf, vcf_gds, method="biallelic.only",ignore.chr.prefix="Chr")

genofile <- snpgdsOpen(vcf_gds, readonly=FALSE)
genofile 

#summarize data
snpgdsSummary(vcf_gds)

## Already LD pruned 
pca <- snpgdsPCA(genofile, autosome.only=FALSE, num.thread=2)
show(pca)
summary(pca)
tw <- tracy.widom(pca)

# variance proportion (%)
pca
pc.percent <- pca$varprop*100
eig.vect.values <- head(round(pc.percent, 2))
eig.vect1 <- eig.vect.values[1]
eig.vect2 <- eig.vect.values[2]
eig.vect3 <- eig.vect.values[3]
eig.vect4 <- eig.vect.values[4]

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
################################################################################
## Match dataframe for pca with samples from pops metadata ##

# Get sample id
# Makes df from matching the sample.id from tabs and the Samples from pops
pca$sample.id
pops
tab
matches <- inner_join(tab, pops, c("sample.id" = "Sample"))
#matches <- inner_join(tab, q_mat, c("sample.id" = "Sample"))
head(matches)
names(matches)

morph_sp <- matches$Morph_Species
sample.id <- matches$sample.id
Ad_cluster <- matches$Ad_cluster
sNMF_sp <- matches$sNMF_cluster
mont.unit <- matches$Monitoring.Unit
################################################################################
## Set out directory ##
setwd(out_dir)
################################################################################
## PCA by Morphological Species ##

morph_sp <- matches$Morph_Species
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(morph_sp)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
pop <- tab$pop
pop

## Plot Using ggplot
pdf("PCA_by_morph_species.pdf")
png(file="PCA_by_morph_species.png", width = 800, height = 800) # open the pdf plotting device
tab %>% ggplot() + 
  geom_point(aes(x = EV1, y = EV2, fill = pop), size = 5, shape = 21) +
  #geom_text(aes(x = EV1, y = EV2, label = sample.id)) + # adds individual labels
  #scale_fill_manual(values = myPalette[c(3:1)], na.value = "grey80") +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"), 
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       fill = "Morphological Groups") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
dev.off()

################################################################################
## PCA by Monitoring Unit##

head(matches)
mont.unit <- matches$Monitoring.Unit
class(mont.unit)
mont.unit <- as.factor(mont.unit)
levels(mont.unit)
mont.unit

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(mont.unit)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)

# Plot with ggplot: genetic PC1 vs. genetic PC2
pdf("PCA_by_Monitoring_Unit.pdf")
png(file="PCA_by_Monitoring_Unit.png", width = 800, height = 800) # open the pdf plotting device
gen_gen_plot = tab %>% ggplot() + geom_point(aes(x = EV1, y = EV2, fill = pop), size = 5, shape = 21) +
  #geom_text(aes(x = EV1, y = EV2, label = sample.id)) + # adds individual labels
  #scale_fill_manual(values = myPalette[c(3:1)], na.value = "grey80") +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       fill = "Monitoring Unit") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
gen_gen_plot
dev.off()

################################################################################
## PCA by Morph and Monitoring Unit##

tab <- data.frame(sample.id = pca$sample.id,
                  pop.mu = factor(mont.unit)[match(pca$sample.id, sample.id)],
                  pop.msp = factor(morph_sp)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)

# Plot with ggplot: genetic PC1 vs. genetic PC2
pdf("PCA_by_MU_Msp.pdf")
png(file="PCA_by_MU_Msp.png", width = 800, height = 800) # open the pdf plotting device

gen_gen_plot <- ggplot() + 
  geom_point(data = tab, 
             aes(x = EV1, y = EV2, color = factor(pop.mu), 
                 shape = factor(pop.msp)), size = 5) +
  scale_shape_manual(values = c(1, 16, 17, 4)) +
  #scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+ # assigns individual colors
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       color = "Monitoring Unit", shape = "Morphological Species") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
gen_gen_plot
dev.off()

################################################################################
## PCA by Admixture Groups ##

head(matches)
ad_clust <- matches$Ad_cluster
class(ad_clust)
ad_clust <- as.factor(ad_clust)
levels(ad_clust)
ad_clust

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(ad_clust)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)

# Plot with ggplot: genetic PC1 vs. genetic PC2
pdf("PCA_by_Ad_cluster.pdf")
png(file="PCA_by_Ad_cluster.png", width = 800, height = 800) # open the pdf plotting device
gen_gen_plot = tab %>% ggplot() + geom_point(aes(x = EV1, y = EV2, fill = pop), size = 5, shape = 21) +
  #geom_text(aes(x = EV1, y = EV2, label = sample.id)) + # adds individual labels
  #scale_fill_manual(values = myPalette[c(3:1)], na.value = "grey80") +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       fill = "Population") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
gen_gen_plot
dev.off()
################################################################################
## PCA by Morph and Admixture Groups ##

tab <- data.frame(sample.id = pca$sample.id,
                  ad_clust = factor(ad_clust)[match(pca$sample.id, sample.id)],
                  pop.msp = factor(morph_sp)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)


# Plot with ggplot: genetic PC1 vs. genetic PC2
pdf("PCA_by_AdCluster_Msp.pdf")
png(file="PCA_by_AdCluster_Msp.png", width = 800, height = 800) # open the pdf plotting device

gen_gen_plot <- ggplot() + 
  geom_point(data = tab, 
             aes(x = EV1, y = EV2, color = factor(ad_clust), 
                 shape = factor(pop.msp)), size = 5) +
  scale_shape_manual(values = c(1, 16, 17, 4)) +
  #scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+ # assigns individual colors
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       color = "Admixture Populations", shape = "Morphological Species") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
gen_gen_plot
dev.off()
################################################################################
## PCA by MU and Admixture Groups ##

tab <- data.frame(sample.id = pca$sample.id,
                  ad_clust = factor(ad_clust)[match(pca$sample.id, sample.id)],
                  pop.mu = factor(mont.unit)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)


# Plot with ggplot: genetic PC1 vs. genetic PC2
pdf("PCA_by_AdCluster_MU.pdf")
png(file="PCA_by_AdCluster_MU.png", width = 800, height = 800) # open the pdf plotting device

gen_gen_plot <- ggplot() + 
  geom_point(data = tab, 
             aes(x = EV1, y = EV2, color = factor(ad_clust), 
                 shape = factor(morph_sp)), size = 5) +
  scale_shape_manual(values = c(1, 16, 17, 4)) +
  #scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+ # assigns individual colors
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       color = "Admixture Populations", shape = "Morphological Species") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
gen_gen_plot
dev.off()
################################################################################
## PCA by sNMF Groups ##

head(matches)
sNMF <- matches$sNMF_cluster
class(sNMF)
sNMF <- as.factor(sNMF)
levels(sNMF)
sNMF

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(sNMF)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
head(tab)

# Plot with ggplot: genetic PC1 vs. genetic PC2
pdf("PCA_by_sNMF_cluster.pdf")
png(file="PCA_by_sNMF_cluster.png", width = 800, height = 800) # open the pdf plotting device
gen_gen_plot = tab %>% ggplot() + geom_point(aes(x = EV1, y = EV2, fill = pop), size = 5, shape = 21) +
  #geom_text(aes(x = EV1, y = EV2, label = sample.id)) + # adds individual labels
  #scale_fill_manual(values = myPalette[c(3:1)], na.value = "grey80") +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       fill = "Population") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
gen_gen_plot
dev.off()

