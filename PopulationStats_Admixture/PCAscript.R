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
library(scales)

###########################################################################################################################
## noreponly_v2: All samples with no repeats - 273 individuals ##
# For filtered through snpFiltR data 
data_dir <- "/PopulationStats_Admixture"

## Set working directory ##
setwd(data_dir)

################################################################################
## Get in MetaData ## 

## add population data ##
# Population data
md <- read_csv("/admixture_k3_pops/popmap.csv")
md$Ad_cluster
md

## Change names for Morphological Species 
md$Morph_Species <- gsub("MCA", "Mainland C. acutus", md$Morph_Species)
md$Morph_Species <- gsub("CA", "Cayes C. acutus", md$Morph_Species)
md$Morph_Species <- gsub("HY", "Putative Hybrids", md$Morph_Species)
md$Morph_Species <- gsub("CM", "C. moreletii", md$Morph_Species)
md$Morph_Species
as.factor(md$Morph_Species)

# Replace NA values in a specific column
# pops <- pops %>%
#mutate(Monitoring.Unit = ifelse(is.na(Monitoring.Unit), "unknown", Monitoring.Unit))
# pops$Monitoring.Unit
md$Monitoring.Unit <- str_replace_na(md$Monitoring.Unit, "unknown")

md$Monitoring.Unit <- gsub("N/S Lagoon Watershed", "North/South Lagoon Watershed", md$Monitoring.Unit)
md$Monitoring.Unit <- as.factor(md$Monitoring.Unit)
md$Monitoring.Unit

## Adjust Population names ##
pops<-md
pops$Ad_cluster
pops$Monitoring.Unit
pops <- pops %>% select(Sample, Ad_cluster, Morph_Species, sNMF_cluster, Monitoring.Unit)

pops$Ad_cluster <- gsub("CM_90", "C. moreletii", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("CM_99", "C. moreletii", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("MCA_90", "acutus_A (Mainland lineage)", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("MCA_99", "acutus_A (Mainland lineage)", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("CA_90", "acutus_B (Cayes lineage)", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("CA_99", "acutus_B (Cayes lineage)", pops$Ad_cluster) 
#pops$Ad_cluster <- gsub(".*Hybrid.*", "Hybrid", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("Acutus_backcross", "Hybrid_acutus", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("F1_Hybrid", "Hybrid", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("Hybrid_MCA_backcross", "Hybrid_acutus_A", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("Hybrid_CA_backcross", "Hybrid_acutus_B", pops$Ad_cluster) 
pops$Ad_cluster <- gsub("Hybrid_CM_backcross", "Hybrid_moreletii", pops$Ad_cluster) 
as.factor(pops$Ad_cluster)

pops$sNMF_cluster
pops$sNMF_cluster <- gsub("CM_90", "moreletii", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("CM_99", "moreletii", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("MCA_90", "acutus_A (Mainland lineage)", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("MCA_99", "acutus_A (Mainland lineage)", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("CA_90", "acutus_B (Cayes lineage)", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("CA_99", "acutus_B (Cayes lineage)", pops$sNMF_cluster) 
#pops$sNMF_cluster <- gsub(".*Hybrid.*", "Hybrid", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("Acutus_backcross", "Hybrid_acutus", pops$sNMF_cluster)
pops$sNMF_cluster <- gsub("F1_Hybrid", "Hybrid", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("Hybrid_MCA_backcross", "Hybrid_acutus_A", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("Hybrid_CA_backcross", "Hybrid_acutus_B", pops$sNMF_cluster) 
pops$sNMF_cluster <- gsub("Hybrid_CM_backcross", "Hybrid_moreletii", pops$sNMF_cluster) 
as.factor(pops$sNMF_cluster)

################################################################################
################################################################################
basefile <- "noreponly_v2.75.renamed.LDpruned"

vcf_gds <-paste0(out_dir,"/", basefile,".gds")
vcf_gds

## noreponly_v2 PCA
# genofile used for admixture 
path_vcf <- "/Filtering/filteredVCF/noreponly_v2.75.renamed.LDpruned.vcf.gz"

## convert to GDS format for SNPRelate analysis
#snpgdsVCF2GDS(path_vcf, vcf_gds, method="biallelic.only",ignore.chr.prefix="Chr")

genofile <- snpgdsOpen(vcf_gds, readonly=FALSE)
genofile 

#summarize data
snpgdsSummary(vcf_gds)

## Already LD pruned data
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
out_dir <- paste0(data_dir, "/PCA/")

if(!dir.exists(out_dir)){ # check if the directory exists
  dir.create(out_dir)   # and create it if it does not
}

setwd(out_dir)
list.files(out_dir)
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

#hue_pal()(4) # see colors 

## Plot Using ggplot
PCA.msp <- tab %>% ggplot() + 
  geom_point(aes(x = EV1, y = EV2, fill = pop), size = 5, shape = 21) +
  #geom_text(aes(x = EV1, y = EV2, label = sample.id)) + # adds individual labels
  #scale_fill_manual(values = myPalette[c(3:1)], na.value = "grey80") +
  scale_fill_manual("Morphological Species Groups", values = c("Cayes C. acutus" = "#00BFC4", 
                                        "Mainland C. acutus" = "#7CAE00", 
                                        "Putative Hybrids" = "#C77CFF", 
                                        "C. moreletii" = "#F8766D")) + 
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"), 
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       fill = "Morphological Species Groups") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
PCA.msp

ggsave(
  filename = "PCA_msp.png",
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", PCA.msp
)

ggsave(
  filename = "PCA_msp.pdf",
  width = 8, height = 7, dpi = 600,
  device = "pdf", bg = "white", PCA.msp
)

pdf("PCA_by_morph_species.pdf")
PCA.msp
dev.off()

png(file="PCA_by_morph_species.png", width = 800, height = 800) 
PCA.msp
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
PCA.mu <- tab %>% ggplot() + geom_point(aes(x = EV1, y = EV2, fill = pop), size = 5, shape = 21) +
  #geom_text(aes(x = EV1, y = EV2, label = sample.id)) + # adds individual labels
  #scale_fill_manual(values = myPalette[c(3:1)], na.value = "grey80") +
  scale_fill_manual("Sampling Localities", values = c(
    "Belize River Watershed" = "#F8766D",
    "Chiquibul Forest" = "#ff7f0e",
    "Coastal Cayes" = "#A3A500",
    "Cockscomb Basin" = "#39B600",
    "North/South Lagoon Watershed" = "#8c564b",
    "New River Watershed" = "#00BFC4",
    "Northern Cayes" = "#00B0F6",
    "Northern Toledo Watershed" = "#9590FF",
    "Rio Hondo Watershed" = "#E76BF3",
    "Southern Toledo Watershed" = "#FF62BC",
    "unknown" = "grey"
  )) +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       fill = "Sampling Localities") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
PCA.mu

ggsave(
  filename = "PCA_mu.png",
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", PCA.mu
)

ggsave(
  filename = "PCA_mu.pdf",
  width = 8, height = 7, dpi = 600,
  device = "pdf", bg = "white", PCA.mu
)


pdf("PCA_by_Monitoring_Unit.pdf")
PCA.mu
dev.off()

png(file="PCA_by_Monitoring_Unit.png", width = 800, height = 800)
PCA.mu
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
PCA.msp.mu <- ggplot() + 
  geom_point(data = tab, 
             aes(x = EV1, y = EV2, color = factor(pop.mu), 
                 shape = factor(pop.msp)), size = 5) +
  scale_shape_manual("Morphological Species Groups", 
                     values = c("Cayes C. acutus" = 1, 
                                "Mainland C. acutus" = 16, 
                                "Putative Hybrids" = 17, 
                                "C. moreletii" = 4)) +
  scale_color_manual("Sampling Localities", values = c(
    "Belize River Watershed" = "#F8766D",
    "Chiquibul Forest" = "#ff7f0e",
    "Coastal Cayes" = "#A3A500",
    "Cockscomb Basin" = "#39B600",
    "North/South Lagoon Watershed" = "#8c564b",
    "New River Watershed" = "#00BFC4",
    "Northern Cayes" = "#00B0F6",
    "Northern Toledo Watershed" = "#9590FF",
    "Rio Hondo Watershed" = "#E76BF3",
    "Southern Toledo Watershed" = "#FF62BC",
    "unknown" = "grey"
  )) +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       color = "Sampling Localities", shape = "Morphological Species") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
PCA.msp.mu


ggsave(
  filename = "PCA_msp.mu.png",
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", PCA.msp.mu
)

ggsave(
  filename = "PCA_msp.mu.pdf",
  width = 8, height = 7, dpi = 600,
  device = "pdf", bg = "white", PCA.msp.mu
)

pdf("PCA_by_MU_Msp.pdf")
PCA.msp.mu
dev.off()

png(file="PCA_by_MU_Msp.png", width = 800, height = 800) # open the pdf plotting device
PCA.msp.mu
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
PCA.ad = tab %>% ggplot() + geom_point(aes(x = EV1, y = EV2, fill = pop), size = 5, shape = 21) +
  #geom_text(aes(x = EV1, y = EV2, label = sample.id)) + # adds individual labels
  #scale_fill_manual(values = myPalette[c(3:1)], na.value = "grey80") +
  scale_fill_manual("Admixture Group Assignments", 
                     breaks = c(
                       "acutus_A (Mainland lineage)",
                       "acutus_B (Cayes lineage)",
                       "Hybrid",
                       "Hybrid_acutus",
                       "Hybrid_acutus_A",
                       "Hybrid_acutus_B",
                       "Hybrid_moreletii",
                       "C. moreletii"),
                     values = c(
                       "acutus_A (Mainland lineage)" = "#F8766D",
                       "acutus_B (Cayes lineage)" = "#CD9600",
                       "Hybrid"=  "#7CAE00",
                       "Hybrid_acutus" = "#00BE67",
                       "Hybrid_acutus_A" = "#00BFC4",
                       "Hybrid_acutus_B" = "#00A9FF",
                       "Hybrid_moreletii" = "#C77CFF",
                       "C. moreletii" = "#FF61CC")) +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       fill = "Admixture Group Assignments") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
PCA.ad

ggsave(
  filename = "PCA_ad.png",
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", PCA.ad
)

ggsave(
  filename = "PCA_ad.pdf",
  width = 8, height = 7, dpi = 600,
  device = "pdf", bg = "white", PCA.ad
)

pdf("PCA_by_Ad_cluster.pdf")
PCA.ad
dev.off()

png(file="PCA_by_Ad_cluster.png", width = 800, height = 800) # open the pdf plotting device
PCA.ad
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
unique(tab$ad_clust)

# Plot with ggplot: genetic PC1 vs. genetic PC2
PCA.ad.msp <- ggplot() + 
  geom_point(data = tab, 
             aes(x = EV1, y = EV2, color = factor(ad_clust), 
                 shape = factor(pop.msp)), size = 5) +
  scale_shape_manual("Morphological Species Groups", 
                     values = c("Cayes C. acutus" = 1, 
                                "Mainland C. acutus" = 16, 
                                "Putative Hybrids" = 17, 
                                "C. moreletii" = 4)) +
  scale_color_manual("Admixture Group Assignments", 
                     breaks = c(
                       "acutus_A (Mainland lineage)",
                       "acutus_B (Cayes lineage)",
                       "Hybrid","Hybrid_acutus",
                       "Hybrid_acutus_A",
                       "Hybrid_acutus_B",
                       "Hybrid_moreletii",
                       "C. moreletii"),
                     values = c(
                       "acutus_A (Mainland lineage)" = "#F8766D",
                       "acutus_B (Cayes lineage)" = "#CD9600",
                       "Hybrid"=  "#7CAE00",
                       "Hybrid_acutus" = "#00BE67",
                       "Hybrid_acutus_A" = "#00BFC4",
                       "Hybrid_acutus_B" = "#00A9FF",
                       "Hybrid_moreletii" = "#C77CFF",
                       "C. moreletii" = "#FF61CC")) +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       color = "Admixture Group Assignments", shape = "Morphological Species") 
#theme(plot.title = element_text(size = 30, hjust = 0.5),
#axis.text = element_text(size = 24),
#axis.title = element_text(size = 22),
#panel.grid = element_blank()
PCA.ad.msp
dev.off()

ggsave(
  filename = "PCA_ad.msp.png",
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", PCA.ad.msp
)

ggsave(
  filename = "PCA_ad.msp.pdf",
  width = 8, height = 7, dpi = 600,
  device = "pdf", bg = "white", PCA.ad.msp
)

ggsave(
  filename = "PCA_ad.msp.tiff",
  width = 8.5, height = 7, dpi = 600,
  device = "tiff", bg = "white", PCA.ad.msp
)

pdf("PCA_by_AdCluster_Msp.pdf")
PCA.ad.msp
dev.off()

png(file="PCA_by_AdCluster_Msp.png", width = 800, height = 800) # open the pdf plotting device
PCA.ad.msp
dev.off()

tiff(file="PCA_by_AdCluster_Msp.tiff", width = 800, height = 800) # open the pdf plotting device
PCA.ad.msp
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
PCA.ad.mu <- ggplot() + 
  geom_point(data = tab, 
             aes(x = EV1, y = EV2, color = factor(pop.mu), 
                 shape = factor(ad_clust)), size = 5) +
  scale_shape_manual("Admixture Populations", 
                     breaks = c(
                       "acutus_A (Mainland lineage)",
                       "acutus_B (Cayes lineage)",
                       "Hybrid","Hybrid_acutus",
                       "Hybrid_acutus_A",
                       "Hybrid_acutus_B",
                       "Hybrid_moreletii",
                       "C. moreletii"),
                     values = c(
                       "acutus_A (Mainland lineage)" = 16,
                       "acutus_B (Cayes lineage)" = 1,
                       "Hybrid"=  13,
                       "Hybrid_acutus" = 8,
                       "Hybrid_acutus_A" = 17,
                       "Hybrid_acutus_B" = 2,
                       "Hybrid_moreletii" = 7,
                       "C. moreletii" = 4)) +
  scale_color_manual("Sampling Localities", values = c(
    "Belize River Watershed" = "#F8766D",
    "Chiquibul Forest" = "#ff7f0e",
    "Coastal Cayes" = "#A3A500",
    "Cockscomb Basin" = "#39B600",
    "North/South Lagoon Watershed" = "#8c564b",
    "New River Watershed" = "#00BFC4",
    "Northern Cayes" = "#00B0F6",
    "Northern Toledo Watershed" = "#9590FF",
    "Rio Hondo Watershed" = "#E76BF3",
    "Southern Toledo Watershed" = "#FF62BC",
    "unknown" = "grey"
  )) +
  theme_light() +
  labs(y = paste0("PC2", " ", "(", eig.vect2, "%",")"),
       x = paste0("PC1", " ", "(", eig.vect1, "%",")"),
       shape = "Admixture Group Assignment", color = "Sampling Locality") 
PCA.ad.mu

ggsave(
  filename = "PCA_ad.mu.png",
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", PCA.ad.mu
)

ggsave(
  filename = "PCA_ad.mu.pdf",
  width = 8, height = 7, dpi = 600,
  device = "pdf", bg = "white", PCA.ad.mu
)

pdf("PCA_by_AdCluster_MU.pdf")
PCA.ad.mu
dev.off()

png(file="PCA_by_AdCluster_MU.png", width = 800, height = 800) # open the pdf plotting device
PCA.ad.mu
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

