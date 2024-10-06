################################################################################
####################### Genomic FIle Converstion/Manipulation ##################
################################################################################
## Load libraries ##
library(SeqArray)
library(tidyverse)
library(tidyr)
library(dplyr)
library(vcfR)
library(SNPRelate)
library(gdsfmt) # tutorial for formating GDS:  https://bioconductor.org/packages/release/bioc/vignettes/gdsfmt/inst/doc/gdsfmt.html
library(SeqVarTools)
library(dartR)
library(ape)

################################################################################
########################### Set working directory ##############################
################################################################################
setwd("/Filtering/filteredVCF") 

################################################################################
############################## Load data files #################################
################################################################################
# VCF file 
vcffile <- "noreponly_v2.filtered.85.vcf.gz"   # input 
vcffile <- "noreponly_v2.filtered.75.vcf.gz"   # input 

vcf <- read.vcfR(vcffile)
vcf

# GDS file
gdsfile <- "noreponly_v2.75.renamed.gds" # what to write it as 
gdsfile <- "noreponly_v2.85.renamed.gds" # what to write it as 

# SNPrelate VCF to GDS file 
snpgdsVCF2GDS(vcffile, gdsfile, method="biallelic.only")

# Summary
snpgdsSummary(gdsfile)

################################################################################
########################### Open/Close GDS file data ###########################
################################################################################
## NOTES: https://uw-gac.github.io/SISG_2021/gds-format.html ##

# open a GDS file
genofile <- seqOpen(gdsfile)

## Read in GDS using SNPrelate
snpgdsSummary(gdsfile)

# Close gds file
seqClose(gds)

## Close all open gds files 
gdsfmt::showfile.gds(closeall=TRUE) # make sure file is not already open

# SNPrelate - Close the GDS file
snpgdsClose(genofile)

################################################################################
############################# GDS file Summary data ############################
################################################################################
## the unique sample identifier comes from the VCF header
sample.id <- seqGetData(genofile, "sample.id") # can read whatever field into object you define 
length(sample.id)
head(sample.id)

# Get the attributes of chromosome coding - SNPrelate
get.attr.gdsn(index.gdsn(genofile, "snp.chromosome"))

# Chromosome is stored in the `seqnames` column. The `ranges` column has variant position, 
# which can be a single base pair or a range
gr <- granges(gds)
gr

## Annotation - "VCF metadata"
# Depth
annot <- seqGetData(gds, "annotation/info/DP")  # when you want nested structure use "/"
head(annot)

# SNP ID
id <- seqGetData(gds, "annotation/id")  # when you want nested structure use "/"
head(id)

# SNPrelate - Read population information
pop <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))
table(pop)

# look at reference and alternate alleles
refChar(gds)  # extract reference allele are for each variant
altChar(gds)
################################################################################
## Genotype Information ##

# Read in Genotypes -  3 dimensions but not most intuitve way to look at genotypes 
geno <- seqGetData(gds, "genotype") 
length(geno)
head(geno)
dim(geno)

# print the first two variants
geno[,,1:2]
geno[,1:5,1:5]

# SeqVarTools Package - Also can use for getting Genotypes
geno <- getGenotype(gds)
dim(geno)
head(geno)

# SNPrelate - Take out genotype data for the first 3 samples and the first 5 SNPs
(g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(5,3)))

# SNPrelate - take out genotype data with sample and SNP IDs, and four possible values are returned 0, 1, 2 and NA (3 is replaced by NA):
g <- snpgdsGetGeno(genofile, sample.id=..., snp.id=...)
g <- snpgdsGetGeno(genofile, sample.id=sample.id, snp.id=snp.id)
head(g)

# SNPrelate - Get the attribute of genotype
get.attr.gdsn(index.gdsn(genofile, "genotype"))

# Genotype alleles
geno.allele <- getGenotypeAlleles(gds)  # returns the letters on the alleles 
head(geno.allele)

# count of how many copies of the reference alleles (good for association testing)
refdos <- refDosage(gds) 
head(refdos)

dos <- alleleDosage(gds, n=0) 
head(dos)

# count of *any* alternate alleles
altdos <- altDosage(gds)
head(altdos)

# count of the first alternate allele 
dos <- alleleDosage(gds, n=1)
head(dos)

# count of the third alternate allele
dos <- alleleDosage(gds, n=3)
head(dos)

# # count of *each* of the alternate alleles separately
# returns multiple columns per variant
dos <- expandedAltDosage(gds)
head(dos)
################################################################################
## Variant Information ##

# a unique integer ID is assigned to each variant
variant.id <- seqGetData(gds, "variant.id") 
length(variant.id)
head(variant.id)

# Variant position
pos <- seqGetData(gds, "position")
head(pos)

# look at reference and alternate alleles
refChar(gds)  # extract reference allele are for each variant
altChar(gds) # alternate allele for each variant

# data.frame of variant information
variantInfo(gds) 

# how many alleles for each variant?
n <- seqNumAllele(gds)
table(n)

# some variants have more than one alternate allele
multi.allelic <- which(n > 2)
altChar(gds)[multi.allelic]

# select snps from GDS file
snpset <- snpgdsSelectSNP(gds, autosome.only=FALSE, maf=0.05, missing.rate=0.95, remove.monosnp=FALSE)
length(snpset)
snpset

# SNPrelate - Take out snp.id
head(read.gdsn(index.gdsn(genofile, "snp.id")))
################################################################################
## Minor allele frequency of each variant ##
maf <- seqAlleleFreq(gds, minor = TRUE)   # from gds file, it will compute minor allele freq for all variants 
head(maf)
summary(maf)
hist(maf, breaks=50)

################################################################################
############################### Data Filters  ##################################
################################################################################
## Set Filters - subsequent reads from the `gds` object are restricted to the selected subset of data,
seqSetFilter(gds, variant.id=, sample.id=)  

# set sample and variant filters
seqSetFilter(gds, variant.id=pruned, ret.idx=TRUE)  
# when you don't want to read in all data at once you use filter and can filter
# on variants and samples only returns on that filter 
# use filter to read in which slice of data you want to read 

# only returns data for the filtered variants
seqGetData(gds, "variant.id")

# Reset Filters for GDS
seqResetFilter(gds) 


################################################################################
################################ LD Pruning ####################################
################################################################################
## SNPrelate 
# LD pruning has a random element; so make this reproducible
set.seed(1000)

snpset <- snpgdsLDpruning(gds, 
                          method="corr", 
                          autosome.only = F,
                          remove.monosnp = T,
                          #slide.max.bp=50000, # 50 kb window
                          ld.threshold=sqrt(0.2))
str(snpset)
names(snpset)

# how many variants on each chr?
sapply(snpset, length)

# get the full list of LD-pruned variants 
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)
pruned

## Create new subsetted gds including only LD-pruned variants 
gds

set.seed(100)

# set sample and variant filters
seqSetFilter(gds, variant.id=pruned, ret.idx=TRUE)  


# the unique sample identifier comes from the VCF header
sample.id <- seqGetData(gds, "sample.id") # can read whatever field into object you define 
length(sample.id)
head(sample.id)

# a unique integer ID is assigned to each variant
variant.id <- seqGetData(gds, "variant.id") 
length(variant.id)
head(variant.id)

geno <- seqGetData(gds, "genotype") # 3 dimensions but not most intuitve way to look at genotypes 
length(geno)
head(geno)


## export 
LDpruned_subset <- "noreponly_v2.filtered.85.renamed.LDpruned.gds"
seqExport(gds, LDpruned_subset)

seqResetFilter(gds) 
seqClose(gds)

# show file
(LDpruned.gds <- seqOpen(LDpruned_subset))

# convert to vcf
seqGDS2VCF(LDpruned.gds, "noreponly_v2.filtered.85.renamed.LDpruned.vcf.gz")

vcffile <- "noreponly_v2.85.renamed.LDpruned.vcf.gz"  # input 
vcf <- read.vcfR(vcffile)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
str(snpset)
names(snpset)

# Get all selected snp id
snpset.id <- unlist(unname(snpset))
head(snpset.id)

################################################################################
################## Converting between Other Genomic filetypes  #################
################################################################################
### Extra steps if you want individual data files based on subsets
################################################################################
# read vcffile 
vcf <- read.vcfR(vcffile)
gendata_gl<-vcfR2genlight(vcf) 

## dartR - Genlight to genofile (.geno) ##
df.geno <- gl2geno(gendata_gl, outfile = "filename", outpath = "data/outpath/")

################################################################################
################# Subset metadata by Popmap - Admixture Groups #################
################################################################################
# Metadata file:
### NOTE: Sampling Localities are under "Monitoring.Unit" in scripts and files: 
md <- read_csv("/PopulationStats_Admixture/admixture_k3_pops/popmap.csv")
md$Ad_cluster

## Subset by population ##
# Morelets
CM_99 <- md %>% 
  filter(grepl("CM_99", Ad_cluster, ignore.case = TRUE))
CM_90 <- md %>% 
  filter(grepl("CM_90", Ad_cluster, ignore.case = TRUE))
pure_morelets <- md %>% 
  filter(Ad_cluster %in% c("CM_90", "CM_99"))
All_CM <- md %>% 
  filter(grepl("CM", Ad_cluster, ignore.case = TRUE))
CM_NRW <- pure_morelets %>% 
  filter(grepl("NRW", Abbrv_Monitoring.Unit, ignore.case = TRUE))
CM_notNRW <- pure_morelets %>% 
  filter(!grepl("NRW", Abbrv_Monitoring.Unit, ignore.case = TRUE))

## Change the "Ad_cluster" pops for the CM
CM_NRW <- CM_NRW %>%
  mutate(Ad_cluster = "CM_NRW")
CM_notNRW <- CM_notNRW %>%
  mutate(Ad_cluster = "CM_notNRW")

# CA
CA_99 <- md %>%
  filter(grepl("CA_99", Ad_cluster, ignore.case = TRUE))
CA_90 <- md %>%
  filter(grepl("CA_90", Ad_cluster, ignore.case = TRUE))
CA <- md %>% 
  filter(Ad_cluster %in% c("CA_90", "CA_99") & !Ad_cluster %in% c("MCA_90", "MCA_99"))
All_CA <- md %>% 
  filter(Ad_cluster %in% c("CA_90", "CA_99","Hybrid_CA_backcross") & !Ad_cluster %in% c("MCA_90", "MCA_99"))

# MCA
MCA_99 <- md %>%
  filter(grepl("MCA_99", Ad_cluster, ignore.case = TRUE))
MCA_90 <- md %>%
  filter(grepl("MCA_90", Ad_cluster, ignore.case = TRUE))
MCA <- md %>% 
  filter(Ad_cluster %in% c("MCA_90", "MCA_99"))
All_MCA <- md %>% 
  filter(grepl("MCA", Ad_cluster, ignore.case = TRUE))

# Acutus
All_acutus <- md %>% 
  filter(grepl("CA", Ad_cluster, ignore.case = TRUE))
All_acutus$Ad_cluster
Acutus_90 <- md %>%
  filter(Ad_acutus >= 0.90)

# Hybrids
F1_Hybrid <- md %>%
  filter(grepl("F1_Hybrid", Ad_cluster, ignore.case = TRUE))
Hybrid_CM_backcross <- md %>%
  filter(grepl("Hybrid_CM_backcross", Ad_cluster, ignore.case = TRUE))
Hybrid_CA_backcross <- md %>%
  filter(grepl("Hybrid_CA_backcross", Ad_cluster, ignore.case = TRUE))
Hybrid_MCA_backcross <- md %>%
  filter(grepl("Hybrid_MCA_backcross", Ad_cluster, ignore.case = TRUE))
All_hybrid <- md %>%
  filter(grepl("Hybrid", Ad_cluster, ignore.case = TRUE))
Hybrid_25.75 <- md %>%
  filter(Ad_CM >= 0.25,Ad_CM <= 0.75)
Acutus_backcross <- md %>%
  filter(grepl("Acutus_backcross", Ad_cluster, ignore.case = TRUE))

# Subset population file to only include retained individuals

################################################################################
################# Subset metadata by Monitoring Unit/Morph Sp ##################
################################################################################
### Set working metadata files ### 
md <- read_csv("noreponly_metadata_fixed.csv")
md

## Data Organization ##

# Monitoring Units 
mu <- md$Monitoring.Unit
mu <- as.factor(mu)
unique(mu)

ab.mu <- md$Abbrv_Monitoring.Unit
ab.mu <- as.factor(ab.mu)
ab.mu
unique(ab.mu)

msp <- md$Morph_Species
msp <- as.factor(msp)
unique(msp)

## Subset populations by ab.MU
(BRW <- md %>% 
    filter(grepl("BRW", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(CB <- md %>% 
    filter(grepl("CB", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(CC <- md %>% 
    filter(grepl("CC", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(CF <- md %>% 
    filter(grepl("CF", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(NBDW <- md %>% 
    filter(grepl("NBDW", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(NC <- md %>% 
    filter(grepl("NC", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(NRW <- md %>% 
    filter(grepl("NRW", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(NSLW <- md %>% 
    filter(grepl("NSLW", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(NTW <- md %>% 
    filter(grepl("NTW", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(RHW <- md %>% 
    filter(grepl("RHW", Abbrv_Monitoring.Unit, ignore.case = TRUE)))
(STW <- md %>% 
    filter(grepl("STW", Abbrv_Monitoring.Unit, ignore.case = TRUE)))

NA.mu <- md %>%
  filter(is.na(Abbrv_Monitoring.Unit))

tibble_list <- list(
  BRW = BRW,
  CB = CB,
  CC = CC,
  CF = CF,
  NBDW = NBDW,
  NC = NC,
  NRW = NRW,
  NSLW = NSLW,
  NTW = NTW,
  RHW = RHW,
  STW = STW,
  NA.mu = NA.mu)
#tibble_list <- list(BRW,CB,CC, CF,NBDW,NC, NRW,NSLW,NTW,RHW,STW,NA.mu)
################################################################################
######################### Make subset gds and vcf data #########################
################################################################################
gdsfile <- "noreponly_v2.75.renamed.gds"

## open a connection to the GDS file - SeqArray Package 
gds <- seqOpen(gdsfile, readonly = FALSE)
gds

################################################################################
## Create each subset file in a new folder named "subset_data" 

CM.gds <- "subset_data/noreponly_v2.75.renamed.pure_CM.gds"
CM.gds <- seqSetFilter(gds, sample.id=pure_morelets$Sample, ret.idx=TRUE)  
seqExport(gds, CM.gds)
## open a connection to the GDS file - SeqArray Package 
gds <- seqOpen(CM.gds, readonly = FALSE)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.pure_CM.vcf.gz")
seqResetFilter(gds) 

all.CM.gds <- "subset_data/noreponly_v2.75.renamed.all_CM.gds"
seqSetFilter(gds, sample.id=All_CM$Sample)  
seqExport(gds, all.CM.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.All_CM.vcf.gz")
seqResetFilter(gds) 

CA.gds <- "subset_data/noreponly_v2.75.renamed.pure_CA.gds"
seqSetFilter(gds, sample.id=CA$Sample)  
seqExport(gds, CA.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.pure_CA.vcf.gz")
seqResetFilter(gds) 
vcf <- read.vcfR("subset_data/noreponly_v2.75.renamed.pure_CA.vcf.gz")
vcf 

all.CA.gds <- "subset_data/noreponly_v2.75.renamed.all_CA.gds"
seqSetFilter(gds, sample.id=All_CA$Sample)  
seqExport(gds, all.CA.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.All_CA.vcf.gz")
seqResetFilter(gds) 

pure.acutus.gds <- "subset_data/noreponly_v2.75.renamed.pure_acutus.gds"
seqSetFilter(gds, sample.id=Acutus_90$Sample)  
seqExport(gds, pure.acutus.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.pure_acutus.vcf.gz")
seqResetFilter(gds) 

MCA.gds <- "subset_data/noreponly_v2.75.renamed.pure_MCA.gds"
seqSetFilter(gds, sample.id=MCA$Sample)  
seqExport(gds, MCA.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.pure_MCA.vcf.gz")
seqResetFilter(gds) 

all.MCA.gds <- "subset_data/noreponly_v2.75.renamed.all_MCA.gds"
seqSetFilter(gds, sample.id=All_MCA$Sample)  
seqExport(gds, all.MCA.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.All_MCA.vcf.gz")
seqResetFilter(gds) 

acutus.gds <- "subset_data/noreponly_v2.75.renamed.all_Acutus.gds"
seqSetFilter(gds, sample.id=All_acutus$Sample)  
seqExport(gds, acutus.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.All_Acutus.vcf.gz")
seqResetFilter(gds) 

All_Hybrid.gds <- "subset_data/noreponly_v2.75.renamed.All_Hybrid.gds"
seqSetFilter(gds, sample.id=All_hybrid$Sample)  
seqExport(gds, All_Hybrid.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.All_Hybrid.vcf.gz")
seqResetFilter(gds) 

F1_Hybrid.gds <- "subset_data/noreponly_v2.75.renamed.F1_Hybrid.gds"
seqSetFilter(gds, sample.id=F1_Hybrid$Sample)  
seqExport(gds, F1_Hybrid.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.F1_Hybrid.vcf.gz")
seqResetFilter(gds) 

Hybrid_CM_backcross.gds <- "subset_data/noreponly_v2.75.renamed.Hybrid_CM_backcross.gds"
seqSetFilter(gds, sample.id=Hybrid_CM_backcross$Sample)  
seqExport(gds, Hybrid_CM_backcross.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.Hybrid_CM_backcross.vcf.gz")
seqResetFilter(gds) 

Hybrid_CA_backcross.gds <- "subset_data/noreponly_v2.75.renamed.Hybrid_CA_backcross.gds"
seqSetFilter(gds, sample.id=Hybrid_CA_backcross$Sample)  
seqExport(gds, Hybrid_CA_backcross.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.Hybrid_CA_backcross.vcf.gz")
seqResetFilter(gds) 

Hybrid_MCA_backcross.gds <- "subset_data/noreponly_v2.75.renamed.Hybrid_MCA_backcross.gds"
seqSetFilter(gds, sample.id=Hybrid_MCA_backcross$Sample)  
seqExport(gds, Hybrid_MCA_backcross.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.Hybrid_MCA_backcross.vcf.gz")
seqResetFilter(gds) 

Hybrid_25.75.gds <- "subset_data/noreponly_v2.75.renamed.Hybrid_25.75.gds"
seqSetFilter(gds, sample.id=Hybrid_25.75$Sample)  
seqExport(gds, Hybrid_25.75.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/noreponly_v2.75.renamed.Hybrid_25.75.vcf.gz")
seqResetFilter(gds) 

vcf <- read.vcfR("subset_data/noreponly_v2.75.renamed.pure_MCA.vcf.gz")
vcf 
################################################################################
## FROM LD pruned data ##
################################################################################
gdsfile <- "noreponly_v2.75.renamed.LDpruned.gds"

## open a connection to the GDS file - SeqArray Package 
gds <- seqOpen(gdsfile, readonly = FALSE)
gds
################################################################################
## Make subset gds and vcf data ##
################################################################################
CM.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.pure_CM.gds"
seqSetFilter(gds, sample.id=pure_morelets$Sample, ret.idx=TRUE)  
seqExport(gds, CM.gds)
## open a connection to the GDS file - SeqArray Package 
gds <- seqOpen(CM.gds, readonly = FALSE)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.pure_CM.vcf.gz")
seqResetFilter(gds) 

all.CM.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.all_CM.gds"
seqSetFilter(gds, sample.id=All_CM$Sample)  
seqExport(gds, all.CM.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.All_CM.vcf.gz")
seqResetFilter(gds) 

CA.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.pure_CA.gds"
seqSetFilter(gds, sample.id=CA$Sample)  
seqExport(gds, CA.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.pure_CA.vcf.gz")
seqResetFilter(gds) 

all.CA.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.all_CA.gds"
seqSetFilter(gds, sample.id=All_CA$Sample)  
seqExport(gds, all.CA.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.All_CA.vcf.gz")
seqResetFilter(gds) 

pure.acutus.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.pure_acutus.gds"
seqSetFilter(gds, sample.id=Acutus_90$Sample)  
seqExport(gds, pure.acutus.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.pure_acutus.vcf.gz")
seqResetFilter(gds) 

MCA.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.pure_MCA.gds"
seqSetFilter(gds, sample.id=MCA$Sample)  
seqExport(gds, MCA.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.pure_MCA.vcf.gz")
seqResetFilter(gds) 

all.MCA.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.all_MCA.gds"
seqSetFilter(gds, sample.id=All_MCA$Sample)  
seqExport(gds, all.MCA.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.All_MCA.vcf.gz")
seqResetFilter(gds) 

acutus.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.all_Acutus.gds"
seqSetFilter(gds, sample.id=All_acutus$Sample)  
seqExport(gds, acutus.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.All_Acutus.vcf.gz")
seqResetFilter(gds) 

All_Hybrid.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.All_Hybrid.gds"
seqSetFilter(gds, sample.id=All_hybrid$Sample)  
seqExport(gds, All_Hybrid.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.All_Hybrid.vcf.gz")
seqResetFilter(gds) 

F1_Hybrid.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.F1_Hybrid.gds"
seqSetFilter(gds, sample.id=F1_Hybrid$Sample)  
seqExport(gds, F1_Hybrid.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.F1_Hybrid.vcf.gz")
seqResetFilter(gds) 

Hybrid_CM_backcross.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.Hybrid_CM_backcross.gds"
seqSetFilter(gds, sample.id=Hybrid_CM_backcross$Sample)  
seqExport(gds, Hybrid_CM_backcross.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.Hybrid_CM_backcross.vcf.gz")
seqResetFilter(gds) 

Hybrid_CA_backcross.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.Hybrid_CA_backcross.gds"
seqSetFilter(gds, sample.id=Hybrid_CA_backcross$Sample)  
seqExport(gds, Hybrid_CA_backcross.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.Hybrid_CA_backcross.vcf.gz")
seqResetFilter(gds) 

Hybrid_MCA_backcross.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.Hybrid_MCA_backcross.gds"
seqSetFilter(gds, sample.id=Hybrid_MCA_backcross$Sample)  
seqExport(gds, Hybrid_MCA_backcross.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.Hybrid_MCA_backcross.vcf.gz")
seqResetFilter(gds) 

Hybrid_25.75.gds <- "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.Hybrid_25.75.gds"
seqSetFilter(gds, sample.id=Hybrid_25.75$Sample)  
seqExport(gds, Hybrid_25.75.gds)
# convert to vcf
seqGDS2VCF(gds, "subset_data/LDpruned/noreponly_v2.75.renamed.LDpruned.Hybrid_25.75.vcf.gz")
seqResetFilter(gds) 
