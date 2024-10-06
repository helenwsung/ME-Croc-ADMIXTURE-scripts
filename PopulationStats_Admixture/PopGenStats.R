################################################################################
####################### Population genetic summary statistics ##################
################################################################################
## Load libraries 
library(adegenet)
library(hierfstat)
library(vcfR)
library(tidyverse)
library(tidyr)
library(dplyr)
################################################################################
# set working directory
data_dir <- "/your/data/dir/"

setwd(data_dir)

# Create out directory 
out_dir <- paste0(data_dir, "/PopulationStats_Admixture/diversity_stats") 

if(!dir.exists(out_dir)){ # check if the directory exists
  dir.create(out_dir)   # and create it if it does not
}

################################################################################
# Read vcf into genlight format
vcffile <- "/Filtering/filteredVCF/noreponly_v2.filtered.85.vcf.gz"   
vcf <- read.vcfR(vcffile)

# convert vcf to genind
genind_data <- vcfR2genind(vcf)
genind_data
indNames(genind_data)

################################################################################
####################### Add in population information ##########################
################################################################################
# Metadata file
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

# MCA
MCA_99 <- md %>%
  filter(grepl("MCA_99", Ad_cluster, ignore.case = TRUE))
MCA_90 <- md %>%
  filter(grepl("MCA_90", Ad_cluster, ignore.case = TRUE))
MCA <- md %>% 
  filter(Ad_cluster %in% c("MCA_90", "MCA_99"))
All_MCA <- md %>% 
  filter(grepl("MCA", Ad_cluster, ignore.case = TRUE))

# CA
CA_99 <- md %>%
  filter(grepl("CA_99", Ad_cluster, ignore.case = TRUE))
CA_90 <- md %>%
  filter(grepl("CA_90", Ad_cluster, ignore.case = TRUE))
CA <- md %>% 
  filter(Ad_cluster %in% c("CA_90", "CA_99") & !Ad_cluster %in% c("MCA_90", "MCA_99"))
All_CA <- md %>% 
  filter(Ad_cluster %in% c("CA_90", "CA_99","Hybrid_CA_backcross") & !Ad_cluster %in% c("MCA_90", "MCA_99"))

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

################################################################################
######## Subset genind object to retain only "pure" species groups #############
################################################################################
## Make a vector of Samples to keep - CA, MCA, Pure morelets ##
P1 <- pure_morelets
P2 <- MCA
P3 <- CA
pop <- bind_rows(P1, P2, P3)

# Change pop names to get rid of "_90" and "_99"
pop$Ad_cluster <- gsub("_90", "", pop$Ad_cluster) 
pop$Ad_cluster <- gsub("_99", "", pop$Ad_cluster) 
pop$Ad_cluster

# get sample names
desired_samples <- pop$Sample

# Subset genind object by sample names
subset_genind <- genind_data[indNames(genind_data) %in% desired_samples,]
subset_genind 

# Assuming all samples matched and there are no NAs, you can now update the 
#     population information in your genind object
subset_genind$pop <- factor(pop$Ad_cluster)
population <- pop(subset_genind)
table(population)

#############################################
## Make a separate dataset for pure acutus ##
#############################################
P1 <- pure_morelets
P2 <- Acutus_90
pop.2 <- bind_rows(P1, P2)

# Change pop names to get rid of "_90" and "_99"
pop.2$Ad_cluster <- gsub(".*CA.*", "acutus", pop.2$Ad_cluster) 
pop.2$Ad_cluster <- gsub("_90", "", pop.2$Ad_cluster) 
pop.2$Ad_cluster <- gsub("_99", "", pop.2$Ad_cluster) 
pop.2$Ad_cluster

# get sample names
desired_samples <- pop.2$Sample

# Subset genind object by sample names
genind2 <- genind_data[indNames(genind_data) %in% desired_samples,]
genind2

# Assuming all samples matched and there are no NAs, you can now update the population information in your genind object
genind2$pop <- factor(pop.2$Ad_cluster)
pop_2 <- pop(genind2)
table(pop_2)

################################################################################
# make subsetted genind working data 
RAD <- subset_genind  # CA, MCA, CM
RAD2 <- genind2       # pure acutus, pure moreletii


# convert to hierfstat format
RADhf <- genind2hierfstat(RAD, pop=population)
RADhf2 <- genind2hierfstat(RAD2, pop=pop_2)
RADhf

################################################################################
################ Calculate observed and expected heterozygosity ################
################################################################################
# Get basic stats 
hfstats <- basic.stats(RADhf)
hfstats2 <- basic.stats(RADhf2)

## calculate observed heterozygosity population-specific values
Hobs <- as.data.frame(hfstats$Ho)
Hobs2 <- as.data.frame(hfstats2$Ho)

# across K=3 
CA <- mean(Hobs$CA) 
MCA <- mean(Hobs$MCA)
CM <- mean(Hobs$CM)
mean_Hobs_K3 <- cbind.data.frame(CA, MCA, CM)
mean_Hobs_K3

# across K=2
acutus <- mean(Hobs2$acutus)
morelet <- mean(Hobs2$CM)
mean_Hobs_K2 <-cbind.data.frame(acutus, morelet)
mean_Hobs_K2

# calculate expected heterozygosity population-specific values 
Hess <- as.data.frame(hfstats$Hs)
Hess2 <- as.data.frame(hfstats2$Hs)

# across K=3 
CA <- mean(Hess$CA) 
MCA <- mean(Hess$MCA)
CM <- mean(Hess$CM)
mean_Hess_K3 <- cbind.data.frame(CA, MCA, CM)
mean_Hess_K3

# across K=2
acutus <- mean(Hess2$acutus)
morelet <- mean(Hess2$CM)
mean_Hess_K2 <-cbind.data.frame(acutus, morelet)
mean_Hess_K2

# across K=3 
t.test(Hess$CA, Hobs$CA, paired=TRUE, var.equal = TRUE)
t.test(Hess$MCA, Hobs$MCA, paired=TRUE, var.equal = TRUE)
t.test(Hess$CM, Hobs$CM, paired=TRUE, var.equal = TRUE)
t.test(hfstats$Hs, hfstats$Ho, paired=TRUE, var.equal = TRUE)

# across K=2
t.test(Hess2$acutus, Hobs2$acutus, paired=TRUE, var.equal = TRUE)
t.test(Hess2$CM, Hobs2$CM, paired=TRUE, var.equal = TRUE)
t.test(hfstats2$Hs, hfstats2$Ho, paired=TRUE, var.equal = TRUE)

## For 3 pop groups 
# calculate allelic richness
AR <- allelic.richness(RADhf)
ARdf <- as.data.frame(AR$Ar)

CA <- mean(ARdf$CA)
MCA <- mean(ARdf$MCA)
CM <- mean(ARdf$CM)
mean_AR_K3 <- cbind.data.frame(CA, MCA, CM)
mean_AR_K3

# calculate inbreeding coefficients
Fis <- as.data.frame(hfstats$Fis)

# calculate population-specific values 
CA <- mean(Fis$CA, na.rm = TRUE) 
MCA <- mean(Fis$MCA, na.rm = TRUE)
CM <- mean(Fis$CM, na.rm = TRUE)
mean_Fis_K3 <- cbind.data.frame(CA,MCA,CM)
mean_Fis_K3

## For 2 pop groups 
# calculate allelic richness
AR2 <- allelic.richness(RADhf2)
ARdf2 <- as.data.frame(AR2$Ar)

acutus <- mean(ARdf2$acutus)
morelet <- mean(ARdf2$CM)
mean_AR_K2 <- cbind.data.frame(acutus, morelet)
mean_AR_K2 

# calculate inbreeding coefficients
Fis2 <- as.data.frame(hfstats2$Fis)

# calculate population-specific values
acutus <- mean(Fis2$acutus, na.rm = TRUE)
morelet <- mean(Fis2$CM, na.rm = TRUE)
mean_Fis_K2 <- cbind.data.frame(acutus, morelet)
mean_Fis_K2 

## Make combined datatables and save as .csv ##
## for K=3 
# Set row names for each data frame
rownames(mean_Hobs_K3) <- "Hobs"
rownames(mean_Hess_K3) <- "Hess"
rownames(mean_AR_K3) <- "AR"
rownames(mean_Fis_K3) <- "Fis"

# Combine the data frames with rbind
combined_df <- rbind(mean_Hobs_K3, mean_Hess_K3, mean_AR_K3, mean_Fis_K3)

# Now combined_df is a combined data frame where columns are "CA", "MCA", and "CM", and rows are the mean values for Hobs, Hess, AR, and Fis
print(combined_df)
write.csv(combined_df, file = paste0(out_dir,"/divstats_k3.csv"))

## for K=2
# Set row names for each data frame
rownames(mean_Hobs_K2) <- "Hobs"
rownames(mean_Hess_K2) <- "Hess"
rownames(mean_AR_K2) <- "AR"
rownames(mean_Fis_K2) <- "Fis"

# Combine the data frames with rbind
combined_df <- rbind(mean_Hobs_K2, mean_Hess_K2, mean_AR_K2, mean_Fis_K2)

# Now combined_df is a combined data frame where columns are "CA", "MCA", and "CM", and rows are the mean values for Hobs, Hess, AR, and Fis
print(combined_df)
write.csv(combined_df, file = paste0(out_dir,"/divstats_K2.csv"))

################################################################################
####################### Nucleotide diversity from VCFtools #####################
################################################################################
## calculate nucleotide diversity for each population across all sites, then 
## average across sites (I was lazy and did this in excel)

# vcftools \
# --vcf cuba-19-cmu.recode.vcf
# --out cuba-19-cmu-pi
# --site-pi
################################################################################
################################################################################
## From https://eacooper400.github.io/gen8900/exercises/tajd.html
################################################################################
## Load Libraries and Functions ##
library(devtools)
################################################################################
### Functions from:
## From https://eacooper400.github.io/gen8900/exercises/ld-2.html ##
################################################################################
### Calculate Allele Frequencies from Genotype Counts
allele.freq <- function(genotypeCounts) {
  n = sum(genotypeCounts) - genotypeCounts["NN"]
  p = ((2*genotypeCounts["AA"]) + genotypeCounts["Aa"])/(2*n)
  q = 1-p
  freqs = c(p,q)
  names(freqs) = c("p", "q")
  return(freqs)
}
### Calculate the LD parameter R-squared from 2 rows of a VCF file
calc_r2 <- function(row1, row2) {
  g1 = get.field(row1[10:length(row1)], row1[9], "GT")
  g2 = get.field(row2[10:length(row2)], row2[9], "GT")
  pA = unname(allele.freq(count.genotypes(g1))["p"])
  pB = unname(allele.freq(count.genotypes(g2))["p"])
  h = get.haplotypes(g1, g2)
  pAB = (length(h[h=="00"]))/(length(h))
  D = pAB - (pA * pB)
  rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))
  return(rsq)
}
### Count the AA, Aa, and aa genotypes in a sample
count.genotypes <- function(genotypes) {
  genotypes = gsub("(\\||/)", "", genotypes) 
  gen.patterns = c("00", "01", "10", "11", "..") 
  my.counts=table(factor(genotypes, levels=gen.patterns)) 
  final.counts = c(my.counts[1], (my.counts[2]+my.counts[3]), my.counts[4:5]) 
  names(final.counts) = c("AA", "Aa", "aa", "NN") 
  return(final.counts)
}
### Count the number of derived alleles for a VCF row
derivedCount <- function(row) {
  row=as.vector(row, mode="character")
  x=count.genotypes(get.field(row[10:length(row)], row[9], "GT"))
  dc=(2*x["aa"])+x["Aa"]
  return(unname(dc))
}
### Calculate Nucleotide Divergence (Dxy)
dxy <- function(vcf1, vcf2, perBP=TRUE) {
  g1=t(apply(vcf1, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  g2=t(apply(vcf2, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  af1=apply(g1, 1, function(x) allele.freq(count.genotypes(x))["p"])
  af2=apply(g2, 1, function(x) allele.freq(count.genotypes(x))["p"])
  ## Let x be the allele frequency (p) in pop1 * (1-p) in Pop2
  x = af1 * (1-af2)
  ## Let y be the allele frequency (p) in pop2 * (1-p) in Pop1
  y = af2 * (1-af1)
  dxy=sum((x+y))
  if (perBP) {
    c = unique(vcf1$CHROM)
    s = sapply(c, function(x,y) min(y[which(y$CHROM==x),2]), y=vcf1)
    e = sapply(c, function(x,y) max(y[which(y$CHROM==x),2]), y=vcf1)
    bp=sum(e-s)
    return(dxy/bp)
  } else { 
    return(dxy) 
  }
}
### Calculate Expected Het. and return, p, n, and H
expected.het <- function(genotypes) {
  obs.counts = count.genotypes(genotypes)
  n = sum(obs.counts) - obs.counts["NN"]
  freqs = allele.freq(obs.counts)
  Hexp = 2 * freqs[1] * freqs[2]
  res = c(freqs["p"], n, Hexp)
  res=as.numeric(unname(res))
  return(res)
}
### Determine which SNPs are Polymorphic vs Fixed in 2 species
fixed.poly <- function(vcf1, vcf2) {
  g1=t(apply(vcf1, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  g2=t(apply(vcf2, 1, function(x) get.field(x[10:length(x)], x[9], "GT")))
  af1=apply(g1, 1, function(x) allele.freq(count.genotypes(x))["p"])
  af2=apply(g2, 1, function(x) allele.freq(count.genotypes(x))["p"])
  res = rep("Polymorphic", nrow(g1))
  res[which(abs(af1-af2)==1)] = "Fixed"
  return(res)
}
### Get a Specified Field From a VCF Sample/Genotype String
get.field <- function(samples, format, fieldName) {
  x=strsplit(samples, split=":")
  fields=unlist(strsplit(format, split=":")) 
  i=which(fields==fieldName)
  if (!(fieldName %in% fields)) stop('fieldName not found in format fields') 
  return(sapply(x, `[[`, i)) 
}
### Get the HAPLOTYPES for a pair of genotype strings
get.haplotypes <- function(genotypes1, genotypes2) {
  a1 = gsub("\\|", "", genotypes1) 
  a2 = gsub("\\|", "", genotypes2)
  a1=unlist(strsplit(paste0(a1, collapse=""), split="")) 
  a2=unlist(strsplit(paste0(a2, collapse=""), split=""))
  haps = paste0(a1,a2)
  return(haps)
}
### Calculate minor allele frequency
maf <- function(vcf.row) {
  temp=as.vector(vcf.row, mode="character")
  af=allele.freq(count.genotypes(get.field(temp[10:length(temp)], temp[9], "GT")))
  maf=min(unname(af))
  return(maf)
}
### Calculate Nucleotide Diversity (Pi)
pi.diversity <- function(vcf, perBP=TRUE) {
  J=apply(vcf, 1, derivedCount)
  N=apply(vcf, 1, function(x) sum(count.genotypes(get.field(x[10:length(x)], x[9], "GT"))))
  C=2*N
  pi = sum((2*J*(C-J))/(C*(C-1)))
  if (perBP) {
    c = unique(vcf$CHROM)
    s = sapply(c, function(x,y) min(y[which(y$CHROM==x),2]), y=vcf)
    e = sapply(c, function(x,y) max(y[which(y$CHROM==x),2]), y=vcf)
    bp=sum(e-s)
    return(pi/bp)
  } else { return(pi) }
}	      
### Read in a VCF file as a table
read.vcf <- function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"))
}
### Calculate the variance for Tajima's D
variance.d <- function(n,S) {
  a1=sum(1/(seq(from=1, to=(n-1), by=1)))
  a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
  b1=(n+1)/(3*(n-1))
  b2=(2*((n**2)+n+3))/((9*n)*(n-1))
  c1=b1 - (1/a1)
  c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
  e1=c1/a1
  e2=c2/((a1**2)+a2)
  var=(e1*S) + (e2*S*(S-1))
  return(var)
}
### Calculate Waterson's theta
waterson.theta <- function(data, perBP=TRUE) {
  num.bp=nrow(data)
  check=apply(data, 1, FUN=maf)
  filter=data[which(check>0),]
  Sn=nrow(filter)
  n.ind=ncol(filter)-9
  i=seq(1, ((2*n.ind)-1))
  theta=Sn/sum(1/i)
  if (perBP) { return(theta/num.bp) }
  else { return(theta) }
}

################################################################################
## get individually subsetted VCF files w/ LDpruning_Dataconversion.R ##
setwd("/Filtering/filteredVCF/subset_data/")

## Load VCF files ## 
pure_acutus_vcf <- "noreponly_v2.75.renamed.pure_acutus.vcf.gz"
pure_morelet_vcf <- "noreponly_v2.75.renamed.pure_CM.vcf.gz"
pure_CA_vcf <- "noreponly_v2.75.renamed.pure_CA.vcf.gz"
pure_MCA_vcf <- "noreponly_v2.75.renamed.pure_MCA.vcf.gz"
################################################################################
## Acutus ##
my.data=read.vcf(pure_acutus_vcf, header=TRUE, stringsAsFactors=FALSE)
head(my.data, n=5)
dim(my.data)

## Morelets ##
my.data=read.vcf(pure_morelet_vcf, header=TRUE, stringsAsFactors=FALSE)
head(my.data, n=5)
dim(my.data)

## CA ##
my.data=read.vcf(pure_CA_vcf, header=TRUE, stringsAsFactors=FALSE)
head(my.data, n=5)
dim(my.data)

## MCA ## 
my.data=read.vcf(pure_MCA_vcf, header=TRUE, stringsAsFactors=FALSE)
head(my.data, n=5)
dim(my.data)

################################################################################
################################ Fst Stats #####################################
################################################################################
## https://popgen.nescent.org/DifferentiationSNP.html#observed-and-expected-heterozygosity-f_st ##
################################################################################
setwd(out_dir)

## Set up Genind objects ## 

############################
## K = 3 ;  CA, MCA, CM ##
############################
RAD1 <- RAD     # CA, MCA, CM

pop_info <- pop %>% dplyr::select(Sample, Abbrv_Monitoring.Unit, Morph_Species)
# Assuming RAD is your genind object and pop_info is your dataframe
# Extract the individual names from the genind object
ind_names_RAD <- indNames(RAD)

# Reorder pop_info to match the order of ind_names_RAD
pop_info_ordered <- pop_info %>% 
  mutate(Sample = factor(Sample, levels = ind_names_RAD)) %>% 
  arrange(Sample)

pop_info$Hierarchy <- with(pop_info, paste(MU, MSp, sep = "_"))

RAD@pop <- factor(pop_info$Hierarchy)

# Convert indNames to a data frame for easy manipulation
RAD_df <- data.frame(Sample = indNames(RAD), stringsAsFactors = FALSE)

# Join the additional information from pop based on matching Sample names
RAD_df <- left_join(RAD_df, pop, by = "Sample")

# Check for any unmatched samples if necessary
unmatched_samples <- RAD_df[is.na(RAD_df$MU) | is.na(RAD_df$MSp), ]

# Update the RAD object with the new hierarchy information
# Assuming RAD@other is a list, if not, you might need to initialize it
if(is.null(RAD@other)) RAD@other <- list()
RAD@other$MU <- RAD_df$Abbrv_Monitoring.Unit 
RAD@other$MSp <- RAD_df$Morph_Species
RAD
names(other(RAD))
strata(RAD) <- data.frame(other(RAD))
RAD

ind <- as.character(pop$Sample) # use later with adegenet (individual labels)
population <- as.character(population) # use later with adegenet (population labels)
MU <- pop$Abbrv_Monitoring.Unit 
MSp <- pop$Morph_Species

## Create new Hierfstat objects ##
RADhf <- genind2hierfstat(RAD) 

##########################################
## K = 2 ;  pure acutus, pure moreletii ##
##########################################
RAD2  

# Convert indNames to a data frame for easy manipulation
RAD_df2 <- data.frame(Sample = indNames(RAD2), stringsAsFactors = FALSE)

# Join the additional information from pop based on matching Sample names
RAD_df2 <- left_join(RAD_df2, pop_2, by = "Sample")

# Check for any unmatched samples if necessary
unmatched_samples <- RAD_df2[is.na(RAD_df2$Abbrv_Monitoring.Unit) | is.na(RAD_df2$Morph_Species), ]

# Update the RAD object with the new hierarchy information
# Assuming RAD@other is a list, if not, you might need to initialize it
if(is.null(RAD2@other)) RAD2@other <- list()
RAD2@other$MU <- RAD_df2$Abbrv_Monitoring.Unit 
RAD2@other$MSp <- RAD_df2$Morph_Species
RAD2
names(other(RAD2))
strata(RAD2) <- data.frame(other(RAD2))
RAD2

ind <- as.character(pop_2$Sample) # use later with adegenet (individual labels)
pop_2 <- as.character(pop_2) # use later with adegenet (population labels)
MU2 <- pop_2$Abbrv_Monitoring.Unit 
MSp2 <- pop_2$Morph_Species

## Create new Hierfstat objects ##
RADhf2 <- genind2hierfstat(RAD2) 


################################################################################
## Observed and expected heterozygosity: Fst
# Fst following Nei (1987) on genind object
wc(RAD) # Weir and Cockerham's estimate
  # $FST 0.09242777
  # $FIS 0.6121305

wc(RAD2)
  # $FST 0.08872248
  # $FIS 0.6102203

RAD
################################################################################
## Hierarchical Fst tests (=AMOVA for SNP dataset) ##

## Set up Pop data
loci <- RADhf[, -1] # Remove the population column

## produces a Hierarchical Fst (=AMOVA for SNPs or bi-allelic markers) 
##    It is possible to make permutations on the different levels
varcomp.glob(levels = data.frame(population), loci, diploid = TRUE) 

## tests the effect of the population on genetic differentiation. 
##    Individuals are randomly permuted among states. The states influence genetic 
##    differentiation at a 5% level. 
test.g(loci, level = population) 

## the MU are permuted among states. The states influence significantly genetic structuring.
test.between(loci, test.lev = population, rand.unit = MU, nperm = 100) 

## Pairwise Fst
genet.dist(RAD, method = "WC84")

################################################################################
################################################################################
################# Other Methods for Calculating F-statistics ###################
################################################################################
################################################################################
## From Physalia workshop script 
####################################################################################
## want to have data that is NOT LDpruned 
setwd(data_dir)
list.files()

library(ggplot2)
vcf <- read.vcfR(vcffile)
gl<-vcfR2genlight(vcf)
gl
gendata_names <- indNames(gl) # get the sample names
gendata_names
gl@pop <- as.factor(md$Morph_Species)
gl@pop


pops<-pops[pops$Sample %in% gendata_names,]
dim(pops)

# Obtained through vcftools - PopStats.sh #
CA_CM_Fst<-read.table("FST/fst_CACM_morph_win.windowed.weir.fst", header=T)
CA_HY_Fst<-read.table("FST/fst_CAHY_morph_win.windowed.weir.fst", header=T)
CA_MCA_Fst<-read.table("FST/fst_CAMCA_morph_win.windowed.weir.fst", header=T)
CM_HY_Fst<-read.table("FST/fst_CMHY_morph_win.windowed.weir.fst", header=T)
MCA_CM_Fst<-read.table("FST/fst_MCACM_morph_win.windowed.weir.fst", header=T)
MCA_HY_Fst<-read.table("FST/fst_MCAHY_morph_win.windowed.weir.fst", header=T)

fst <- CA_CM_Fst

# Visualising per SNP Fst
library(tidyverse)

# read in the file
Fst <- CA_CM_Fst
Fst

fst <- Fst
Fst$CHROM <- sub("^MDVP010000", "", Fst$CHROM)
Fst$CHROM

names(Fst)[1] <- "SCAFF"
Fst
# capdown the headers
names(fst) <- tolower(fst)
names(fst)[grep("fst")] <- "fst"
fst

ggplot(fst, aes(POS, WEIR_AND_COCKERHAM_FST)) + geom_point()

# Genome window
Fst$midPos<-(Fst$BIN_START+Fst$BIN_END)/2
head(Fst)

ggplot(Fst, aes(x=midPos/1000000, y=WEIGHTED_FST, col=SCAFF))+
  geom_point()+ geom_smooth()+ theme_classic()+
  facet_grid(cols = vars(SCAFF), scales = "free_x", space="free_x")+
  labs(  x = "position (in MB)")


######################################################
## 1.1 Pairwise differentiation between populations ##
######################################################
# Libraries
library(dplyr)
library(magrittr)
library(tibble)
library(gplots)
library(RColorBrewer)
library(corrplot)
# --------------

makeSymm <- function(m, position) {
  # Add symetrical triangle matrix (upper or lower)
  if (position == 'upper'){
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  if (position == 'lower'){
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
}

#-------------- fst matrix for all SNPs --------------------
fst.mat <- read.table("FST/FILE_fst_matrix.txt")
#use the function to fill the full matrix
fst.all.mat<- fst.mat %>%
  as.matrix(.) %>%
  makeSymm(., 'upper')
fst.all.mat[is.na(fst.all.mat)] <-  0 #replace NAs by 0 (NAs unaccepted for the heatmap function)
fst.all.mat[1:10,1:10] #check the fst_matrix

#visualise values
corrplot(fst.all.mat, is.corr = FALSE, method="number", addgrid.col = FALSE, diag=F, type="lower", number.digits = 3, number.cex=0.7)

#Visualize pairwise FST through a heatmap
gplots::heatmap.2(fst.all.mat, trace = 'none',
                  col= colorRampPalette(brewer.pal(9, "Reds"))(15),
                  key.xlab='FST')

## 1.2 Isolation by distance ##

library(reshape2)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)

#We switched SoDA to geosphere package.
if(!require(geosphere)){
  install.packages("geosphere")
  library(geosphere)
}


#import information about populations
info_pop <- read.table("documents/info_pop_geo_eco.txt", header=T)
head(info_pop)

#calculate geogrpahic (euclidian) distances between all pairs of populations
distance <- geosphere::distm(info_pop[,c(3,4)], fun=distGeo) %>%
  as.matrix(.)
#change it from meters to km
distance <- distance/1000

#set the colnames and rownames of the distance matrix
dimnames(distance) <- list(info_pop$pop,info_pop$pop)
distance

#prepare datasets
#linearize the distance matrix
dist.melt <- reshape2::melt(distance) %>%
  set_colnames(., c('pop1', 'pop2','distance'))
head(dist.melt)

#linearize the fst matrix
fst.melt <- reshape2::melt(fst.all.mat) %>%
  set_colnames(., c('pop1', 'pop2','FST'))

#join the distance and fst
IBD.df <- left_join(dist.melt, fst.melt, by=c('pop1','pop2')) %>%
  filter(., distance > 0)
head(IBD.df)

#test association with FST
cor.test(log(IBD.df$distance), IBD.df$FST/(1-IBD.df$FST))

#plot IBD
ggplot(IBD.df) + aes(x=log(distance), y=FST/(1-FST)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()

info_ind<-read.table("documents/info_samples.txt", header=T)
head(info_ind)
table (info_ind$pop,info_ind$sex)

################################################################################
####################### Basic genetic diversity stats ##########################
################################################################################
## From Alves et al 2023 - Data and R scripts ##
################################################################################
# needed packages
library("dartR")
library("adegenet")
library("hierfstat")
library("pegas")

options(scipen=999) # disable the scientific notations

# load filtered data
gl.pard <- readRDS("data/gl.filt_no_related_ind.rds")

# Convert genlight object to {adegenet} genind (genind objets store genetic data at an individual level)
gi.pard <- gl2gi(gl.pard, v=0)

# Create hierfstat object that will be needed for some analysis
df.pard <- genind2hierfstat(gi.pard) 

# Basic statistics with hierfstat. The function basic.stats() provides the observed heterozygosity (Ho), mean gene diversities within populations (Hs), Fis, and Fst. The function boot.ppfis() provides confidence interval for Fis. 

# calculate overall genetic diversity metrics 
basicstat <- basic.stats(df.pard, diploid = TRUE, digits = 2) 
names(basicstat)
basicstat$overall
summary(basicstat$Fis)# Calculate mean Fis for each populations

# calculate fis CIs
fis_ci <- boot.ppfis(df.pard, nboots = 1000)

# calculate genetic diversity with dartr for populations (created a table of it for the paper)
pop.div.stats <-gl.report.heterozygosity(gl.pard, method = "pop", verbose=3)
pop.div.stats


# calculates pairwise fst values based on the implementation in the StAMPP package. It allows to run bootstrap to estimate probability of fst values to be different from zero. 
# The fixation index (FST) is a measure of population differentiation due to genetic structure. Here I'll calculate for each pop (table in the manuscript)
fst.pop <- gl.fst.pop(gl.pard, nboots = 1000, percent = 95, nclusters = 1)
fst.pop$Fsts
fst.pop$Pvalues
?gl.fst.pop
fst.pop

