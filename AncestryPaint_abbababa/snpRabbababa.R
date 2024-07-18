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
library(snpR)
library(poppr)

################################################################################
########################### Set working directory ##############################
################################################################################
setwd("/datadir/filteredVCF_ab") 

list.files()

################################################################################
############################## Load data files #################################
################################################################################
# VCF file 
vcffile <- "data/noreponly_v2.filtered.85.renamed.vcf.gz"   # input 
vcf<-read.vcfR(vcffile) 
head(vcf@fix)
vcf@fix[,"REF"]
vcf@fix[,"ALT"]

genlight<-vcfR2genlight(vcf) 
genlight
ploidy(genlight) <- 2

## import genlight to snpR data file ##
dat1 <- import.snpR.data(genlight)
dat_vcf <- import.snpR.data(vcf)
dat_vcf
sample.meta(dat_vcf)
snp.meta(dat_vcf)["ID"]
snp.meta(dat_vcf)$snpID <- snp.meta(dat_vcf)["ID"]
snp.meta(dat_vcf)

x<-dat1
x
################################################################################
## Add snp metadata from genlight object ##
# Note that ref and anc columns are suggested in the SNP metadata, containing
# the derived (ref) and ancestral (anc) character states, respectively. 
sample.meta(x)
snp.meta(x)
maf <- get.snpR.stats(x, ".base", stats = "maf")$single
head(maf)
snp.meta(x)$ref <- maf$major
snp.meta(x)$anc <- maf$minor

## NOTE: try using the reference and alternate alleles from gds variant info for
## reference and ancestral alleles 
snp.meta(x)$ref <- vcf@fix[,"ALT"]
snp.meta(x)$anc <- vcf@fix[,"REF"]
#ref.count <- tabulate_allele_frequency_matrix(x, facets = "ref")

## Get chromosome info
snp.meta(x)$chr <- genlight@chromosome
snp.meta(x)$chr <- gsub("\\.", "_", snp.meta(x)$chr)
snp.meta(x)$chr <- gsub("MDVP010000", "", snp.meta(x)$chr)
snp.meta(x)$chr <- gsub("_1", "", snp.meta(x)$chr)

snp.meta(dat_vcf)$chr <- snp.meta(dat_vcf)["CHROM"]
snp.meta(dat_vcf)$chr <- gsub("\\.", "_", snp.meta(dat_vcf)$chr)
snp.meta(dat_vcf)$chr <- gsub("MDVP010000", "", snp.meta(dat_vcf)$chr)
snp.meta(dat_vcf)$chr <- gsub("_1", "", snp.meta(dat_vcf)$chr)
head(snp.meta(dat_vcf)$chr)

snp.meta(x)$position <- genlight@position
snp.meta(x)

################################################################################
## add population data ##
# Population data
md <- read_csv("admixture/popdata/popmap.csv")
md$Ad_cluster

## Subset populations 
# Morelets
CM_99 <- md %>% 
  filter(grepl("CM_99", Ad_cluster, ignore.case = TRUE))
CM_90 <- md %>% 
  filter(grepl("CM_90", Ad_cluster, ignore.case = TRUE))
Morelets <- md %>% 
  filter(Ad_cluster %in% c("CM_90", "CM_99"))
All_CM <- md %>% 
  filter(grepl("CM", Ad_cluster, ignore.case = TRUE))

# CA
CA_99 <- md %>%
  filter(grepl("CA_99", Ad_cluster, ignore.case = TRUE) & !grepl("MCA_99", Ad_cluster, ignore.case = TRUE))
CA_90 <- md %>%
  filter(grepl("CA_90", Ad_cluster, ignore.case = TRUE) & !grepl("MCA_90", Ad_cluster, ignore.case = TRUE))
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
Acutus_90 <- md %>%
  filter(Ad_acutus > 0.90)

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
Acutus_backcross <- md %>%
  filter(grepl("Acutus_backcross", Ad_cluster, ignore.case = TRUE))
Hybrid_25.75 <- md %>%
  filter(Ad_CM >= 0.25,Ad_CM <= 0.75)
################################################################################
## Create popmap for each test ##

## ACUTUS -> MORELET ##
# Test 1: CA -> CM 
P1 <- Morelets
P2 <- bind_rows(Hybrid_CM_backcross, F1_Hybrid)
P3 <- CA
test1.pop <- bind_rows(P1, P2, P3)
test1.pop$Sample
test1.pop <- test1.pop %>% select(Sample, Ad_cluster)
names(test1.pop) <- c("sampID", "pop")
test1.pop$pop <- gsub("_90", "", test1.pop$pop) 
test1.pop$pop <- gsub("_99", "", test1.pop$pop) 
test1.pop$pop <- gsub(".*Hybrid.*", "Hybrid", test1.pop$pop) 
test1.pop$pop
test1.ref <- tibble(sampID = "Reference", pop = "reference")
test1.pop <- add_row(test1.pop, test1.ref)
test1.pop

# Test 2: MCA -> CM 
P1 <- Morelets
P2 <- bind_rows(Hybrid_CM_backcross, F1_Hybrid)
P3 <- MCA
test2.pop <- bind_rows(P1, P2, P3)
test2.pop <- test2.pop %>% select(Sample, Ad_cluster)
names(test2.pop) <- c("sampID", "pop")
test2.pop$pop <- gsub("_90", "", test2.pop$pop) 
test2.pop$pop <- gsub("_99", "", test2.pop$pop) 
test2.pop$pop <- gsub(".*Hybrid.*", "Hybrid", test2.pop$pop) 
test2.pop$pop

# Test 3: Acutus -> CM 
P1 <- Morelets
P2 <- bind_rows(Hybrid_CM_backcross, F1_Hybrid)
P3 <- Acutus_90
test3.pop <- bind_rows(P1, P2, P3)
test3.pop <- test3.pop %>% select(Sample, Ad_cluster)
names(test3.pop) <- c("sampID", "pop")
test3.pop$pop <- gsub("_90", "", test3.pop$pop) 
test3.pop$pop <- gsub("_99", "", test3.pop$pop) 
test3.pop$pop <- gsub(".*Hybrid.*", "Hybrid", test3.pop$pop) 
test3.pop$pop

## MORELET -> ACUTUS ##
# Test 4: CM -> CA
P1 <- CA
P2 <- bind_rows(Hybrid_CA_backcross, F1_Hybrid)
P3 <- Morelets
test4.pop <- bind_rows(P1, P2, P3)
test4.pop <- test4.pop %>% select(Sample, Ad_cluster)
names(test4.pop) <- c("sampID", "pop")
test4.pop$pop <- gsub("_90", "", test4.pop$pop) 
test4.pop$pop <- gsub("_99", "", test4.pop$pop) 
test4.pop$pop <- gsub(".*Hybrid.*", "Hybrid", test4.pop$pop) 
test4.pop$pop

# Test 5: CM -> MCA
P1 <- MCA
P2 <- bind_rows(Hybrid_MCA_backcross, F1_Hybrid)
P3 <- Morelets
test5.pop <- bind_rows(P1, P2, P3)
test5.pop <- test5.pop %>% select(Sample, Ad_cluster)
names(test5.pop) <- c("sampID", "pop")
test5.pop$pop <- gsub("_90", "", test5.pop$pop) 
test5.pop$pop <- gsub("_99", "", test5.pop$pop) 
test5.pop$pop <- gsub(".*Hybrid.*", "Hybrid", test5.pop$pop) 
test5.pop 

# Test 6: CM -> ACUTUS
P1 <- Acutus_90
P2 <- bind_rows(Hybrid_CA_backcross, Hybrid_MCA_backcross, Acutus_backcross, F1_Hybrid)
P3 <- Morelets
test6.pop <- bind_rows(P1, P2, P3)
test6.pop <- test6.pop %>% select(Sample, Ad_cluster)
names(test6.pop) <- c("sampID", "pop")
test6.pop$pop <- gsub("_90", "", test6.pop$pop) 
test6.pop$pop <- gsub("_99", "", test6.pop$pop) 
test6.pop$pop <- gsub(".*Hybrid.*", "Hybrid", test6.pop$pop) 
test6.pop$pop <- gsub(".*backcross.*", "Hybrid", test6.pop$pop) #find things "_backcross" at end of string and replace with empty
test6.pop$pop

## Non-admixed population geneflow ##

# Test 7: CM -> MCA
P1 <- CA
P2 <- MCA
P3 <- Morelets
test7.pop <- bind_rows(P1, P2, P3)
test7.pop <- test7.pop %>% select(Sample, Ad_cluster)
names(test7.pop) <- c("sampID", "pop")
test7.pop$pop <- gsub("_90", "", test7.pop$pop) 
test7.pop$pop <- gsub("_99", "", test7.pop$pop) 
test7.pop$pop

# Test 8: CA -> MCA 
P1 <- MCA
P2 <- Acutus_backcross
P3 <- CA
test8.pop <- bind_rows(P1, P2, P3)
test8.pop <- test8.pop %>% select(Sample, Ad_cluster)
names(test8.pop) <- c("sampID", "pop")
test8.pop$pop <- gsub("_90", "", test8.pop$pop) 
test8.pop$pop <- gsub("_99", "", test8.pop$pop) 
test8.pop$pop

# Test 9: MCA -> CM
P1 <- CM_99
P2 <- CM_90
P3 <- MCA
test9.pop <- bind_rows(P1, P2, P3)
test9.pop <- test9.pop %>% select(Sample, Ad_cluster)
names(test9.pop) <- c("sampID", "pop")
test9.pop$pop <- gsub("MCA_90", "MCA", test9.pop$pop) 
test9.pop$pop <- gsub("MCA_99", "MCA", test9.pop$pop) 
#test9.pop$pop <- gsub("_90", "", test9.pop$pop) 
#test9.pop$pop <- gsub("_99", "", test9.pop$pop) 
test9.pop$pop

# Test 9b: MCA -> CM: CM_NRW vs CM_notNRW
## Subset CM populations by Monitoring Unit
#CM_NRW <- All_CM %>% 
#  filter(grepl("NRW", Abbrv_Monitoring.Unit, ignore.case = TRUE))
#CM_notNRW <- All_CM %>% 
#  filter(!grepl("NRW", Abbrv_Monitoring.Unit, ignore.case = TRUE))
CM_NRW <- Morelets %>% 
  filter(grepl("NRW", Abbrv_Monitoring.Unit, ignore.case = TRUE))
CM_notNRW <- Morelets %>% 
  filter(!grepl("NRW", Abbrv_Monitoring.Unit, ignore.case = TRUE))

## Change the "Ad_cluster" pops for the CM
CM_NRW <- CM_NRW %>%
  mutate(Ad_cluster = "CM_NRW")
CM_notNRW <- CM_notNRW %>%
  mutate(Ad_cluster = "CM_notNRW")

#P1 <- CM_NRW 
#P2 <- CM_notNRW
#P3 <- MCA
#test9b.pop <- bind_rows(P1, P2, P3)
#test9b.pop <- test9b.pop %>% select(Sample, Ad_cluster)
#names(test9b.pop) <- c("sampID", "pop")
## Change the MCA pops to be uniform
#test9b.pop$pop <- gsub("MCA_90", "MCA", test9b.pop$pop) 
#test9b.pop$pop <- gsub("MCA_99", "MCA", test9b.pop$pop) 
#test9b.pop$pop

## changed temporarily to test9
P1 <- CM_NRW 
P2 <- CM_notNRW
P3 <- MCA
test9.pop <- bind_rows(P1, P2, P3)
test9.pop <- test9.pop %>% select(Sample, Ad_cluster)
names(test9.pop) <- c("sampID", "pop")
## Change the MCA pops to be uniform
test9.pop$pop <- gsub("MCA_90", "MCA", test9.pop$pop) 
test9.pop$pop <- gsub("MCA_99", "MCA", test9.pop$pop) 
test9.pop$pop

# Test 10: MCA -> CA
P1 <- CA
P2 <- Acutus_backcross
P3 <- MCA
test10.pop <- bind_rows(P1, P2, P3)
test10.pop <- test10.pop %>% select(Sample, Ad_cluster)
names(test10.pop) <- c("sampID", "pop")
test10.pop$pop <- gsub("_90", "", test10.pop$pop) 
test10.pop$pop <- gsub("_99", "", test10.pop$pop) 
test10.pop$pop

# Test 11: CA -> CM
#P1 <- CM_99
#P2 <- CM_90
#P3 <- CA
#test11.pop <- bind_rows(P1, P2, P3)
#test11.pop <- test11.pop %>% select(Sample, Ad_cluster)
#names(test11.pop) <- c("sampID", "pop")
#test11.pop$pop <- gsub("CA_90", "CA", test11.pop$pop) 
#test11.pop$pop <- gsub("CA_99", "CA", test11.pop$pop) 
#test11.pop$pop <- gsub("_90", "", test11.pop$pop) 
#test11.pop$pop <- gsub("_99", "", test11.pop$pop) 
#test11.pop$pop

P1 <- CM_NRW 
P2 <- CM_notNRW
P3 <- CA
test11.pop <- bind_rows(P1, P2, P3)
test11.pop <- test11.pop %>% select(Sample, Ad_cluster)
names(test11.pop) <- c("sampID", "pop")
test11.pop$pop <- gsub("CA_90", "CA", test11.pop$pop) 
test11.pop$pop <- gsub("CA_99", "CA", test11.pop$pop) 
#test11.pop$pop <- gsub("_90", "", test11.pop$pop) 
#test11.pop$pop <- gsub("_99", "", test11.pop$pop) 
test11.pop$pop

# Test 12: CM -> CA
P1 <- MCA
P2 <- CA
P3 <- Morelets
test12.pop <- bind_rows(P1, P2, P3)
test12.pop <- test12.pop %>% select(Sample, Ad_cluster)
names(test12.pop) <- c("sampID", "pop")
test12.pop$pop <- gsub("_90", "", test12.pop$pop) 
test12.pop$pop <- gsub("_99", "", test12.pop$pop) 
test12.pop$pop


# Test 13: MCA -> CM: CM_CF vs CM_notCF
## Subset CM populations by Monitoring Unit
CM_CF <- All_CM %>% 
  filter(grepl("CF", Abbrv_Monitoring.Unit, ignore.case = TRUE))
CM_notCF <- All_CM %>% 
  filter(!grepl("CF", Abbrv_Monitoring.Unit, ignore.case = TRUE))

## Change the "Ad_cluster" pops for the CM
CM_CF <- CM_CF %>%
  mutate(Ad_cluster = "CM_CF")
CM_notCF <- CM_notCF %>%
  mutate(Ad_cluster = "CM_notCF")

P1 <- CM_CF 
P2 <- CM_notCF
P3 <- MCA
test13.pop <- bind_rows(P1, P2, P3)
test13.pop <- test13.pop %>% select(Sample, Ad_cluster)
names(test13.pop) <- c("sampID", "pop")
## Change the MCA pops to be uniform
test13.pop$pop <- gsub("MCA_90", "MCA", test13.pop$pop) 
test13.pop$pop <- gsub("MCA_99", "MCA", test13.pop$pop) 
test13.pop$pop


# Test 14: CA -> CM: CM_CF vs CM_notCF
## Subset CM populations by Monitoring Unit
P1 <- CM_CF 
P2 <- CM_notCF
P3 <- CA
test14.pop <- bind_rows(P1, P2, P3)
test14.pop <- test14.pop %>% select(Sample, Ad_cluster)
names(test14.pop) <- c("sampID", "pop")
## Change the MCA pops to be uniform
test14.pop$pop <- gsub("CA_90", "CA", test14.pop$pop) 
test14.pop$pop <- gsub("CA_99", "CA", test14.pop$pop) 
test14.pop$pop

# Test 15: acutus -> CM: CM_CF vs CM_notCF
## Subset CM populations by Monitoring Unit
P1 <- CM_CF 
P2 <- CM_notCF
P3 <- Acutus_90
test15.pop <- bind_rows(P1, P2, P3)
test15.pop <- test15.pop %>% select(Sample, Ad_cluster)
names(test15.pop) <- c("sampID", "pop")
## Change the MCA pops to be uniform
test15.pop$pop <- gsub("CA_90", "CA", test15.pop$pop) 
test15.pop$pop <- gsub("CA_99", "CA", test15.pop$pop) 
test15.pop$pop <- gsub(".*CA.*", "acutus", test15.pop$pop) 
test15.pop$pop
################################################################################
####################
## Other Popmaps ##
####################
## All samples including backcross and hybrids
#pop <- read_tsv("admixture/popdata/ad_cluster.pop1.txt",col_names = F) 
#names(pop) <- c("sampID", "pop")
#pop$sampID
#pop$pop

## pure (90%) samples only
#pop_pure <- read_tsv("admixture/popdata/unadmixedpops.txt",col_names = F) 
#names(pop_pure) <- c("sampID", "pop")
#pop_pure$sampID
#pop_pure$pop

# replace string to rewrite names 
#pop_pure$pop <- gsub("_90", "", pop_pure$pop) 
#pop_pure$pop <- gsub("_99", "", pop_pure$pop) 

################################################################################
## For Loop that subsets snpR object by populations 

# Create a list to store subsetted snpR objects
subsetted_objects <- list()

# Loop through each test
for (i in 1:12) {
  dat <- x
  # view sample metadata
  head(sample.meta(dat))
  sampID <- sample.meta(dat)$sampID
  sampID
  # Match and subset based on sampID
  # Subset df1 based on matching values in df2
  pop_name <- paste0("test",i ,".pop")
  pop <- get(pop_name)
  subset_df <- sample.meta(dat)[sample.meta(dat)$sampID %in% pop$sampID, ]
  row_numbers <- as.integer(rownames(subset_df))
  
  # New subsetted snpR object by test
  subsetted_objects[[i]] <- dat[, row_numbers]
  
  # Print or save the subsetted object as needed
  # You can customize this part based on your requirements
  print(paste("Test", i, "subsetted object created."))

}

# Access the subsetted snpR objects using subsetted_objects[[1]], subsetted_objects[[2]], etc.
# make separate snpR objects per subetted object 
(test1<-subsetted_objects[[1]])
(test2<-subsetted_objects[[2]])
(test3<-subsetted_objects[[3]])
(test4<-subsetted_objects[[4]])
(test5<-subsetted_objects[[5]])
(test6<-subsetted_objects[[6]])
(test7<-subsetted_objects[[7]])
(test8<-subsetted_objects[[8]])
(test9<-subsetted_objects[[9]])
(test10<-subsetted_objects[[10]])
(test11<-subsetted_objects[[11]])
(test12<-subsetted_objects[[12]])

## Doing separately for the test9b ##
dat <- x
# view sample metadata
head(sample.meta(dat))
sampID <- sample.meta(dat)$sampID
sampID
# Match and subset based on sampID
# Subset df1 based on matching values in df2
pop_name <- paste0("test","9b" ,".pop")
pop <- get(pop_name)
subset_df <- sample.meta(dat)[sample.meta(dat)$sampID %in% pop$sampID, ]
row_numbers <- as.integer(rownames(subset_df))

# New subsetted snpR object by test
test9b <- dat[, row_numbers]
test9b

## Test 13-15
# Create a list to store subsetted snpR objects
subsetted_objects <- list()

# Loop through each test
for (i in 13:15) {
  dat <- x
  # view sample metadata
  head(sample.meta(dat))
  sampID <- sample.meta(dat)$sampID
  sampID
  # Match and subset based on sampID
  # Subset df1 based on matching values in df2
  pop_name <- paste0("test",i ,".pop")
  pop <- get(pop_name)
  subset_df <- sample.meta(dat)[sample.meta(dat)$sampID %in% pop$sampID, ]
  row_numbers <- as.integer(rownames(subset_df))
  
  # New subsetted snpR object by test
  subsetted_objects[[i]] <- dat[, row_numbers]
  
  # Print or save the subsetted object as needed
  # You can customize this part based on your requirements
  print(paste("Test", i, "subsetted object created."))
  
}
subsetted_objects

(test13<-subsetted_objects[[13]])
(test14<-subsetted_objects[[14]])
(test15<-subsetted_objects[[15]])
################################################################################
############################### abbababa test ##################################
################################################################################
## An excess of either ABBA or BABA, resulting in a D-statistic that is significantly 
## different from zero, is indicative of gene flow between two taxa.
## A positive D-statistic (i.e. an excess of ABBA) points to introgression between P2 and P3, 
## whereas a negative D-statistic (i.e. an excess of BABA) points to introgression between P1 and P3.
################################################################################
## input population data from each pop and input into respective snpR object ##
################################################################################
# Test 1: CA -> CM 
## input pop info from tibble
(sample.meta(test1) <- sample.meta(test1) %>% left_join(test1.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test1)$sampID <- gsub("\\.", "_", sample.meta(test1)$sampID)
summarize_facets(test1)
sample.meta(test1)$pop

## Run Abba-baba test on pops ##
# run with jackknifing, 1000kb windows!
# Test if CA or MCA have more geneflow with CM...
#x <- calc_abba_baba(x, "pop", "CA", "MCA","Morelets", TRUE, sigma = 1000)
#test1 <- calc_abba_baba(test1, "pop", "CM", "Hybrid","CA", "reference", jackknife = TRUE, sigma = 1000)
test1.pop <- calc_abba_baba(test1, "pop", "CM", "Hybrid","CA", jackknife = TRUE, sigma = 1000)

get.snpR.stats(test1.pop, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test1.pop, "pop", "abba_baba") # gets the overall pop results

# smoothed windowed averages
test1.pop <- calc_smoothed_averages(test1.pop, "pop", sigma = 200, step = 200, 
                            nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test1.pop, "pop", "abba_baba"))

test1.pairwise <- get.snpR.stats(test1.pop, "pop", "abba_baba")$pairwise
test1.pairwise.window <- get.snpR.stats(test1.pop, "pop", "abba_baba")$pairwise.window
test1.weighted.means <- get.snpR.stats(test1.pop, "pop", "abba_baba")$weighted.means

get.snpR.stats(test1.pop, "pop.chr", "abba_baba")$pairwise
get.snpR.stats(test1.pop, "pop.chr", "abba_baba")$pairwise.window
get.snpR.stats(test1.pop, "pop.chr", "abba_baba")$weighted.means

## Save pop-based abba-baba results as .csv 
(abba_baba_pairwise.stats <- get.snpR.stats(test1, "pop", "abba_baba")$pairwise)
filename <- paste0("abbababa/out/", "test1_", "abba_baba_pairwise",".csv") 
write_csv(test1.pairwise, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test1, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test1_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(test1.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test1, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test1_", "pairwise.window", ".csv") 
write_csv(test1.pairwise.window , file = filename)

## Run Abba-baba test on pops.chr ##
test1 <- calc_abba_baba(test1, "pop.chr", "CM", "Hybrid","CA", jackknife = TRUE, sigma = 1000) # sigma = 1000 uses 1000kb windows

get.snpR.stats(test1, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test1, "pop", "abba_baba") # gets the overall pop results

# smoothed windowed averages for pop.chr
test1 <- calc_smoothed_averages(test1, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
get.snpR.stats(test1, "pop.chr", "abba_baba")
names(get.snpR.stats(test1, "pop.chr", "abba_baba"))

## Save per chr abba-baba results as .csv 
(abba_baba_pairwise.stats <- get.snpR.stats(test1, "pop.chr", "abba_baba")$pairwise)
filename <- paste0("abbababa/out/", "test1_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test1, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test1_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test1, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test1_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
# Test 2: MCA -> CM 
# input pop info from tibble
(sample.meta(test2) <- sample.meta(test2) %>% left_join(test2.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test2)$sampID <- gsub("\\.", "_", sample.meta(test2)$sampID)
summarize_facets(test2)
sample.meta(test2)$pop

## Run Abba-baba test on pops ##
# run with jackknifing, 1000kb windows!
# Test if CA or MCA have more geneflow with CM...
#x <- calc_abba_baba(x, "pop", "CA", "MCA","Morelets", TRUE, sigma = 1000)
test2 <- calc_abba_baba(test2, "pop", "CM", "Hybrid","MCA", TRUE, sigma = 1000)

get.snpR.stats(test2, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test2, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test2 <- calc_smoothed_averages(test2, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test2, "pop", "abba_baba"))

(abba_baba_pairwise.stats <- get.snpR.stats(test2, "pop", "abba_baba")$pairwise)
filename <- paste0("abbababa/out/", "test2_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test2, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test2_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test2, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test2_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

weighted.means.alltests <- rbind(test1.weighted.means, abba_baba_pairwise.weighted.means)

## Run Abba-baba test on pops.chr ##
test2 <- calc_abba_baba(test2, "pop.chr", "CM", "Hybrid","MCA", TRUE, sigma = 1000)

get.snpR.stats(test2, "pop.chr", "abba_baba") # gets the per chr results
#get.snpR.stats(test2, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test2 <- calc_smoothed_averages(test2, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test2, "pop", "abba_baba"))
names(get.snpR.stats(test2, "pop.chr", "abba_baba"))

(abba_baba_pairwise.stats <- get.snpR.stats(test2, "pop.chr", "abba_baba")$pairwise)
filename <- paste0("abbababa/out/", "test2_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test2, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test2_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test2, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test2_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
# Test 3: Acutus -> CM 
test3.pop$pop 
test3.pop$pop <- gsub(".*CA.*", "Acutus", test3.pop$pop) 

# input pop info from tibble
(sample.meta(test3) <- sample.meta(test3) %>% left_join(test3.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test3)$sampID <- gsub("\\.", "_", sample.meta(test3)$sampID)
summarize_facets(test3)
sample.meta(test3)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
# Test if CA or MCA have more geneflow with CM...
#x <- calc_abba_baba(x, "pop", "CA", "MCA","Morelets", TRUE, sigma = 1000)
test3.pop <- calc_abba_baba(test3, "pop", "CM", "Hybrid","Acutus", TRUE, sigma = 1000)

get.snpR.stats(test3.pop, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test3.pop, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test3.pop <- calc_smoothed_averages(test3.pop, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test3.pop, "pop", "abba_baba"))

(abba_baba_pairwise.stats <- get.snpR.stats(test3.pop, "pop", "abba_baba")$pairwise)
filename <- paste0("abbababa/out/", "test3_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test3.pop, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test3_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test3.pop, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test3_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.weighted.means 

## for chrm ##
test3.chr <- calc_abba_baba(test3, "pop.chr", "CM", "Hybrid","Acutus", TRUE, sigma = 1000)

get.snpR.stats(test3.chr, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test3.chr, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test3.chr <- calc_smoothed_averages(test3.chr, "pop", sigma = 200, step = 200, 
                                    nk = TRUE, stats.type = "pairwise")
test3.chr <- calc_smoothed_averages(test3.chr, "pop.chr", sigma = 200, step = 200, 
                                    nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test3.chr, "pop.chr", "abba_baba"))

get.snpR.stats(test3.chr, "pop", "abba_baba")$weighted.means
abba_baba_pairwise.weighted.means 

(abba_baba_pairwise.stats <- get.snpR.stats(test3.chr, "pop.chr", "abba_baba")$pairwise)
filename <- paste0("abbababa/out/", "test3_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test3.chr, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test3_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test3.chr, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test3_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
# Test 4: CM -> CA
# input pop info from tibble
(sample.meta(test4) <- sample.meta(test4) %>% left_join(test4.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test4)$sampID <- gsub("\\.", "_", sample.meta(test4)$sampID)
summarize_facets(test4)
sample.meta(test4)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
# Test if CA or MCA have more geneflow with CM...
#x <- calc_abba_baba(x, "pop", "CA", "MCA","Morelets", TRUE, sigma = 1000)
test4 <- calc_abba_baba(test4, "pop.chr", "CA", "Hybrid","CM", TRUE, sigma = 1000)

get.snpR.stats(test4, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test4, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test4 <- calc_smoothed_averages(test4, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test4 <- calc_smoothed_averages(test4, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test4, "pop", "abba_baba"))
names(get.snpR.stats(test4, "pop.chr", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test4, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test4_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test4, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test4_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test4, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test4_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test4, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test4_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test4, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test4_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test4, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test4_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
# Test 5: CM -> MCA
# input pop info from tibble
(sample.meta(test5) <- sample.meta(test5) %>% left_join(test5.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test5)$sampID <- gsub("\\.", "_", sample.meta(test5)$sampID)
summarize_facets(test5)
sample.meta(test5)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
# Test if CA or MCA have more geneflow with CM...
#x <- calc_abba_baba(x, "pop", "CA", "MCA","Morelets", TRUE, sigma = 1000)
test5 <- calc_abba_baba(test5, "pop.chr", "MCA", "Hybrid","CM", TRUE, sigma = 1000)

get.snpR.stats(test5, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test5, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test5 <- calc_smoothed_averages(test5, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test5 <- calc_smoothed_averages(test5, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test5, "pop", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test5, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test5_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test5, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test5_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test5, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test5_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test5, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test5_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test5, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test5_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test5, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test5_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
# Test 6: CM -> ACUTUS
test6.pop$pop <- gsub(".*CA.*", "Acutus", test6.pop$pop) 
# Remove duplicated rows
test6.pop <- distinct(test6.pop)
(table(test6.pop$pop))

#write_delim(test6.pop, file = "admixture/popdata/acutuspopmap.txt")
# input pop info from tibble
(sample.meta(test6) <- sample.meta(test6) %>% left_join(test6.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test6)$sampID <- gsub("\\.", "_", sample.meta(test6)$sampID)
summarize_facets(test6)
sample.meta(test6)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
# Test if CA or MCA have more geneflow with CM...
#x <- calc_abba_baba(x, "pop", "CA", "MCA","Morelets", TRUE, sigma = 1000)
test6 <- calc_abba_baba(test6, "pop.chr", "Acutus", "Hybrid","CM", TRUE, sigma = 1000)

get.snpR.stats(test6, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test6, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test6 <- calc_smoothed_averages(test6, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test6 <- calc_smoothed_averages(test6, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test6, "pop", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test6, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test6_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test6, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test6_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test6, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test6_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test6, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test6_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test6, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test6_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test6, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test6_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
######################## For non-admixed populations ###########################
################################################################################
# Test 7: CM -> MCA
test7.pop$pop 
# Remove duplicated rows
test7.pop <- distinct(test7.pop)
(table(test7.pop$pop))
#write_tsv(test7.pop, file="data/pops/test7pop.txt")

# input pop info from tibble
(sample.meta(test7) <- sample.meta(test7) %>% left_join(test7.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test7)$sampID <- gsub("\\.", "_", sample.meta(test7)$sampID)
summarize_facets(test7)
sample.meta(test7)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
# Test if CA or MCA have more geneflow with CM...
#x <- calc_abba_baba(x, "pop", "CA", "MCA","Morelets", TRUE, sigma = 1000)
test7 <- calc_abba_baba(test7, "pop.chr", "CA", "MCA","CM", TRUE, sigma = 1000)

get.snpR.stats(test7, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test7, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test7 <- calc_smoothed_averages(test7, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test7 <- calc_smoothed_averages(test7, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test7, "pop", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test7, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test7_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test7, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test7_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test7, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test7_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test7, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test7_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test7, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test7_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test7, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test7_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

################################################################################
# Test 8: CA -> MCA 
test8.pop$pop 
# Remove duplicated rows
test8.pop <- distinct(test8.pop)
(table(test8.pop$pop))
#write_tsv(test8.pop, file="data/pops/test8pop.txt")

# input pop info from tibble
(sample.meta(test8) <- sample.meta(test8) %>% left_join(test8.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test8)$sampID <- gsub("\\.", "_", sample.meta(test8)$sampID)
summarize_facets(test8)
sample.meta(test8)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
# Test if CA or MCA have more geneflow with CM...
test8 <- calc_abba_baba(test8, "pop.chr", "MCA", "Acutus_backcross","CA", TRUE, sigma = 1000)

get.snpR.stats(test8, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test8, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test8 <- calc_smoothed_averages(test8, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test8 <- calc_smoothed_averages(test8, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test8, "pop", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test8, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test8_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test8, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test8_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test8, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test8_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test8, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test8_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test8, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test8_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test8, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test8_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

################################################################################
# Test 9: MCA -> CM
test9.pop$pop 
# Remove duplicated rows
test9.pop <- distinct(test9.pop)
(table(test9.pop$pop))
#write_tsv(test9.pop, file="data/pops/test9pop.txt")

# input pop info from tibble
(sample.meta(test9) <- sample.meta(test9) %>% left_join(test9.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test9)$sampID <- gsub("\\.", "_", sample.meta(test9)$sampID)
summarize_facets(test9)
sample.meta(test9)$pop
## changed temporarily to test9

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
#test9 <- calc_abba_baba(test9, "pop.chr", "CM_99","CM_90","MCA", jackknife=TRUE, sigma = 1000)
test9 <- calc_abba_baba(test9, "pop.chr", "CM_NRW","CM_notNRW","MCA", jackknife=TRUE, sigma = 1000)

get.snpR.stats(test9, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test9, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test9 <- calc_smoothed_averages(test9, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test9 <- calc_smoothed_averages(test9, "pop.chr", sigma = 200, step = 200, 
                                 nk = TRUE, stats.type = "pairwise")

names(get.snpR.stats(test9, "pop", "abba_baba"))
names(get.snpR.stats(test9, "pop.chr", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test9, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test9_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test9, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test9_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test9, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test9_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test9, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test9_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test9, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test9_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test9, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test9_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

################################################################################
# Test 10: MCA -> CA
test10.pop$pop 
# Remove duplicated rows
test10.pop <- distinct(test10.pop)
(table(test10.pop$pop))
#write_tsv(test10.pop, file="data/pops/test10pop.txt")

# input pop info from tibble
(sample.meta(test10) <- sample.meta(test10) %>% left_join(test10.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test10)$sampID <- gsub("\\.", "_", sample.meta(test10)$sampID)
summarize_facets(test10)
sample.meta(test10)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
test10 <- calc_abba_baba(test10, "pop.chr", "CA", "Acutus_backcross", "MCA", TRUE, sigma = 1000)

get.snpR.stats(test10, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test10, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test10 <- calc_smoothed_averages(test10, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test10 <- calc_smoothed_averages(test10, "pop.chr", sigma = 200, step = 200, 
                                 nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test10, "pop", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test10, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test10_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test10, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test10_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test10, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test10_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test10, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test10_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test10, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test10_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test10, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test10_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
# Test 11: CA -> CM
test11.pop$pop

# Remove duplicated rows
test11.pop <- distinct(test11.pop)
(table(test11.pop$pop))
#write_tsv(test11.pop, file="data/pops/test11pop.txt")

# input pop info from tibble
(sample.meta(test11) <- sample.meta(test11) %>% left_join(test11.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test11)$sampID <- gsub("\\.", "_", sample.meta(test11)$sampID)
summarize_facets(test11)
sample.meta(test11)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
#test11 <- calc_abba_baba(test11, "pop.chr", "CM_99","CM_90","CA", TRUE, sigma = 1000)
test11 <- calc_abba_baba(test11, "pop.chr", "CM_NRW","CM_notNRW","CA", TRUE, sigma = 1000)

get.snpR.stats(test11, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test11, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test11 <- calc_smoothed_averages(test11, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test11 <- calc_smoothed_averages(test11, "pop.chr", sigma = 200, step = 200, 
                                 nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test11, "pop", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test11, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test11_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test11, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test11_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test11, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test11_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test11, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test11_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test11, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test11_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test11, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test11_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
# Test 12: CM -> CA
test12.pop$pop 
# Remove duplicated rows
test12.pop <- distinct(test12.pop)
(table(test12.pop$pop))
#write_tsv(test12.pop, file="abbababa/data/pops/test12pop.txt")

# input pop info from tibble
(sample.meta(test12) <- sample.meta(test12) %>% left_join(test12.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test12)$sampID <- gsub("\\.", "_", sample.meta(test12)$sampID)
summarize_facets(test12)
sample.meta(test12)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
test12 <- calc_abba_baba(test12, "pop.chr", "MCA", "CA","CM", TRUE, sigma = 1000)

get.snpR.stats(test12, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test12, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test12 <- calc_smoothed_averages(test12, "pop", sigma = 200, step = 200, 
                                     nk = TRUE, stats.type = "pairwise")

test12 <- calc_smoothed_averages(test12, "pop.chr", sigma = 200, step = 200, 
                                     nk = TRUE, stats.type = "pairwise")

test12.pop <- calc_abba_baba(test12, "pop", "MCA", "CA","CM", TRUE, sigma = 1000)
get.snpR.stats(test12.pop, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test12.pop, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test12.pop <- calc_smoothed_averages(test12.pop, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")

test12.pop <- calc_smoothed_averages(test12.pop, "pop.chr", sigma = 200, step = 200, 
                                 nk = TRUE, stats.type = "pairwise")

names(get.snpR.stats(test12, "pop", "abba_baba"))
get.snpR.stats(test12.pop, "pop", "abba_baba")$weighted.means



abba_baba_pairwise.stats <- get.snpR.stats(test12, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test12_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test12, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test12_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test12, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test12_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test12, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test12_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test12, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test12_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test12, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test12_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)


################################################################################
# Test 13: MCA -> CM_CF
test13.pop$pop 
# Remove duplicated rows
test13.pop <- distinct(test13.pop)
(table(test13.pop$pop))
#write_tsv(test13.pop, file="data/pops/test13pop.txt")

# input pop info from tibble
(sample.meta(test13) <- sample.meta(test13) %>% left_join(test13.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test13)$sampID <- gsub("\\.", "_", sample.meta(test13)$sampID)
summarize_facets(test13)
sample.meta(test13)$pop

## Run Abba-baba test and save to file ##
test13 <- calc_abba_baba(test13, "pop.chr", "CM_CF","CM_notCF","MCA", jackknife=TRUE, sigma = 200)

get.snpR.stats(test13, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test13, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test13 <- calc_smoothed_averages(test13, "pop", sigma = 100, step = 100, 
                                nk = TRUE, stats.type = "pairwise")
test13 <- calc_smoothed_averages(test13, "pop.chr", sigma = 100, step = 100, 
                                nk = TRUE, stats.type = "pairwise")

names(get.snpR.stats(test13, "pop", "abba_baba"))
names(get.snpR.stats(test13, "pop.chr", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test13, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test13_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test13, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test13_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test13, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test13_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test13, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test13_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test13, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test13_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test13, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test13_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

################################################################################
# Test 14: CA -> CM_CF
test14.pop$pop 
# Remove duplicated rows
test14.pop <- distinct(test14.pop)
(table(test14.pop$pop))
#write_tsv(test14.pop, file="data/pops/test14pop.txt")

# input pop info from tibble
(sample.meta(test14) <- sample.meta(test14) %>% left_join(test14.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test14)$sampID <- gsub("\\.", "_", sample.meta(test14)$sampID)
summarize_facets(test14)
sample.meta(test14)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
test14 <- calc_abba_baba(test14, "pop.chr", "CM_CF","CM_notCF","CA", jackknife=TRUE, sigma = 1000)

get.snpR.stats(test14, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test14, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test14 <- calc_smoothed_averages(test14, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test14 <- calc_smoothed_averages(test14, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")

names(get.snpR.stats(test14, "pop", "abba_baba"))
names(get.snpR.stats(test14, "pop.chr", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test14, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test14_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test14, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test14_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test14, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test14_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test14, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test14_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test14, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test14_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test14, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test14_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)


################################################################################
# Test 15: acutus -> CM_CF
test15.pop$pop 
# Remove duplicated rows
test15.pop <- distinct(test15.pop)
(table(test15.pop$pop))
#write_tsv(test15.pop, file="data/pops/test15pop.txt")

# input pop info from tibble
(sample.meta(test15) <- sample.meta(test15) %>% left_join(test15.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test15)$sampID <- gsub("\\.", "_", sample.meta(test15)$sampID)
summarize_facets(test15)
sample.meta(test15)$pop
## changed temporarily to test15

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
#test15 <- calc_abba_baba(test15, "pop.chr", "CM_99","CM_90","MCA", jackknife=TRUE, sigma = 1000)
test15 <- calc_abba_baba(test15, "pop.chr", "CM_CF","CM_notCF","acutus", jackknife=TRUE, sigma = 1000)

get.snpR.stats(test15, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test15, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test15 <- calc_smoothed_averages(test15, "pop", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")
test15 <- calc_smoothed_averages(test15, "pop.chr", sigma = 200, step = 200, 
                                nk = TRUE, stats.type = "pairwise")

names(get.snpR.stats(test15, "pop", "abba_baba"))
names(get.snpR.stats(test15, "pop.chr", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test15, "pop", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test15_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test15, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test15_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test15, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test15_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test15, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test15_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test15, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test15_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test15, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test15_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

################################################################################

#########################################
# Test 9b: MCA -> CM: CM_NRW vs CM_notNRW
#########################################
test9b.pop$pop 
# Remove duplicated rows
test9b.pop <- distinct(test9b.pop)
(table(test9b.pop$pop))
#write_tsv(test9b.pop, file="data/pops/test9bpop.txt")

# input pop info from tibble
(sample.meta(test9b) <- sample.meta(test9b) %>% left_join(test9b.pop, by = "sampID"))

# Replace sample names "." with "_"
sample.meta(test9b)$sampID <- gsub("\\.", "_", sample.meta(test9b)$sampID)
summarize_facets(test9b)
sample.meta(test9b)$pop

## Run Abba-baba test and save to file ##
# run with jackknifing, 1000kb windows!
test9b <- calc_abba_baba(test9b, "pop.chr", "CM_NRW","CM_notNRW","MCA", jackknife=TRUE, sigma = 1000)

get.snpR.stats(test9b, "pop.chr", "abba_baba") # gets the per chr results
get.snpR.stats(test9b, "pop", "abba_baba") # gets the overall results

# smoothed windowed averages
test9b <- calc_smoothed_averages(test9b, "pop", sigma = 200, step = 200, 
                                 nk = TRUE, stats.type = "pairwise")
test9b <- calc_smoothed_averages(test9b, "pop.chr", sigma = 200, step = 200, 
                                 nk = TRUE, stats.type = "pairwise")
names(get.snpR.stats(test9b, "pop", "abba_baba"))

abba_baba_pairwise.stats <- get.snpR.stats(test9b, "pop", "abba_baba")$pairwise
abba_baba_pairwise.stats 
filename <- paste0("abbababa/out/", "test9b_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test9b, "pop", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test9b_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test9b, "pop", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test9b_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)

abba_baba_pairwise.stats <- get.snpR.stats(test9b, "pop.chr", "abba_baba")$pairwise
filename <- paste0("abbababa/out/", "test9b_chr_", "abba_baba_pairwise",".csv") 
write_csv(abba_baba_pairwise.stats, file = filename)

(abba_baba_pairwise.weighted.means <- get.snpR.stats(test9b, "pop.chr", "abba_baba")$weighted.means)
filename <- paste0("abbababa/out/", "test9b_chr_", "abba_baba_pairwise.weighted.means", ".csv") 
write_csv(abba_baba_pairwise.weighted.means, file = filename)

(abba_baba_pairwise.window <- get.snpR.stats(test9b, "pop.chr", "abba_baba")$pairwise.window)
filename <- paste0("abbababa/out/", "test9b_chr_", "pairwise.window", ".csv") 
write_csv(abba_baba_pairwise.window , file = filename)
################################################################################
## make a data table of all the pairwise weighted means ##
################################################################################
# Load necessary library
library(fs) # for file manipulation

# Define the source and destination folders
source_folder <- "abbababa/out"
destination_folder <- "abbababa/out/chr_results"

# Create the destination folder if it doesn't exist
if (!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

# List all files in the source folder that contain "_chr_" in the filename
file_list <- list.files(path = source_folder, pattern = "_chr_", full.names = TRUE)
file_list 

# Move each file to the destination folder
for (file_path in file_list) {
  # Extract the filename from the path
  file_name <- basename(file_path)
  
  # Define the new path in the destination folder
  new_path <- file.path(destination_folder, file_name)
  
  # Move the file
  file_move(file_path, new_path)
}

setwd(source_folder)
list.files()

# Function to read a single file, add a filename column, and return its content as a data frame
read_file <- function(file) {
  # Read the file
  df <- read_csv(file)
  
  # Extract the filename, remove path and extension, and clean it up
  filename <- gsub("_abba_baba_pairwise.weighted.means.csv$", "", basename(file))
  
  # Add the cleaned filename as a new column in the data frame
  df <- df %>% mutate(Filename = filename)
  
  return(df)
}

# Use lapply to read each file, add the filename column, and bind them into a single data frame
combined_data <- bind_rows(lapply(file_list, read_file), .id = "File_ID")

# View the combined data table
print(combined_data)

# Optionally, you can write the combined data table to a new CSV file
write_csv(combined_data, "abba_baba_combined_pairwise.weighted.means.csv")


################################################################################