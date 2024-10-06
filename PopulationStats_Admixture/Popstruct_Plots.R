################################################################################
#################### Custom Script to create Admixture Plots ###################
################################################################################
## Load packages ## 
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
library("RColorBrewer")
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")
library("sf")
library("tools")
library("rgdal")  # readOGR() spTransform()
library("raster")  # intersect()
library("ggsn")  # north2() scalebar()
library("paletteer")
library("hierfstat")
library("plotrix")
library("mapdata")
library("rworldmap")
library("fossil")
library("MASS")
library("Hmisc") 
library("cowplot")
library("viridis") # Use as color palette for colorblind friendly 

################################################################################
############################## DATA DIRECTORIES ################################
################################################################################
## set path to working directory ##
data_dir <- "/path/to/files"

## Set path to the ADMIXTURE results ##
ad_dir <- paste0(data_dir, "/PopulationStats_Admixture/")

## set working directory ##
setwd(data_dir)

################################################################################
############## Make popmap files with ADMIXTURE Q-scores and metadata ##########
################################################################################

## Subset metadata to retain only inidividuals from vcf file ##
vcfR <- read.vcfR("/Filtering/filteredVCF/noreponly_v2.75.renamed.LDpruned.vcf.gz")

### get metadata files ### 
md <- read_csv("noreponly_metadata_fixed.csv")
md

## subset via matching names to create popfiles ##
subset_metadata <- md %>%
  filter(Sample %in% colnames(vcfR@gt)) %>%
  slice(match(colnames(vcfR@gt), Sample))
dim(subset_metadata)
subset_metadata
glimpse(subset_metadata)

## set directory to admixture results
setwd(ad_dir)

## Read in Q3 data from Admixture analyses ## 
admixture_q3 <- read_delim("/admixture/K1-10/noreponly_v2.75.renamed.LDpruned.3.Q", col_names = F)
names(admixture_q3) <- c("Ad_P1", "Ad_P2", "Ad_P3")

## make another column for assignments based on highest value in respective cluster
clusterK3 <- apply(admixture_q3 , 1, which.max) #this corresponds with the 1:24 order 
clusterK3
table(clusterK3)

qmatrix_K3 <- cbind.data.frame(admixture_q3, clusterK3)

q_mat <-  cbind.data.frame(qmatrix_K3, subset_metadata)
## NOTE: I edited the q_md manually to see which cluster corresponded with
## which morph_sp group and relabed Ad_P1-P3 to its respective group for ease. I 
## also created a new group called Ad_acutus which was the sum of Ad_CA and Ad_MCA.

## Create column Ad_cluster by the Q score proportions defined in the manuscript.
# Assuming your tibble is named 'your_tibble', replace this with your actual tibble name
q_mat <- q_mat %>%
  mutate(Ad_cluster = case_when(
    Ad_CM > 0.98 ~ "CM_99",
    Ad_CM > 0.90 ~ "CM_90",
    Ad_CA > 0.98 ~ "CA_99",
    Ad_MCA > 0.98 ~ "MCA_99",
    Ad_CA > 0.90 ~ "CA_90",
    Ad_MCA > 0.90 ~ "MCA_90",
    Ad_CM >= 0.4 & Ad_CM <= 0.6 ~ "F1_Hybrid",
    Ad_CM < 0.90 & Ad_CM > 0.6 ~ "Hybrid_CM_backcross",
    Ad_CA < 0.90 & Ad_CA > 0.6 ~ "Hybrid_CA_backcross",
    Ad_MCA < 0.90 & Ad_MCA > 0.6 ~ "Hybrid_MCA_backcross",
    TRUE ~ "Acutus_backcross"
  ))
glimpse(dat)
dat$Ad_cluster

write.csv(q_md, "popmap.csv", row.names = FALSE)

## Write out separate .txt files for each unique Ad_cluster value ##

# Assuming your tibble is named 'q_mat', replace this with your actual tibble name
unique_clusters <- unique(q_mat$Ad_cluster)

# Create separate .txt files for each unique Ad_cluster value
for (cluster in unique_clusters) {
  # Filter rows for the current Ad_cluster value
  cluster_data <- q_mat[q_mat$Ad_cluster == cluster, "Sample"]
  
  # Convert the cluster_data to a character vector and format as desired
  cluster_data_str <- paste0('"', cluster_data$ID, '"', collapse = ", ")
  # Define the file name
  file_name <- paste0(cluster, ".txt")
  
  # Write the cluster_data_str to the .txt file
  writeLines(cluster_data_str, con = file_name)
}
################################################################################
# Assign the first argument to prefix
prefix="noreponly_v2.75.renamed.LDpruned"

# Get .list file from ADMIXTURE.md script for Monitoring Unit and Morph Species
# File should contain individual names in the correct order 
infofile <- "/admixture/popdata/noreponly_v2.75.renamed.LDpruned.MU.list.txt"
infofile <- "/admixture/popdata/noreponly_v2.75.renamed.LDpruned.Msp.list.txt"

#labels<-read.table(opt$infofile)
labels<-read.table(infofile)
# Name the columns
names(labels)<-c("ind","pop")
labels

# replace all NA values in pop to "Unk"
labels$pop[is.na(labels$pop)] <- "N/A"
labels

# Max K value
maxK <- 6

# Min K value 
minK <- 2

# Population order 
## abbreviated MU 
populations <- c("CC","NC","STW","NTW","NSLW","CB","N/A","BRW","CF","RHW","NRW")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
#labels$n<-factor(labels$pop,levels=unlist(strsplit(opt$populations,",")))
labels$n<-factor(labels$pop,levels=unlist(strsplit(populations,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))
labels

# read in the different admixture output files
#minK=opt$minK
#maxK=opt$maxK
tbl<-lapply(minK:maxK, function(x) read.table(paste0(prefix,".",x,".Q")))

# Prepare spaces to separate the populations/species
rep<-as.vector(table(labels$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),1)}
spaces<-spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual for each K as a separate row
#tiff(file=paste0(opt$outPrefix,".tiff"),width = 2000, height = 1200,res=200)
# par(mfrow=c(maxK-1,1),mar=c(0,1,0,0),oma=c(2,1,9,1),mgp=c(0,0.2,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)
 # Plot minK
# bp<-barplot(t(as.matrix(tbl[[1]][order(labels$n),])), col=rainbow(n=minK),xaxt="n", border=NA,ylab=paste0("K=",minK),yaxt="n",space=spaces)
# axis(3,at=bp,labels=labels$ind[order(labels$n)],las=2,tick=F,cex=0.6)
 # Plot higher K values
# if(maxK>minK)lapply(2:(maxK-1), function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=rainbow(n=x+1),xaxt="n", border=NA,ylab=paste0("K=",x+1),yaxt="n",space=spaces))
# axis(1,at=c(which(spaces==0.5),bp[length(bp)])-diff(c(1,which(spaces==0.5),bp[length(bp)]))/2,
#     labels=unlist(strsplit(opt$populations,",")))
#dev.off()

#tiff(file=paste0(opt$outPrefix,".tiff"),width = 2000, height = 1200,res=200)
#par(mfrow=c(maxK-1,1),mar=c(0,1,0,0),oma=c(2,1,9,1),mgp=c(0,0.2,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)
# Plot minK
#bp<-barplot(t(as.matrix(tbl[[1]][order(labels$n),])), col=rainbow(n=minK),xaxt="n", border=NA,ylab=paste0("K=",minK),yaxt="n",space=spaces)
#axis(3,at=bp,labels=labels$ind[order(labels$n)],las=2,tick=F,cex=0.6)
# Plot higher K values
#if(maxK>minK)lapply(2:(maxK-1), function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=rainbow(n=x+1),xaxt="n", border=NA,ylab=paste0("K=",x+1),yaxt="n",space=spaces))
#axis(1,at=c(which(spaces==0.5),bp[length(bp)])-diff(c(1,which(spaces==0.5),bp[length(bp)]))/2,
#     labels=unlist(strsplit(populations,",")))
#dev.off()

################################################################################
# Assuming you have the rest of the script setup for reading files and preparing the data

# Define your color schemes here
# Example color schemes for K=2, K=3, and so on. Add more as needed.
# Each list element is a vector of colors for each ancestry component at a given K
color_schemes <- list(
  c("#E41A1CFF", "#377EB8FF"),             # Colors for K=2
  c("#E41A1CFF","#4DAF4AFF","#377EB8FF"),  # Colors for K=3
  c("#E41A1CFF","#984EA3FF", "#377EB8FF","#4DAF4AFF"), # Example for K=4, etc.
  c("#4DAF4AFF","#984EA3FF","#377EB8FF","#FF7F00FF", "#E41A1CFF"),
  c( "#E41A1CFF","#4DAF4AFF","#984EA3FF","#FF7F00FF","#FFFF33FF","#377EB8FF"))

color_schemes

figure<- "/admixture/figures"
setwd(figure)

tiff(file=paste0(prefix, "_ADMIXTURE.MU",".tiff"),width = 2500, height = 1200,res=200)
png(file=paste0(prefix, "_ADMIXTURE.MU",".png"),width = 2500, height = 1100, res=300)
png("ADMIXTURE.MU.png",width = 2500, height = 1100, res=300)
svg(file=paste0(prefix, "_ADMIXTURE.MU",".svg"), width = 12, height = 6)
pdf(file=paste0(prefix, "_ADMIXTURE.MU",".pdf"))

par(mfrow=c(maxK-1,1),mar=c(0,1,0,0),oma=c(2,1,2,1),mgp=c(0,0.3,0),xaxs="i",cex.lab=1.2,cex.axis=0.9)

# Initial population plot for minK
bp <- barplot(t(as.matrix(tbl[[1]][order(labels$n),])), col=color_schemes[[1]], xaxt="n", border=NA, ylab=paste0("K=",minK), yaxt="n", space=spaces)
#axis(3, at=bp, labels=labels$ind[order(labels$n)], las=2, tick=F, cex=0.6)

# Population labels for minK (Adjust the 'axis' function here if needed)

# Calculate midpoints for population labels more accurately
pop_names <- unlist(strsplit(populations, ","))

# Calculate midpoints by adding half the group width to the start positions

## MU ##
#midpoints <- c(0, cumsum(group_widths))[-length(group_widths)+1] + group_widths / 2
midpoints <- c(19.0,  52.0,  73.5,  82.0,  90, 124.0, 159.5, 173.0, 183, 186, 217) # changed the last position manually

# Plot for higher K values, using the color schemes defined above
if (maxK > minK) {
  lapply(2:(maxK-1), function(x) {
    barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=color_schemes[[x]], xaxt="n", border=NA, ylab=paste0("K=",x+1), yaxt="n", space=spaces)
  })
}

# Adjusted labels for populations, centered accurately
axis(1, at=bp[midpoints], labels=pop_names, cex.axis=0.9)
dev.off()

################################################################################
## For Morph_Sp Groups ##
################################################################################
# Adjust plotting parameters for displaying all populations
par(mfrow=c(maxK-1,1),mar=c(0,1,0,0),oma=c(2,1,2,1),mgp=c(0,0.3,0),xaxs="i",cex.lab=1.2,cex.axis=0.9)

# Plot for minK
bp <- barplot(t(as.matrix(tbl[[1]][order(labels$n),])), col=color_schemes[[1]], xaxt="n", border=NA, ylab=paste0("K=",minK), yaxt="n", space=spaces)
#axis(3, at=bp, labels=labels$ind[order(labels$n)], las=2, tick=F, cex=0.6)

## MSp
populations <- c("Cayes C. acutus", "Mainland C. acutus", "Putative Hybrids", "C. moreletii")

# Plot for higher K values, using the color schemes defined above
if (maxK > minK) {
  lapply(2:(maxK-1), function(x) {
    barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=color_schemes[[x]], xaxt="n", border=NA, ylab=paste0("K=",x+1), yaxt="n", space=spaces)
  })
}
axis(1, at=c(which(spaces==1), bp[length(bp)]) - diff(c(1, which(spaces==1), bp[length(bp)]))/2,
     labels=unlist(strsplit(populations, ",")), cex.axis=1)
dev.off()

tiff(file=paste0(prefix, "_ADMIXTURE.Msp",".tiff"),width = 2500, height = 1200,res=300)
png(file=paste0(prefix, "_ADMIXTURE.Msp",".png"),width = 2500, height = 1100, res=300)
png("ADMIXTURE.Msp.png",width = 2500, height = 1100, res=300)
svg(file=paste0(prefix, "_ADMIXTURE.Msp",".svg"), width = 12, height = 6)
pdf(file=paste0(prefix, "_ADMIXTURE.Msp",".pdf"))

################################################################################
## Barchart + admixture maps ##
################################################################################
# From q_md created above:
admix_data <- read_csv("/admixture/popdata/popmap.csv")
admix_data

world <- getMap(resolution = "low")
world_belize <- world[world@data$ADMIN == "Belize", ]
world_belize
#bz <- CRS("+init=EPSG:2028")  #projected to Belize UTM
Belize1<-getData("GADM", country="BZ", level=1)
Belize1
class(Belize1)

par(mar = c(0, 0, 0, 0))
plot(Belize1) # better map of Belize with district lines 


admix_data <- admix_pops
plot(Belize1)
# Plot pies at each locality
# Plot pies at each locality
#for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
#                                        c(admix_K3$Ad_P1[x],admix_K3$Ad_P2[x],admix_K3$Ad_P3[x]), radius=0.05, col=c("#E41A1CFF","#4DAF4AFF","#377EB8FF"))}

# Define custom colors
color_CM = "#E41A1CFF"     # Replace with your desired color or hexadecimal color code
color_CA = "#377EB8FF"     # Replace with your desired color or hexadecimal color code
color_MCA = "#4DAF4AFF"  # Replace with your desired color or hexadecimal color code


plot(Belize1)
# Loop over each row in the admix_data to plot pies
for (x in 1:nrow(admix_data)) {
  floating.pie(admix_data$Longitude[x], admix_data$Latitude[x],
               c(admix_data$Ad_CM[x], admix_data$Ad_CA[x], admix_data$Ad_MCA[x]),
               radius=0.05,
               col=c(color_CM, color_CA, color_MCA))
}
dev.off()


## K = 2 ##
Ad_Map_K3 <-paste0("/admixture/figures/", "Ad_Map_K3_", basefile,".png")
png(file=Ad_Map_K3, width = 700, height = 700, res = 600) # open the pdf plotting device
png(file=paste0("K3_ADMIX_MAP",".png"),width = 700, height = 700)
#par(mar=c(1,0.25,1, 0.25))
par(mar = c(0, 0, 0, 0))
plot(Belize1)
# Loop over each row in the admix_data to plot pies
for (x in 1:nrow(admix_data)) {
  floating.pie(admix_data$Longitude[x], admix_data$Latitude[x],
               c(admix_data$Ad_CM[x], admix_data$Ad_CA[x], admix_data$Ad_MCA[x]),
               radius=0.05,
               col=c(color_CM, color_CA, color_MCA))
}
dev.off()

par(mfrow = c(2, 2)) # Create a 2 x 2 plotting matrix
# Plot 1 

# Plot 2 

################################################################################
# Adjust plotting parameters for displaying all populations
par(mfrow=c(maxK-1,1), mar=c(0,1,0,0), oma=c(2,1,9,1), mgp=c(0,0.2,0), xaxs="i", cex.lab=1.2, cex.axis=0.8)

# Plot for minK
bp <- barplot(t(as.matrix(tbl[[1]][order(labels$n),])), col=color_schemes[[1]], xaxt="n", border=NA, ylab=paste0("K=",minK), yaxt="n", space=spaces)
axis(3, at=bp, labels=labels$ind[order(labels$n)], las=2, tick=F, cex=0.6)


# Plot for higher K values, using the color schemes defined above
if (maxK > minK) {
  lapply(2:(maxK-1), function(x) {
    barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=color_schemes[[x]], xaxt="n", border=NA, ylab=paste0("K=",x+1), yaxt="n", space=spaces)
  })
}
axis(1, at=c(which(spaces==1), bp[length(bp)]) - diff(c(1, which(spaces==1), bp[length(bp)]))/2,
     labels=unlist(strsplit(populations, ",")), cex.axis=1)
dev.off()

## Adjust plotting parameters for displaying all populations
figure<- "/admixture/figures"
setwd(figure)

tiff(file=paste0(prefix, ".MU_v3",".tiff"),width = 2500, height = 1200,res=200)
par(mfrow=c(maxK-1,1), mar=c(0,1,0,0), oma=c(2,1,9,1), mgp=c(0,0.2,0), xaxs="i", cex.lab=1.2, cex.axis=0.8)

# Initial population plot for minK
bp <- barplot(t(as.matrix(tbl[[1]][order(labels$n),])), col=color_schemes[[1]], xaxt="n", border=NA, ylab=paste0("K=",minK), yaxt="n", space=spaces)
axis(3, at=bp, labels=labels$ind[order(labels$n)], las=2, tick=F, cex=0.6)

# Population labels for minK (Adjust the 'axis' function here if needed)

# Calculate midpoints for population labels more accurately
pop_names <- unlist(strsplit(populations, ","))
# Assuming 'spaces' marks the start of a new population, find positions where spaces are placed
space_positions <- which(spaces == 1)
space_positions
#as.character(space_positions)

# Calculate the number of bars in each group to find group widths
group_widths <- diff(c(0, space_positions, length(bp)+1)) - 1  # Minus one because spaces count as extra here
group_widths
#group_widths <- c(36,28,13, 6, 5, 60, 9, 16,4,1, 54)
group_widths/2

# Calculate midpoints by adding half the group width to the start positions
#midpoints <- c(0, cumsum(group_widths))[-length(group_widths)+1] + group_widths / 2
midpoints <- space_positions-(group_widths/2)
midpoints <- c(19.0,  52.0,  73.5,  82.0,  90, 124.0, 159.5, 173.0, 183, 186.5, 217) # changed the last position manually

# Plot for higher K values, using the color schemes defined above
if (maxK > minK) {
  lapply(2:(maxK-1), function(x) {
    barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=color_schemes[[x]], xaxt="n", border=NA, ylab=paste0("K=",x+1), yaxt="n", space=spaces)
  })
}

# Adjusted labels for populations, centered accurately
axis(1, at=bp[midpoints], labels=pop_names, cex.axis=0.8)
dev.off()
