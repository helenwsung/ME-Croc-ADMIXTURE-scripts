################################################################################
## Script to run Population assignment in R using DAPC, IBD, sNMF, Admixture ##
################################################################################
## https://wyoibc.github.io/popgen_workshop/week4/DAPC_snmf/dapc_ibd_snmf.html 
## NOTE: Make sure you have all relevant outfiles after running Ipyrad 
################################################################################
################################################################################
## Load packages ## 

## For editing genid files
#install.packages("dartR")
#BiocManager::install("SNPRelate")
library("SNPRelate")
library("dartR")
library("snpR")

########################
## For running PCA, sNMF
library("adegenet")
library("RColorBrewer")
library("dplyr")
library("magrittr")
library("tibble")
library("ggplot2")
library("reshape2")
library("vegan")
library("gplots")
library("edgeR")
library("vcfR")
library("LEA")

################################################
## For running DAPC and IBD and mapping
library("hierfstat")
library("plotrix")
library("mapdata")
library("rworldmap")
library("fossil")
library("MASS")
library("plyr")
library("tidyverse")
library("Hmisc") # for imputing data for IBD
library("cowplot")
library("tidyr")
library("viridis") # Use as color palette for colorblind friendly 

################################################
## Packages for making my Belize map
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")
library("sf")
library("tools")
library("rgdal")  # readOGR() spTransform()
library("raster")  # intersect()
library("ggsn")  # north2() scalebar()
library("paletteer")

################################################################################
############################## DATA DIRECTORIES ################################
################################################################################
################################################################################
## Then we will specify a number of file paths and read in a few files so that 
## we don’t have to repeatedly hardcode file paths farther down in the script. 
## This makes it easier to reuse the script on different datasets or the same 
## data with different filtering schemes without having to search through the 
## script for every time an absolute file path is specified.
################################################################################

## Set up an object containing the path to the data ##

data_dir <- "/path/to/files"
data_dir <- getwd()

################################################################################
## For filtered data through snpFiltR and LD-pruned via SeqArray & SNPrelate 
filteredVCF <- paste0(data_dir, "/filteredVCF_ab")

#data_dir <- filteredVCF
#data_dir
################################################################################

## Set our working directory ##
setwd(data_dir)
getwd()
list.files()
################################################################################
## Set up object for Output Directory ##
################################################################################

out_dir <- paste0(filteredVCF, "/Pop_structr_out/") 

if(!dir.exists(out_dir)){ # check if the directory exists
  dir.create(out_dir)   # and create it if it does not
}

################################################################################
## Set up an object that contains the base file name of files in the output directory. 
## Data files are all this basename with varying extensions
## we won't call this 'basename' because that is a function in R

# noreponly_v2: All 3RADmerged samples with no repeats from ipyrad and filltered
# 242 individuals
# 12862 variants 
basefile <- "noreponly_v2"

######################################################################################
## Set up paths to input files using the base file name specified above ##
path_ugeno<-paste0(data_dir,"/", basefile,".ugeno")
path_ustr<-paste0(data_dir,"/", basefile,".ustr")
path_vcf<-paste0(data_dir,"/", basefile,".vcf")
path_coords<-paste0(data_dir,"/", basefile,"_coords.csv")

## Path for SNPfiltR filtered data 
filtered_thinned_vcf <-paste0(filteredVCF,"/data/", basefile,".75.renamed.LDpruned.vcf.gz")
path_vcf <- filtered_thinned_vcf

filtered_thinned_geno <-paste0(filteredVCF,"/data/", basefile,".75.renamed.LDpruned.geno")
path_ugeno <- filtered_thinned_geno
path_ugeno

######################################################################################
## Get in MetaData ## 
md_csv <- read_csv("/datadir/noreponly_metadata_fixed.csv")
md <- read_csv(md_csv)
md
names(md)
glimpse(md)

# make coords dataframe for IBD and mapping
coords.tbl <- md %>% dplyr::select(Sample, Longitude, Latitude)
coords.tbl

coords <- as.data.frame(coords.tbl)
coords
################################################################################
# make dataframe for 'population' later
glimpse(md)

pops <- md %>% dplyr::select(Sample, Longitude, Latitude, Morph_Species, Monitoring.Unit, Abbrv_Monitoring.Unit, Route, Site, Subdivision, Ad_cluster)
pops
names(pops)

################################################################################
################ Discriminant analysis of principal components #################
################################################################################
## Discriminant analysis of principal components, or DAPC, is a method of 
## classifying individuals into clusters that does not include any explicit 
## population genetic model. We'll use the unlinked Structure-formatted file for 
## this, but because Structure-formatted files can come in a variety of different 
## configurations, we need to tell the function how many individuals and loci 
## are present in the file. Read in the 'ugeno' file and use the dimensions of 
## that file to get the number of individuals and snps.
################################################################################
## read in vcf file and convert to genlight object
path_vcf
gendata_vcf<-read.vcfR(path_vcf) # read in all of the genetic data
gendata_vcf 

gendata_gl<-vcfR2genlight(gendata_vcf) # make the genetic data a biallelic matrix of alleles in genlight format
#quick summary of genlight to see what info we have
gendata_gl

gendata_names <- indNames(gendata_gl) # get the sample names
gendata_names

#subset popmap to only include retained individuals
pops
pops<-pops[pops$Sample %in% gendata_names,]
dim(pops)

################################################################################
## read in the geno file to get the number of individuals and snps for this assembly
geno_txt<-readLines(path_ugeno)

## The number of lines is the number of loci
nums_snps<-length(geno_txt) 
nums_snps # 12866 loci 

## the number of columns is the number of individuals
num_ind<-length(strsplit(geno_txt[[1]], "")[[1]]) # here we split apart the first line into individual characters and count the length
num_ind 


## convert genlight to genind for LDthinned data
df_genind <- gl2gi(gendata_gl)
DAPC_ustr <- df_genind
DAPC_ustr

# fill in pop with morph_species
DAPC_ustr@pop <- as.factor(md$Morph_Species)

# fill in pop with monitoring unit
pops$Abbrv_Monitoring.Unit[is.na(pops$Abbrv_Monitoring.Unit)] <- "UK"
DAPC_ustr@pop <- as.factor(pops$Abbrv_Monitoring.Unit)
#md$Abbrv_Monitoring.Unit <- md$Abbrv_Monitoring.Unit %>% replace_na("unknown")
#DAPC_ustr@pop <- as.factor(md$Abbrv_Monitoring.Unit)
DAPC_ustr@pop

## Get the individual names in the order that they show up in the various files 
## this is important farther down for getting coordinates into the right order for plotting
ind_names<-rownames(DAPC_ustr@tab)
ind_names

## Run DAPC ##
dapc1<-dapc(DAPC_ustr)

dapc1<-dapc(DAPC_ustr, n.pca=num_ind, n.da=1) 
dapc2<-dapc(DAPC_ustr, n.pca=num_ind, n.da=2) 
dapc3<-dapc(DAPC_ustr, n.pca=num_ind, n.da=3)

# we can assess the trade-off between power of discrimination and over-fitting 
# by calculating the alpha-score, which is the difference between the proportion 
# of successful reassignment of the analysis (observed discrimination) and values
# obtained using random groups (random discrimination).
temp1 <- optim.a.score(dapc1) # 207 for MU; 39 for morph sp
temp1

temp2 <- optim.a.score(dapc2) # 165 for MU; 22 for morph sp
temp2

temp3 <- optim.a.score(dapc3) # 142 for MU; 20 for morph sp
temp3

## start by determining how many clusters we want to use, if best fit >= max.n.clust, 
## expand max number of clusters to ensure that we aren't artificially limiting 
## number of clusters that our data can fall into.
grp <- find.clusters(DAPC_ustr, max.n.clust=8) # test up to 8 clusters, 
## spits out plot of variance by each PC in the PCA transformation, 
## we want to retain all PCs up to point at which the cumulative variance plateaus, 
## if no clear plateau then retail all PCs 
# input number of pc found = # of individuals 
## input number of clusters we want to retain, with the fit of each number of 
## clusters estimated by BIC. lowest BIC = best fit or look for inflection point 
### w/ sharpest decrease in BIC (the elbow in the BIC curve as function of K)
# input = 2

grp <- find.clusters(DAPC_ustr, max.n.clust=8, n.pca = 242, n.clust = 2) # test up to 8 clusters, 
grp3 <- find.clusters(DAPC_ustr, max.n.clust=8, n.pca = 242, n.clust = 3) # group by 3 instead of 2
grp4 <- find.clusters(DAPC_ustr, max.n.clust=8, n.pca = 242, n.clust = 4) # group by 4 instead of 2


## Output from find.clusters is a list
names(grp)
grp$Kstat  # What is the lowest kstat?
grp$stat
head(grp$grp, 10)
grp$size
table(pop(DAPC_ustr), grp$grp)

# Rows correspond to our morph group, columns to inferred groups 
table.value(table(pop(DAPC_ustr), grp$grp), col.lab=paste("inf_groups", 1:length(levels(grp$grp))),
            row.lab=levels(pop(DAPC_ustr)))

## grp object now contains the groupings of these individuals into the 2 clusters we decided on, and we can use the function dapc to describe these groups.
dapc1 <- dapc(DAPC_ustr, grp$grp) # run DAPC
# input number of pc found = 275; 230; 256; best alpha score
# We will then be asked how many discriminant functions to retain. 
# With only 2 groups, only 1 is possible --> input number = 1
dapc2 <- dapc(DAPC_ustr, grp3$grp) # run DAPC with 2 discriminant functions from 3 groups
dapc3 <- dapc(DAPC_ustr, grp4$grp) # run DAPC with 3 discriminant functions from 4 groups

DAPC_scatterplot <-paste0(out_dir, "/DAPC/", "DAPC_scatterplot_", basefile,".pdf")
pdf(DAPC_scatterplot)
scatter(dapc1)
scatter(dapc2)
scatter(dapc3) # why don't my points show behind clusters? - lower PC = more visual
## membership probabilities based on the retained discriminant functions 
# stored in dapc objects in the slot "posterior" 
dev.off()

dapc1
class(dapc1$posterior)
dim(dapc1$posterior)
summary(dapc1)
dapc_groups2 <- round(dapc1$posterior,2) # Each row corresponds to an individual, each column to a group.
dapc_groups2

class(dapc2$posterior)
dim(dapc2$posterior)
round(head(dapc2$posterior),3) # Each row corresponds to an individual, each column to a group.
summary(dapc2)
dapc_groups3 <- round(dapc2$posterior,3)

dim(dapc3$posterior)
round(head(dapc3$posterior),3) # Each row corresponds to an individual, each column to a group.
summary(dapc3)
dapc_groups4 <- round(dapc3$posterior,3)
dapc_groups4

dapc_groups <- list(
  dapc.g1 = round(dapc1$posterior, 2),
  dapc.g2 = round(dapc2$posterior, 3),
  dapc.g3 = round(dapc3$posterior, 3)
)
DAPC_scatterplot <-paste0(out_dir, "/DAPC/", obj_name, "_",basefile,".csv")

# Loop through the list of objects and write them to separate CSV files
for (obj_name in names(dapc_groups)) {
  data <- dapc_groups[[obj_name]]
  DAPC_scatterplot <-paste0(out_dir, "/DAPC/", obj_name, "_",basefile,".csv")
  write.csv(data, DAPC_scatterplot, row.names = TRUE)
}

## Do all at once ##
my_k <- 2:4
grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))
dapc_l

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(DAPC_ustr, n.pca = 200, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(DAPC_ustr, pop = grp_l[[i]]$grp, n.pca = 200, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

grp <- grp_l[[1]]
grp
grp3 <- grp_l[[2]]
grp4 <- grp_l[[3]]
grp_l
##########################################################
################# Plotting out everything ################
##########################################################

setwd(out_dir) # set to out directory 
setwd("/datadir/filteredVCF_ab/Pop_structr_out/DAPC/")

## Plot out our confidence in assigning individuals to each of the groups w/dapc1 ##
scatter(dapc_l[[1]], posi.da="bottomleft")
dapc_l[[1]]
## plot the DAPC the ugly way
scatter(dapc1, col=colors_2,  bg="white",
        legend=FALSE, posi.da = "bottomright",
        solid=.5)
scatter(dapc_l[[2]], posi.da="bottomright",  bg="white",
        pch=17:22, cstar=0, col=colors_3, scree.pca=TRUE,
        posi.pca="bottomleft")
scatter(dapc_l[[3]], posi.da = "bottomleft")
scatter(dapc3,1,1, col=colors_4, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

## Plot from $assign.per.pop indicating proportions of successful reassignment (based on the discriminant functions) of individuals to their original clusters. 
# Large values indicate clear-cut clusters, while low values suggest admixed groups.
# Heat colors represent membership probabilities (red=1, white=0); blue crosses represent the prior cluster provided to DAPC.
# DAPC classification is consistent with the original clusters (blue crosses are on red rectangles),
# if blue cross on yellow, then there is discrepancy where its classified in red group but DAPC would assign to yellow group
# useful when prior biological groups are used, as one may infer admixed or misclassified individuals.
assignplot(dapc1)
assignplot(dapc2)
assignplot(dapc3)

## To find most "admixed" individuals, consider admixed individual having no more than 90% prob of membership in single cluster
temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.9)))
temp

## Plot in STRUCTURE-like way with with bar chart
compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:length(grp$size)), lab=rownames(dapc1),
          n.col=1, xlab="individuals", show.lab = TRUE, col = colors_2)

compoplot(dapc2, posi="bottomright",
          txt.leg=paste("Cluster", 1:length(grp3$size)), lab=rownames(dapc2),
          n.col=1, xlab="individuals", show.lab = TRUE, col = colors_3)

compoplot(dapc3, posi="bottomright",
          txt.leg=paste("Cluster", 1:length(grp4$size)), lab=rownames(dapc3),
          n.col=1, xlab="individuals", show.lab = TRUE, col = colors_4)

compoplot(dapc3, posi="bottomright",
          txt.leg=paste("Cluster", 1:length(grp4$size)), lab=rownames(dapc3),
          n.col=1, xlab="individuals", show.lab = TRUE, col = colors_4)

my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")
for(i in 1:3){
  compoplot(dapc_l[[i]], posi="bottomright",
            txt.leg=paste("Cluster", 1:length(grp_l[[i]]$size)), lab=rownames(dapc_l[[i]]),
            n.col=1, xlab="individuals", show.lab = TRUE, col = my_pal)
}
## Plot using spatial pattern on map ##
## first process the coordinates (remember that we read these into R earlier) 
## to make sure that the coordinates are in the right order and that every 
## individual we want to plot has locality data.

## We'll first do a little processing of the coordinates (remember that we read these into R earlier) 
##    to make sure that the coordinates are in the right order and that every individual we want to plot 
##    has locality data.

## make sure there aren't any individuals that don't have coordinates
ind_names[which(!ind_names %in% coords[,"Sample"])]

## match up the coordinates to the order of the individuals from genetic data
match_coords<-match(ind_names, coords[,"Sample"])
coords<-coords[match_coords,]
coords

## Plot out the weird map I don't like
map("worldHires", "Belize", xlim=c(-89.5,-87.5), ylim=c(15.8,18.5),col="gray90", fill=TRUE)

## Add pie charts showing the probability of membership in each cluster for each sample
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc1$posterior[x,1], dapc1$posterior[x,2]), radius=0.1, col=my_pal)}

###############################################
## Better Map from my own script ##
###############################################

world <- getMap(resolution = "low")
world_belize <- world[world@data$ADMIN == "Belize", ]
world_belize
#bz <- CRS("+init=EPSG:2028")  #projected to Belize UTM
Belize1<-getData("GADM", country="BZ", level=1)
Belize1
class(Belize1)

par(mar = c(1, 1, 1, 1))
plot(Belize1) # better map of Belize with district lines 
dev.off()

###############################################
## Even better map with ggplot ## 
###############################################

## Identify individuals missing geographic data 
missing.ind <- coords[rowSums(is.na(coords)) > 0, ]   
missing.ind <- missing.ind$Sample 
missing.ind

# Create new dataframe without individuals missing geographic data
coords_nomissing<- filter(coords.tbl, !(Sample %in% missing.ind))
coords_nomissing <- as.data.frame(coords_nomissing)
coords_nomissing

## Convert data frame point plots into sf coordinate objects
(coord_sites <- st_as_sf(coords_nomissing, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant"))

## Convert districts to sf data and use to get names for map
districts <- st_as_sf(Belize1, plot = FALSE, fill = TRUE)
districts
class(districts)
head(districts)

districts <- cbind(districts, st_coordinates(st_centroid(districts)))
head(districts)

##ggplot map with districts and legend ##
Belize <- ggplot() +
  #geom_point(data = saltycroc,  # Add and plot species data
  #aes(x = long, y = lat, color = factor(Species, labels = c("C. acutus", "C. moreletii", "Hybrid")))) +
  #labs(color = "Putative Species") +
  geom_sf(data = districts, fill = NA) +
  #geom_text(data = districts, aes(X , Y , label = NAME_1), size = 2.5, fontface = "bold") +
  coord_sf(xlim = c(-89.5, -87.5), ylim = c(15.8, 18.6), expand = FALSE) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotation_scale(location = "br", width_hint = 0.4) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_classic()  # Remove ugly grey background 

plot(Belize)

paletteer_d("ggsci::default_ucscgb")
paletteer_d("colorBlindness::paletteMartin")

##############################################################################################
par(mar=c(1,0.5,1, 0.5))
plot(Belize1)

# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc1$posterior[x,1], dapc1$posterior[x,2]), radius=0.05, col=colors_2)}
# DAPC2
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc2$posterior[x,1], dapc2$posterior[x,2], dapc2$posterior[x,3]), radius=0.05, col=colors_3)}
# DAPC3
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc3$posterior[x,1], dapc3$posterior[x,2], dapc3$posterior[x,3], dapc3$posterior[x,4]), radius=0.05, col=colors_4)}

# DAPC3
dapc2<- dapc_l[[2]]
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc_l[[3]]$posterior[x,1], dapc_l[[3]]$posterior[x,2], dapc_l[[3]]$posterior[x,3], dapc_l[[3]]$posterior[x,4]), radius=0.05, col=my_pal)}



## MAKE PDF image of maps
DAPC_Map <-paste0(out_dir,"DAPC","/DAPC_Map_", basefile,".png")
DAPC_Map
#pdf(file=DAPC_Map, width=10, height=10) # open the pdf plotting device
png(DAPC_Map)
par(mar=c(1,0.25,1, 0.25))
plot(Belize1) # Plot out the map
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x], 
                                        c(dapc1$posterior[x,1], dapc1$posterior[x,2]), radius=0.05, col=c("#E41A1CFF", "#377EB8FF"))}
dev.off() # close the pdf plotting device

DAPC_Map <-paste0(out_dir,"DAPC/","DAPC2_Map_", basefile,".pdf")
DAPC_Map
pdf(file=DAPC_Map, width=10, height=10) # open the pdf plotting device
par(mar=c(1,0.25,1, 0.25))
plot(Belize1) # Plot out the map
# DAPC2
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc2$posterior[x,1], dapc2$posterior[x,2], dapc2$posterior[x,3]), col=c("#4DAF4AFF","#E41A1CFF", "#377EB8FF"),radius=0.05)}
dev.off() # close the pdf plotting device

DAPC_Map <-paste0(out_dir,"DAPC/","DAPC3_Map_", basefile,".pdf")
DAPC_Map
pdf(file=DAPC_Map, width=10, height=10) # open the pdf plotting device
par(mar=c(1,0.25,1, 0.25))
plot(Belize1) # Plot out the map
# DAPC3
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc3$posterior[x,1], dapc3$posterior[x,2], dapc3$posterior[x,3], dapc3$posterior[x,4]), radius=0.05, col=c( "#984EA3FF","#E41A1CFF","#4DAF4AFF", "#377EB8FF"))}
dev.off() # close the pdf plotting device

## because this is a classification problem and isn’t explicitly modeling admixture, 
## we only see the confidence with which each sample is assigned to each cluster. 
## I.e., this does not indicate that there is no evidence of admixture among populations. 

################################################################################
## Plots for model-based blustering of ancestral pops ## 
# Plotting structure, DAPC, or Admixture results with ggplot2 
# https://luisdva.github.io/rstats/model-cluster-plots/
library(dplyr)
library(tibble)
library(purrr)

## DAPC Plot with population assignments ##

#https://luisdva.github.io/rstats/dapc-plot/
# Load libraries. Install first if needeed
library(readxl)    # CRAN v1.3.1
library(janitor)   # CRAN v2.1.0
library(dplyr)     # CRAN v1.0.7
library(tidyr)     # CRAN v1.1.3
library(ggplot2)   # CRAN v3.3.5
library(forcats)   # CRAN v0.5.1
library(stringr)   # CRAN v1.4.0
library(ggh4x)     # [github::teunbrand/ggh4x] v0.2.0.9000
library(paletteer) # CRAN v1.4.0
library(extrafont) # CRAN v0.17


# create an object with membership probabilities based on the retained discriminant functions
postprobs_dapc1 <- as.data.frame(round(dapc1$posterior, 3))  # Each row corresponds to an individual, each column to a group.
postprobs_dapc1

postprobs_dapc2 <- as.data.frame(round(dapc2$posterior, 3))  # Each row corresponds to an individual, each column to a group.
postprobs_dapc2

postprobs_dapc3 <- as.data.frame(round(dapc3$posterior, 3))  # Each row corresponds to an individual, each column to a group.
postprobs_dapc3

pops<- arrange(pops, Sample) # rearrange to match gendata
pops

# put probabilities in a tibble with IDS and labels for sites
clusters <- tibble::rownames_to_column(postprobs_dapc2, var = "Sample") %>%
  left_join(pops, by = "Sample")
clusters 

head(clusters)
head(coords)
head(postprobs_dapc2)

# melt into long format
croc_long <- clusters %>% pivot_longer(2:4, names_to = "cluster", values_to = "prob")
head(croc_long)
croc_long$cluster
# manual relevel of the sampling sites (to avoid alphabetical ordering)
croc_long$Morph_Species <- fct_relevel(as.factor(croc_long$Morph_Species), "CA", "CM", "HY", "MCA")
croc_long

# set up custom facet strips
facetstrips <- strip_nested(
  text_x = elem_list_text(size = c(12, 4)),
  by_layer_x = TRUE, clip = "off"
)

setwd(out_dir) # set to out directory 
DAPC_barchart <-paste0(out_dir,"DAPC/","DAPC2_barchart_MorphSp", basefile,".pdf")
pdf(DAPC_barchart, width=10, height=10) # open the pdf plotting device
ggplot(croc_long, aes(factor(Sample), prob, fill = factor(cluster))) +
  geom_col(color = "gray", size = 0.01) +
  facet_nested(~ Morph_Species,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips,
  ) +
  labs(x = "Individuals", y = "Membership probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  #scale_fill_paletteer_d("ghibli::PonyoMedium", guide = "none") +
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )
dev.off() # close the pdf plotting device
#######################################
## Different way of plotting##

#tidy the data for plotting. 
dapc_data_df <-
  # as_tibble() converts the ind.coord matrix to a special data frame called a tibble. 
  as_tibble(dapc1$ind.coord, rownames = "Sample") %>%
  # mutate changes or adds columns to a data frame. Here, we're adding the population and group assignment columns to the data frame
  mutate(Morph_Species = pops$Morph_Species,
         group = dapc1$grp)

dapc_data_df

#plot the data. Can color the points according to your pre-defined populations and the dapc groups to see if it conforms to your hypothesis.
dapc_plot <-
  ggplot(dapc_data_df, aes(
    x = LD1,
    #y = LD2,
    fill = Morph_Species
  )) +
  geom_point(shape = 21, size = 3) 

dapc_plot
################################################################################
############################ Isolation by distance #############################
################################################################################
## if isolation by distance (IBD) exists in our dataset, programs that seek to 
## cluster individuals but do not model continuous spatial structure can be 
## positively misled by IBD. This can result in overestimating the number of 
## population clusters, potentially identifying discrete population structure 
## when no such structure exists. The method conStruct can explicitly model both 
## processes, but can be finicky to run and have long run times. Do quick test 
## and visualization of isolation by distance to try to determine if IBD is
## misleading our population structure analyses.
################################################################################

#genlight_gl<-read.vcfR(path_vcf) # read in all of the genetic data
#gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
gendata <- gendata_gl
gendata_names <- indNames(gendata) # get the sample names
gendata_names

## make sure there aren't any individuals that don't have coordinates
gendata_names[which(!gendata_names %in% coords[,"Sample"])]
coords

## match up the coordinates to the order of the individuals from genetic data
match_coords<-match(gendata_names, coords[,"Sample"])
coords<-coords[match_coords,]
coords

################################################################################
## Identify individuals missing geographic data 
missing.ind <- coords[rowSums(is.na(coords)) > 0, ]   
missing.ind <- missing.ind$Sample 
missing.ind

# Create new dataframe without individuals missing geographic data
coords_nomissing<- filter(coords.tbl, !(Sample %in% missing.ind))
coords_nomissing <- as.data.frame(coords_nomissing)
coords_nomissing

# Drop individuals with missing geographic data
gendata@ind.names
gendata
missing.ind.list <- as.list(missing.ind)
gl2 <- gl.drop.ind(gendata, ind.list=missing.ind.list,recalc=TRUE)

################################################################################
## Run Mantel test (correlates two different distance matrices); need to convert DNA and geographic data into pairwise distances for this test
Dgeo<-earth.dist(coords[,c("Longitude", "Latitude")]) # get the geographic distances
Dgeo<-earth.dist(coords_nomissing[,c("Longitude", "Latitude")]) # get the geographic distances with no missing individuals

Dgen<-dist(gendata) # get the genetic distances
Dgen<-dist(gl2) # get the genetic distances without missing individuals

ibd<-mantel.randtest(Dgen,Dgeo) # run the mantel test
ibd

## Impute data for missing data if error message: Error in testmantel(nrepet, col, as.matrix(m1), as.matrix(m2)) : NA/NaN/Inf in foreign function call (arg 3)
Dgen1 <- impute(Dgen, "random")
ibd<-mantel.randtest(Dgen1,Dgeo) # run the mantel test with imputed data 

## We can then visualize the result of our empirical estimate of isolation by distance 
## compared to a permuted null distribution to see how significant our result is. 
## We'll also plot out a kernel density plot of genetic vs. geographic distances to 
## visualize how these distances are associated.

## PDF of mantel output and a kernel density plot of genetic vs. geograophic distance
Mantel_KD <-paste0("Mantel_KD_", basefile,".pdf")
pdf(file=Mantel_KD, width=8, height=8)
plot(ibd, main=paste0("mantel p = ", ibd$pvalue)) # plot out the IBD significance
## make kernel density plot of genetic and geographic distances
dens <- kde2d(as.numeric(Dgeo),as.numeric(Dgen), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5, xlab="Geographic Distance", ylab="Genetic Distance")
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.numeric(Dgen)~as.numeric(Dgeo)))
lines(loess.smooth(Dgeo, Dgen), col="red")
title("IBD plot")
dev.off()

## Plot With imputation
Mantel_KD <-paste0("Mantel_KD_withimputation_", basefile,".pdf")
pdf(file=Mantel_KD, width=8, height=8)
plot(ibd, main=paste0("mantel p = ", ibd$pvalue)) # plot out the IBD significance
## make kernel density plot of genetic and geographic distances
dens <- kde2d(as.numeric(Dgeo),as.numeric(Dgen1), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen1, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.numeric(Dgen1)~as.numeric(Dgeo)))
title("IBD plot")
dev.off()
## Result --> 2-page pdf of these plots ##
## Page 1: we can see that we have highly significant isolation by distance, 
## with the lowest possible significance value given our number of permutations 
## in the Mantel test (999 by default). The diamond-shaped point with the line
## indicates our empirical estimate, with the permutations shown as the gray histogram.
## Page 2: If we look at the second page, we can see how the geographic and genetic 
## distances are correlated. What we see is a positive relationship, but with 
## a major disjunction. This type of disjunction indicates the presence of some 
## level of not fully continuous spatial genetic structure.

################################################################################
################################## sNMF ########################################
################################################################################
## The snmf function requires a geno file as input, and requires that it has the 
## extension .geno. We want to use only unlinked SNPs here (i.e., 1 SNP per RAD 
## locus, assumed to be unlinked), and the geno file of unlinked snps has the 
## extension ugeno, so we’ll copy the file and give it a new extension
################################################################################

##########################
# Manage an snmf project #
##########################

# All the runs of snmf for a given file are 
# automatically saved into an snmf project directory and a file.
# The name of the snmfProject file is the same name as 
# the name of the input file with a .snmfProject extension 
# ("genotypes.snmfProject").
# The name of the snmfProject directory is the same name as
# the name of the input file with a .snmf extension ("genotypes.snmf/")
# There is only one snmf Project for each input file including all the runs.

# An snmfProject can be load in a different session.
project_path <- paste0(basefile, ".u.snmfProject")
project_path <- "/Users/hwsung/Library/CloudStorage/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly_v2/filteredVCF_ab/data/noreponly_v2.75.renamed.LDpruned.snmfProject"

obj.at = load.snmfProject(project_path) 
obj.at
################################################################################

## Use a regular expression substitution to generate the new file name
gendata_gl
geno_file <- gl2geno(gendata_gl, outfile = "testgenofile", outpath = "test/")
path_ugeno <- "/Users/hwsung/Library/CloudStorage/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly_v2/filteredVCF_ab/data/noreponly_v2.75.renamed.LDpruned.geno"
#file.copy(path_ugeno, path_geno) # do the copying with the new name

## Now we're ready to run sNMF. We'll run this using 1 to 10 possible ancestral 
## populations and evaluate the fit of these different numbers of populations 
## (referred to as k values) to the data using the cross entropy criterion.
obj.at <- snmf(input.file = path_ugeno,  # input file is the .geno format file
               K = 1:8, # we will test for k=1 through 8
               ploidy = 2, 
               entropy = T, # use the cross entropy criterion for assessing the best k 
               repetitions = 10, # Run 10 independent replicate analyses
               CPU = 4, 
               project = "new", tolerance = 0.00001, iterations = 10000)

## Let's make a pdf of the cross-entropy plot ##
snmf_cross_ent <-paste0(out_dir, "/LEA/", "snmf_cross_ent_", basefile,".pdf")
pdf(snmf_cross_ent, width = 8, height=5)
plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)
dev.off()

## As for DAPC, the best fit model is one with 2 populations, as shown by the 
## lowest cross-entropy score. We can also look at a numeric summary of this result:
outstats <- summary(obj.at)
outstats # best fit is 2 populations shown by lowest cross-entropy score
# the value for which the function plateaus or increases is our estimate of K 

## We can also confirm cross entropy values for K are consistent across runs and get the single best run for K=2
(ce <- cross.entropy(obj.at, K = 2))
ce ## best ce is run 6 at 0.3860456
(best.run <- which.min(ce)) # find the run with the lowest cross validation error

# Entropy value corresponding to best run for best K
e = round(ce[best.run], digits = 4) # rounding to four digits
e

(ce_3 <- cross.entropy(obj.at, K = 3))
(best.run_3 <- which.min(ce_3)) # find the run with the lowest cross validation error
# Entropy value corresponding to best run for best K
e = round(ce_3[best.run_3], digits = 4) # rounding to four digits
e

(ce_4 <- cross.entropy(obj.at, K = 4))
(best.run_4 <- which.min(ce_4)) # find the run with the lowest cross validation error

## Then we can get the snmf Q matrix from the best run at the best k, which is 
## a matrix of the proportion of ancestry that each sample derives from each population.
qmatrix <- Q(obj.at, K = 2, run = best.run)
qmatrix_K3 <- Q(obj.at, K = 3, run = best.run_3)
qmatrix_K4 <- Q(obj.at, K = 4, run = best.run_4)

# cluster assignment for each individual
cluster<- apply(qmatrix, 1, which.max) #this corresponds with the 1:24 order 
cluster
table(cluster)

## Ignore this - from admixture stuff #
#qmatrix_K3 <- read_csv("/Users/hwsung/Library/CloudStorage/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly_v2/filteredVCF_ab/admixture/admix_snmf_q3.csv")
#qmatrix_K3 
#qmatrix_K3 <- qmatrix_K3 %>% dplyr::select(Ad_P1,   Ad_P2,   Ad_P3)

clusterK3<- apply(qmatrix_K3, 1, which.max) #this corresponds with the 1:24 order 
clusterK3
table(clusterK3)

clusterK4<- apply(qmatrix_K4, 1, which.max) #this corresponds with the 1:24 order 
table(clusterK4)

colnames(qmatrix) <- paste0("P", 1:2)
colnames(qmatrix_K3) <- paste0("P", 1:3)
colnames(qmatrix_K4) <- paste0("P", 1:4)
qmatrix
qmatrix_K3
qmatrix_K3 <- cbind.data.frame(qmatrix_K3, clusterK3)
#names(qmatrix_K3)[9] <- "clusterK3_ad"
qmatrix_K3
#write.csv(qmatrix_K3, file = "admix_snmf_q3_clust.csv", quote = FALSE, row.names = FALSE)

#non_matching_rows <- qmatrix_K3$clusterK3 != qmatrix_K3$clusterK3_ad
#non_matching_rows
#non_matching_df <- qmatrix_K3[non_matching_rows, ]
#non_matching_df$Sample
#####################################################################
## Convert the qmatrix into a dataframe to order based on ancestry ##
## IGNORE ##
# Meh way
admix<-as.data.frame(qmatrix)
admix_K3<-as.data.frame(qmatrix_K3)
admix_K4<-as.data.frame(qmatrix_K4)

#sort by ascending order for P1
admix_ordered <- admix[order(admix$P1),]
admix_ordered3 <- admix_K3[order(admix_K3$P1, admix_K3$P2),]
#admix_ordered3 <- admix_K3[order(admix_K3[,1], -admix_K3[,3] ),]
admix_ordered4 <- admix_K4[order(admix_K4[,1], -admix_K4[,4]),]
#####################################################################
## Read in Admixture analyses ## 
# Q2
admixture_q2 <- read_delim("/filteredVCF_ab/admixture/K2-10/noreponly_v2.75.renamed.LDpruned.2.Q", col_names = F)
names(admixture_q2) <- c("Ad_P2", "Ad_P1")
admixture_q2
#admix_snmf_q2<- cbind.data.frame(q_df, admixture_q2)
admix_snmf_q2<- cbind.data.frame(pops, admixture_q2)

# Q3
admixture_q3 <- read_delim("/filteredVCF_ab/admixture/K2-10/noreponly_v2.75.renamed.LDpruned.3.Q", col_names = F)
names(admixture_q3) <- c("Ad_P1", "Ad_P2", "Ad_P3")
admixture_q3 
#admix_snmf_q3 <- cbind.data.frame(q_df3, admixture_q3)
admix_snmf_q3 <- cbind.data.frame(pops, admixture_q3)

# Q4 
admixture_q4 <- read_delim("/filteredVCF_ab/admixture/K2-10/noreponly_v2.75.renamed.LDpruned.4.Q", col_names = F)
names(admixture_q4) <- c("Ad_P4", "Ad_P3", "Ad_P1", "Ad_P2")
admixture_q4 
#admix_snmf_q4 <- cbind.data.frame(q_df4, admixture_q4)
admix_snmf_q4 <- cbind.data.frame(pops, admixture_q4)

setwd("/filteredVCF_ab/Admixture")
write.csv(admix_snmf_q2, file = "admix_snmf_q2.csv", row.names = FALSE)
write.csv(admix_snmf_q3, file = "admix_snmf_q3.csv", row.names = FALSE)
write.csv(admix_snmf_q4, file = "admix_snmf_q4.csv", row.names = FALSE)

# Order data by decending P1 and then ascending P2 for graphing
q_df <- arrange(q_df, desc(P1)) 
q_df <- mutate(q_df, order = 1:num_ind)
q_df

## Create csv file for qmatrix and metadata
q_mat <- q_df %>% left_join(pops, by = "Sample")
q_mat

qmatrix_csv<-paste0(out_dir,"/", basefile,"_qmatrix.csv")
qmatrix_csv
write.csv(q_mat, file = qmatrix_csv, row.names = FALSE)

## Create qmatrix dataframe for other K values ##
# convert the q matrix to a data frame
pops<- arrange(pops, Sample)

q_df3 <- qmatrix_K3 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Sample = pops$Sample,) 
q_df3

# Order data by decending P1 and then ascending P2 for graphing
q_df3 <- arrange(q_df3, desc(P1)) 
q_df3 <- mutate(q_df3, order = 1:num_ind)
q_df3

# Create csv file for qmatrix and metadata
q_mat3 <- q_df3 %>% left_join(pops, by = "Sample")
q_mat3

qmatrix3_csv<-paste0(out_dir,"/", basefile,"_qmatrix3.csv")
write.csv(q_mat3, file = qmatrix3_csv, row.names = FALSE)

q_df4 <- qmatrix_K4 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Sample = pops$Sample,) 
q_df4

# Order data by decending P1 and then ascending P2 for graphing
q_df4 <- arrange(q_df4, desc(P1)) 
q_df4 <- mutate(q_df4, order = 1:num_ind)
q_df4

# Create csv file for qmatrix and metadata
q_mat4 <- q_df4 %>% left_join(pops, by = "Sample")
q_mat4

qmatrix4_csv<-paste0(out_dir,"/", basefile,"_qmatrix4.csv")
write.csv(q_mat4, file = qmatrix4_csv, row.names = FALSE)

##########################################
## PLOTTING RESULTS: Map ##
##########################################
coords
pops

colors_4 <-c("#984EA3FF", "#4DAF4AFF", "#377EB8FF","#E41A1CFF")
colors_2 <- c("#E41A1CFF","#4DAF4AFF")
colors_3 <- c("#E41A1CFF", "#377EB8FF","#4DAF4AFF")

# For results from admixture 
admix.md <- read_csv("/filteredVCF_ab/admixture/popdata/popmap.csv")
admix.md

admix_k2 <- read_csv("/filteredVCF_ab/admixture/admix_snmf_q2.csv")

# cbind coord data with admixture data
admix <- admix_k2
admix.coord <- admix %>% left_join(coords.tbl, by = "Sample")
glimpse(admix.coord)

## Plot results onto a map like we did for our DAPC results ##
sNMF_Map_K2 <-paste0(out_dir, "/LEA/", "sNMF_Map_K2_", basefile,".pdf")
pdf(file=sNMF_Map_K2, width=8, height=10) # open the pdf plotting device
# Plot out the map
# map("worldHires", "Belize", xlim=c(-89.5,-87.5), ylim=c(15.8,18.5),col="gray90", fill=TRUE)
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix$P1[x],admix$P2[x]), radius=0.05, col=c("#377EB8FF","#E41A1CFF"))}
dev.off() # close the pdf plotting device


## Plot results onto map for admixture data ##
coords <- admix.coord

Ad_Map_K2 <-paste0(filteredVCF, "/admixture/figures/", "Ad_Map_K2_", basefile,".pdf")
pdf(file=Ad_Map_K2, width=8, height=10) # open the pdf plotting device
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix$Ad_P1[x],admix$Ad_P2[x]), radius=0.05, col=c("#377EB8FF","#E41A1CFF"))}
dev.off() # close the pdf plotting device

Ad_Map_K2 <-paste0(filteredVCF, "/admixture/figures/", "Ad_Map_K2_", basefile,".png")
png(file=Ad_Map_K2, width = 700, height = 700) # open the pdf plotting device
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix$Ad_P1[x],admix$Ad_P2[x]), radius=0.05, col=c("#377EB8FF","#E41A1CFF"))}
dev.off() # close the pdf plotting device
##################################################################################### 
admix_K3 <- read_csv("/filteredVCF_ab/admixture/admix_snmf_q3.csv")

# cbind coord data with admixture data
admix <- admix_K3
admix
coords <- admix.coord

admix.coord <- admix %>% left_join(coords.tbl, by = "Sample")
glimpse(admix.coord)

## Regular way - SNMF plot for K = 3 ##
sNMF_Map_K3 <-paste0(out_dir, "/LEA/", "sNMF_Map_K3_", basefile,".pdf")
pdf(file=sNMF_Map_K3, width=8, height=10) # open the pdf plotting device
# Plot out the map
# map("worldHires", "Belize", xlim=c(-89.5,-87.5), ylim=c(15.8,18.5),col="gray90", fill=TRUE)
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix_K3$P1[x],admix_K3$P2[x],admix_K3$P3[x]), radius=0.05, col=c("#E41A1CFF","#377EB8FF","#4DAF4AFF"))}
dev.off() # close the pdf plotting device

## Admixture plot ##
## Plot results onto map for admixture data ##
coords <- admix.coord

Ad_Map_K3 <-paste0(filteredVCF, "/admixture/figures/", "Ad_Map_K3_", basefile,".pdf")
pdf(file=Ad_Map_K3, width=8, height=10) # open the pdf plotting device
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
# Plot pies at each locality
for (x in 1:nrow(admix_snmf_q3)) {floating.pie(admix_snmf_q3$Longitude[x],admix_snmf_q3$Latitude[x],
                                               c(admix_snmf_q3$Ad_P1[x],admix_snmf_q3$Ad_P2[x],admix_snmf_q3$Ad_P3[x]), radius=0.05, col=c("#E41A1CFF","#4DAF4AFF","#377EB8FF"))}
dev.off() # close the pdf plotting device

Ad_Map_K3 <-paste0(filteredVCF, "/admixture/figures/", "Ad_Map_K3_", basefile,".png")
png(file=Ad_Map_K3, width = 700, height = 700) # open the pdf plotting device
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix_K3$Ad_P1[x],admix_K3$Ad_P2[x],admix_K3$Ad_P3[x]), radius=0.05, col=c("#E41A1CFF","#4DAF4AFF","#377EB8FF"))}
dev.off() # close the pdf plotting device
####################################################################################
admix_K4 <- read_csv("/filteredVCF_ab/admixture/admix_snmf_q4.csv")

# cbind coord data with admixture data
admix <- admix_K4
admix
coords <- admix.coord

admix.coord <- admix %>% left_join(coords.tbl, by = "Sample")
glimpse(admix.coord)

## SNMF plot for K = 4 ##
sNMF_Map_K4 <-paste0(out_dir, "/LEA/", "sNMF_Map_K4_", basefile,".pdf")
pdf(file=sNMF_Map_K4, width=10, height=10) # open the pdf plotting device
# Plot out the map
# map("worldHires", "Belize", xlim=c(-89.5,-87.5), ylim=c(15.8,18.5),col="gray90", fill=TRUE)
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix_K4$P1[x],admix_K4$P2[x],admix_K4$P3[x], admix_K4$P4[x]), radius=0.05,col=c("#E41A1CFF", "#984EA3FF","#4DAF4AFF","#377EB8FF"))}
dev.off() # close the pdf plotting device

## Admixture plot ##
## Plot results onto map for admixture data ##
coords <- admix.coord

Ad_Map_K4 <-paste0(filteredVCF, "/admixture/figures/", "Ad_Map_K4_", basefile,".pdf")
pdf(file=Ad_Map_K4, width=8, height=10) # open the pdf plotting device
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix_K4$Ad_P1[x],admix_K4$Ad_P2[x],admix_K4$Ad_P3[x], admix_K4$Ad_P4[x]), radius=0.05,col=c("#377EB8FF","#4DAF4AFF","#984EA3FF","#E41A1CFF"))}
dev.off() # close the pdf plotting device

Ad_Map_K4 <-paste0(filteredVCF, "/admixture/figures/", "Ad_Map_K4_", basefile,".png")
png(file=Ad_Map_K4, width = 700, height = 700) # open the pdf plotting device
par(mar=c(1,0.25,1, 0.25))
plot(Belize1)
# Plot pies at each locality
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix_K4$Ad_P1[x],admix_K4$Ad_P2[x],admix_K4$Ad_P3[x], admix_K4$Ad_P4[x]), radius=0.05,col=c("#377EB8FF","#4DAF4AFF","#984EA3FF","#E41A1CFF"))}
dev.off() # close the pdf plotting device

##########################################
## PLOTTING RESULTS: stacked barchart ##
##########################################
## Basic ggplot stacked barchart ##
# Example                                                                                                                                    -93L))
tbl # example data

plot_data <- tbl %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  theme_classic()

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  facet_grid(~likely_assignment, scales = 'free', space = 'free')

ggplot() + 
  geom_col(data = q_df, aes(x=order, ))
############################
# Now, you're transforming the data to a "long" format for plotting. The population column names get their own column and the ancestry proportions (q) get their own column.  
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 
q_df_long

which.max(q_df_long$q)

q_df_ordered <- q_df_long %>%
  # assign the population assignment according to the max q value (ancestry proportion) and include the assignment probability of the population assignment
  group_by(Sample) %>%
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>%
  # arrange the data set by the ancestry coefficients
  arrange(likely_assignment, desc(assignment_prob)) %>%
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  ungroup() %>%
  mutate(Sample = forcats::fct_inorder(factor(Sample)))
q_df_ordered

q_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = order, y = q, fill = pop)) +
  scale_fill_manual(values = colors_2, labels = c("CA", "CM")) +
  #scale_fill_viridis_d() +
  labs(fill = "Struc_pop") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  ) 
q_df_ordered

admix_ordered4


## Plot bar chart like standard structure plot 
barplot(t(admix_ordered),col = colors_2, border = NA, space = .2, xlab = "Individuals", ylab = "Admixture coefficients", main = "Ancestry matrix", horiz = FALSE, names.arg = q_df$Sample ,cex.names=0.4, las = 2)
barplot(t(admix_ordered3), col = colors_3, border = NA, space = .2,xlab = "Individuals", ylab = "Admixture coefficients", main = "Ancestry matrix", names.arg = q_df3$Sample,  cex.names=0.4, las = 2)
barplot(t(admix_ordered4), col = colors_4, border = NA, space = .2,xlab = "Individuals", ylab = "Admixture coefficients", main = "Ancestry matrix", names.arg = q_df4$Sample,  cex.names=0.4, las = 2)

colors_4 <-c("#984EA3FF", "#4DAF4AFF", "#377EB8FF","#E41A1CFF")
## create pdf of barchart
sNMF_boxplot_K2 <-paste0(out_dir, "/LEA/", "sNMF_barplot_K2_", basefile,".pdf")
pdf(file=sNMF_boxplot_K2, width=10, height=10) # open the pdf plotting device
bp <- barplot(t(admix_ordered), col = c("#377EB8FF","#E41A1CFF"), border = NA, space = .2, ylab = "Admixture coefficients", main = "Ancestry matrix",cex.names=0.4, las = 2,xaxt="n") 
dev.off() # close the pdf plotting device

sNMF_boxplot_K3 <-paste0(out_dir, "/LEA/", "sNMF_barplot_K3_", basefile,".pdf")
pdf(file=sNMF_boxplot_K3, width=10, height=10) # open the pdf plotting device
barplot(t(admix_ordered3), col = c("#E41A1CFF","#377EB8FF","#4DAF4AFF"), border = NA, space = .2,ylab = "Admixture coefficients", main = "Ancestry matrix", cex.names=0.4, las = 2, xaxt="n")
dev.off() # close the pdf plotting device

sNMF_boxplot_K4 <-paste0(out_dir, "/LEA/", "sNMF_barplot_K4_", basefile,".pdf")
pdf(file=sNMF_boxplot_K4, width=10, height=10) # open the pdf plotting device
barplot(t(admix_ordered4), col = colors_4, border = NA, space = .2,ylab = "Admixture coefficients", main = "Ancestry matrix", cex.names=0.4, las = 2, xaxt="n")
dev.off() # close the pdf plotting device

#####################################################################################
##Some formatting crap## 
q_df <- q_df %>% dplyr::select(P1, P2, P3, ID, Sample, order, Longitude, Latitude, Monitoring.Unit, Route, Site, Subdivision, Species)
head(q_df)
as_tibble(q_df)
glimpse(q_df)


q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

glimpse(q_df_long)

q_df_prates <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(ID = forcats::fct_inorder(factor(ID)))

q_df_prates

q_df_ordered <- q_df_long %>% 
  # assign the population assignment according to the max q value (ancestry proportion) and include the assignment probability of the population assignment
  group_by(ID) %>%
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>%
  # arrange the data set by the ancestry coefficients
  arrange(likely_assignment, assignment_prob) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  ungroup() %>% 
  mutate(ID = forcats::fct_inorder(factor(ID)))

q_df_ordered

# a custom palette for plotting
q_palette <- c("#fde725", "#35b779", "#440154")

q_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = ID, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette, labels = c("CA", "CM", "HY")) +
  #scale_fill_viridis_d() +
  labs(fill = "Species") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  )
#####################################################################################
#####################################################################################
## Organized barplot by ancestral coeff ##
# Adding ID names for plotting
individuals = pops$CrocID
q_mat$ID = individuals
q_mat$Sample = pops$Sample
head(q_mat)
q_mat %>% add_column(plot_order = 1:273)
pops

# Arrange order or samples in plot
df <- arrange(coords_md, desc(CrocID))
df <- df %>% add_column(plot_order = 1:273)
df$plot_order
df
admix$Samples
df$Sample

## match up the coordinates to the order of the individuals from genetic data
match_names<-match(admix$Samples, df$Sample)
df_reordered <- df %>% slice(match(admix$Samples, df$Sample))
df_reordered$plot_order

admix$plot_order = df_reordered$plot_order
admix_ord = plyr::arrange(admix, plot_order)

admix = as_tibble(admix) # dplyr likes the tibble format
admix <- q_df
admix_ord$ID = factor(admix_ord$ID, levels = unique(admix_ord$ID)) # transforming ID into factor to keep order in plot
admix_ord = as.data.frame(admix_ord)

# Save qmatrix
write.csv(admix, file = paste0("qmatrix_", basefile, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(q_df, file = paste0("qmatrix_", basefile, ".csv"), quote = FALSE, row.names = FALSE)

# "Melt" dataframe using the gather function, as required by ggplot
admix_melt = gather(admix, key = cluster, value = coeff, 1:2)

# Plot with ggplot
SNMF_plot = admix_melt %>% ggplot(aes(x= order())) +
  geom_bar(aes(x= order, y = coeff, fill = cluster), stat = "identity", position = "fill") +
  #scale_fill_manual("Population", values = myPalette[c(1:bestK)], labels = c("Atlantic Forest", "Amazonia 1", "Amazonia 2")) + 
  #scale_fill_manual("Population", values = myPalette[c(1:bestK)], labels = c("Amazonia", "Atlantic Forest")) + 
  scale_color_manual(values = colors_2, guide = "none") +
  theme_minimal(base_size = 16) +
  labs(y = "Ancestrality coefficient", x = "") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 5),
        panel.grid.major.x = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0))

SNMF_plot

#####################
## Diff method ## 
#####################

# plot the ancestry coefficients for the best run and K = 9
library(MetBrewer)

pdf(file=out_dir, "/LEA/", "sNMF_barchart.pdf", width=10, height=10) # open the pdf plotting device
bp <- barchart(obj.at, K = 2, run = best.run, sort.by.Q = TRUE,
               border = NA, space = 0, col = colors_2, 
               xlab = "Individuals", ylab = "Ancestry proportions", 
               main = "Ancestry matrix")

axis(1, at = 1:length(bp$order), 
     labels = bp$order, las = 3, 
     cex.axis = .4)
dev.off() # close the pdf plotting device
bp

b <- c(1:length(gendata_names))
b <- as.numeric(b)
df <- cbind(gendata_names, b)
df

match_df<-match(bp$order, df[,b])
match_df
df <- df[match_df,]
df
coords<-coords[match_coords,]
coords
coords1 <- na.omit(coords)
coords1

# An snmfProject can be exported to be imported in another directory
# or in another computer

export.snmfProject("dataset3.u.snmfProject")
