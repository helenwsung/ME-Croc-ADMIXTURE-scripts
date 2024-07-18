################################################################################
################# visualize LD decay results from plink in R ###################
################################################################################
## got data from genomics bootcamp tutorial and then used scripts from 
## https://eacooper400.github.io/gen8900/exercises/ld-2.html for analyses 
################################################################################
# my plink run:
# Create Plink format files
system("plink --vcf vcffile --double-id --allow-extra-chr --chr-set 69 --out outfile")
# Calculate R-squared values
system("plink --bfile bfile --double-id --allow-extra-chr --chr-set 69 --r2 gz --ld-window 10000 --ld-window-kb 5000 --ld-window-r2 0 --out bfile_r2")
# Calculate D-prime values 
system("plink --bfile bfile --double-id --allow-extra-chr --chr-set 69 --r2 dprime --out")
################################################################################

## Load libraries 
library("poppr")
library("magrittr")
library(dplyr)
library(tidyverse)

library(adegenet)
#library(gstudio)
library(tibble)
library(here)
library(vcfR)
library(pinfsc50)
library(utils)


################################################################################
# read in LD results file
data_dir <- "/datadir/filteredVCF_ab/LD/data_plink/plinkfiles_10Mb"
setwd(data_dir)

## noreponly_v2
out_dir <- "/datadir/filteredVCF_ab/LD/out_plink_10Mb"
################################################################################
## Read in plink output
LdValues <- read_table("noreponly_v2.75.renamed.pure_CA_r2.ld.gz")
basefile <- "pure_CA"
n <- 61 # samples

LdValues <- read_table("noreponly_v2.75.renamed.pure_CM_r2.ld.gz")
basefile <- "pure_CM"
n <- 90 # samples 

LdValues <- read_table("noreponly_v2.75.renamed.pure_MCA_r2.ld.gz")
basefile <- "pure_MCA"
n <- 17 

LdValues <- read_table("noreponly_v2.75.renamed.pure_acutus_r2.ld.gz")
basefile <- "pure_acutus"
n <- 88 # samples

LdValues <- read_table("noreponly_v2.75.renamed.All_Hybrid_r2.ld.gz")
basefile <- "All_Hybrid"
n <- 66 # samples 

LdValues <- read_table("noreponly_v2.75.renamed.F1_Hybrid_r2.ld.gz")
basefile <- "F1_Hybrid"
n <- 12 # 12 samples 

LdValues <- read_table("noreponly_v2.75.renamed.Hybrid_25.75_r2.ld.gz")
basefile <- "Hybrid_25.75"
n <- 23 # samples 


################################################################################
head(LdValues)

# calculate LD in 20 kb bins to display the trendline
averageLD <- LdValues %>%
  mutate(markerDistance = abs(BP_B - BP_A)/1000) %>%
  dplyr::filter(markerDistance < 10000) %>%
  mutate(intervals = cut_width(markerDistance, 20, boundary = 0)) %>%
  group_by(intervals) %>%
  summarise(averageR2 = mean(R2, na.rm = TRUE))

# calculate inter marker distances in 20kb bins
fullLD <- LdValues %>%
  mutate(markerDistance = abs(BP_B - BP_A)/1000) %>%
  dplyr::filter(markerDistance < 10000) %>%
  mutate(intervals = cut_width(markerDistance, 20, boundary = 0))

#merge the two data sets (full LD info and average per bin)
mergedLD <- full_join(fullLD,averageLD, by = "intervals")
names(mergedLD)

# visualize LD decay
LDdecayplot <- ggplot(mergedLD) +
  geom_point(aes(x=markerDistance, y=R2)) +
  geom_line(aes(x=markerDistance, y=averageR2), color="red", linewidth=2) +
  xlab("marker distance (kb)") + ylab("r-squared")+
  labs(title = "LD Decay (20kb bins)") +
  theme(plot.title = element_text(hjust = 0.5)) 
#LDdecayplot
#ggsave(filename = "LDdecay.pdf", plot = LDdecayplot)


# calculate LD in 2.5 kb bins to display the trendline
averageLD.2.5 <- LdValues %>%
  mutate(markerDistance = abs(BP_B - BP_A)/1000) %>%
  dplyr::filter(markerDistance < 5000) %>%
  mutate(intervals = cut_width(markerDistance, 2.5, boundary = 0)) %>%
  group_by(intervals) %>%
  summarise(across(R2, ~mean(., na.rm = TRUE), .names = "average{.col}"))

# calculate inter marker distances in 2.5 kb bins
fullLD.2.5 <- LdValues %>%
  mutate(markerDistance = abs(BP_B - BP_A)/1000) %>%
  dplyr::filter(markerDistance < 5000) %>%
  mutate(intervals = cut_width(markerDistance, 2.5, boundary = 0))

#merge the two data sets (full LD info and average per bin)
mergedLD.2.5 <- full_join(fullLD,averageLD, by = "intervals")
mergedLD.2.5

# visualize LD decay
LDdecayplot.2.5 <- ggplot(mergedLD.2.5) +
  geom_point(aes(x=markerDistance, y=R2)) +
  geom_line(aes(x=markerDistance, y=averageR2), color="red", linewidth=2) +
  xlab("marker distance (kb)") + ylab("r-squared")+
  labs(title = "LD Decay (20kb bins)") +
  theme(plot.title = element_text(hjust = 0.5)) 
#LDdecayplot

################################################################################
################### Find the Average Linkage Block Size ########################
################################################################################
my.means <- mergedLD$averageR2
bins <- mergedLD$markerDistance
LD.averages=data.frame(bins, my.means) # Create the results table (to hold calculations)
LD.averages
LD.averages$bins        # marker distances 
LD.averages$my.means    # average r2 

max(LD.averages$my.means)   # maximum r2 value (theoretically at distance 0)

# Sorting the data frame from decreasing to increasing by "bins"
LD.averages.sorted <- LD.averages[order(LD.averages$bins), ]
# Print the sorted data frame to get the min r2 at min value 
print(head(LD.averages.sorted))     


my.results <- data.frame(mergedLD$markerDistance, mergedLD$R2)
names(my.results) <- c("Distance", "Rsquared")
my.results


### Find the point when LD drops below 0.1
### The which statement tells me which rows of the table have means
### greater than 0.1,
### the brackets are giving me the subset of the "bins" with averages
### that correspond to these rows, and the max statement tells me which bin has
### the greatest distance value
LD.drop = max(LD.averages$bins[which(LD.averages$my.means>0.1)])
LD.drop
LD.drop <- round(LD.drop, 2)
LD.drop

#LDdecayplot + 
#  geom_vline(xintercept=LD.drop, color="blue", linetype="dashed")  # Add a vertical line corresponding to the drop off point

### Find the LD half life
LD.half = (max(LD.averages$my.means))/2 # Calculate half of the maximum average LD value
LD.half
LD.half.point = max(LD.averages$bins[which(LD.averages$my.means>LD.half)]) # Same strategy as when I found the 0.1 drop off point
LD.half.point
LD.half.point <- round(LD.half.point, 2)

#LDdecayplot + 
#  geom_vline(xintercept=LD.half.point, color="green", lty=2)  # Add a vertical line corresponding to LD half life


### Calculate rho based on Hill Weir equation:
# n <- number of samples
hill.weir.eq = (Rsquared~(((10+(rho*Distance))/((2+(rho*Distance))*(11+(rho*Distance))))*(1+(((3+(rho*Distance))*(12+(12*(rho*Distance))+((rho*Distance)**2)))/(n*(2+(rho*Distance))*(11+(rho*Distance))))))) # Define the formula for the Hill-Weir equation

rho.start=0.1 # Set an arbitrary starting value for rho
m=nls(formula=hill.weir.eq, data=my.results, start=list(rho=rho.start)) # Run the function to test different rho values and find the best fit to the data
results.m=summary(m) # Get a summary of results from the regression model
results.m

### Plot Expected R-squared
rho.estimate=results.m$parameters[1] # extract the rho estimate that produced the best fit
rho.estimate
# Round rho.estimate to 3 decimal places
rounded_rho <- round(rho.estimate, 3)

Distance=sort(my.results$Distance) # sort the distance results from lowest to highest

exp.rsquared=(((10+(rho.estimate*Distance))/((2+(rho.estimate*Distance))*(11+(rho.estimate*Distance))))*(1+(((3+(rho.estimate*Distance))*(12+(12*(rho.estimate*Distance))+((rho.estimate*Distance)**2)))/(n*(2+(rho.estimate*Distance))*(11+(rho.estimate*Distance)))))) # Use the best rho estimate inside of the Hill-Weir formula to calculate the expected r-squared at each distance value
exp.rsquared

print(results.m)

LDdecayplot <- ggplot(mergedLD) +
  geom_point(aes(x=markerDistance, y=R2)) +
  geom_line(aes(x=markerDistance, y=averageR2, color="LD Decay"), size=2) +  # Include color in aes() for legend
  geom_vline(aes(xintercept=LD.drop, color="LD drops below 0.1"), show.legend = TRUE) +  # Vertical line for drop-off point with legend
  geom_vline(aes(xintercept=LD.half.point, color="LD half point"), linetype="dashed", show.legend = TRUE) + # Dashed line for LD half-life with legend
  geom_line(aes(x=Distance, y=exp.rsquared, color="Population recombination rate"), size=2, show.legend = TRUE) + # Line for expected values with legend
  scale_color_manual(
    name = "",  # No need for a color legend title
    values = c("LD Decay" = "red", "LD drops below 0.1" = "dodgerblue", "LD half point" = "green", "Population recombination rate" = "purple"),
    labels = c(`LD Decay` = "LD Decay (20kb bins)",
               `LD drops below 0.1` = paste("LD drops below 0.1 =", LD.drop), 
               `LD half point` = paste("LD half point =", LD.half.point),
               `Population recombination rate` = paste("Population recombination rate (rho) =", rounded_rho))
  ) +
  labs(x="marker distance (kb)", y="r-squared") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position="bottom"  # Position legend at the bottom
  )
LDdecayplot

filename <- paste0(out_dir,"/", basefile, ".", "LDdecay_1MB_r2_", "20kb", ".png")
png(file=filename, width = 800, height = 800) # open the pdf plotting device
LDdecayplot
dev.off()

filename <- paste0(out_dir,"/", basefile, ".", "LDdecay_1MB_r2_", "20kb", ".pdf")
pdf(filename, width = 10, height = 10)
LDdecayplot
dev.off()
