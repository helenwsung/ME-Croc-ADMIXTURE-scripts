## This script takes the output of bootstrap reps and converts parameter estimates, and summarizes across the boot reps

# just modifications to Par_conv_FSC_getula.R
## Some comments from that script copied here 

#### This script takes the parameter estimates from the best runs of the best models for FSC runs on
####    datasets and converts the units and makes a table of them
####    this *should* work with any number of populations with divergence times, secondary contact times
####    and migration rates: requires that migration rates are specified such that the parameter names include the numbers of 
####    the populations between which migration is occurring and that these same numbers are included in the 
####    population names in a consistent way: e.g., MIG10 is migration bwtween POP1 and POP0; MIG102 is migration between 
####    an ancestral POP10 (or POP10anc - just so long as the numbers an only those numbers are there) amd POP2
####    (in an alternate tree structure, that could be migration between POP1 and POP02 instead - I beleive I have this correctly handled to parse out the correct pops)


#### This should be effective for any number of populations (of course, assuming they're structured exactly the way I 
####      name things in FSC), including migration among ancestral populations



library(dplyr)
library(stringr)
library(ggplot2)


setwd("/datadir/best_L_allMods_boot_16.2yr")

plots_out_dir <- "/datadir/best_L_allMods_boot_16.2yr/Boot_out"

## List all of the .bestlhoods files
files_lnl<-list.files(recursive=TRUE, pattern=".bestlhoods")
files_lnl

## Read in these files
output<-lapply(files_lnl, read.table, header=TRUE)
names(output)<-gsub("\\.bestlhoods", "", files_lnl) # name each element

gen_time <- 16.2
## Let's just get estimates for all of the models - this code will not universally convert, but will handle divergence times, population sizes, and present migration rates
##     have now added in functionality to do migration rates among ancestral populations - this works only if the migration rates are structure such that ancestral migration rates 
##     between population 1 and ancestor of pops 2&3 is written as MIG123$ or MIG132$ and the ancestral pop size of 
mods_convert<-names(output)
converted_ests<-vector("list", length(mods_convert)) # empty list to dump coverted output into 
for(j in 1:length(mods_convert)){ # loop over the models doing the conversions we want
  raw_est<-output[mods_convert[j]][[1]] # where the raw estimates are
  converted<-rep(NA, length(raw_est)-2) # Make an empty object that the converted numbers will end up in 
  names(converted)<-gsub("\\.", "", names(raw_est)[1:(length(raw_est)-2)]) # add on the parameter names
  names(raw_est)<-gsub("\\.", "", names(raw_est)) ## Strip out the . in the original names, too
  
  ### Go through and populate the new matrix with values converted from the raw values
  ind_npop<-grep("NPOP", names(converted)) # indices of elements that are population sizes
  converted[ind_npop]<-round(raw_est[ind_npop]/2, -2) # pop sizes get divided by 2 to be diploid and round to hundreds place
  ind_times<-grep("TDIV|TCONT|TSEP|TRESIZE", names(converted))  # index of elements that are times (divergences or secondary contact)
  converted[ind_times]<-round(raw_est[ind_times]*gen_time, -2) ## times get multiplied by 3 years for generation time - round to hundreds place also
  # migration rates are multipled by 2 and by population size of the population that the migration is into (e.g., MIG20 is migration 2->0 and is multiplied by pop size 0 NPOP0)
  # this is slightly more complex, so handle with a loop
  ind_mig<-grep("MIG", names(converted)) # indices of elements that are migration rates
  ind_migA<-grep("MIGA", names(converted))# after handling standard migration rates, handle migration rates among the same pops at an earlier time - I designate these MIGA
  if(length(ind_migA)>0){ # if we have MIGA in here, then these need to be removed from ind_mig and handled separately
    ind_mig <- ind_mig[!ind_mig %in% ind_migA]
  }
  # for handling any migration among ancestral populations, need to get out the number of digits in each population name to identify ancestral populations - do this outside of the loop below to avoid repeating calculations
  pop_names<-names(converted)[ind_npop] # get population names
  nums_pops<-str_split(str_extract(pop_names, "[[:digit:]]+"), pattern="")
  num_digits_pops<-sapply(nums_pops, function(x) sum(!is.na(x))) # get the number of numerical digits that are not NA
  ind_pop_names_anc<-which(num_digits_pops>1) # indices of pop names that are ancestral pops
  dig_anc<-nums_pops[ind_pop_names_anc]  # digits in ancestral pops
  names(dig_anc)<-pop_names[ind_pop_names_anc] # name those
  for(i in ind_mig){
    if(exists("pop1_anc")){ # if we've run this function before or a previous iteration of the loop defined pop1_anc and pop2_anc, clear them out
      rm(pop1_anc)
    }
    if(exists("pop2_anc")){
      rm(pop2_anc) 
    }
    ## Need to get and count the number of digits in the migration rate name, if it's more than 2, then ancestral populations are involved
    ##     this works for 4 popuations where the tree is symmetric right now... gets even more complicated past that
    nums_mig<-str_split(str_extract(names(converted)[i], "[[:digit:]]+"), pattern="")[[1]]
    if(length(nums_mig)==2){
      pop_mult<-gsub("MIG.", "NPOP", names(converted)[i]) # for each index, replace MIG. (i.e.,. MIG and the following character) with NPOP resulting in e.g., changing MIG20 to NPOP0, the population we want to multiply by
      converted[i]<-round(raw_est[i]*2*raw_est[pop_mult], 2) # multiply the migration rate by 2 and the pop size of interest - round to 2 decimal places
    }
    if(length(nums_mig)>2){ 
      ## we don't know where the split is among the two ancestral populations, so we'll try all possible groupings
      for(splt in 2:length(nums_mig)){
        ptpop1<-nums_mig[1:splt-1] #potential ancestral population 1
        ptpop2<-nums_mig[splt:length(nums_mig)] # potetential ancestral population 2
        ## Now we need to test if each of these populations exists
        if(length(ptpop1)>1){ # if there is more than 1 number, check if all of them are in any of the population names
          for(m in 1:length(dig_anc)){
            anc<-unlist(dig_anc[m])
            # if(length(dig_anc)==1){  ## Think this was a failed attempt at handling what I handle with unlist in the previous line now
            #   anc<-anc[[1]]
            # }
            if(sum(ptpop1 %in% anc)==length(ptpop1) && length(ptpop1)==length(anc)){ # check that all numbers in the potential population based on splitting out migration numbers are in the numbers of an ancestral population, and that the lengths of these are the same (e.g., without checking the lengths, we could find 01 in an ancestral populaion 012, but that's not what we're looking for)
              pop1_anc<-names(dig_anc)[m] # if we've found the right ancestral population, assign it to pop1_anc object
            }
          }
        }else{ # if there is not more than 1 number (i.e., there is only 1) then that number is a single present population and necessarily exists
          pop1_anc<-paste0("NPOP", ptpop1)
        }
        if(length(ptpop2)>1){
          for(m in 1:length(dig_anc)){
            anc<-unlist(dig_anc[m])
            # if(length(dig_anc)==1){  ## Think this was a failed attempt at handling what I handle with unlist in the previous line now
            #   anc<-anc[[1]]
            # }
            if(sum(ptpop2 %in% anc)==length(ptpop2) && length(ptpop2)==length(anc)){ # check that all numbers in the potential population based on splitting out migration numbers are in the numbers of an ancestral population, and that the lengths of these are the same (e.g., without checking the lengths, we could find 01 in an ancestral populaion 012, but that's not what we're looking for)
              pop2_anc<-names(dig_anc)[m] # if we've found the right ancestral population, assign it to pop1_anc object
            }
          }
        }else{
          pop2_anc<-paste0("NPOP", ptpop2)
        }
        if(exists("pop1_anc") && exists("pop2_anc")){ # for a given split, if both pops exist, then break and we're done
          break
        }else{ # if not, remove these objects (if either exists) - in any case where we split out the first or last digit, that population definitely does exist, but we only want to keep COMBINATIONS where both populations exist
          if(exists("pop1_anc")){
            rm(pop1_anc)
          }
          if(exists("pop2_anc")){
            rm(pop2_anc) 
          }
        }
      }
      converted[i]<-round(raw_est[i]*2*raw_est[pop2_anc], 2) # multiply the migration rate by 2 and the pop size of interest - round to 2 decimal places
    }
  }
  if(length(ind_migA)>0){
    for(i in ind_migA){
      if(exists("pop1_anc")){ # if we've run this function before or a previous iteration of the loop defined pop1_anc and pop2_anc, clear them out
        rm(pop1_anc)
      }
      if(exists("pop2_anc")){
        rm(pop2_anc) 
      }
      ## Need to get and count the number of digits in the migration rate name, if it's more than 2, then ancestral populations are involved
      ##     this works for 4 popuations where the tree is symmetric right now... gets even more complicated past that
      nums_mig<-str_split(str_extract(names(converted)[i], "[[:digit:]]+"), pattern="")[[1]]
      if(length(nums_mig)==2){
        pop_mult<-gsub("MIGA.", "NPOP", names(converted)[i]) # for each index, replace MIGA. (i.e.,. MIGA and the following character) with NPOP resulting in e.g., changing MIGA20 to NPOP0, the population we want to multiply by
        converted[i]<-round(raw_est[i]*2*raw_est[pop_mult], 2) # multiply the migration rate by 2 and the pop size of interest - round to 2 decimal places
      }
      if(length(nums_mig)>2){ 
        ## we don't know where the split is among the two ancestral populations, so we'll try all possible groupings
        for(splt in 2:length(nums_mig)){
          ptpop1<-nums_mig[1:splt-1] #potential ancestral population 1
          ptpop2<-nums_mig[splt:length(nums_mig)] # potetential ancestral population 2
          ## Now we need to test if each of these populations exists
          if(length(ptpop1)>1){ # if there is more than 1 number, check if all of them are in any of the population names
            for(m in 1:length(dig_anc)){
              anc<-unlist(dig_anc[m])
              if(length(dig_anc)==1){
                anc<-anc[[1]]
              }
              if(sum(ptpop1 %in% anc)==length(ptpop1) && length(ptpop1)==length(anc)){ # check that all numbers in the potential population based on splitting out migration numbers are in the numbers of an ancestral population, and that the lengths of these are the same (e.g., without checking the lengths, we could find 01 in an ancestral populaion 012, but that's not what we're looking for)
                pop1_anc<-names(dig_anc)[m] # if we've found the right ancestral population, assign it to pop1_anc object
              }
            }
          }else{ # if there is not more than 1 number (i.e., there is only 1) then that number is a single present population and necessarily exists
            pop1_anc<-paste0("NPOP", ptpop1)
          }
          if(length(ptpop2)>1){
            for(m in 1:length(dig_anc)){
              anc<-unlist(dig_anc[m])
              if(length(dig_anc)==1){
                anc<-anc[[1]]
              }
              if(sum(ptpop2 %in% anc)==length(ptpop2) && length(ptpop2)==length(anc)){ # check that all numbers in the potential population based on splitting out migration numbers are in the numbers of an ancestral population, and that the lengths of these are the same (e.g., without checking the lengths, we could find 01 in an ancestral populaion 012, but that's not what we're looking for)
                pop2_anc<-names(dig_anc)[m] # if we've found the right ancestral population, assign it to pop1_anc object
              }
            }
          }else{
            pop2_anc<-paste0("NPOP", ptpop2)
          }
          if(exists("pop1_anc") && exists("pop2_anc")){ # for a given split, if both pops exist, then break and we're done
            break
          }else{ # if not, remove these objects (if either exists) - in any case where we split out the first or last digit, that population definitely does exist, but we only want to keep COMBINATIONS where both populations exist
            if(exists("pop1_anc")){
              rm(pop1_anc)
            }
            if(exists("pop2_anc")){
              rm(pop2_anc) 
            }
          }
        }
        converted[i]<-round(raw_est[i]*2*raw_est[pop2_anc], 2) # multiply the migration rate by 2 and the pop size of interest - round to 2 decimal places
      }
    }
  }
  if("MUTRATE" %in% names(converted)){ # if MUTRATE is a parameter (i.e., if we let FSC estimate the mutation rate rather than fixing it)
    converted["MUTRATE"] <- raw_est["MUTRATE"]   ## then just add in the mutation rate - note that this is still mutation rate per generation
  }  
  converted_ests[[j]]<-unlist(converted)
  names(converted_ests)[[j]]<-mods_convert[j]
}


converted_table<-as.data.frame(do.call(bind_rows, converted_ests)) # combine into 1 big table
rownames(converted_table)<-names(converted_ests) # add rownames onto it

## get the upper and lower ranges of each par
max<-apply(converted_table, 2, max)
min<-apply(converted_table, 2, min)
par_ranges<-rbind(max, min)
rownames(par_ranges)<-c("max", "min")
write.csv(par_ranges, "boot_ranges_16.2.csv")




## Plot some of these out as distributions
colnames(converted_table)


## Read in the maximum pseudolikelihood estimates for each parameter
#point_ests_all <- read.csv("/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/FSC/Clust_out_best_runs/fix_root_time/K3/Lgetula_Kk3_fix_root_time_FSC_estimates.csv")
point_ests_all <- read.csv("/datadir/filteredVCF_ab/fastsimcoal/best_L_allMods_16.2yr/noreponly_v2_16.2yr_pars_conv.csv")

ests <- point_ests_all[1,]  # reduce this down to just the best fit model

# Set up times
  time_ests <- ests[c("TDIV1", "TRESIZE")]
  tests_2 <- as.data.frame(t(time_ests))
  tests_3 <- cbind(rownames(tests_2), tests_2)
  colnames(tests_3) <- c("Event", "Time")
  #tests_3$Event <- gsub("TCONT2", "TCont2", tests_3$Event)
  tests_3$Event <- gsub("TDIV1", "TDiv1", tests_3$Event)
  #tests_3$Event <- gsub("TCONT1", "TCont1", tests_3$Event)
  tests_3$Event <- gsub("TRESIZE", "TResize", tests_3$Event)
  
# set up migration rates
  mig_ests <- ests[c("MIG02", "MIG20", "MIG12", "MIG21", "MIG10","MIG012", "MIG120")]
  mests_2 <- as.data.frame(t(mig_ests))
  mests_3 <- cbind(rownames(mests_2), mests_2)
  colnames(mests_3) <- c("Rate", "IndperGen")
  mests_3$Rate <- c(paste0("moreletii", "\U2192", "acutus_A"), 
                    paste0("acutus_A", "\U2192", "moreletii"), 
                    paste0("acutus_B", "\U2192", "acutus_A"), 
                    paste0("acutus_A", "\U2192", "acutus_B"), 
                    paste0("acutus_B", "\U2192", "moreletii") , 
                    paste0("moreletii", "\U2192", "acutus"),
                    paste0("acutus", "\U2192", "moreletii"))
  
  

# Set up pop sizes
  pop_ests <- ests[c("NPOProotanc", "NPOP12anc", "NPOP0PB", "NPOP1PB", "NPOP2PB", "NPOP0", "NPOP1", "NPOP2")]
  pests_2 <- as.data.frame(t(pop_ests))
  pests_3 <- cbind(rownames(pests_2), pests_2)
  colnames(pests_3) <- c("Population", "Individuals")
  pests_3$Population <- c("RootAnc", "acutus_Anc", "moreletii PreResize", "acutus_B PreResize", "acutus_A PreResize", "moreletii", "acutus_B","acutus_A")
  
################################################################################
################################################################################
### Plot out times as violin plots
#tcont2<-cbind(converted_table$TCONT2, "TCont2")
tdiv1<-cbind(converted_table$TDIV1, "TDiv1")
#tcont1<-cbind(converted_table$TCONT1, "TCont1")
tresize<-cbind(converted_table$TRESIZE, "TResize")

#times <- rbind(tcont2, tdiv1, tcont1, tresize)
times <- rbind(tdiv1, tresize)
colnames(times) <- c("Time", "Event")
times <- as.data.frame(times)
times$Time <- as.numeric(times$Time)
times$Event <- factor(times$Event , levels=c("TDiv1", "TResize"))

ggplot(times, aes(x=Event, y=Time)) + 
  geom_violin(trim=FALSE, fill="blue") +
  theme_minimal()

##############################################
## add in the best estimates - VIOLIN PLOT ##
##############################################
times.plot<-ggplot(times, aes(x=Event, y=Time)) + 
  geom_violin(trim=FALSE, fill="lightblue") +
  geom_point(data = tests_3, aes(x = Event, y = Time), size = 3, shape = 23, fill = "black") +
  scale_x_discrete(labels = c("TDiv1" = "Time of Divergence",
                              "TResize" = "Time of Population Resize")) +
  xlab("") +
  ylab("Time (years)") +
  theme_minimal()
times.plot

plots_out_dir

ggsave(
  filename = paste0(plots_out_dir,"/FSCboot_times.plot.png"),
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", times.plot
)

ggsave(
  filename = paste0(plots_out_dir,"/FSCboot_times.plot.pdf"),
  width = 8.5, height = 7, dpi = 600,
  device = "pdf", bg = "white", times.plot
)

pdf(file="Boot_out/noreponly_v2_16.2yr_FSCboot_times.pdf", width=5.5, height=5.5)
times.plot
dev.off()

png(file="Boot_out/noreponly_v2_16.2yr_FSCboot_times.png")
times.plot
dev.off()

################################################################################
################################################################################
### Plot out migration as violin plots
MIG02<-cbind(converted_table$MIG02, paste0("moreletii", "\U2192", "acutus_A"))    # MIG02
MIG20<-cbind(converted_table$MIG20, paste0("acutus_A", "\U2192", "moreletii"))    # MIG20
MIG12<-cbind(converted_table$MIG12, paste0("acutus_B", "\U2192", "acutus_A"))    # MIG12
MIG21<-cbind(converted_table$MIG21, paste0("acutus_A", "\U2192", "acutus_B"))    # MIG21
MIG10<-cbind(converted_table$MIG10, paste0("acutus_B", "\U2192", "moreletii"))    # MIG10
MIG012<-cbind(converted_table$MIG012, paste0("moreletii", "\U2192", "acutus")) # MIG012
MIG120<-cbind(converted_table$MIG120, paste0("acutus", "\U2192", "moreletii")) # MIG120

migs <- rbind(MIG02, MIG20, MIG12, MIG21, MIG10,MIG012, MIG120)
colnames(migs) <- c("IndperGen", "Rate")
migs <- as.data.frame(migs)
migs$IndperGen <- as.numeric(migs$IndperGen)
#migs$Rate <- as.factor(migs$Rate)
migs$Rate <- factor(migs$Rate , levels=c(paste0("acutus_B", "\U2192", "acutus_A"),
                                         paste0("acutus_A", "\U2192", "acutus_B"),
                                         paste0("acutus_A", "\U2192", "moreletii"),
                                         paste0("moreletii","\U2192", "acutus_A"),
                                         paste0("acutus_B", "\U2192", "moreletii"),
                                         paste0("acutus", "\U2192", "moreletii"),
                                         paste0("moreletii", "\U2192", "acutus")))

## set the ylim lower - there are some extreme outliers for mig120
migs.plot <- ggplot(migs, aes(x=Rate, y=IndperGen)) + 
  geom_violin(trim=FALSE, fill="green") +
  xlab("") +
  #ylab("Migration Rate (Individuals per Generation)") +
  labs(y="Migration Rate \n (Individuals per generation)") +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", angle = 325, vjust=0.3)) +
  ylim(0,15)

migs.plot
##############################################
## add in the best estimates - VIOLIN PLOT ##
##############################################
pdf(file="Boot_out/croc_FSCboot_mig.pdf", width = 10, height = 10)
png(filename = "Boot_out/croc_FSCboot_mig.png",width = 800, height = 800)

mig.plot <- ggplot(migs, aes(x=Rate, y=IndperGen)) + 
  geom_violin(trim=FALSE, fill="green") +
  # ylim(0,30) +
  geom_point(data = mests_3, aes(x = Rate, y = IndperGen), size = 3, shape = 23, fill = "black") +
  scale_y_continuous(breaks = seq(0, 15, by = 5), minor_breaks = seq(0 , 15, 1), limits = c(0, 15)) +
  labs(y="Migration Rate \n (Individuals per generation)") +
  xlab(NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", angle = 325, vjust=0.3)) 
mig.plot

dev.off()

pdf(file="Boot_out/croc_FSCboot_mig.pdf", width = 10, height = 10)
mig.plot
dev.off()

png(filename = "Boot_out/croc_FSCboot_mig.png",width = 800, height = 800)
mig.plot
dev.off()

ggsave(
  filename = paste0(plots_out_dir,"/FSCboot_mig.plot.png"),
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", mig.plot
)

ggsave(
  filename = paste0(plots_out_dir,"/FSCboot_mig.plot.pdf"),
  width = 8.5, height = 7, dpi = 600,
  device = "pdf", bg = "white", mig.plot
)
################################################################################
################################################################################
### Plot out pop sizes as violin plots
NPOProotanc<-cbind(converted_table$NPOProotanc, "RootAnc")
NPOP12anc<-cbind(converted_table$NPOP12anc, "acutus_Anc")
NPOP0PB<-cbind(converted_table$NPOP0PB, "moreletii PreResize")
NPOP1PB<-cbind(converted_table$NPOP1PB, "acutus_B PreResize")
NPOP2PB<-cbind(converted_table$NPOP2PB, "acutus_A PreResize")
NPOP0<-cbind(converted_table$NPOP0, "moreletii")
NPOP1<-cbind(converted_table$NPOP1, "acutus_B")
NPOP2<-cbind(converted_table$NPOP2, "acutus_A")

popsizes <- rbind(NPOProotanc, NPOP12anc, NPOP0PB, NPOP1PB, NPOP2PB, NPOP0, NPOP1, NPOP2)
colnames(popsizes) <- c("Individuals", "Population")
popsizes <- as.data.frame(popsizes)
popsizes$Individuals <- as.numeric(popsizes$Individuals)
popsizes$Population <- factor(popsizes$Population , levels=c("acutus_B", "acutus_A", "moreletii", "acutus_B PreResize", "acutus_A PreResize", "moreletii PreResize", "acutus_Anc", "RootAnc"))
#popsizes$Population <- as.factor(popsizes$Population)

ggplot(popsizes, aes(x=Population, y=Individuals)) + 
  geom_violin(trim=FALSE, fill="pink") +
  theme_minimal()

##############################################
## add in the best estimates - VIOLIN PLOT ##
##############################################

ggplot(popsizes, aes(x=Population, y=Individuals)) + 
  geom_violin(trim=FALSE, fill="pink") +
  geom_point(data = pests_3, aes(x = Population, y = Individuals), size = 3, shape = 23, fill = "black") +
  theme_minimal()

## edit this so that units are thousands
popsizes_thou <- popsizes
popsizes_thou$Individuals <- popsizes_thou$Individuals/1000
pests_4 <- pests_3
pests_4$Individuals <- pests_4$Individuals/1000

pop.plot <- ggplot(popsizes_thou, aes(x=Population, y=Individuals)) + 
  geom_violin(trim=FALSE, fill="pink") +
  geom_point(data = pests_4, aes(x = Population, y = Individuals), size = 3, shape = 23, fill = "black") +
  labs(y="Effective Population size \n (1,000 individuals)") +
  scale_y_continuous(breaks = seq(0, 550, by = 100), minor_breaks = seq(0 , 550, 50)) +
  xlab(NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", angle = 325, vjust=0.3)) 
pop.plot

dev.off()


ggsave(
  filename = paste0(plots_out_dir,"/FSCboot_pop.plot.png"),
  width = 8.5, height = 7, dpi = 600,
  device = "png", bg = "white", pop.plot
)

ggsave(
  filename = paste0(plots_out_dir,"/FSCboot_pop.plot.pdf"),
  width = 8.5, height = 7, dpi = 600,
  device = "pdf", bg = "white", pop.plot
)

pdf(file="Boot_out/croc_FSCboot_pop.pdf", width = 10, height = 10)
pop.plot
dev.off()

png(filename = "Boot_out/croc_FSCboot_pop.png",width = 800, height = 800)
pop.plot
dev.off()

##############################
## Put all 3 plots together ##
##############################
library(ggplot2)
library(gridExtra)

p1 <- pop.plot
p2 <- mig.plot
p3 <- times.plot

# Stack the plots
grid.arrange(pop.plot, mig.plot, times.plot, ncol=1)

library(patchwork)

# Labels for each plot
pop.plot + mig.plot + times.plot + plot_annotation(tag_levels = "A") 

# Stack the plots
(p1 | p2 | p3) & plot_layout(ncol = 1)
dev.off()

library(ggplot2)
library(cowplot)
plot_grid(p1, p2, p3) 
# Create an empty plot for spacing
spacer <- ggplot() + theme_void()

pdf(file="Boot_out/croc_FSCboot_all.pdf", width = 10, height = 10)
png(filename = "Boot_out/croc_FSCboot_all.png",width = 800, height = 800)
fsc.plots <- plot_grid(p1, spacer, p2, spacer, p3, ncol = 1, labels = c("A", "", "B", "", "C"), rel_heights = c(1, 0.1, 1, 0.1, 1))
dev.off()

pdf(file="Boot_out/croc_FSCboot_all.pdf", width = 10, height = 10)
fsc.plots
dev.off()

png(filename = "Boot_out/croc_FSCboot_all.png",width = 800, height = 800)
fsc.plots
dev.off()

# Use plot_grid with spacers between plots
plot_grid(p1, spacer, p2, spacer, p3, ncol = 1, labels = c("A", "", "B", "", "C"), rel_heights = c(1, 0.1, 1, 0.1, 1))

ggsave(
  filename = paste0(plots_out_dir,"/FSCboot.plot.png"),
  width = 10, height = 10, dpi = 300,
  device = "png", bg = "white", fsc.plots
)

ggsave(
  filename = paste0(plots_out_dir,"/FSCboot.plot.pdf"),
  width = 10, height = 10, dpi = 600,
  device = "pdf", bg = "white", fsc.plots
)



