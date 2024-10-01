##########################################################################################
########################## Demographic modeling: Fastsimcoal2 ############################
##########################################################################################
module load lang/Anaconda3
module load lang/R
module load lang/Python/2.7.18-GCCcore-11.3.0-bare

## easySFS for fastsimcoal2 ##
$ source activate easySFS
	$ conda install -c conda-forge numpy pandas scipy -y		# easySFS dependencies

source activate easySFS

## clone entire repository 
cd /yourdir/analyses/fastsimcoal/easySFS
git clone https://github.com/isaacovercast/easySFS.git

## make executable 
chmod +x ./easySFS.py
ls -ltrh 

## Your Data ##
VCF=/data/noreponly_v2.filtered.75.renamed.vcf.gz
POPMAP=popmap.txt

##########################################################################################
############################## Preparing the input files #################################
##########################################################################################
## To run fastsimcoal2 we need three input files which are all named in a consistent way. 
# They are all just plain text files:
	# observed SFS - ${PREFIX}_jointDAFpop1_0.obs
	# template file defining the demographic model - ${PREFIX}.tpl
	# estimation file defining the parameters - ${PREFIX}.est

##########################################################################################
#################### easySFS: multidimensional SFS: '_MSFS.obs' ##########################
##########################################################################################
# For easySFS, to identify the best projection value - Run Preview to see if works - 
# Each column is the number of samples in the projection and the number of segregating sites
# at that projection value, the dadi manual recommends maximizing the number of segregating 
# sites, but at the same time if you have lots of missing data then you might have to 
# balance # of segregating sites against # of samples to avoid downsampling too far this 
# line will calculate the necessary SFS files - the SFS cannot be computed with sites that 
# contain missing data. By projecting down the number of individuals per population, more 
# sites can be recovered

##################
## noreponly_v2 ##
##################

## Ad_pop_90:
./easySFS.py -i ../data/noreponly_v2/$VCF -p ../data/noreponly_v2/$POPMAP --preview -a
			# pop0 = CM
			# pop1 = CA
			# pop2 = MCA

# picked the largest projection segregating sites
	./easySFS.py -i ../data/noreponly_v2/$VCF -p ../data/noreponly_v2/$POPMAP --proj=130,84,24 -a -o ./noreponly_v2.filtered.75_16.4y_SFS/ -f

## K=3 Clusters:
./easySFS.py -i ../data/noreponly_v2/$VCF -p ../data/noreponly_v2/$POPMAP --preview -a
			# pop0 = CM => (164, 30022)
			# pop1 = CA => (100, 24748)
			# pop2 = MCA => (70, 26996)

##########################################################################################
########################## Preparing Models: .tpl & .est files ###########################
##########################################################################################

######################
## 16 Tested Models ##
######################

curMK3 = Current Migration, no acutus_B->M        
curMK3_CACM = Current Migration
curM_RecResize_CACM = Current Migration, recent resize
curM_Resize_CACM = Current Migration, old resize

MallK3 = All migration, no acutus_B->M     
MallAsymK3_CACM = All migration
MallAsym_RecResize_CACM = All migration, recent resize
MallAsymResize_CACM = All migration, old resize

noMK3 = No migration
HistSec_CACM = Historic Secondary Contact
HistSecRecResize_CACM = Historic Secondary Contact, recent resize
HistSecResize_CACM = Historic Secondary Contact, old resize

SecK3 = Recent Secondary Contact, no acutus_B->M         
Sec_CACM = Recent Secondary Contact
Sec_RecResize_CACM = Recent Secondary Contact, recent resize
Sec_Resize_CACM = Recent Secondary Contact, old resize


################################################################################################################################
## NOTE: YEARS PER GENERATION = C. acutus --> average generation time = 16.2 years (Briggs-Gonzalez et al 2017)
## CHANGE TO MAKE YEARS PER GENERATION = C. acutus --> average generation time = 16.2 years (Briggs-Gonzalez et al 2017)
# [RULES]
# TDIV1$ < 246914
246,914 generations = (4 my/16.2 yo) 

## Other potential dates: ##
# Upper bound:
500000 generations (4 million years divergence/8 years per generation) 
	# More than 5% of female population os C. acutus were predicted to be reproductive by age of 8yo (Briggs-Gonzalez et al. 2017)
# Lower bound: 
160,000 generations = (4 my/25 yo)

# Can also try 18yo b/c more than 95% of female population were predicted to be reproductive by 18yo
222,222 generations = (4 my/18 yo)
################################################################################################################################

## Run script from inside Models folder
cd /datadir/fastsimcoal/Models # make folder if haven't already

# copy easysfs file into Models folder
cp /datadir/fastsimcoal/easySFS/noreponly_v2.filtered.75_16.2y_SFS/fastsimcoal2/noreponly_v2_MSFS.obs ./

sh Prep_FSC_reps.sh MSFS 
# uses script: fsc-selectbestrun.sh

## edit .slurm script to array --> 50 runs * 16 models = 800
sbatch FS_FixRootTime_K3_Mods.slurm 
	# come back tomorrow when slurm script ran 
	
## Added fsc26 executable to home directory --> edit .bash_profile to add executable to bash profile --> chmod +x /home/hwsung/fsc26 --> source ./bash_profile 
cd /datadir/fastsimcoal/noreponly_v2/Reps

# run script 
sh Get_best_FSCacross_mods.sh # makes folder best_L_allMods in /Reps/

cd /datadir/fastsimcoal/noreponly_v2/Reps/best_L_allMods

Rscript ~/scripts/Get_AIC_across_mods.R 

################################################################################################################################
## Get_AIC_across_mods.R SCRIPT ##
################################################################################################################################
#### This script reads in multiple bestlhoods files and calculates the AIC and deltaAIC for each
################################################################################################################################
####  NOTE & WARNING: This will only work if you output 1 and exactly 1 parameter to this file for each estimated parameter
####     i.e., if you output a complex parameter AND the estimated parameters, your number of parameters will
####     be wrong. Or if you estimated a parameter but didn't output it, your number of parameters will also be 
####    wrong -- e.g., if you don't care about population sizes and hid these from the output for some reason
################################################################################################################################
model<-list.files(pattern="bestlhoods") # list out the bestlhoods files
output<-lapply(model, read.table, header=TRUE) # read in the actual values of parameters and likelihoods
lhoods<-sapply(output, function(x) x$MaxEstLhood) # get just the estimated likelihood
npars<-sapply(output, function(x) length(x)-2)  # get the number of parameters, again assuming the note at the top of this script
AIC<-(-2*lhoods)+(2*npars) # calculate the AIC score from number of parameters and the LnL
fit_table<-cbind(lhoods, npars, AIC)  # Combine this information together into a single object
rownames(fit_table)<-model
fit_table<-fit_table[order(fit_table[,"AIC"]),]  # sort by AIC
deltaAIC<-fit_table[,"AIC"]-fit_table[1,"AIC"]  ## calculate deltaAIC
rel_lhood<-exp(-0.5*deltaAIC)  # calculate the relative likelihood - step in getting to AIC
sum_rel_lhood<-sum(rel_lhood)  # sum of relative likelihoods, another step in getting to AICw
AICw<-rel_lhood/sum_rel_lhood  #calculate AICw
models_summ_fit<-cbind(round(fit_table, 1), round(deltaAIC,1), round(AICw, 3)) # combine everything relevant together
colnames(models_summ_fit)<-c("LnL", "N_pars", "AIC", "deltaAIC", "AICw") # rename the columns
##Write out a csv with the model fit data
write.csv(models_summ_fit, file="FSC_model_Fits.csv")

# In Reps/best_L_allMods see .csv file generated and look at best model (lowest deltaAIC)

################################################################################################################################
######################################### on computer run: FSC2.R to get data tables ###########################################
################################################################################################################################
## Obtains FSC estimates ##
FSC2.R 

# generates output: noreponly_v2_pars_conv.csv

## Get the boot inputs from "Set up bootstrapping" from FSC2.R section which gives you a folder "boot_input" containing files
## from your bestfit model:
# 		.tpl 
#		.pv 
# 		.est
#  		_maxL.par

## Copy "boot_input" folder into working directory in cluster to run FSC boots for next step
################################################################################################################################

#######################################
## Bootstrapping from best fit model ## 
#######################################

# Step 8: Parametric bootstrapping for confidence intervals 
# after getting boot_input from R script copy the folder into hpc 
cd boot_input 

## Running this:
fsc26 -i {bestfitmodel}_maxL.par -n100 -j -m -s0 -x –I -q -u

## submit as a slurm
sbatch Parametric_boot_for_CI.slurm

##### Use script Prep_FSC_reps_of_bootreps.sh after running Parametric_boot_for_CI.slurm is done 
### 	This script generates multiple replicates of fastsimcoal models within bootstrap replicates
###     for each bootstrap rep, we want to make make 50 replicate FSC runs

cd ../boot_input/
cp  *.pv *.est *.tpl ./{bestfitmodel}_maxL

# change all files to incorporate "_maxL" into each .est .tpl .pv --> name_maxL.est
cd .../boot_input/{bestfitmodel}_maxL

mv {bestfitmodel}.est {bestfitmodel}_maxL.est
mv {bestfitmodel}.tpl {bestfitmodel}_maxL.tpl
mv {bestfitmodel}.pv {bestfitmodel}_maxL.pv

# Run this from inside the directory in which each directory is a bootstrap rep - should also contain the est, tpl, and pv files 
sh Prep_FSC_reps_of_bootreps.sh

## FSC_croc_boot.slurm will then run parameter estimation on each bootstrap replicate. - array = 50 runs * 100 boot = 5000

sbatch ~/scripts/FSC_croc_boot.slurm 	# changed to shared partition 

# run script 

sh Get_best_FSCacross_boots.sh # makes folder best_L_allMods in /Reps/

################################################################################################################################
############################# on computer run: Get_pars_across_bootreps_crocs.R to get bootstrapped plots ######################
################################################################################################################################
























