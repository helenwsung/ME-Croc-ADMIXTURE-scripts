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
git clone https://github.com/isaacovercast/easySFS.git

cd /yourdir/analyses/fastsimcoal/easySFS

## make executable 
chmod +x ./easySFS.py
ls -ltrh 

VCF=noreponly_v2.filtered.75.renamed.vcf.gz
POPMAP=popmap.txt

## noreponly
# or easySFS, to identify the best projection value - Run Preview to see if works - 
# Each column is the number of samples in the projection and the number of segregating sites at that projection value
# The dadi manual recommends maximizing the number of segregating sites, but at the same time if you have lots of missing data then you might have to balance # of segregating sites against # of samples to avoid downsampling too far
# this line will calculate the necessary SFS files - he SFS cannot be computed with sites that contain missing data. By projecting down the number of individuals per population, more sites can be recovered
./easySFS.py -i ../data/noreponly_v2/$VCF -p ../data/noreponly_v2/$POPMAP --preview -a
			# pop0 = CM
			# pop1 = CA
			# pop2 = MCA

# picked the largest projection segregating sites
./easySFS.py -i ../data/noreponly_v2/$VCF -p ../data/noreponly_v2/$POPMAP --proj=130,84,24 -a -o ./noreponly_v2.filtered.75_16.4y_SFS/ -f

################################################################################################################################
## NOTE:
500000 generations (4 million years divergence/8 years per generation) 
	# More than 5% of female population os C. acutus were predicted to be reproductie by age of 8yo (Briggs-Gonzalez et al. 2017)

## NOTE UPDATED: CHANGE TO MAKE YEARS PER GENERATION = C. acutus --> average generation time = 16.2 years (Briggs-Gonzalez et al 2017)
246,914 generations = (4 my/16.2 yo)

# also try 18yo b/c more than 95% of female population were predicted to be reproductive by 18yo
222,222 generations = (4 my/18 yo)
## live at least 25 years - lower bound
160,000 generations = (4 my/25 yo)
################################################################################################################################

## Run script from inside Models folder
cd /datadir/fastsimcoal/Models

# copy easysfs file into Models folder
cp /datadir/fastsimcoal/easySFS/noreponly_v2.filtered.75_16.2y_SFS/fastsimcoal2/noreponly_v2_MSFS.obs ./

sh fsc-selectbestrun.sh
sh Prep_FSC_reps.sh MSFS 

## edit .slurm script to array --> 50 runs * 16 models = 800; 50*25 = 1250
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
# MallAsymK3_CACM.bestlhoods was best so add model to everything
	## New FSC_model_Fits.csv --> best model = curM_RecResize_CACM.bestlhoods

# on computer run: FSC2.R to get data tables 

#######################################################
## Bootstrapping ## 
########################################################

# Step 8: Parametric bootstrapping for confidence intervals 
# after getting boot_input from R script copy the folder into hpc 
cd boot_input 
fsc26 -i MallAsym_RecResize_CACM_maxL.par -n100 -j -m -s0 -x –I -q -u


### 	This script generates multiple replicates of fastsimcoal models within bootstrap replicates
###     for each bootstrap rep, we want to make make 50 replicate FSC runs
cd /datadir/fastsimcoal/noreponly_v2/boot_input/MallAsym_RecResize_CACM_maxL

# Run this from inside the directory in which each directory is a bootstrap rep - should also contain the est, tpl, and pv files 
# change all files to incorporate "_maxL" into each .est .tpl .pv --> name_maxL.est
sh Prep_FSC_reps_of_bootreps.sh

## FSC_croc_boot.slurm will then run parameter estimation on each bootstrap replicate. - array = 50 runs * 100 boot = 5000

sbatch ~/scripts/FSC_croc_boot.slurm 	# changed to shared partition 

# run script 

sh Get_best_FSCacross_boots.sh # makes folder best_L_allMods in /Reps/

##########################################################################################
## METHOD 2? ##
##########################################################################################
## Bootstrapping with https://speciationgenomics.github.io/fastsimcoal2/ tutorial ##
PREFIX="MallAsym_RecResize_CACM"
cd /datadir/fastsimcoal/noreponly_v2/Reps_16.2y/MallAsym_RecResize_CACM/bestrun
VCF=/datadir/fastsimcoal/data/noreponly_v2/noreponly_v2.filtered.75.renamed.vcf.gz

# Get all lines with genomic data
zgrep -v "^#" $VCF > $PREFIX.allSites

# Get the header
zgrep "^#" $VCF > header

# get 100 files with 4338 sites each (number 101 removed due to only 90 sites)
split -l 36297 $PREFIX.allSites $PREFIX.sites.

# Generate 50 files each with randomly concatenated blocks and compute the SFS for each:
for i in {1..50}
do
  # Make a new folder for each bootstrapping iteration:
  mkdir bs$i
  cd bs$i

  # Add the header to our new bootstrapped vcf file
  cat ../header > $PREFIX.bs.$i.vcf
  # Randomly add 100 blocks
  for r in {1..100}
  do
    cat `shuf -n1 -e ../$PREFIX.sites.*` >> ${PREFIX}.bs.$i.vcf
  done
  # Compress the vcf file again
  gzip ${PREFIX}.bs.$i.vcf

  # Make an SFS from the new bootstrapped file
  /datadir/fastsimcoal/easySFS/easySFS.py -i ${PREFIX}.bs.$i.vcf.gz -p /datadir/fastsimcoal/data/noreponly_v2/noreponly_v2.75.renamed.LDpruned.popmap90.txt -a -f --proj=130,84,24

  # Copy the observed SFS file into this folder renaming it to match the .tpl prefix
  cp ../${PREFIX}_MSFS.obs  ${PREFIX}.bs.${i}_MSFS.obs

  # Say that it is finished with iteration $i
  echo bs$i" ready"

  cd ..
done

# Copy .est .par .tpl files into each bootstrap folder 
for i in {1..50}
do
cp ./*.est ./bs$i/${PREFIX}.bs.${i}.est
cp ./*.tpl ./bs$i/${PREFIX}.bs.${i}.tpl
done



## Now we would run the parameter estimation under the best model 100 times with each of these boostrapped SFS. This would take very long.
for bs in {1..50}
do
  cd bs$bs
  # Run fastsimcoal 100 times:
  for i in {1..100}
  do
    mkdir run$i
    cd run$i
    cp ${PREFIX}.bs.$bs.* ./
    fastsimcoal2 -t ${PREFIX}.bs.$bs.tpl -e ${PREFIX}.bs.$bs.est -m -0 -C 10 -n 10000 -L 40 -s0 -M -q
    cd ..
  done
  # Find the best run:
  fsc-selectbestrun.sh

  cd ..
done