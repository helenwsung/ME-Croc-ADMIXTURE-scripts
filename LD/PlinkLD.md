##########################################################################################
#################### CALCULATING LD DECAY with PLINK #####################################
##########################################################################################
module load bio/BCFtools/1.15.1-GCC-11.3.0
module load lang/Anaconda3
source activate plink2

##########################################################################################

# move to your data directory
cd /data

##########################################################################################
## Create Plink format files from filtered VCF files for computing LD ##
plink --vcf $VCF --double-id --allow-extra-chr --chr-set 69 --make-pgen --sort-vars --out $OUT
	# --vcf - specified the location of our VCF file.
	# --double-id - allows for "_" in sample name e.g. HS_PC1_A10; told plink to duplicate the id of our samples (this is because plink typically expects a family and individual id - i.e. for pedigree data - this is not necessary for us.
	# --allow-extra-chr - allow additional chromosomes beyond the human chromosome set. This is necessary as otherwise plink expects chromosomes 1-22 and the human X chromosome
	# --chr-set - required if you have something other than human. You provide the number of “chromosomes”, scaffolds in our case = 69
	# --make-pgen - creates all our plink files 
	# --out - output file name header
plink --pfile $OUT --double-id --allow-extra-chr --chr-set 69 --max-alleles 2 --make-bed --out $OUT
	# --pfile - input ped files from previous
	# --max-alleles - set for only biallelic snps
	# --make-bed - output bed files for plink 1.9 to calculate LD
	 
## Get R-squared values for calculating LD
plink --bfile $OUT --double-id --allow-extra-chr --chr-set 69 --r2 gz --ld-window 100 --ld-window-kb 2000 --ld-window-r2 0 --out $(echo $OUT)_r2
	# --bfile points to binary format (PLINK2) variant files (via. VCF, generated above)
	# --r2 gz indicates that you want to run variant association and save as .gz file
	# --ld-window - this allows us to set the size of the lower end of the LD window. In this case, we set it to 100 bp - i.e. any sites with < 100 sites between them are ignored.
	# --ld-window-kb - window size for distance between snps you accept; this is the upper end of the LD window. Here we set it to 2000, meaning that we ignore any two sites more than 2 Mb apart in the genome.
	# --ld-window-r2 - is the threshold for reporting r2.; this sets a filter on the final output but we want all values of LD to be written out, so we set it to 0
	# --ld-window-cm - It is the maximum distance for association testing, by defaullt only 1000bp. I really should test up to 100,000bp apart (this will be huge data file I believe

## Get D-prime values for calculating LD 
plink --bfile plink2_out --double-id --allow-extra-chr --chr-set 69 --r2 dprime --out $(echo $OUT)_r2_dprime
	# --r2 dprime indicates that you want to run variant association for dprime; also provides r2

## Get Rsquared values as a matrix to create heatmaps ##
# THIS FILE IS HUGE!!!! DON'T OPEN ON MY COMPUTER 
plink --bfile plink2_out --double-id --allow-extra-chr --chr-set 69 --r2 square --out $(echo $OUT)_r2_square
gzip $(echo $OUT)_r2_square

######################################################################################################
## Loop through all files in folder ##
# Directory where the .vcf.gz files are located

# Change this to your directory if it's different - includes vcf files for each pop group created by LDpruning_Dataconversion.R script
VCF_DIR=/data/filtered.75 
							
# Directory to store the output files
OUT_DIR=/data/filtered.75/plinkfiles_10Mb

# Create the output directory if it doesn't exist
mkdir -p "$OUT_DIR"

## For Plink1 
# Loop through each .vcf.gz file in the directory - created .vcf.gz files per population group from ADMIXTURE
for VCF in "$VCF_DIR"/*.vcf.gz; do

    # Extract the base name of the file for output
    BASE_NAME=$(basename "$VCF" .vcf.gz)
    OUT="$OUT_DIR/$BASE_NAME"
    
    # Create Plink format files
    plink --vcf "$VCF" --double-id --allow-extra-chr --chr-set 69 --out "$OUT"
   
    # Calculate R-squared values
    plink --bfile "$OUT" --double-id --allow-extra-chr --chr-set 69 --r2 gz --ld-window 10000 --ld-window-kb 10000 --ld-window-r2 0 --out "${OUT}_r2"

    # Calculate D-prime values
    plink --bfile "$OUT" --double-id --allow-extra-chr --chr-set 69 --r2 dprime --out "${OUT}_r2_dprime"
done


################################################################################
## visualize LD decay results from plink in R - LDscript ##
################################################################################
