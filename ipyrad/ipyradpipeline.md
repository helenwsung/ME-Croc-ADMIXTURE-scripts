##########################################################################################
############################### Ipyrad SNP Calling Pipeline ##############################
##########################################################################################
## Load programs ##
module load bio/BCFtools/1.15.1-GCC-11.3.0
module load lang/Anaconda3
module load lang/R

## Quality Control: multiqc ## 
$ source activate multiqc
	# FastQC (also in trimmomatic_hws) 
		$ conda install bioconda::fastqc					 	# fastqc v0.12.1
	# MultiQC 
		$ conda install bioconda::multiqc 						# multiqc v1.22.2
		
## Trimming: trimmomatic ##
$ source activate trimmomatic
	# Trimmomatic 
		$ conda install -c bioconda/label/main trimmomatic		# trimmomatic v0.39
		
## ipyrad ##
$ source activate ipyrad
	$ conda install bioconda::ipyrad							# v.0.9.96

##########################################################################################
##################################### QC Reads ###########################################
##########################################################################################
source ~/.bash_profile

# Load Anaconda
module load lang/Anaconda3

# Activate environment
source activate multiqc

mkdir rawFastQC

# Directory paths
data_dir="./rawfastq"
output_dir="./multiqc/"

# Move to the data directory
cd $data_dir

# Get the list of fastq.gz files
files=(*.fq.gz)

# Get the current file to process from the job array index
current_file=${files[$SLURM_ARRAY_TASK_ID - 1]}

# Run FastQC on the current file
fastqc "$current_file" -o $output_dir

cd $output_dir

# Generate MultiQC report
multiqc *_fastqc.zip

##########################################################################################
##################################### Trimmomatic ########################################
##########################################################################################

# Activate environment
source activate trimmomatic

data_dir="./rawfastq"
out_dir="./trimmmed_fastq"

cd $data_dir

# Make Trimmomatic command 
for R1 in *R1*
do
R2=${R1//R1_.fq.gz/R2_.fq.gz}
R1paired=${R1//.fq.gz/paired.fq.gz}
R1unpaired=${R1//.fq.gz/unpaired.fq.gz}
R2paired=${R2//.fq.gz/paired.fq.gz}
R2unpaired=${R2//.fq.gz/unpaired.fq.gz}	
echo "trimmomatic PE -threads 6 -phred33 $R1 $R2 $R1paired $R1unpaired $R2paired $R2unpaired \
ILLUMINACLIP:/home/hwsung/.conda/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10:8:keepBothReads \
HEADCROP:10 LEADING:7 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:50" >> trimmomatic.cmds
done
cat trimmomatic.cmds

## Run trimmomatic.cmds in Sbatch 
## RUN QC again on trimmed reads to double check 

##########################################################################################
################################## IPYRAD PIPELINE  ######################################
##########################################################################################

module load lang/Anaconda3
source activate ipyrad # activates created conda environment

## Get working data from trimmed files
fq_dir="./trimmmed_fastq"
data_dir="./ipyrad"

cd $data_dir
cp -r $fq_dir ./	# copying trimmed files into working directory for ipyrad

## changing file names for loop (changed from *paired.fq.gz to *fq.gz)
for i in *paired.fq.gz
do
  mv -- "$i" "${i/%paired.fq.gz/.fq.gz}"
done

# Create list of fastq files used in analysis
for fastq in *R1*.fq.gz
do
echo "${fastq%*_R1*}" >> samplenamelist.txt
done


# creates param file - params-noreponly_v2.txt 
ipyrad -n noreponly_v2 # creates param file 

## ipyrad steps 1-2 ## 
ipyrad -p params-noreponly_v2.txt  -s 12


# -r fetches informative results from currently executed steps
ipyrad -p params-noreponly_v2.txt -r

## ipyrad steps 3-7 on each new branch param file ##
ipyrad -p params-noreponly_v2.txt -s 34567 -c 8 		# repeated samples taken out

##########################################################################################
## EXAMPLE SLURM SCRIPT ## 

#!/bin/bash
#SBATCH --job-name=ipyrad_3RAD_trimmed_complete
#SBATCH --partition=kill-shared

## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=03-00:00:00 ## time format is DD-HH:MM:SS

#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
##SBATCH --ntasks=4
#SBATCH --mem=200G
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node

#SBATCH --error=job%A.err ## %A - filled with jobid
#SBATCH --output=job%A.out ## %A - filled with jobid


module load lang/Anaconda3
source activate ipyrad

cd ~/ipyrad/

## ipyrad steps 3-7 on each new branch param file
ipyrad -p params-noreponly_v2.txt -s 34567 -c 8
































