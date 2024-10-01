##########################################################################################
################################### ADMIXTURE ############################################
##########################################################################################
## Load programs ##
module load bio/BCFtools/1.15.1-GCC-11.3.0
module load lang/Anaconda3
module load lang/R

## Conda: Plink 2 ##
$ source activate plink2
	$ conda install bioconda::plink2						# plink2 v. 2.00a5	
	$ conda install bioconda/label/cf201901::plink 			# plink v. 1.90b4

## ADMIXTURE ##
$ source activate admixture 
	$ conda install bioconda::admixture						# v. 1.3.0
	
cd <admixturedirectory>

## Output of filtered vcf file 
VCF=noreponly_v2.75.renamed.LDpruned.vcf.gz
echo $VCF

VCFFILE=$(basename $VCF)
BASENAME=$(echo ${VCFFILE%%.*})
FILE=${VCFFILE%.vcf.gz} 
echo "VCF = $VCFFILE File = $FILE Basename = $BASENAME"

OUT_file=$(echo ${VCFFILE%.vcf.gz})


OUT=/youroutdirectory/$FILE
OUT=/youroutputdirectory/plink_output/$FILE

# Generate the input file in plink format
source activate plink2

plink --vcf $VCF --make-bed --out $OUT --allow-extra-chr --double-id

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
FILE=noreponly_v2.75.renamed.LDpruned
awk '{$1="0";print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

source activate admixture 

## Run ADMIXTURE. We will run it with cross-validation (the default is 5-fold CV, for higher, choose e.g. cv=10) and K=2
admixture --cv=10 ./plink_output/$FILE.bed 2 > ./admixture_output/log2.out
	# ADMIXTURE produced 2 files: 
	#	.Q which contains cluster assignments for each individual
	#	.P which contains for each SNP the population allele frequencies.
	# --cv=10 cross validation at 10-fold

## Run it in a for loop with K=2 to K=10 and direct the output into log files
for i in {1..8}
do
admixture --cv=10 ./plink_output/$FILE.bed $i > ./admixture_output/logfiles/log${i}.out
done

## To identify the best value of k clusters which is the value with lowest cross-validation error, we need to collect the cv errors. 
## to extract the number of K and the CV error for each corresponding K
# Method 1:
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > $FILE.cv.error
	# 2 0.22044
	# 3 0.21377
	# 4 0.21105
	# 5 0.20811

# Method 2:
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > $FILE.cv.error
	# 2 0.22044
	# 3 0.21377
	# 4 0.21105
	# 5 0.20811
	# 6 0.20743
	# 7 0.21003
	# 8 0.21086
	# 9 0.20791
	# 10 0.20868

# Method 3:
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > $FILE.cv.error

# Method 4: for plotting CV errors
grep -h CV log*.out>cross_validation.txt

## Extract the right order for individual id from your vcf file, using the tfam file information. To do so, you need to use VCFTOOLS with the following command:
cut -f 1 nameofyourfile.tfam > id_admixture.txt
##########################################################################################
## Find optimal number of clusters - plot
R

library(stringr)
library(ggplot2)
library(dplyr)

# Download the cross-validation results you have previously created via bash command.
cv <- read.table("cross_validation.txt")
# cv <- read.table("noreponly_v2.75.renamed.LDpruned.cv.error")

# Analyze the cross-validation results Then, add a K-cluster column indicating the number of K you test and select only two columns of interest, CV and K.
cv$K <-gsub("[\\(\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3)))
names(cv) <- c("cv", "error", "Kvalues", "CVerrors", "K")
CV <- select(cv, CVerrors,K)

# Rename your two columns CV and K-cluster
colnames(CV) <- c("CV","K")

## Change order of columns so that K=2 is first and K=10 is last
# Step 1: Extract numeric part and convert to numeric
CV$K <- as.numeric(gsub("K=", "", CV$K))

# Step 2: Order the dataframe by the numeric part of K
CV <- CV[order(CV$K), ]

# Do a graph showing the cross validation results. Then select the optimal number of clusters regarding:
# the lowest cross validation error or when the cross-validation error decrease the most
graph_title="Cross-Validation plot"
x_title="K"
y_title="Cross-validation error"
graph_1<-ggplot(CV,aes(x=K,y=CV))
graph_1+geom_line()+scale_x_continuous(breaks=c(2,3,4,5,6,7,8,9,10))+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))
        
ggsave("Admixture_cross-validation.pdf",width=7,height=5,dpi=600)
dev.off()
##########################################################################################

## To make plotting easier, we can make a file with the individual names in one column and the species names in the second column. 
## As the species name is in the individual name, it is very easy to extract the species name from the individual name:

awk '{split($1,name,"."); print $1,name[3]}' ./plink_output/$FILE.nosex > $FILE.MU.list # get monitoring unit from name
awk '{split($1,name,"."); print $1,name[1]}' ./plink_output/$FILE.nosex > $FILE.MSp.list # get morph species from name


bcftools query -l ./data/$FILE.vcf.gz > $FILE.samples.txt # save sample ID as txt file
	# manually added morph_sp and saved as $FILE.list
	
## Plot the results in R using script plotADMIXTURE.r
wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r

module load lang/R

# For Morph_Sp
Rscript plotADMIXTURE.r -p $FILE -i $FILE.MSp.list -k 6 -l CA,MCA,HY,CM
	# requires four arguments:
	# (1) the prefix for the ADMIXTURE output files (-p )
	# (2) the file with the species information (-i )
	# (3) the maximum number of K to be plotted (-k 5)
	# (4) a list with the populations or species separated by commas (-l <pop1,pop2...>) 
	# 		The list of populations provided with -l gives the order in which the populations or species shall be plotted

# For monitoring unit: 
Rscript plotADMIXTURE.r -p $FILE -i $FILE.montunit.list -k 6 -l CC,NC,BC,BRW,NSLW,CF,RHW,NRW,CB,NTW,STW

## To evaluate if the ADMIXTURE plot is a good fit, you can use evaladmix and we highly recommend to also use other methods to help infer the demographic history 
## and evidence of hybridisation such as Dstatistics, demographic modeling etc.
https://github.com/GenisGE/evalAdmix
http://www.popgen.dk/software/index.php/EvalAdmix

# Run evalAdmix for all admixture K values
./evalAdmix -plink inputPlinkPrefix -fname inputPlinkPrefix.K.P -qname inputPlinkPrefix.K.Q -P 10 -o output.corres.txt

# For noreponly_v2
# http://127.0.0.1:30905/graphics/plot_zoom_png?width=1728&height=1376
for i in {2..8}
do
./evalAdmix -plink ./plink/noreponly_v2.75.renamed.LDpruned -fname ./admixture/noreponly_v2.75.renamed.LDpruned.$i.P -qname ./admixture/noreponly_v2.75.renamed.LDpruned.$i.Q -P 10 -o noreponly_v2.75.renamed.LDpruned.Q$i.corres.txt
done

# Plot results in R for visualization 
source("visFuns.R")

pop <- as.vector(read.table("inputPlinkPrefix.fam")$V1) # N length character vector with each individual population assignment
q <- as.matrix(read.table("inputPlinkPrefix.K.Q")) # admixture porpotions q is optional for visualization but if used for ordering plot might look better
r <- as.matrix(read.table("output.corres.txt"))

ord <- orderInds(pop=pop, q=q) # ord is optional but this make it easy that admixture and correlation of residuals plots will have individuals in same order

plotAdmix(q=q, pop=pop, ord=ord)
plotCorRes(cor_mat = r, pop = pop, ord=ord, title = "Admixture evaluation as correlation of residuals", max_z=0.25, min_z=-0.25)

