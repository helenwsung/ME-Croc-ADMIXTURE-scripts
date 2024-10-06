######################################################################################################################################################
########################################### Estimating pi, Tajima’s D and Fst along a genome #########################################################
######################################################################################################################################################
# https://www.york.ac.uk/res/dasmahapatra/teaching/MBiol_sequence_analysis/workshop4_2019.html#plotting_fst
cd /datadir/diversity_stats

module load bio/BCFtools/1.15.1-GCC-11.3.0
source activate vcftools
modules

VCF=noreponly_v2.filtered.85.renamed.vcf.gz
OUT=/datadir/populations


## Nucleotide Diversity 
vcftools --gzvcf ./data/$VCF --window-pi 10000 --out ./out/noreponly_v2.filtered.85.renamed_10kb

## Tajima's D
vcftools --gzvcf ./data/$VCF --TajimaD 10000 --out ./out/noreponly_v2.filtered.85.renamed_10kb

## Estimating population divergence with Fst 

# Acutus samples determined by ADMIXTURE K=3
CAsamples=CA90.txt		# 90% cayes acutus
MCAsamples=MCA90.txt	# 90% mainland acutus 
acutus=acutus90.txt		# 90% combined acutus

# Morelets samples determined by ADMIXTURE K=3
CMsamples=morelets90.txt		# 90% moreletii

# Admixed samples determined by ADMIXTURE K=3
HYsamples=admixed.txt		# all hybrids
F1samples=F1hybrids.txt		# 40-60% hybrids

# CA90-MCA90
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$CAsamples --weir-fst-pop ./pops/admixture_k3/$MCAsamples  --out ./out/fst_CA90-MCA90

# CA90-Morelets90
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$CAsamples --weir-fst-pop ./pops/admixture_k3/$CMsamples  --fst-window-size 10000 --out ./out/fst_CA90-CM90

# Acutus90-Morelets90
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$acutus --weir-fst-pop ./pops/admixture_k3/$CMsamples  --fst-window-size 10000 --out ./out/fst_acutus90-CM90

# MCA90-Morelets90
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$MCAsamples --weir-fst-pop ./pops/admixture_k3/$CMsamples  --fst-window-size 10000 --out ./out/fst_MCA90-CM90


## Admixed Samples
# CA90-HY
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$CAsamples --weir-fst-pop ./pops/admixture_k3/$HYsamples  --fst-window-size 10000 --out ./out/fst_CA90-HY

# MCA90-HY
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$MCAsamples --weir-fst-pop ./pops/admixture_k3/$HYsamples  --fst-window-size 10000 --out ./out/fst_MCA90-HY

# Morelets90-HY
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$HYsamples --weir-fst-pop ./pops/admixture_k3/$CMsamples  --fst-window-size 10000 --out ./out/fst_CM90-HY

# Acutus90-HY
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$acutus --weir-fst-pop ./pops/admixture_k3/$HYsamples  --fst-window-size 10000 --out ./out/fst_acutus90-HY

## Clusters 
cluster1=CMcluster1_noHY.txt
cluster2=MCAcluster2_noHY.txt
cluster3=CAcluster3_noHY.txt

# CM-MCA
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$cluster1 --weir-fst-pop ./pops/admixture_k3/$cluster2  --fst-window-size 10000 --out ./out/fst_cluster_CM-MCA

# CM-CA
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$cluster1 --weir-fst-pop ./pops/admixture_k3/$cluster3 --fst-window-size 10000 --out ./out/fst_cluster_CM-CA

# MCA-CA
vcftools --gzvcf ./data/$VCF --weir-fst-pop ./pops/admixture_k3/$cluster2 --weir-fst-pop ./pops/admixture_k3/$cluster3 --fst-window-size 10000 --out ./out/fst_cluster_MCA-CA


## OPEN R - Plotting Results - PopGenStats.R
module load lang/R
R