# Download two ruby scripts for ancestry painting from github 
if [ ! -f get_fixed_site_gts.rb ]
then
    wget https://raw.githubusercontent.com/mmatschiner/tutorials/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rb
fi
if [ ! -f plot_fixed_site_gts.rb ]
then
    wget https://raw.githubusercontent.com/mmatschiner/tutorials/master/analysis_of_introgression_with_snp_data/src/plot_fixed_site_gts.rb
fi


######################################################################################################################################################
################################################################# ANCESTRY PAINTING ##################################################################
######################################################################################################################################################
# Load the ruby and bcftools modules.
module load lang/Anaconda3
module load bio/BCFtools/1.15.1-GCC-11.3.0
module load lang/Ruby/3.2.2-GCCcore-12.2.0
module load lang/Python/3.9.6-GCCcore-11.2.0 

cd /datadirectory/

VCF=/datadirectory/data/noreponly_v2.filtered.85.renamed.vcf.gz
VCFFILE=$(basename $VCF)
BASENAME=noreponly_v2.filtered.85.renamed

gunzip $VCF > noreponly_v2.filtered.85.renamed.vcf

OUT=/datadirectory/out/$BASENAME
echo $OUT 

# Count the number of variables (words) in the string
word_count=$(echo "$CA_99_string" | wc -w)

# Print the word count
echo "Number of variables in the string: $word_count"

######################################################################################################################################################
## Running script: ##
## Example test script on subset individuals ##
  ruby get_fixed_site_gts.rb $VCF $OUTPUT.txt $PUTATIVEPARENTSPECIESSTRING $HYBRIDSTRING $PARENT2STRING 0.8 
        # the name of the uncompressed VCF input file, [NC_031969.f5.sub1.vcf,
        # the name of an output file, which will be a tab-delimited table,
        # a string of comma-separated IDs of samples for the first putative parent species,
        # a string of comma-separated IDs of samples for the putative hybrid species,
        # another string of comma-separated IDs of samples for the second putative parent species,
        # a threshold value for the required completeness of parental genotype information so that sites with too much missing data are discarded.

## SCRIPT 2:
  ruby plot_fixed_site_gts.rb pops1.fixed.txt pops1.fixed.svg 1.0 1000
    # the name of the file written by script get_fixed_site_gts.rb,
        # the name of an output file which will be a plot in SVG format,
        # a threshold value for the required completeness, which now applies not only to the parental species but also to the putative hybrid species,
        # the minimum chromosomal distance in bp between SNPs included in the plot. This last argument aims to avoid that the ancestry painting is overly dominated by high-divergence regions.
######################################################################################################################################################
##########################################################################################
## Noreponly_v2 ##

cd /home/hwsung/kantar_koastore/helen/analyses/datafiles/3RADmerged_noreponly/filteredVCF/ancestrypaint

# Set the name of the vcf file.
# noreponly_v2
VCF=/datadirectory/noreponly_v2.filtered.85.renamed.vcf
FILE=noreponly_v2.filtered.85
OUT=/datadirectory/out/$FILE

## Example test script on subset individuals ##
  ruby get_fixed_site_gts.rb $VCF ./data/$FILE.test1.txt ${CA_99_string} ${F1_Hybrid_string} ${CM_99_string} 0.8
    # CA_99_string: "MCA.P4G2.CB","CA.P3G10.CC","CA.P3F10.CC","CA.P2E8.CC","CA.P2E3.CC","CA.P2A3.CC","CA.P1H12.NC","CA.P2C7.CC","CA.P3C9.NC","CA.PC11.NC","CA.P1A12.NC","CA.P1C12.NC","CA.PC10.NC","CA.P1F12.NC","CA.P1G11.NC","CA.P3D9.NC","CA.P2C3.CC","CA.P2B3.CC","CA.P2G3.CC","CA.P2B8.CC","CA.P1E11.NC","CA.P2G5.CC","CA.P2A5.CC","CA.P2B4.CC","CA.P2B7.CC","CA.P2C4.CC","CA.P2D2.NC","CA.P2E4.CC","CA.P2E5.CC","CA.P2E7.CC","CA.P2F4.CC","CA.P2F7.CC","CA.P2G4.CC","CA.P2G8.CC","CA.P2H5.NA","CA.P2H8.CC","CA.PC12.NC","CA.PC14.CC","CA.PC15.CC","CA.PC16.CC","CA.PC9.NC","MCA.P1B11.NRW"
    # F1_Hybrid_string: "HY.P2H6.CB","CM.P4D1.NSLW","HY.P2B11.CB","HY.P1B4.NA","MCA.P1E7.CB","MCA.P1B7.CB","HY.P1B5.CB","MCA.P1A8.CB","MCA.P1D7.CB","MCA.P1A7.CB","MCA.P1C7.CB","MCA.P1C6.NTW"
    # CM_99_string: "CM.PC6.NRW", "CM.PC8.NRW", "CM.P1D3.NRW", "CM.P3A3.NRW", "CM.P3B3.NRW", "CM.P3C1.NRW", "CM.P3C3.NRW", "CM.P3D3.NRW", "CM.P3G3.NRW", "CM.PC7.NRW", "HY.P3A2.NRW", "HY.P3B5.NRW", "HY.P3C5.NRW", "HY.P3D5.NRW", "HY.P3H4.NRW", "CM.P1F10.NRW", "CM.P3F3.NRW", "CM.P1G10.NRW", "HY.P3B4.NRW", "CM.P1A10.NRW", "CM.P3C8.CB", "HY.P3A5.NRW", "CM.P3E3.NRW", "CM.P3B1.NRW", "HY.P3F4.BRW", "CM.P3D6.RHW", "CM.P3G4.NRW", "HY.P3B2.NRW", "CM.P3A4.NRW", "CM.P2D10.NRW", "HY.P3E4.BRW", "CM.P1D10.NRW", "CM.P2C6.NRW", "CM.P2H12.NRW", "CM.P2F10.NRW", "HY.P2D11.RHW", "CM.P1E10.NRW", "CM.P2A12.NRW", "CM.P1D9.NRW", "CM.P2C10.NRW", "CM.P3C6.NRW", "CM.P3H3.NRW", "CM.P2C12.NRW", "CM.P2B12.NRW", "HY.P3E8.NA", "CM.P3A1.NRW", "CM.P2D12.NRW", "CM.P3F6.NTW", "CM.P2E12.NRW", "CM.P3A7.CB", "CM.P2G12.NRW", "HY.P1G3.NRW", "CM.P1H10.NRW"
    # Found 27264 sites with sufficient completeness (at least 0.8).
    # Found 804 sites with complete fixation of different alleles in the parental populations.
  ruby get_fixed_site_gts.rb $VCF ./data/$FILE.test2.txt ${CA_99_string} ${F1_Hybrid_string} ${CM_99_string} 0.9
    # Found 19374 sites with sufficient completeness (at least 0.9).
    # Found 613 sites with complete fixation of different alleles in the parental populations.
  
# Get ancestry painting figure 
  ruby plot_fixed_site_gts.rb ./data/$FILE.test1.txt $OUT.test1.svg 0.8 1000 > $OUT.test1.log
    # at minimum distance of 100 bp to each other 
    ## IT WORKED!!!!
  ruby plot_fixed_site_gts.rb ./data/$FILE.test2.txt $OUT.test2.svg 0.9 1000 > $OUT.test2.log

##########################################################################################
##################
## noreponly_v2 ##
##################

# Set the name of the sample table based on admixture K3.
table=/datadirectory/data/ad_cluster.pop.txt

################
## CA strings ## 
################

# Define a list of 99% CA individuals (>99% CA, <.01 MCA, <.02 CM)
table=/datadirectory/data/ad_cluster.pop.txt

## Define a list of 99% CA individuals (>99% CA, <.01 MCA, <.02 CM)
CA_99_string=()
# Use awk to extract the objects without quotes and with no spaces between commas
CA_99_string=$(awk -F'\t' '$2 ~ /CA_99/ && $2 !~ /MCA_99/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
CA_99_string=${CA_99_string%,}
echo ${CA_99_string}
    # CA_99_string="MCA.P4G2.CB,CA.P3G10.CC,CA.P3F10.CC,CA.P2E8.CC,CA.P2E3.CC,CA.P2A3.CC,CA.P1H12.NC,CA.P2C7.CC,CA.P3C9.NC,CA.PC11.NC,CA.P1A12.NC,CA.P1C12.NC,CA.PC10.NC,CA.P1F12.NC,CA.P1G11.NC,CA.P3D9.NC,CA.P2C3.CC,CA.P2B3.CC,CA.P2G3.CC,CA.P2B8.CC,CA.P1E11.NC,CA.P2G5.CC,CA.P2A5.CC,CA.P2B4.CC,CA.P2B7.CC,CA.P2C4.CC,CA.P2D2.NC,CA.P2E4.CC,CA.P2E5.CC,CA.P2E7.CC,CA.P2F4.CC,CA.P2F7.CC,CA.P2G4.CC,CA.P2G8.CC,CA.P2H5.NA,CA.P2H8.CC,CA.PC12.NC,CA.PC14.CC,CA.PC15.CC,CA.PC16.CC,CA.PC9.NC,MCA.P1B11.NRW"
    # removed: CA.P2A4.CC

## Define a list of >90% CA individuals (>90% CA, <1% MCA, <1% CM)
CA_90_string=()
CA_90_string=$(awk -F'\t' '$2 ~ /CA_(90|99)/ && $2 !~ /MCA_(90|99)/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
CA_90_string=${CA_90_string%,}
echo ${CA_90_string}
    # CA_90_string="CA.P3E10.CC,CA.P3B10.NC,MCA.P4G2.CB,CA.P3A11.CC,CA.P3D10.CC,CA.P3G10.CC,CA.P3F10.CC,CA.P2E8.CC,CA.P2F3.CC,CA.P1E12.NC,CA.P2E3.CC,CA.P2A3.CC,CA.P1H12.NC,CA.P1B12.NC,CA.P2C7.CC,CA.P3C9.NC,CA.PC11.NC,CA.P1A12.NC,CA.P1C12.NC,CA.PC10.NC,CA.P3G9.NC,CA.P1F12.NC,CA.P1G11.NC,CA.P3D9.NC,CA.P2C3.CC,CA.P2B3.CC,CA.P2G3.CC,CA.P2D7.CC,CA.P2B8.CC,CA.P1E11.NC,CA.P2G5.CC,CA.P3B9.NC,CA.P2G7.CC,CA.P2F8.NC,CA.P2F2.NC,CA.P2C5.CC,CA.P3H9.NC,CA.P2B2.NC,CA.P2A5.CC,CA.P2B4.CC,CA.P2B7.CC,CA.P2C4.CC,CA.P2D2.NC,CA.P2E4.CC,CA.P2E5.CC,CA.P2E7.CC,CA.P2F4.CC,CA.P2F7.CC,CA.P2G4.CC,CA.P2G8.CC,CA.P2H5.NA,CA.P2H8.CC,CA.PC12.NC,CA.PC14.CC,CA.PC15.CC,CA.PC16.CC,CA.PC9.NC,MCA.P1B11.NRW"
    # removed: CA.P3H10.CC,CA.P4B1.CC,CA.P2A4.CC
 
## Define a list of CA backcrossed (mostly acutus) individuals (<90% CA, <2% MCA, <.1% CM)
CA_backcross_string=()
CA_backcross_string=$(awk -F'\t' '$2 ~ /CA_backcross/ && $2 !~ /MCA_backcross/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
CA_backcross_string=${CA_backcross_string%,}
echo ${CA_backcross_string}
    # CA_backcross_string="MCA.P1D5.BRW,CA.P1D12.NC,MCA.P1F7.BC,CA.P2A2.NC,CA.P2C2.NC,CA.P1D11.NC,CA.P1G12.NC,CA.P3E9.NC"
    # removed: CA.P3C10.NC
    
# Define a list of 99% MCA individuals (>99% MCA, <.01 CA, <.02 CM)
MCA_99_string=()
MCA_99_string=$(awk '$2 ~ /MCA_99/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
MCA_99_string=${MCA_99_string%,}
echo ${MCA_99_string}
    # MCA_99_string="MCA.P1F2.STW,MCA.P4C1.STW,MCA.P3D11.STW,MCA.P1H2.STW,MCA.P1C2.STW,MCA.P1D2.STW,MCA.P1E2.STW,MCA.P1E6.STW,MCA.P1F6.STW,MCA.P1G2.STW,MCA.P2H10.NA"

## Define a list of >90% MCA individuals (>90% MCA, <1% CA, <1% CM)
MCA_90_string=()
MCA_90_string=$(awk '$2 ~ /MCA_9/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
MCA_90_string=${MCA_90_string%,}
echo ${MCA_90_string}
    # MCA_90_string="MCA.P3C11.STW,MCA.P1F2.STW,MCA.P4C1.STW,MCA.P3D11.STW,MCA.P1H2.STW,MCA.P1G8.CB,MCA.P1C2.STW,MCA.P1D2.STW,MCA.P1E2.STW,MCA.P1E6.STW,MCA.P1F6.STW,MCA.P1G2.STW,MCA.P2H10.NA,MCA.P1G6.CB,MCA.P2B9.CB,HY.P1C4.CB"
    # removed: MCA.P1E5.CB,
    
## Define a list of MCA backcrossed (mostly acutus) individuals (<90% MCA, <2% CA, <.1% CM)
MCA_backcross_string=()
MCA_backcross_string=$(awk -F'\t' '$2 ~ /MCA_backcross/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
MCA_backcross_string=${MCA_backcross_string%,}
echo ${MCA_backcross_string}
    # MCA_backcross_string="MCA.P1H7.CB,MCA.P3E11.STW,MCA.P1G5.CB,MCA.P2G10.NA,MCA.P3B11.CB,MCA.P1B2.STW,HY.P1A2.NTW,MCA.P1B8.CB,MCA.P1H8.BRW,HY.P1A5.CB,MCA.P1C5.CB,MCA.P1D6.STW,HY.P1D4.NA,MCA.P1D8.CB,HY.P1A4.CB,MCA.P1B6.CB,HY.P2C9.CB,MCA.P1E8.CB,MCA.P1F8.NTW"
    # removed: MCA.P2A9.CB,MCA.P1H6.CB,MCA.P1H5.CB
    
####################
## Acutus strings ##
####################

# Define a list of >90% C.acutus individuals (>90% CA, >90% MCA, <1% CM)
table=/datadirectory/data/acutuspopmap.txt

acutus_90_string=()
acutus_90_string=$(awk '$2 ~ /Acutus/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
acutus_90_string=${acutus_90_string%,}
echo ${acutus_90_string}
    # acutus_90_string="CA.P3E10.CC,CA.P3B10.NC,MCA.P4G2.CB,CA.P3A11.CC,MCA.P3C11.STW,CA.P3D10.CC,CA.P3G10.CC,CA.P3F10.CC,MCA.P1F2.STW,CA.P2E8.CC,CA.P2F3.CC,CA.P1E12.NC,CA.P2E3.CC,MCA.P4C1.STW,CA.P2A3.CC,CA.P1H12.NC,CA.P1B12.NC,CA.P2C7.CC,CA.P3C9.NC,CA.PC11.NC,CA.P1A12.NC,CA.P1C12.NC,MCA.P3D11.STW,CA.PC10.NC,MCA.P1H2.STW,CA.P3G9.NC,CA.P1F12.NC,MCA.P1G8.CB,CA.P1G11.NC,CA.P3D9.NC,CA.P2C3.CC,CA.P2B3.CC,CA.P2G3.CC,CA.P2D7.CC,CA.P2B8.CC,CA.P1E11.NC,CA.P2G5.CC,CA.P3B9.NC,CA.P2G7.CC,MCA.P1C2.STW,MCA.P1D2.STW,MCA.P1E2.STW,MCA.P1E6.STW,MCA.P1F6.STW,MCA.P1G2.STW,MCA.P2H10.NA,MCA.P1G6.CB,MCA.P2B9.CB,HY.P1C4.CB,CA.P2F8.NC,CA.P2F2.NC,CA.P2C5.CC,CA.P3H9.NC,CA.P2B2.NC,CA.P2A5.CC,CA.P2B4.CC,CA.P2B7.CC,CA.P2C4.CC,CA.P2D2.NC,CA.P2E4.CC,CA.P2E5.CC,CA.P2E7.CC,CA.P2F4.CC,CA.P2F7.CC,CA.P2G4.CC,CA.P2G8.CC,CA.P2H5.NA,CA.P2H8.CC,CA.PC12.NC,CA.PC14.CC,CA.PC15.CC,CA.PC16.CC,CA.PC9.NC,MCA.P1B11.NRW"
    # removed=CA.P4B1.CC,CA.P2A4.CC,CA.P3H10.CC
    
table=/datadirectory/data/ad_cluster.pop.txt
acutus_backcross_string=()
acutus_backcross_string=$(awk '$2 ~ /Acutus_backcross/  { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
acutus_backcross_string=${acutus_backcross_string%,}
echo ${acutus_backcross_string}
    # acutus_backcross_string="HY.P3H11.NSLW,HY.P2E9.NA,HY.P2F9.NA,HY.P3G11.NSLW,CM.P1B9.NSLW,MCA.P1A6.NTW,CM.P1A9.NSLW,MCA.P1G7.CB"

################
## CM strings ##
################

# Define a list of 99% CM individuals (>99% CM)
CM_99_string=()
CM_99_string=$(awk '$2 ~ /CM_99/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
CM_99_string=${CM_99_string%,}
echo ${CM_99_string}
    # CM_99_string="CM.PC6.NRW,CM.PC8.NRW,CM.P1D3.NRW,CM.P3A3.NRW,CM.P3B3.NRW,CM.P3C1.NRW,CM.P3C3.NRW,CM.P3D3.NRW,CM.P3G3.NRW,CM.PC7.NRW,HY.P3A2.NRW,HY.P3B5.NRW,HY.P3C5.NRW,HY.P3D5.NRW,HY.P3H4.NRW,CM.P1F10.NRW,CM.P3F3.NRW,CM.P1G10.NRW,HY.P3B4.NRW,CM.P1A10.NRW,CM.P3C8.CB,HY.P3A5.NRW,CM.P3E3.NRW,CM.P3B1.NRW,HY.P3F4.BRW,CM.P3D6.RHW,CM.P3G4.NRW,HY.P3B2.NRW,CM.P3A4.NRW,CM.P2D10.NRW,HY.P3E4.BRW,CM.P1D10.NRW,CM.P2C6.NRW,CM.P2H12.NRW,CM.P2F10.NRW,HY.P2D11.RHW,CM.P1E10.NRW,CM.P2A12.NRW,CM.P1D9.NRW,CM.P2C10.NRW,CM.P3C6.NRW,CM.P3H3.NRW,CM.P2C12.NRW,CM.P2B12.NRW,HY.P3E8.NA,CM.P3A1.NRW,CM.P2D12.NRW,CM.P3F6.NTW,CM.P2E12.NRW,CM.P3A7.CB,CM.P2G12.NRW,HY.P1G3.NRW"
    # removed: CM.P2H11.NRW
    
## Define a list of >90% CM individuals (>90% CM)
CM_90_string=()
CM_90_string=$(awk -F'\t' '$2 ~ /CM_9/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
CM_90_string=${CM_90_string%,}
echo ${CM_90_string}
    # CM_90_string="CM.PC6.NRW,CM.PC8.NRW,CM.P1D3.NRW,CM.P3A3.NRW,CM.P3B3.NRW,CM.P3C1.NRW,CM.P3C3.NRW,CM.P3D3.NRW,CM.P3G3.NRW,CM.PC7.NRW,HY.P3A2.NRW,HY.P3B5.NRW,HY.P3C5.NRW,HY.P3D5.NRW,HY.P3H4.NRW,CM.P1F10.NRW,CM.P3F3.NRW,CM.P1G10.NRW,HY.P3B4.NRW,CM.P1A10.NRW,CM.P3C8.CB,HY.P3A5.NRW,CM.P3E3.NRW,CM.P3B1.NRW,HY.P3F4.BRW,CM.P3D6.RHW,CM.P3G4.NRW,CM.P1B10.NRW,HY.P3B2.NRW,CM.P3A4.NRW,CM.P2D10.NRW,HY.P3E4.BRW,CM.P1D10.NRW,CM.P2C6.NRW,CM.P2H12.NRW,CM.P2F10.NRW,HY.P2D11.RHW,CM.P1E10.NRW,CM.P2A12.NRW,CM.P1D9.NRW,CM.P2C10.NRW,CM.P3C6.NRW,CM.P3H3.NRW,CM.P2C12.NRW,CM.P2B12.NRW,CM.P3B6.NRW,CM.P1C10.NRW,HY.P3E8.NA,CM.P3A1.NRW,HY.P3E5.CB,CM.P2D12.NRW,CM.P3F6.NTW,CM.P2E12.NRW,CM.P2G6.NRW,CM.P3A7.CB,HY.P3D2.BRW,HY.P3D4.BRW,CM.P3G7.CB,CM.PC5.CF,HY.P3C2.BRW,HY.P3H5.BRW,CM.PC2.CF,CM.P2G9.NRW,CM.P2E10.NRW,HY.P3E2.BRW,CM.P2G12.NRW,HY.P3B8.CB,CM.P2H9.NRW,CM.P1B3.NRW,CM.P2A10.NRW,CM.P1H9.NRW,HY.P3G5.CB,CM.P3G8.CB,HY.P1G4.CB,HY.P1G3.NRW,CM.P1C9.NSLW,HY.P2D9.CB,CM.P2F11.CB,HY.P3F5.CB,CM.P3F8.CB,CM.P4E1.CB,HY.P1E4.NRW,HY.P1F4.CB,CM.P3F7.CB,CM.P3F2.NA,CM.P3A6.CF,HY.P1F5.NRW,CM.P1H10.NRW,HY.P4E2.CB"
    # removed: CM.P2H11.NRW
    
# Define a list of CM backcrossed individuals (90% > CM > 10%)
CM_backcross_string=()
CM_backcross_string=$(awk -F'\t' '$2 ~ /CM_backcross/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
CM_backcross_string=${CM_backcross_string%,}
echo ${CM_backcross_string}
    # CM_backcross_string="CM.PC6.NRW,CM.PC8.NRW,CM.P1D3.NRW,CM.P3A3.NRW,CM.P3B3.NRW,CM.P3C1.NRW,CM.P3C3.NRW,CM.P3D3.NRW,CM.P3G3.NRW,CM.PC7.NRW,HY.P3A2.NRW,HY.P3B5.NRW,HY.P3C5.NRW,HY.P3D5.NRW,HY.P3H4.NRW,CM.P1F10.NRW,CM.P3F3.NRW,CM.P1G10.NRW,HY.P3B4.NRW,CM.P1A10.NRW,CM.P3C8.CB,HY.P3A5.NRW,CM.P3E3.NRW,CM.P3B1.NRW,HY.P3F4.BRW,CM.P3D6.RHW,CM.P3G4.NRW,CM.P1B10.NRW,HY.P3B2.NRW,CM.P3A4.NRW,CM.P2D10.NRW,HY.P3E4.BRW,CM.P1D10.NRW,CM.P2C6.NRW,CM.P2H12.NRW,CM.P2F10.NRW,HY.P2D11.RHW,CM.P1E10.NRW,CM.P2A12.NRW,CM.P1D9.NRW,CM.P2C10.NRW,CM.P3C6.NRW,CM.P3H3.NRW,CM.P2C12.NRW,CM.P2B12.NRW,CM.P3B6.NRW,CM.P1C10.NRW,HY.P3E8.NA,CM.P3A1.NRW,HY.P3E5.CB,CM.P2D12.NRW,CM.P3F6.NTW,CM.P2E12.NRW,CM.P2G6.NRW,CM.P3A7.CB,HY.P3D2.BRW,HY.P3D4.BRW,CM.P3G7.CB,CM.PC5.CF,HY.P3C2.BRW,HY.P3H5.BRW,CM.PC2.CF,CM.P2G9.NRW,CM.P2E10.NRW,HY.P3E2.BRW,CM.P2G12.NRW,HY.P3B8.CB,CM.P2H9.NRW,CM.P1B3.NRW,CM.P2A10.NRW,CM.P1H9.NRW,HY.P3G5.CB,CM.P3G8.CB,HY.P1G4.CB,HY.P1G3.NRW,CM.P1C9.NSLW,HY.P2D9.CB,CM.P2F11.CB,HY.P3F5.CB,CM.P3F8.CB,CM.P4E1.CB,HY.P1E4.NRW,HY.P1F4.CB,CM.P3F7.CB,CM.P3F2.NA,HY.P1F5.NRW"
    # removed: "CM.P1H10.NRW", ,CM.P3A6.CF,CM.P2H11.NRW

################################
## Consolidated Hybrid string ##
################################

# All hybrids and backcrosses 
HY_string=()
HY_string=$(awk -F'\t' '$2 ~ /Hybrid|backcross/ { objects = objects (objects ? "," : "") $1 } END { print objects }' "$table")
echo ${HY_string}
    # HY_string="CM.PC4.CF,CM.PC3.CF,CM.P3A8.CB,HY.P2A6.CB,CM.P3G6.CB,HY.P3G1.CB,HY.P3D8.NA,CM.P3E7.CB,CM.P3D7.CB,CM.P4B2.NA,CM.P3H2.NRW,CM.P3C7.CB,HY.P2B6.CB,HY.P3C4.CB,CM.P3E6.NA,HY.P3H6.CB,HY.P3B7.CB,CM.P3G2.CB,CM.P2B10.NRW,HY.P2H6.CB,HY.P3H8.CB,CM.P4D1.NSLW,HY.P2B11.CB,HY.P1B4.NA,MCA.P1E7.CB,MCA.P1B7.CB,HY.P1B5.CB,MCA.P1A8.CB,MCA.P1D7.CB,MCA.P1A7.CB,HY.P3H11.NSLW,,MCA.P1C7.CB,MCA.P1H7.CB, HY.P2E9.NA,HY.P2F9.NA,MCA.P1C6.NTW,HY.P3G11.NSLW,MCA.P3E11.STW,MCA.P1G5.CB,CM.P1B9.NSLW,MCA.P1D5.BRW,MCA.P2G10.NA,MCA.P3B11.CB,CA.P1D12.NC,MCA.P1A6.NTW,MCA.P1B2.STW,HY.P1A2.NTW,CM.P1A9.NSLW,MCA.P1F7.BC,MCA.P1B8.CB,CA.P2A2.NC,MCA.P1H8.BRW,CA.P2C2.NC,CA.P1D11.NC,CA.P1G12.NC,CA.P3E9.NC,HY.P1A5.CB,MCA.P1C5.CB,MCA.P1D6.STW,HY.P1D4.NA,MCA.P1D8.CB,HY.P1A4.CB,MCA.P1B6.CB,HY.P2C9.CB,MCA.P1E8.CB,MCA.P1F8.NTW,MCA.P1G7.CB"
    # removed:  MCA.P1H6.CB,CM.P4G1.CB,HY.P4A2.CB,MCA.P1H5.CB,MCA.P1C11.NRW,CA.P3C10.NC
    
# F1 hybrids as determined as anything between 40-60%
F1_Hybrid_string=()
F1_Hybrid_string=$(awk -F'\t' '$2 ~ /F1_Hybrid/ { printf("%s,", $1) } END { ORS=""; print "" }' "$table")
F1_Hybrid_string=${F1_Hybrid_string%,}
echo ${F1_Hybrid_string}
    # F1_Hybrid_string="HY.P2H6.CB,CM.P4D1.NSLW,HY.P2B11.CB,HY.P1B4.NA,MCA.P1E7.CB,MCA.P1B7.CB,HY.P1B5.CB,MCA.P1A8.CB,MCA.P1D7.CB,MCA.P1A7.CB,MCA.P1C7.CB,MCA.P1C6.NTW"

######################################################################################################################################################
## Run each test comparison and put all results in "/output/ folder"    
VCF=/datadirectory/data/noreponly_v2.filtered.85.renamed.vcf
FILE=noreponly_v2.filtered.85.renamed
OUT=/datadirectory/out/$FILE

## Figure 6: ACUTUS_99 -> CM_99
# P1 <- "CM_99"
# P2 <- Hybrid_CA_backcross, Hybrid_MCA_backcross, Acutus_backcross, F1_Hybrid, CM_backcross
Hybrid_strings=()
    Hybrid_strings="$CM_backcross_string,$F1_Hybrid_string,$acutus_backcross_string,$MCA_backcross_string,$CA_backcross_string"
    # Remove the last comma and any trailing spaces
    Hybrid_strings=${Hybrid_strings%,}
    echo ${Hybrid_strings}
        # Hybrid_strings="CM.PC6.NRW,CM.PC8.NRW,CM.P1D3.NRW,CM.P3A3.NRW,CM.P3B3.NRW,CM.P3C1.NRW,CM.P3C3.NRW,CM.P3D3.NRW,CM.P3G3.NRW,CM.PC7.NRW,HY.P3A2.NRW,HY.P3B5.NRW,HY.P3C5.NRW,HY.P3D5.NRW,HY.P3H4.NRW,CM.P1F10.NRW,CM.P3F3.NRW,CM.P1G10.NRW,HY.P3B4.NRW,CM.P1A10.NRW,CM.P3C8.CB,HY.P3A5.NRW,CM.P3E3.NRW,CM.P3B1.NRW,HY.P3F4.BRW,CM.P3D6.RHW,CM.P3G4.NRW,CM.P1B10.NRW,HY.P3B2.NRW,CM.P3A4.NRW,CM.P2D10.NRW,HY.P3E4.BRW,CM.P1D10.NRW,CM.P2C6.NRW,CM.P2H12.NRW,CM.P2F10.NRW,HY.P2D11.RHW,CM.P1E10.NRW,CM.P2A12.NRW,CM.P1D9.NRW,CM.P2C10.NRW,CM.P3C6.NRW,CM.P3H3.NRW,CM.P2C12.NRW,CM.P2B12.NRW,CM.P3B6.NRW,CM.P1C10.NRW,HY.P3E8.NA,CM.P3A1.NRW,HY.P3E5.CB,CM.P2D12.NRW,CM.P3F6.NTW,CM.P2E12.NRW,CM.P2G6.NRW,CM.P3A7.CB,HY.P3D2.BRW,HY.P3D4.BRW,CM.P3G7.CB,CM.PC5.CF,HY.P3C2.BRW,HY.P3H5.BRW,CM.PC2.CF,CM.P2G9.NRW,CM.P2E10.NRW,HY.P3E2.BRW,CM.P2G12.NRW,HY.P3B8.CB,CM.P2H9.NRW,CM.P1B3.NRW,CM.P2A10.NRW,CM.P1H9.NRW,HY.P3G5.CB,CM.P3G8.CB,HY.P1G4.CB,HY.P1G3.NRW,CM.P1C9.NSLW,HY.P2D9.CB,CM.P2F11.CB,HY.P3F5.CB,CM.P3F8.CB,CM.P4E1.CB,HY.P1E4.NRW,HY.P1F4.CB,CM.P3F7.CB,CM.P3F2.NA,HY.P1F5.NRW,HY.P2H6.CB,CM.P4D1.NSLW,HY.P2B11.CB,HY.P1B4.NA,MCA.P1E7.CB,MCA.P1B7.CB,HY.P1B5.CB,MCA.P1A8.CB,MCA.P1D7.CB,MCA.P1A7.CB,MCA.P1C7.CB,MCA.P1C6.NTW,HY.P3H11.NSLW,HY.P2E9.NA,HY.P2F9.NA,HY.P3G11.NSLW,CM.P1B9.NSLW,MCA.P1A6.NTW,CM.P1A9.NSLW,MCA.P1G7.CB,MCA.P1H7.CB,MCA.P3E11.STW,MCA.P1G5.CB,MCA.P2G10.NA,MCA.P3B11.CB,MCA.P1B2.STW,HY.P1A2.NTW,MCA.P1B8.CB,MCA.P1H8.BRW,HY.P1A5.CB,MCA.P1C5.CB,MCA.P1D6.STW,HY.P1D4.NA,MCA.P1D8.CB,HY.P1A4.CB,MCA.P1B6.CB,HY.P2C9.CB,MCA.P1E8.CB,MCA.P1F8.NTW,MCA.P1D5.BRW,CA.P1D12.NC,MCA.P1F7.BC,CA.P2A2.NC,CA.P2C2.NC,CA.P1D11.NC,CA.P1G12.NC,CA.P3E9.NC"
# P3 <- CA_99, MCA_99
    Acutus_99_strings=()
    Acutus_99_strings="$MCA_99_string,$CA_99_string"
    # Remove the last comma and any trailing spaces
    Acutus_99_strings=${Acutus_99_strings%,}
    echo ${Acutus_99_strings}
        # Acutus_99_strings="MCA.P1F2.STW,MCA.P4C1.STW,MCA.P3D11.STW,MCA.P1H2.STW,MCA.P1C2.STW,MCA.P1D2.STW,MCA.P1E2.STW,MCA.P1E6.STW,MCA.P1F6.STW,MCA.P1G2.STW,MCA.P2H10.NA,MCA.P4G2.CB,CA.P3G10.CC,CA.P3F10.CC,CA.P2E8.CC,CA.P2E3.CC,CA.P2A3.CC,CA.P1H12.NC,CA.P2C7.CC,CA.P3C9.NC,CA.PC11.NC,CA.P1A12.NC,CA.P1C12.NC,CA.PC10.NC,CA.P1F12.NC,CA.P1G11.NC,CA.P3D9.NC,CA.P2C3.CC,CA.P2B3.CC,CA.P2G3.CC,CA.P2B8.CC,CA.P1E11.NC,CA.P2G5.CC,CA.P2A5.CC,CA.P2B4.CC,CA.P2B7.CC,CA.P2C4.CC,CA.P2D2.NC,CA.P2E4.CC,CA.P2E5.CC,CA.P2E7.CC,CA.P2F4.CC,CA.P2F7.CC,CA.P2G4.CC,CA.P2G8.CC,CA.P2H5.NA,CA.P2H8.CC,CA.PC12.NC,CA.PC14.CC,CA.PC15.CC,CA.PC16.CC,CA.PC9.NC,MCA.P1B11.NRW"
ruby get_fixed_site_gts.rb $VCF ./data/$FILE.test15b.txt ${CM_99_string} ${Hybrid_strings} ${Acutus_99_strings} 0.98
    # Found 5130 sites with sufficient completeness (at least 0.98).
    # Found 114 sites with complete fixation of different alleles in the parental populations.
ruby plot_fixed_site_gts.rb ./data/$FILE.test15b.txt $OUT.test15b.svg 0.98 > $OUT.test15b.log

## SUPP FIGURE S9: ACUTUS_99 -> CM_99
# P1 <- "CM_99"
# P2 <- F1_Hybrid
    F1_strings="$F1_Hybrid_string"
    F1_strings=${F1_strings%,}
# P3 <- CA_99, MCA_99
    Acutus_99_strings=()
    Acutus_99_strings="$MCA_99_string,$CA_99_string"
ruby get_fixed_site_gts.rb $VCF ./data/$FILE.test16.txt ${CM_99_string} ${F1_strings} ${Acutus_99_strings} 0.95
    # Found 9899 sites with sufficient completeness (at least 0.95).
    # Found 275 sites with complete fixation of different alleles in the parental populations.
ruby plot_fixed_site_gts.rb ./data/$FILE.test16.txt $OUT.test16.svg 0.95 > $OUT.test16.log

## SUPP FIGURE S10: CM_99 --> CA_99
# P1 <- "CM_99"
# P2 <- Hybrid_CA_backcross, Hybrid_MCA_backcross, Acutus_backcross, F1_Hybrid, CM_backcross
Hybrid_strings=()
    Hybrid_strings="$CM_backcross_string,$F1_Hybrid_string,$acutus_backcross_string,$MCA_backcross_string,$CA_backcross_string"
    # Remove the last comma and any trailing spaces
    Hybrid_strings=${Hybrid_strings%,}
    echo ${Hybrid_strings}
        # Hybrid_strings="CM.PC6.NRW,CM.PC8.NRW,CM.P1D3.NRW,CM.P3A3.NRW,CM.P3B3.NRW,CM.P3C1.NRW,CM.P3C3.NRW,CM.P3D3.NRW,CM.P3G3.NRW,CM.PC7.NRW,HY.P3A2.NRW,HY.P3B5.NRW,HY.P3C5.NRW,HY.P3D5.NRW,HY.P3H4.NRW,CM.P1F10.NRW,CM.P3F3.NRW,CM.P1G10.NRW,HY.P3B4.NRW,CM.P1A10.NRW,CM.P3C8.CB,HY.P3A5.NRW,CM.P3E3.NRW,CM.P3B1.NRW,HY.P3F4.BRW,CM.P3D6.RHW,CM.P3G4.NRW,CM.P1B10.NRW,HY.P3B2.NRW,CM.P3A4.NRW,CM.P2D10.NRW,HY.P3E4.BRW,CM.P1D10.NRW,CM.P2C6.NRW,CM.P2H12.NRW,CM.P2F10.NRW,HY.P2D11.RHW,CM.P1E10.NRW,CM.P2A12.NRW,CM.P1D9.NRW,CM.P2C10.NRW,CM.P3C6.NRW,CM.P3H3.NRW,CM.P2C12.NRW,CM.P2B12.NRW,CM.P3B6.NRW,CM.P1C10.NRW,HY.P3E8.NA,CM.P3A1.NRW,HY.P3E5.CB,CM.P2D12.NRW,CM.P3F6.NTW,CM.P2E12.NRW,CM.P2G6.NRW,CM.P3A7.CB,HY.P3D2.BRW,HY.P3D4.BRW,CM.P3G7.CB,CM.PC5.CF,HY.P3C2.BRW,HY.P3H5.BRW,CM.PC2.CF,CM.P2G9.NRW,CM.P2E10.NRW,HY.P3E2.BRW,CM.P2G12.NRW,HY.P3B8.CB,CM.P2H9.NRW,CM.P1B3.NRW,CM.P2A10.NRW,CM.P1H9.NRW,HY.P3G5.CB,CM.P3G8.CB,HY.P1G4.CB,HY.P1G3.NRW,CM.P1C9.NSLW,HY.P2D9.CB,CM.P2F11.CB,HY.P3F5.CB,CM.P3F8.CB,CM.P4E1.CB,HY.P1E4.NRW,HY.P1F4.CB,CM.P3F7.CB,CM.P3F2.NA,HY.P1F5.NRW,HY.P2H6.CB,CM.P4D1.NSLW,HY.P2B11.CB,HY.P1B4.NA,MCA.P1E7.CB,MCA.P1B7.CB,HY.P1B5.CB,MCA.P1A8.CB,MCA.P1D7.CB,MCA.P1A7.CB,MCA.P1C7.CB,MCA.P1C6.NTW,HY.P3H11.NSLW,HY.P2E9.NA,HY.P2F9.NA,HY.P3G11.NSLW,CM.P1B9.NSLW,MCA.P1A6.NTW,CM.P1A9.NSLW,MCA.P1G7.CB,MCA.P1H7.CB,MCA.P3E11.STW,MCA.P1G5.CB,MCA.P2G10.NA,MCA.P3B11.CB,MCA.P1B2.STW,HY.P1A2.NTW,MCA.P1B8.CB,MCA.P1H8.BRW,HY.P1A5.CB,MCA.P1C5.CB,MCA.P1D6.STW,HY.P1D4.NA,MCA.P1D8.CB,HY.P1A4.CB,MCA.P1B6.CB,HY.P2C9.CB,MCA.P1E8.CB,MCA.P1F8.NTW,MCA.P1D5.BRW,CA.P1D12.NC,MCA.P1F7.BC,CA.P2A2.NC,CA.P2C2.NC,CA.P1D11.NC,CA.P1G12.NC,CA.P3E9.NC"
# P3 <- CA_99
ruby get_fixed_site_gts.rb $VCF ./data/$FILE.test17.txt ${CM_99_string} ${Hybrid_strings} ${CA_99_string} 0.95
    # Found 12383 sites with sufficient completeness (at least 0.95).
    # Found 402 sites with complete fixation of different alleles in the parental populations.
ruby plot_fixed_site_gts.rb ./data/$FILE.test17.txt $OUT.test17.svg 0.95 > $OUT.test17.log

## SUPP FIGURE S11: CM_99 --> CA_99
# P1 <- "CM_99"
# P2 <- F1_Hybrid
# P3 <- CA_99
ruby get_fixed_site_gts.rb $VCF ./data/$FILE.test17b.txt ${CM_99_string} ${F1_strings} ${CA_99_string} 0.95
# Found 12383 sites with sufficient completeness (at least 0.95).
# Found 402 sites with complete fixation of different alleles in the parental populations.
ruby plot_fixed_site_gts.rb ./data/$FILE.test17b.txt $OUT.test17b.svg 0.95 > $OUT.test17b.log

## SUPP FIGURE S12: : CM_99 --> MCA_99
# P1 <- "CM_99"
# P2 <- Hybrid_CA_backcross, Hybrid_MCA_backcross, Acutus_backcross, F1_Hybrid, CM_backcross
Hybrid_strings=()
    Hybrid_strings="$CM_backcross_string,$F1_Hybrid_string,$acutus_backcross_string,$MCA_backcross_string,$CA_backcross_string"
    # Remove the last comma and any trailing spaces
    Hybrid_strings=${Hybrid_strings%,}
    echo ${Hybrid_strings}
        # Hybrid_strings="CM.PC6.NRW,CM.PC8.NRW,CM.P1D3.NRW,CM.P3A3.NRW,CM.P3B3.NRW,CM.P3C1.NRW,CM.P3C3.NRW,CM.P3D3.NRW,CM.P3G3.NRW,CM.PC7.NRW,HY.P3A2.NRW,HY.P3B5.NRW,HY.P3C5.NRW,HY.P3D5.NRW,HY.P3H4.NRW,CM.P1F10.NRW,CM.P3F3.NRW,CM.P1G10.NRW,HY.P3B4.NRW,CM.P1A10.NRW,CM.P3C8.CB,HY.P3A5.NRW,CM.P3E3.NRW,CM.P3B1.NRW,HY.P3F4.BRW,CM.P3D6.RHW,CM.P3G4.NRW,CM.P1B10.NRW,HY.P3B2.NRW,CM.P3A4.NRW,CM.P2D10.NRW,HY.P3E4.BRW,CM.P1D10.NRW,CM.P2C6.NRW,CM.P2H12.NRW,CM.P2F10.NRW,HY.P2D11.RHW,CM.P1E10.NRW,CM.P2A12.NRW,CM.P1D9.NRW,CM.P2C10.NRW,CM.P3C6.NRW,CM.P3H3.NRW,CM.P2C12.NRW,CM.P2B12.NRW,CM.P3B6.NRW,CM.P1C10.NRW,HY.P3E8.NA,CM.P3A1.NRW,HY.P3E5.CB,CM.P2D12.NRW,CM.P3F6.NTW,CM.P2E12.NRW,CM.P2G6.NRW,CM.P3A7.CB,HY.P3D2.BRW,HY.P3D4.BRW,CM.P3G7.CB,CM.PC5.CF,HY.P3C2.BRW,HY.P3H5.BRW,CM.PC2.CF,CM.P2G9.NRW,CM.P2E10.NRW,HY.P3E2.BRW,CM.P2G12.NRW,HY.P3B8.CB,CM.P2H9.NRW,CM.P1B3.NRW,CM.P2A10.NRW,CM.P1H9.NRW,HY.P3G5.CB,CM.P3G8.CB,HY.P1G4.CB,HY.P1G3.NRW,CM.P1C9.NSLW,HY.P2D9.CB,CM.P2F11.CB,HY.P3F5.CB,CM.P3F8.CB,CM.P4E1.CB,HY.P1E4.NRW,HY.P1F4.CB,CM.P3F7.CB,CM.P3F2.NA,HY.P1F5.NRW,HY.P2H6.CB,CM.P4D1.NSLW,HY.P2B11.CB,HY.P1B4.NA,MCA.P1E7.CB,MCA.P1B7.CB,HY.P1B5.CB,MCA.P1A8.CB,MCA.P1D7.CB,MCA.P1A7.CB,MCA.P1C7.CB,MCA.P1C6.NTW,HY.P3H11.NSLW,HY.P2E9.NA,HY.P2F9.NA,HY.P3G11.NSLW,CM.P1B9.NSLW,MCA.P1A6.NTW,CM.P1A9.NSLW,MCA.P1G7.CB,MCA.P1H7.CB,MCA.P3E11.STW,MCA.P1G5.CB,MCA.P2G10.NA,MCA.P3B11.CB,MCA.P1B2.STW,HY.P1A2.NTW,MCA.P1B8.CB,MCA.P1H8.BRW,HY.P1A5.CB,MCA.P1C5.CB,MCA.P1D6.STW,HY.P1D4.NA,MCA.P1D8.CB,HY.P1A4.CB,MCA.P1B6.CB,HY.P2C9.CB,MCA.P1E8.CB,MCA.P1F8.NTW,MCA.P1D5.BRW,CA.P1D12.NC,MCA.P1F7.BC,CA.P2A2.NC,CA.P2C2.NC,CA.P1D11.NC,CA.P1G12.NC,CA.P3E9.NC"
# P3 <- MCA_99
ruby get_fixed_site_gts.rb $VCF ./data/$FILE.test18.txt ${CM_99_string} ${Hybrid_strings} ${MCA_99_string} 0.95
# Found 6865 sites with sufficient completeness (at least 0.95).
# Found 255 sites with complete fixation of different alleles in the parental populations.
ruby plot_fixed_site_gts.rb ./data/$FILE.test18.txt $OUT.test18.svg 0.95 > $OUT.test18.log

## SUPP FIGURE S13: CM_99 --> MCA_99
# P1 <- "CM_99"
# P2 <- F1_Hybrid
# P3 <- MCA_99
ruby get_fixed_site_gts.rb $VCF ./data/$FILE.test18b.txt ${CM_99_string} ${F1_strings} ${MCA_99_string} 0.95
# Found 6865 sites with sufficient completeness (at least 0.95).
# Found 255 sites with complete fixation of different alleles in the parental populations.
ruby plot_fixed_site_gts.rb ./data/$FILE.test18b.txt $OUT.test18b.svg 0.95 > $OUT.test18b.log

######################################################################################################################################################
