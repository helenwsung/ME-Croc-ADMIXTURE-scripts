# ME-AdmixtureCode
Scripts for ME Manuscript titled: "Out with the old, introgression with the new: Signals of ancient and recent admixture in hybridizing Mesoamerican crocodiles (Crocodylus acutus x Crocodylus moreletii)"

# Files can be found on Data DRYAD: https://doi.org/10.5061/dryad.3bk3j9kt9

Scripts for running analyses and their respective input/output files. All scripts for running the analyses are given as is and need to be edited to match personal data directories and modified to your specific needs).

NOTE: In all scripts and data files, 'Sampling localities' in the main manuscript are listed as "Monitoring unit".

Description of the data and file structure

1. ipyrad

    params-noreponly_v2.txt - input parameter file used for SNP calling and genotyping in ipyrad

    ipyradpipeline.sh - pipeline script to run trimming and ipyrad

2. Filtering 


    SNPfiltR_filtering.R - R script for filtering ipyrad vcf file to working datasets and creating PCA & tSNE plots in Supplementary Figures S3-S4

    LDpruning_Dataconversion.R - R script for LD pruning filtered vcf files and creating subsetted VCF files 

3. PopulationStats_Admixture

    ADMIXTURE.sh - pipeline & script for running ADMIXTURE analyses

    admixture_k3_pops - folder of text files for samples in each population group created from ADMIXTURE K3

    Popstruct_Plots.R - R script for running population structure analyses and generating plots (e.g. ADMIXTURE plots)

    popmap.csv - condensed population map & metadata created for only the 242 working samples used in ADMIXTURE analysis

    PCAscript.R - R script to create PCA plots

    PopGenStats.sh - pipeline script for running VCFtools to generate population genetic summary stats 

    PopGenStats.R - R script for generating population diversity statistics 

4. FSC2 

    FSC2.sh - Pipeline for running Fastsimcoal2 script (need additional scripts to run, descriptions in text)

        FS_FixRootTime_K3_Mods.slurm 

        FSC_croc_boot.slurm

        fsc-selectbestrun.sh

        Get_AIC_across_mods.R

        Get_best_FSCacross_mods.sh

        Get_best_FSCacross_boots.sh

        Get_pars_across_bootreps_crocs.R

        Prep_FSC_reps.sh

    FSC2.R - processing FSC2 results and generating plots 

    Ad_k3_pop90 - Folder containing subfolder for FSC2 input 

          Models_16.2y - folder containing subfolders for each of the 16 tested models and .obs file generated via easySFS. Model folders each contain an .est            and .tpl input parameter file for FSC2 

    noreponly_v2.75.renamed.LDpruned.popmap90.txt - input popmap file for running FSC2

5. LD

    PlinkLD.sh - generating LD files using Plink 

    LDscript.R - Rscript for generating LD decay plots 

6. AncestryPaint_abbababa

    AncestryPainting.sh - script for generating Ancestry painting plots

    snpRabbababa.R - Rscript for ABBA-BABA analyses

noreponly_metadata_fixed.csv - metadata file for all 273 samples 
