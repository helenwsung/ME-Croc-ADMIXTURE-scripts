# ME-Croc-ADMIXTURE-scripts
Scripts for Manuscript titled: "Out with the old, introgression with the new: Signals of ancient and recent admixture in hybridizing Mesoamerican crocodiles (Crocodylus acutus x Crocodylus moreletii)"

Files can be found on Data DRYAD: https://doi.org/10.5061/dryad.3bk3j9kt9

Scripts for running analyses and their respective input/output files. All scripts for running the analyses are given as is and will need to be edited to match personal data directories and modified to your specific needs).

NOTE: The 'Sampling localities' in the main manuscript are listed as "Monitoring units" in all analysis scripts and data files.

Description of the data and folder structure: 

1. ipyrad

        params-noreponly_v2.txt - input parameter file used for SNP calling and genotyping in ipyrad
   
        ipyradpipeline.md - instructions on running ipyrad pipeline 

2. SNPfiltR
   
        SNPfiltR_filtering.R - R script for filtering ipyrad vcf file to working datasets and creating PCA & tSNE plots in Supplementary Figures S3-S4
   
        LDpruning_Dataconversion.R - R script for LD pruning filtered vcf files and creating subsetted VCF files 

3. PopulationStats_Admixture
   
        ADMIXTURE.md - instructions for the pipeline & script for running ADMIXTURE analyses 
   
        admixture_k3_pops - folder of .txt files for the ADMIXTURE population groups assignments containing the respective Sample IDs in each group
   
        Popstruct_Plots.R - R script for generating plots from the ADMIXTURE results.
   
        popmap.csv - condensed population map & metadata for only the 242 working samples used in the ADMIXTURE analysis created from the Popstruct_Plots.R script. 
   
        PCAscript.R - R script to create PCA plots

        PopGenStats.md - instructions for the pipeline & scripts to run VCFtools for generating population genetic summary stats

        PopGenStats.R - R script for generating population diversity statistics

5. FSC2
   
        FSC2.md - Instructions for the pipeline and scripts to run Fastsimcoal2 analyses (needs additional scripts to run, descriptions in text)
   
            FS_FixRootTime_K3_Mods.slurm
   
            FSC_croc_boot.slurm
   
            fsc-selectbestrun.sh
   
            Get_AIC_across_mods.R
   
            Get_best_FSCacross_mods.sh
   
            Get_best_FSCacross_boots.sh
   
            Get_pars_across_bootreps_crocs.R
   
            Prep_FSC_reps.sh
   
        FSC2.R - R script for processing FSC2 results and generating plots
   
        Ad_k3_pop90 - Folder containing subfolder for FSC2 input
   
           Models_16.2y - folder containing subfolders for each of the 16 tested models and .obs file generated via easySFS. Model folders each contain an .est and .tpl input parameter file for FSC2
   
           noreponly_v2.75.renamed.LDpruned.popmap90.txt - input popmap file for running FSC2

6. LD
       PlinkLD.md - instructions for running the pipeline and scripts to generate the LD files using Plink
   
       LDscript.R - Rscript for generating LD decay plots
   
7. AncestryPaint_abbababa
   
        AncestryPainting.md - instructions for the pipeline and scripts used to generate the Ancestry painting plots
    
        snpRabbababa.R - R script for ABBA-BABA analyses

noreponly_metadata_fixed.csv - metadata file for all 273 samples. 
        # Definition of columns for metadata in noreponly_metadata_fixed.csv:
                Sample - Unique Sample ID used for analyses
                Seq_ID - Unique Sequencer ID 
                CrocID - Unique individual ID combining Mark Code and Morph_Species
                Longitude - EDITED Longitude value**	
                Latitude - EDITED Latitude value**		
                Capture_Date - Capture date by date-month-year	
                Monitoring.Unit	- Sampling locality 
                Abbrv_Monitoring.Unit - Abbreviation for Monitoring.Unit 
                Subdivision	- District in Belize
                Water Temp (C) - Water temperature (celcius)
                Air Temp (C) - Air temperature (celcius)	
                Salinity (ppt) - Water Salinity (parts per thousand)	
                pH - water pH
                Mark Code - Unique mark code for sampled individual 	
                Morph_Species - Morphological species group	
                Size Class - Size class (Adult, Subadult, Juvenile, Hatchling) 	
                Sex 
                HL (cm)	- Head length (centimeters)
                SL (cm)	- Snout length (centimeters)
                CW (cm)	- Cranial width (centimeters)
                MAX W (cm) - Maxilary Width (centimeters)
                PMax W (cm)	- Premaxilary width (centimeters)
                POb L (cm) - Preorbital length (centimeters)	
                POb W (cm) - Preorbital width (centimeters)	
                TL (cm)	- Total length (centimeters)
                SVL (cm) - Snout vent length (centimeters)	
                TW (cm)	- Tail width (centimeters)
                HF (cm)	- Hind Right foot lenth (centimeters)
                Weight (kg)	- Weight (kilograms)
                Paratrichosoma - Presence/Absence for Paratrichosoma 
                Tail Condition 
                Skin Condition	
                Musculature Condition	
                Teeth Condition	
                Skeletal Structure	
                Nuchal Photo	
                Post-Occiptal	
                Nuchal Scutes	
                Post-Occipital/Nuchal Pattern	
                Transverse Rows	
                Double Whorls	
                Single Whorls	
                Notes	
                Wounds

        ** NOTE: In order to preserve sensitive occurence data for threatened/at-risk species, we generalized the precision of the geographic coordinates (Lat/Long) by reducing the number of decimal places to 0.1 decimal degrees as recommended by Guide to Best Practices for Generalising Sensitive Species Occurrence Data [Chapman AD (2020) Current Best Practices for Generalizing Sensitive Species Occurrence Data. Copenhagen: GBIF Secretariat. https://doi.org/10.15468/doc-5jp4-5g10.]. Please contact the corresponding author Helen Sung (hwsung@hawaii.edu) if you would like to request this data. 

        
