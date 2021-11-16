#!/bin/sh

### Job Name
#SBATCH --job-name easySFS_preview

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_easySFS_preview
#SBATCH -e stderr_easySFS_preview

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
####SBATCH --mem=0
#SBATCH --mem=250G

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=5-00:00:00



##################################
### 1. SETTING UP THE ANALYSES ###
##################################

###Load libraries
module load gcc
module load htslib
module load miniconda3
module load bcftools
module load r/3.6.1

source activate easySFS

### Paths to directories and/or files ###
#what is the path to the script directory?
script_directory=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts



##################################################################
### 2. RUNNING AND PROCESSING easySFS PREVIEW FOR NORTH REGION ###
##################################################################

### North region paths and names ###
file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/sfs_creation
vcf_file_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/allsites_datasets_6_16_2021/north_allsites_polykaz
vcf_file_name=north_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.vcf.gz
sample_info_file=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/entropy/entropy_6_13_2021/qmodel_fulldataset/processed_results/k2_processed_results.txt
number_of_top_projections=10


vcf_basename=$(basename -s .vcf.gz $vcf_file_name)

### CREATING POP SPEC FILE ###
echo 'creating north region pop spec file'

bcftools query -l $vcf_file_path/$vcf_file_name > $vcf_file_path/SAMPLES_${vcf_basename}.txt
awk 'NR==FNR{a[$1]; next} ($1 in a) {print $1"\t"$4}' $vcf_file_path/SAMPLES_${vcf_basename}.txt $sample_info_file | sort -k 2 > $file_path/pop_spec_${vcf_basename}.txt

if [ "$(wc -l < $file_path/pop_spec_${vcf_basename}.txt)" -ne "$(wc -l < $vcf_file_path/SAMPLES_${vcf_basename}.txt)" ]
then 
  echo "the pop spec file does NOT include the same number of samples as the vcf"
  exit 1
else
  echo "the pop spec file does include the same number of samples as the vcf"
fi


### CONDUCTING easySFS PREVIEW ###
echo 'starting north region preview'

#Step 1: Use the preview function in easySFS to calculate the number of segregating sites for each projection per population
/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $vcf_file_path/$vcf_file_name -p $file_path/pop_spec_${vcf_basename}.txt --preview -a -v > $file_path/PREVIEW_${vcf_basename}.txt

#Step 2: Identify the best projections for each population
Rscript $script_directory/process_easySFS_preview.R -f $file_path/PREVIEW_${vcf_basename}.txt -o $file_path/TOPPROJECTS_PREVIEW_${vcf_basename}.txt -n $number_of_top_projections

echo 'finished north region preview'




################################################################
### 3. RUNNING AND PROCESSING easySFS PREVIEW FOR MID REGION ###
################################################################

### Mid region paths and names ###
file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/sfs_creation
vcf_file_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/allsites_datasets_6_16_2021/mid_allsites_polykaz
vcf_file_name=mid_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.vcf.gz
sample_info_file=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/entropy/entropy_6_13_2021/qmodel_fulldataset/processed_results/k2_processed_results.txt
number_of_top_projections=10


vcf_basename=$(basename -s .vcf.gz $vcf_file_name)

### CREATING POP SPEC FILE ###
echo 'creating mid region pop spec file'

bcftools query -l $vcf_file_path/$vcf_file_name > $vcf_file_path/SAMPLES_${vcf_basename}.txt
awk 'NR==FNR{a[$1]; next} ($1 in a) {print $1"\t"$4}' $vcf_file_path/SAMPLES_${vcf_basename}.txt $sample_info_file | sort -k 2 > $file_path/pop_spec_${vcf_basename}.txt

if [ "$(wc -l < $file_path/pop_spec_${vcf_basename}.txt)" -ne "$(wc -l < $vcf_file_path/SAMPLES_${vcf_basename}.txt)" ]
then
  echo "the pop spec file does NOT include the same number of samples as the vcf"
  exit 1
else
  echo "the pop spec file does include the same number of samples as the vcf"
fi


### CONDUCTING easySFS PREVIEW ###
echo 'starting mid region preview'

#Step 1: Use the preview function in easySFS to calculate the number of segregating sites for each projection per population
/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $vcf_file_path/$vcf_file_name -p $file_path/pop_spec_${vcf_basename}.txt --preview -a -v > $file_path/PREVIEW_${vcf_basename}.txt

#Step 2: Identify the best projections for each population
Rscript $script_directory/process_easySFS_preview.R -f $file_path/PREVIEW_${vcf_basename}.txt -o $file_path/TOPPROJECTS_PREVIEW_${vcf_basename}.txt -n $number_of_top_projections

echo 'finished mid region preview'






#################################
#################################
### CODE NOT CURRENTLY IN USE ###
#################################
#################################

##################################
### 1. SETTING UP THE ANALYSES ###
##################################

### Paths to directories and/or files ###
#what is the path to the script directory?
#script_directory=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts

### North region paths and names ###
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/sfs_creation
##vcf_file_name=NORTH_SUBSET_kazumbe_polyodon_pundamilia_8_19_2020_maf0_maxmissing0.90_minDP5_post_4.recode.vcf
#vcf_file_name=NORTH_SUBSET_kazumbe_polyodon_pundamilia_8_19_2020_maf0_maxmissing0.90_minDP5_post_4_entropysubset.recode.vcf
#pop_spec_name=north_kazumbe_polyodon_popspec_maf0_maxmissing90_1_19_21.txt
#preview_file_name=north_easysfs_preview_maf0_maxmissing90.txt
#output_file_name=north_easysfs_preview_maf0_maxmissing90_processed.txt
#number_of_top_projections=5

### Mid region paths and names ###
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/sfs_creation
#vcf_file_name=MID_SUBSET_kazumbe_polyodon_pundamilia_8_19_2020_maf0_maxmissing0.90_minDP5_post_4.recode.vcf
#pop_spec_name=mid_kazumbe_polyodon_popspec_maf0_maxmissing90_1_19_21.txt
#preview_file_name=mid_easysfs_preview_maf0_maxmissing90.txt
#output_file_name=mid_easysfs_preview_maf0_maxmissing90_processed.txt
#number_of_top_projections=5


############################
### 2. RUNNING  ANALYSES ###
############################

### North region paths and names ###
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/sfs_creation
#vcf_file_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/all_sites_vcf/north_region_vcf/subset_file
#vcf_file_name=NORTH_SUBSET_polyodon_kazumbe_allsites_3_4_2021_subsample_post_4.vcf.gz
#pop_spec_name=pop_spec_NORTH_SUBSET_polyodon_kazumbe_allsites_3_4_2021_subsample_post_4.txt
#preview_file_name=north_easysfs_preview_maxmissing80_monomorph_3_4_2021.txt
#output_file_name=north_easysfs_preview_maxmissing80_monomorph_3_4_2021_topprojects.txt
#number_of_top_projections=5

#awk '{print $1, $2}' /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/final_variant_files/fsc_vcf_processing/subsampled_fsc_dataset_north.txt > $file_path/$pop_spec_name

#echo 'starting north region preview'

#Step 1: Use the preview function in easySFS to calculate the number of segregating sites for each projection per population
#/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $vcf_file_path/$vcf_file_name -p $file_path/$pop_spec_name --preview -a -v > $file_path/$previ$

#Step 2: Identify the best projections for each population
#Rscript $script_directory/process_easySFS_preview.R -f $file_path/$preview_file_name -o $file_path/$output_file_name -n $number_of_top_projections
#R CMD BATCH spp_specific_mods_plot_ALT.R

#echo 'finished north region preview'

### North region paths and names ###
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/sfs_creation
#vcf_file_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/all_sites_vcf/north_region_vcf/file_processing_3_24_2021
#vcf_file_name=NORTH_SUBSET_polyodon_kazumbe_allsites_2_5_2021_post_4_entropycongruent.vcf.gz
#pop_spec_name=pop_spec_NORTH_SUBSET_polyodon_kazumbe_allsites_2_5_2021_post_4_missing75.txt
#preview_file_name=north_easysfs_preview_maxmissing75_monomorph_2_5_2021.txt
#output_file_name=north_easysfs_preview_maxmissing75_monomorph_2_5_2021_topprojects.txt
#number_of_top_projections=5

#bcftools query -l $vcf_file_path/$vcf_file_name        > $vcf_file_path/SAMPLES_NORTH_SUBSET_polyodon_kazumbe_allsites_2_5_2021_post_4_entropycongruent.txt

#awk 'NR==FNR{a[$1]; next} ($1 in a) {print $1, $2}' $vcf_file_path/SAMPLES_NORTH_SUBSET_polyodon_kazumbe_allsites_2_5_2021_post_4_entropycongruent.txt /project/ltcichlidgenomics/alewansk/Petrochrom$

#echo 'starting north region preview'

##Step 1: Use the preview function in easySFS to calculate the number of segregating sites for each projection per population
#/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $vcf_file_path/$vcf_file_name -p $file_path/$pop_spec_name --preview -a -v > $file_path/$previ$

##Step 2: Identify the best projections for each population
#Rscript $script_directory/process_easySFS_preview.R -f $file_path/$preview_file_name -o $file_path/$output_file_name -n $number_of_top_projections
##R CMD BATCH spp_specific_mods_plot_ALT.R

#echo 'finished north region preview'


### Mid region paths and names ###
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/sfs_creation
#vcf_file_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/all_sites_vcf/mid_region_vcf/file_processing_3_24_2021
#vcf_file_name=MID_SUBSET_polyodon_kazumbe_allsites_2_5_2021_post_4_entropycongruent.vcf.gz
#pop_spec_name=pop_spec_MID_SUBSET_polyodon_kazumbe_allsites_2_5_2021_post_4_missing75.txt
#preview_file_name=mid_easysfs_preview_maxmissing75_monomorph_2_5_2021.txt
#output_file_name=mid_easysfs_preview_maxmissing75_monomorph_2_5_2021_topprojects.txt
#number_of_top_projections=5

#bcftools query -l $vcf_file_path/$vcf_file_name > $vcf_file_path/SAMPLES_MID_SUBSET_polyodon_kazumbe_allsites_2_5_2021_post_4_entropycongruent.txt

#awk 'NR==FNR{a[$1]; next} ($1 in a) {print $1, $2}' $vcf_file_path/SAMPLES_MID_SUBSET_polyodon_kazumbe_allsites_2_5_2021_post_4_entropycongruent.txt /project/ltcichlidgenomics/alewansk/Petrochromis$

#echo 'starting mid region preview'

##Step 1: Use the preview function in easySFS to calculate the number of segregating sites for each projection per population
#/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $vcf_file_path/$vcf_file_name -p $file_path/$pop_spec_name --preview -a -v > $file_path/$previ$

##Step 2: Identify the best projections for each population
#Rscript $script_directory/process_easySFS_preview.R -f $file_path/$preview_file_name -o $file_path/$output_file_name -n $number_of_top_projections
##R CMD BATCH spp_specific_mods_plot_ALT.R

#echo 'finished mid region preview'
