#!/bin/sh

### SLURM DIRECTIVES ###
#SBATCH --job-name process_vcf
#SBATCH --account=wagnerlab
#SBATCH -o out_process_vcf
#SBATCH -e err_process_vcf
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=alewansk@uwyo.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0



#####################
### SCRIPT SET-UP ###
#####################

#load modules
module load gcc
module load vcftools
module load bcftools
module load r/3.5.3


### file paths and files ###
home_dir_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/final_variant_files/dataset_processing_6_13_2021
biallelic_dir=biallelic_datasets
biallelic_path=$home_dir_path/$biallelic_dir

allsites_toplevel_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics
allsites_dir=allsites_datasets_6_16_2021
allsites_path=$allsites_toplevel_path/$allsites_dir

processing_script_file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/final_variant_files/dataset_processing_6_13_2021/processing_scripts
cichlid_metadata=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/entropy/monster_23jul20.csv

#biallellic: polyodon, kazumbe (without green, diagramma)
vcf_biallelic_polykaz_unprocessed=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/variants_pundamilia/variants_processing_kazumbe_polyodon/kazumbe_polyodon_pundamilia_8_19_2020_post_1.vcf

# biallellic: polyodon, kazumbe, green, diagramma
vcf_biallelic_polykazout_unprocessed=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/variants_pundamilia/variants_processing_kazumbe_polyodon_green_diagramma/kazumbe_polyodon_green_diagramma_pundamilia_8_13_2020_post_1.vcf

#allsites: polyodon, kazumbe  (without green, diagramma)
vcf_allsites_polykaz_unprocessed=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/all_sites_vcf/polyodon_kazumbe_allsites_2_5_2021_post_1.vcf.gz



#path/directory checks --> making sure that home_dir_path exists but biallelic_dir and allsites_path do not
if [ ! -d $home_dir_path ]
then
  echo "make sure that home_dir_path exists"
  exit 1
fi

if [ -d $biallelic_path ]
then
  echo "make sure that biallelic_dir does not already exist in your home_dir_path directory"
  exit 1
else
  mkdir $biallelic_path
fi

if [ ! -d $allsites_toplevel_path ]
then
  echo "make sure that allsites_path exists"
  exit 1
fi

if [ -d $allsites_path ]
then
  echo "make sure that allsites_path does not already exist in your home_dir_path directory"
  exit 1
else
  mkdir $allsites_path
fi


#script checks --> making sure they exist in the location specified
if [ ! -f $processing_script_file_path/process_vcf.sh ]
then
  echo "make sure that process_vcf.sh is in the directory you specified in processing_script_file_path"
  exit 1
fi

if [ ! -f $processing_script_file_path/subset_region_samples_update.R ]
then
  echo "make sure that subset_region_samples_update.R is in the directory you specified in processing_script_file_path"
  exit 1
fi


#setting up directories and initial filter (removing CEW07F_120)
if [ ! -d  $home_dir_path/sample_files/ ]
then
  mkdir $home_dir_path/sample_files/
fi
#create file to remove CEW07F_120 (creatively named 'remove_CEW07F_120.txt')
echo "CEW07F_120" > $home_dir_path/sample_files/remove_CEW07F_120.txt


#filters consistent across all datasets
min_mean_dp=5
max_mean_dp=75
min_q=20
site_miss=0.70
indiv_miss=0.70



#####################################
### PROCESSING BIALLELIC DATASETS ###
#####################################

cd $biallelic_path

### DATASET 1: Biallelic: polyodon + kazumbe ###
#notes: remove CEW07F_120 before processing
echo "### Biallelic: polyodon + kazumbe ###"

sh $processing_script_file_path/process_vcf.sh -v $vcf_biallelic_polykaz_unprocessed \
-p biallel_polykaz \
-s ^$home_dir_path/sample_files/remove_CEW07F_120.txt \
-q $min_q -i 2 -a 2 -f 0.01 -d $min_mean_dp -b $max_mean_dp -o 4 -g 100 \
-h chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
-x $site_miss -l $indiv_miss \
-n biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss} \
-r -m



### SAMPLES: Processing region-specific sample lists ###
#get samples from the processed file
bcftools query -l biallel_polykaz/biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.recode.vcf > $home_dir_path/sample_files/SAMPLES_biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.txt
filter_sample_reference=$home_dir_path/sample_files/SAMPLES_biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.txt

#north region
Rscript $processing_script_file_path/subset_region_samples_update.R --region north \
-m $cichlid_metadata \
-f $home_dir_path/sample_files/SAMPLES_NORTH_biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.txt \
-v TRUE \
-w $filter_sample_reference \

#mid region
Rscript $processing_script_file_path/subset_region_samples_update.R --region mid \
-m $cichlid_metadata \
-f $home_dir_path/sample_files/SAMPLES_MID_biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.txt \
-v TRUE \
-w $filter_sample_reference \

north_samples_processed=$home_dir_path/sample_files/SAMPLES_NORTH_biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.txt
mid_samples_processed=$home_dir_path/sample_files/SAMPLES_MID_biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.txt



### DATASET 2: NORTH Biallelic: polyodon + kazumbe ###
echo "### NORTH Biallelic: polyodon + kazumbe ###"

mkdir $biallelic_path/north_biallel_polykaz/
bcftools view \
-S $north_samples_processed $biallelic_path/biallel_polykaz/biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.recode.vcf \
-o $biallelic_path/north_biallel_polykaz/north_biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.vcf \



### DATASET 3: MID Biallelic: polyodon + kazumbe ###
echo "### MID Biallelic: polyodon + kazumbe ###"

mkdir $biallelic_path/mid_biallel_polykaz/
bcftools view \
-S $mid_samples_processed $biallelic_path/biallel_polykaz/biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.recode.vcf \
-o $biallelic_path/mid_biallel_polykaz/mid_biallel_polykaz_minq${min_q}_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss}.vcf \



### DATASET 4: Biallelic: polyodon + kazumbe + outgroups ###
echo "### Biallelic: polyodon + kazumbe + outgroups ###"

#creating sample list for subsetting the vcf
#sample list includes polyodon and kazumbe in the processed polyodon/kazumbe vcf + S. diagramma and P. green
bcftools query -l $vcf_biallelic_polykazout_unprocessed > $home_dir_path/sample_files/SAMPLES_polykazoutgroups_unprocessed.txt
awk -F, '$13 == "Simochromis diagramma" || $13 == "Petrochromis green" {print $3}' $cichlid_metadata > $home_dir_path/sample_files/metadata_sdiagramma_pgreen_samples.txt
grep -x -f $home_dir_path/sample_files/metadata_sdiagramma_pgreen_samples.txt $home_dir_path/sample_files/SAMPLES_polykazoutgroups_unprocessed.txt > $home_dir_path/sample_files/sdiagramma_pgreen_invcf.txt
cat $filter_sample_reference $home_dir_path/sample_files/sdiagramma_pgreen_invcf.txt > $home_dir_path/sample_files/biallel_polykazoutgroups_sample_subset.txt


sh $processing_script_file_path/process_vcf.sh -v $vcf_biallelic_polykazout_unprocessed \
-p biallel_polykazoutgroups \
-s $home_dir_path/sample_files/biallel_polykazoutgroups_sample_subset.txt \
-q $min_q -i 2 -a 2 -f 0.01 -d $min_mean_dp -b $max_mean_dp -o 4 -g 100 \
-h chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
-x $site_miss -l $indiv_miss \
-n biallel_polykazoutgroups_minq$min_q_maf1_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss} \
-r -m



####################################
### PROCESSING ALLSITES DATASETS ###
####################################

cd $allsites_path

### DATASET 5: Allsites: polyodon + kazumbe ###
echo "### Allsites: polyodon + kazumbe ###"
sh $processing_script_file_path/process_vcf.sh -v $vcf_allsites_polykaz_unprocessed \
-p allsites_polykaz \
-s $filter_sample_reference \
-q $min_q -i 0 -a 2 -f 0 -d $min_mean_dp -b $max_mean_dp -o 4 -g 100 \
-h chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
-x $site_miss -l $indiv_miss \
-n allsites_polykaz_minq${min_q}_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss} \
-r -z -m



### DATASET 6: NORTH Allsites: polyodon + kazumbe ###
echo "### NORTH Allsites: polyodon + kazumbe ###"
sh $processing_script_file_path/process_vcf.sh -v $vcf_allsites_polykaz_unprocessed \
-p north_allsites_polykaz \
-s $north_samples_processed \
-q $min_q -i 0 -a 2 -f 0 -d $min_mean_dp -b $max_mean_dp -o 4 -g 100 \
-h chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
-x $site_miss -l $indiv_miss \
-n north_allsites_polykaz_minq${min_q}_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss} \
-r -z -m



### DATASET 7: MID Allsites: polyodon + kazumbe ###
echo "### MID Allsites: polyodon + kazumbe ###"
sh $processing_script_file_path/process_vcf.sh -v $vcf_allsites_polykaz_unprocessed \
-p mid_allsites_polykaz \
-s $mid_samples_processed \
-q $min_q -i 0 -a 2 -f 0 -d $min_mean_dp -b $max_mean_dp -o 4 -g 100 \
-h chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
-x $site_miss -l $indiv_miss \
-n mid_allsites_polykaz_minq${min_q}_minDP${min_mean_dp}_maxDP${max_mean_dp}_miss${site_miss} \
-r -z -m
