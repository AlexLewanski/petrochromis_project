#!/bin/sh

### Job Name
#SBATCH --job-name start_pixy

### Declare an account for the job to run under
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_start_pixy
#SBATCH -e stderr_start_pixy

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=0-01:00:00



#####################
### SCRIPT SET-UP ###
#####################

file_name_array=( population_subset_7_1_2021_SUBSAMP_1 population_subset_7_1_2021_SUBSAMP_2 population_subset_7_1_2021_SUBSAMP_3 population_subset_7_1_2021_SUBSAMP_4 population_subset_7_1_2021_SUBSAMP_5 )
vcf=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/allsites_datasets_6_16_2021/allsites_polykaz/allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.vcf.gz
output_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/diversity_calculations/results/pop_subset_kazumbe
pop_file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/diversity_calculations/pop_files/pop_subset_kazumbe_subset_popfile_dir
zarr_path=/gscratch/alewansk/pixy_zarr



###########################
### LAUNCHING PIXY RUNS ###
###########################

[ ! -d $output_path ] && mkdir $output_path

for FILE_NAME in "${file_name_array[@]}"
do
  if [[ ! -d "$zarr_path/zarr_$FILE_NAME" ]]
  then
    mkdir $zarr_path/zarr_$FILE_NAME
  else
    echo 'zarr directory already exists. Using existing directory.'
  fi

  sbatch pixy_launcher.sh $vcf $zarr_path/zarr_$FILE_NAME $pop_file_path/$FILE_NAME.txt $output_path/output_$FILE_NAME

  echo 'launched the following dataset: ' $FILE_NAME
  sleep 1
done
