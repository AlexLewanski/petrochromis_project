#!/bin/bash

### Job Name
#SBATCH --job-name stair

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH -o /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/stairway_plot/out_err_files/stdout_stairyplot_%A_%a
#SBATCH -e /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/stairway_plot/out_err_files/stderr_stairyplot_%A_%a

### The directive below directs that the standard output and error streams are to be merged, \
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=1-00:00:00

#SBATCH --array=1-4



#####################
### SCRIPT SET-UP ###
#####################

###Load libraries
module load swset/2018.05
module load gcc/7.3.0
module load jdk/14.0.1
module load r

#assigning inputs to variables (for clarity)
sfs_name_file=$1
sfs_dir=$2
output_dir=$3

#identifying the SFS
focal_sfs=$(awk -v slurm_arr=${SLURM_ARRAY_TASK_ID} 'NR==slurm_arr' $sfs_name_file)
output_file_base=$(basename -s pop0.obs $focal_sfs)

#echo $output_file_base
#echo $output_dir/$output_file_base


script_dir=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/stairway_plot/scripts
stairwayplot_dir=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/stairway_plot/stairway_plot_v2.1.1

if [ -d "$output_dir/$output_file_base" ]
then
  echo "Output directory already exists. Please specify a different output directory."
  exit 1
else
  mkdir $output_dir/$output_file_base/
fi

#cp -r $stairwayplot_dir/stairway_plot_es/ $output_dir/$output_file_base/



##########################################
### CREATING STAIRWAY PLOT INPUT FILES ###
##########################################

### parse sfs ###

Rscript $script_dir/parse_sfs.R $sfs_dir/$focal_sfs $output_dir/$output_file_base/${output_file_base}_sfs_info.txt

input_indiv=$(grep "^number_individuals" $output_dir/$output_file_base/${output_file_base}_sfs_info.txt | sed 's/number_individuals: //')
input_sites=$(grep "^number_sites: " $output_dir/$output_file_base/${output_file_base}_sfs_info.txt | sed 's/number_sites: //')
input_sfs=$(grep "^processed_sfs: " $output_dir/$output_file_base/${output_file_base}_sfs_info.txt | sed 's/processed_sfs: //')
nrand=$(grep "^nrand: " $output_dir/$output_file_base/${output_file_base}_sfs_info.txt | sed 's/nrand: //')


cd $output_dir/$output_file_base/


### create blueprint file ###
sh $script_dir/create_blueprint.sh -n $output_file_base \
-q "$input_indiv" \
-l "$input_sites" \
-s "$input_sfs" \
-f \
-m 1 \
-p 0.67 \
-r "$nrand" \
-i $output_file_base \
-k $stairwayplot_dir/stairway_plot_es \
-t 1000 \
-u "3.5e-9" \
-g 2 \
-h $output_file_base \
-x "0,0" \
-y "0,0" \
-z 2 \
-a 2 \
-o 12



############################
### LAUNCH STAIRWAY PLOT ###
############################

#create batch file with Stairbuilder.class
java -cp $stairwayplot_dir/stairway_plot_es Stairbuilder ${output_file_base}.blueprint
bash ${output_file_base}.blueprint.sh
