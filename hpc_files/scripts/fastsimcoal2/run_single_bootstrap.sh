#!/bin/sh

### Job Name
#SBATCH --job-name single_bs

### Declare an account for the job to run under - (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH -o /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/prac_fastsimcoal/out_err_files/stdout_fsc_bootstrap
#SBATCH -e /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/prac_fastsimcoal/out_err_files/stderr_fsc_bootstrap

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output.
### mailing options
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=6-00:00:00
module load gcc
module load htslib



#####################
### SCRIPT SET-UP ###
#####################

#assigning inputs to informative variable names
bs_num=$1
file_path=$2
storage_folder=$3
number_of_runs=$4
PREFIX=$5
number_of_sims=$6

#should a file be created that confirms the number of runs? (just to double check that the correct number of runs were completed)
confirm_run_num=YES #YES or NO (anything but YES will result in the number of runs not being recorded)

#should the run directories be deleted after the best run is identified and extracted?
remove_run_dirs=YES #YES or NO (anything but YES will result in no files being deleted)

fsc_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/fsc26_linux64
script_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts



#########################
### RUNNING BOOTSTRAP ###
#########################

#enter $bs bootstrap folder
cd $file_path/$storage_folder/bs$bs_num
  
#complete $number_of_runs runs for the bootstrap
for i in $(eval echo "{1..$number_of_runs}")
do
  mkdir run$i #create folder for run i
  #copy and rename tpl and est files into the run folder
  cp $file_path/$storage_folder/${PREFIX}.tpl run$i/${PREFIX}.bs.${bs_num}.tpl
  cp $file_path/$storage_folder/${PREFIX}.est run$i/${PREFIX}.bs.${bs_num}.est
  cp $file_path/$storage_folder/bs${bs_num}/${PREFIX}.bs.${bs_num}_jointMAFpop1_0.obs run$i"/"

  cd run$i #enter the run folder
        
  #run faststimcoal
  ${fsc_path}/fsc26 -t ${PREFIX}.bs.${bs_num}.tpl -e ${PREFIX}.bs.${bs_num}.est -m -C 10 -n $number_of_sims -L 40 -s0 -M 0.001 -q -c 12
  #/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/fsc26_linux64/fsc26 -t ${PREFIX}.bs.${bs_num}.tpl -e ${PREFIX}.bs.${bs_num}.est -m -0 -C 10 -n $number_of_sims -L 50 -s0 -M -q -c 16
  cd .. #leave the run folder (go up 1 level in the hierarchy)
        
done
  
### identify best run ###
#/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts/bootstrap_selectbestrun.sh
${script_path}/bootstrap_selectbestrun.sh


#Delete all the run directories after the best run directory is identified and extracted (these are never used again)
#helpful source: https://unix.stackexchange.com/questions/23576/how-do-i-recursively-delete-directories-with-wildcard
#find . -type d -name 'run*' -exec rm -r {} +


if [ $confirm_run_num = 'YES' ]
then
     echo 'Number of runs: ' > run_confirmation.txt
     find . -type d -name 'run*' | wc -l >> run_confirmation.txt
fi


#Delete all the run directories after the best run directory is identified and extracted (these are never used again).
#Only delete run directories (if $remove_run_dirs = 'YES') if the bestrun was successfully generated. This second criterion 
#was suggested by Jessi, and I think it is a good idea because it would be good to not delete the runs if an issue occurs,
#which results in the bestrun folder not being created.

#helpful sources:
#helpful source: https://unix.stackexchange.com/questions/23576/how-do-i-recursively-delete-directories-with-wildcard
#https://ryanstutorials.net/bash-scripting-tutorial/bash-if-statements.php
if [ $remove_run_dirs = 'YES' ] && [ -d 'bestrun' ]
then
     find . -type d -name 'run*' -exec rm -r {} +
fi
