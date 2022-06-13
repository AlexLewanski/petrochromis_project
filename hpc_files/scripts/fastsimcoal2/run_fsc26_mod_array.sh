#!/bin/sh

### Job Name ###
#SBATCH --job-name run_fsc26

### Declare an account for the job to run under (check designation with arccquota on Moran) ###
#SBATCH --account=ltcichlidgenomics

#SBATCH --output=out_err_files/fsc_%A_%a.out
#SBATCH --error=out_err_files/fsc_%A_%a.err

###SBATCH -o stdout_run_fsc26_mod
###SBATCH -e stderr_run_fsc26_mod

### The directive below directs that the standard output and error streams are to be merged, ###
### intermixed, as standard output. ###
### mailing options
#SBATCH --mail-type=END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources ###
### 1 nodes, 1 processors (cores) each node ###
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds) ###
#SBATCH --time=0-02:30:00

### Array specification
#SBATCH --array=1-500%30
####SBATCH --array=1-5

### Load packages ###
module load gcc
module load htslib



###########################
### SETTING UP ANALYSES ###
###########################

#the name of the model (passed to script)
PREFIX=${1}

#where is the fsc program located?
fsc_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/fsc26_linux64

#what is the path to the script directory?  (I am currently passing this from start_fsc26_mods.sh)
#script_directory=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts
script_directory=${2}

#path where the analyses will be housed (I am currently passing this from start_fsc26_mods.sh)
#analysis_file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/prac_fastsimcoal
analysis_file_path=${3}

#how many simulations per model run?
number_of_sims=$4



#########################
### RUNNING  ANALYSES ###
#########################

#enter newly-created directory (this is where modelling results will be located)
cd $analysis_file_path/$PREFIX

mkdir run$SLURM_ARRAY_TASK_ID

cp $analysis_file_path/${PREFIX}/${PREFIX}.tpl $analysis_file_path/${PREFIX}/${PREFIX}.est $analysis_file_path/${PREFIX}/${PREFIX}_jointMAFpop1_0.obs run$SLURM_ARRAY_TASK_ID"/"

cd run$SLURM_ARRAY_TASK_ID

${fsc_path}/fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -m -C 10 -n $number_of_sims -L 40 -s0 -M 0.001 -q -c 12



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

#slightly modified for loop header to allow for a variable-defined number of iterations (i.e. number of runs):
#https://www.cyberciti.biz/faq/unix-linux-iterate-over-a-variable-range-of-numbers-in-bash/
#for i in $(eval echo "{1..$number_of_runs}")
#do
#    mkdir run$i
#    cp $analysis_file_path/${PREFIX}/${PREFIX}.tpl $analysis_file_path/${PREFIX}/${PREFIX}.est $analysis_file_path/${PREFIX}/${PREFIX}_jointMAFpop1_0.obs run$i"/"
#    cd run$i
#    ${fsc_path}/fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -m -C 10 -n $number_of_sims -L 40 -s0 -M 0.001 -q -c 12
#    cd ..
#done


#identify run with best likelihood
#${script_directory}/fsc_selectbestrun.sh

#if [ $confirm_run_num = 'YES' ]
#then
#     echo 'Number of runs: ' > run_confirmation.txt
#     find . -type d -name 'run*' | wc -l >> run_confirmation.txt
#fi


##Delete all the run directories after the best run directory is identified and extracted (these are never used again)
##helpful source: https://unix.stackexchange.com/questions/23576/how-do-i-recursively-delete-directories-with-wildcard
#if [ $remove_run_dirs = 'YES' ] && [ -d 'bestrun' ]
#then
#     find . -type d -name 'run*' -exec rm -r {} +
#fi
