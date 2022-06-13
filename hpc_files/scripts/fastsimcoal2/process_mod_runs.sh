#!/bin/sh

### Job Name ###
#SBATCH --job-name run_fsc26

### Declare an account for the job to run under (check designation with arccquota on Moran) ###
#SBATCH --account=ltcichlidgenomics


#SBATCH -o stdout_run_fsc26_mod
#SBATCH -e stderr_run_fsc26_mod

### The directive below directs that the standard output and error streams are to be merged, ###
### intermixed, as standard output. ###
### mailing options
#SBATCH --mail-type=END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources ###
### 1 nodes, 1 processors (cores) each node ###
#SBATCH --nodes=1
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds) ###
#SBATCH --time=0-00:30:00


### Load packages ###
module load gcc
#module load htslib



###########################
### SETTING UP ANALYSES ###
###########################

#analysis_directory_array=( /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec \
#/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec )

#analysis_directory_array=( /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec )
analysis_directory_array=( /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec )


script_directory=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts
#analysis_directory=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_C10_TWO_PHASE_GF_POPGROW_TEST_array
#internal_model_array=( two_cycle_gf_init_isolation two_cycle_gf_init_isolation_symmetric two_cycle_gf_init_isolation_unidirect0 two_cycle_gf_init_isolation_unidirect1 )

#internal_model_array=( early_geneflow_dualpopgrowth_altspec_updated \
#early_geneflow_symmetric_dualpopgrowth_altspec_updated \
#early_geneflow_unidirect0_dualpopgrowth_altspec_updated \
#early_geneflow_unidirect1_dualpopgrowth_altspec_updated \
#recent_geneflow_dualpopgrowth_altspec_updated \
#recent_geneflow_symmetric_dualpopgrowth_altspec_updated \
#recent_geneflow_unidirect0_dualpopgrowth_altspec_updated \
#recent_geneflow_unidirect1_dualpopgrowth_altspec_updated )

internal_model_array=( continuous_geneflow_bottleneck0_grow1_aftermig \
continuous_geneflow_symmetric_bottleneck0_grow1_aftermig \
continuous_geneflow_unidirect0_bottleneck0_grow1_aftermig \
continuous_geneflow_unidirect1_bottleneck0_grow1_aftermig \
early_geneflow_bottleneck0_grow1_aftermig \
early_geneflow_symmetric_bottleneck0_grow1_aftermig \
early_geneflow_unidirect0_bottleneck0_grow1_aftermig \
early_geneflow_unidirect1_bottleneck0_grow1_aftermig \
no_geneflow_bottleneck0_grow1_aftermig \
recent_geneflow_bottleneck0_grow1_aftermig \
recent_geneflow_symmetric_bottleneck0_grow1_aftermig \
recent_geneflow_unidirect0_bottleneck0_grow1_aftermig \
recent_geneflow_unidirect1_bottleneck0_grow1_aftermig \
early_geneflow_dualpopgrowth_altspec_updated \
early_geneflow_symmetric_dualpopgrowth_altspec_updated \
early_geneflow_unidirect0_dualpopgrowth_altspec_updated \
early_geneflow_unidirect1_dualpopgrowth_altspec_updated \
recent_geneflow_dualpopgrowth_altspec_updated \
recent_geneflow_symmetric_dualpopgrowth_altspec_updated \
recent_geneflow_unidirect0_dualpopgrowth_altspec_updated \
recent_geneflow_unidirect1_dualpopgrowth_altspec_updated )


confirm_run_num=YES
remove_run_dirs=YES


#########################
### RUNNING  ANALYSES ###
#########################

for file_path in "${analysis_directory_array[@]}"
do

  for mod in "${internal_model_array[@]}"
  do

    cd $file_path/$mod
    
    #identify run with best likelihood
    ${script_directory}/fsc_selectbestrun.sh

    if [ $confirm_run_num = 'YES' ]
    then
      echo 'Number of runs: ' > run_confirmation.txt
      find . -type d -name 'run*' | wc -l >> run_confirmation.txt
    fi

    ##Delete all the run directories after the best run directory is identified and extracted (these are never used again)
    ##helpful source: https://unix.stackexchange.com/questions/23576/how-do-i-recursively-delete-directories-with-wildcard
    if [ $remove_run_dirs = 'YES' ] && [ -d 'bestrun' ]
    then
      find . -type d -name 'run*' -exec rm -r {} +
    fi

  done

done
