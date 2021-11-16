#!/bin/sh

### Job Name
#SBATCH --job-name fsc_boot

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics


### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=NONE
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=0-01:00:00



#####################
### SCRIPT SET-UP ###
#####################

### PATHS ###
#put the path(s) model directories to the in the this array
file_path_array=(/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec \
/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec)

#what is the name of the directory that contains the bootstrap results?
storage_folder="bootstrap_recent_geneflow_dualpopgrowth_altspec"

PREFIX="recent_geneflow_dualpopgrowth_altspec"



##########################################################################################
### COMPILE THE PARAMETER VALUES OF THE BEST RUN FOR EACH BOOTSTRAP INTO A SINGLE FILE ###
##########################################################################################

for i in "${file_path_array[@]}"
do

  #move to the bootstrap folder
  cd $i/$storage_folder

  #for bs in $(eval echo "{1..$bootstraps}")

  #iterate through each bootstrap folder
  for bootstrap in bs*
  do

    #extract the number from the bootstrap folder name (e.g. bs13 --> 13)
    bs=`echo $bootstrap | grep -o '[0-9]\+'`

    #use the .bestlhoods file from bs1 as the start of the bs_bestrun_param.txt file (in order to include the parameter names in the file)
    if [ $bs == 1 ]
    then
      cp bs$bs/bestrun/${PREFIX}.bs.$bs.bestlhoods bs_bestrun_param.txt
   
    #add the parameter values of the rest of the bs bestruns to bs_bestrun_param.txt
    else
      tail -1 bs$bs/bestrun/${PREFIX}.bs.$bs.bestlhoods >> bs_bestrun_param.txt
    fi
  done
done



#################################
### CODE CURRENTLY NOT IN USE ###
#################################

### NORTH REGION ###
#path to directory where the model directories are located
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec

### MID REGION ###
#path to directory where the model directories are located
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_maxmiss75_monomorph_3_25_2021_C10
