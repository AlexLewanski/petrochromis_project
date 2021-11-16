#!/bin/sh

### Job Name
#SBATCH --job-name fsc_AIC

### Declare an account for the job to run under
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_fsc_AIC
#SBATCH -e stderr_fsc_AIC

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

#load libraries
module load gcc
module load r



###############################
### SETTING UP THE ANALYSES ###
###############################

#where is the calculateAIC_alt.sh script located?
script_directory=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts

#PUT THE NAME(S) OF THE MODEL(S) THAT YOU WANT SUMMARIZED IN THE ARRAY (without quotes, separated by spaces):
#mod_list=(mod1 mod2 mod3)
mod_list=( continuous_geneflow_dualpopgrowth_altspec continuous_geneflow_symmetric_dualpopgrowth_altspec continuous_geneflow_unidirect0_dualpopgrowth_altspec continuous_geneflow_unidirect1_dualpopgrowth_altspec early_geneflow_dualpopgrowth_altspec early_geneflow_symmetric_dualpopgrowth_altspec early_geneflow_unidirect0_dualpopgrowth_altspec early_geneflow_unidirect1_dualpopgrowth_altspec no_geneflow_dualpopgrowth_altspec recent_geneflow_dualpopgrowth_altspec recent_geneflow_symmetric_dualpopgrowth_altspec recent_geneflow_unidirect0_dualpopgrowth_altspec recent_geneflow_unidirect1_dualpopgrowth_altspec )


### MID REGION MODELS ###
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec

### NORTH REGION MODELS ###
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec



########################
### RUNNING ANALYSES ###
########################

#add header to file that stores AIC values
echo -e  model"\t"number_params"\t"deltaL"\t"AIC"\t"MaxObsLhood"\t"MaxEstLhood > $file_path/allmodels.AIC

#loop through each model
for mod in ${mod_list[*]}
do
	cd $file_path/$mod/bestrun/
	$script_directory/calculateAIC_alt.sh $mod
	
	echo -e $mod"\t"`tail -n 1 $mod.AIC` >> $file_path/allmodels.AIC	
done
