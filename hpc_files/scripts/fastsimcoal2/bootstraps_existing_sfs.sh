#!/bin/sh

### Job Name
#SBATCH --job-name fsc_boot

###account to run the job under
#SBATCH --account=ltcichlidgenomics

####SBATCH -o stdout_fsc_bootstrap
####SBATCH -e stderr_fsc_bootstrap

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=1-00:00:00
#module load gcc
#module load htslib



#######################
### 1. SCRIPT SETUP ###
#######################

### 1a. RELEVANT PATHS ###
#file_path --> path to directory where bootstraps will be housed

#path to the tpl and est files
est_tpl_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/template_inputs

#sfs_dir --> path(s) to the bootstrap sfs


### 1b. NAMING DETAILS ###
#what do you want to name the bootstrap results folder?
storage_folder="bootstrap_recent_geneflow_dualpopgrowth_altspec_V2"

#***THIS NEEDS TO MATCH THE BEST MODEL NAME
PREFIX="recent_geneflow_dualpopgrowth_altspec"


### 1c. ADDITIONAL INFO ###
#number of runs per bootstrap dataset
number_of_runs=500

#number of simulations per dataset
number_of_sims=150000

#number of bootstrap datasets
bootstraps=100

change_sample_sizes=YES

### mid dataset setup ###
#sample_size_array=( 16 16 )
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec
#sfs_dir=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/mid_bs_creation_popgrow_altspec/sfs_creation

#north
#sample_size_array=( 50 50 )
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec
#sfs_dir=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/north_bs_creation_popgrow_altspec/sfs_creation


#####################
### 2. BOOTSTRAPS ###
#####################

### 2a. CREATE TOP LEVEL BOOTSTRAP DIRECTORY ###
if [ ! -d  $file_path/$storage_folder ]
then
	mkdir $file_path/$storage_folder
else
	echo "The storage folder for the bootstraps already exists. Please specify a new location."
	exit 1
fi


### 2b. PREPARE FASTSIMCOAL2 INPUT FILES ##

#enter top level bootstrap folder
cd $file_path/$storage_folder

#copy the .est and .tpl files into the main bootstrap folder
cp $est_tpl_path/${PREFIX}.est $est_tpl_path/${PREFIX}.tpl $file_path/$storage_folder/

#optionally change sample sizes in tpl file
if [ $change_sample_sizes = 'YES' ]
then
    for i in "${!sample_size_array[@]}"
    do
        line_num=$(( $i + 7 )) #first sample size is on line 7 in the tpl file (bash loops start indexing at 0)

        #a more portable version of sed in-place editing that should work on sed versions available on macOS and Linux (GNU for many linux distributions)
        #works by creating a temporary file (${mod}.tpl.temp) and then removes it after file is edited
        #see: https://stackoverflow.com/questions/16745988/sed-command-with-i-option-in-place-editing-works-fine-on-ubuntu-but-not-mac
        #info on && (from Stack Overflow): && lets you do something based on whether the previous command completed successfully (1/2)
        #that's why you tend to see it chained as do_something && do_something_else_that_depended_on_something (2/2)
        sed -i.temp "${line_num}s/.*/${sample_size_array[$i]}/" $file_path/$storage_folder/${PREFIX}.tpl && rm $file_path/$storage_folder/${PREFIX}.tpl.temp
    done
fi


### 2c. RUN FASTSIMCOAL2 ON EACH BOOTSTRAP DATASET ###
#submit separate jobs for each bootstrap
for bs in $(eval echo "{1..$bootstraps}")
do
  #make bootstrap directory
  mkdir bs$bs
  
  #copy sfs into the bootstrap folder
  cp ${sfs_dir}/bs_sfs${bs}/fastsimcoal2/${PREFIX}_jointMAFpop1_0.obs  $file_path/$storage_folder/bs$bs/${PREFIX}.bs.${bs}_jointMAFpop1_0.obs
  
  sbatch /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts/run_single_bootstrap.sh $bs $file_path $storage_folder $number_of_runs $PREFIX $number_of_sims
  sleep 1
done
