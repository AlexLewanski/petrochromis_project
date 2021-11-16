#!/bin/sh

### Job Name
#SBATCH --job-name fsc_sim

### Declare an account for the job to run under (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_fsc_sim
#SBATCH -e stderr_fsc_sim

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=4-00:00:00



###############################
### SETTING UP THE ANALYSES ###
###############################

#loading libraries
module load gcc
module load htslib

#where is the fsc program located?
fsc_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/fsc26_linux64

#PUT THE NAME(S) OF THE MODEL(S) THAT YOU WANT SUMMARIZED IN THE ARRAY (without quotes, separated by spaces):
#currently implemented models: mod_list=(geneflow  early_geneflow recent_geneflow no_geneflow)
#mod_list=(geneflow_alt no_geneflow)


#MID REGION MODELS
#mod_list=( continuous_geneflow_dualpopgrowth_altspec continuous_geneflow_symmetric_dualpopgrowth_altspec continuous_geneflow_unidirect0_dualpopgrowth_altspec continuous_geneflow_unidirect1_dualpopgrowth_altspec early_geneflow_dualpopgrowth_altspec early_geneflow_symmetric_dualpopgrowth_altspec early_geneflow_unidirect0_dualpopgrowth_altspec early_geneflow_unidirect1_dualpopgrowth_altspec no_geneflow_dualpopgrowth_altspec recent_geneflow_dualpopgrowth_altspec recent_geneflow_symmetric_dualpopgrowth_altspec recent_geneflow_unidirect0_dualpopgrowth_altspec recent_geneflow_unidirect1_dualpopgrowth_altspec )
#file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec
#storage_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec


#NORTH REGION MODELS
mod_list=( continuous_geneflow_dualpopgrowth_altspec continuous_geneflow_symmetric_dualpopgrowth_altspec continuous_geneflow_unidirect0_dualpopgrowth_altspec continuous_geneflow_unidirect1_dualpopgrowth_altspec early_geneflow_dualpopgrowth_altspec early_geneflow_symmetric_dualpopgrowth_altspec early_geneflow_unidirect0_dualpopgrowth_altspec early_geneflow_unidirect1_dualpopgrowth_altspec no_geneflow_dualpopgrowth_altspec recent_geneflow_dualpopgrowth_altspec recent_geneflow_symmetric_dualpopgrowth_altspec recent_geneflow_unidirect0_dualpopgrowth_altspec recent_geneflow_unidirect1_dualpopgrowth_altspec )
file_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec
storage_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow_altspec


#number of rounds of simulations
#100 is a good number of iterations
number_of_iters=100
likelihood_results_dir=mods_recompute_lhoods

collect_sfs=yes
retain_all_sfs=yes



#########################
### RUNNING  ANALYSES ###
#########################

#make directory to store the likelihood files
mkdir $storage_path/$likelihood_results_dir
mkdir $storage_path/$likelihood_results_dir/likelihood_dir

if [ $collect_sfs = "yes" ]
then
  mkdir $storage_path/$likelihood_results_dir/sfs_dir
fi

#loop through each model
for mod in ${mod_list[*]}
do
	#if [ $collect_sfs = "yes" ]
        #then
        #  mkdir $storage_path/$likelihood_results_dir/sfs_dir/${mod}_sfs
        #fi


	cd $file_path/$mod/bestrun
	
	# create temporary obs file with name _maxL_MSFS.obs
	cp ${mod}_jointMAFpop1_0.obs ${mod}_maxL_jointMAFpop1_0.obs
	

	# Run fastsimcoal $number_of_iters times to get the likelihood of the observed SFS under the best parameter values with 1 million simulated SFS.
	for iter in $(eval echo "{1..$number_of_iters}")
	do
		#/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/fsc26_linux64/fsc26 -i ${mod}_maxL.par -n1000000 -m -q
		$fsc_path/fsc26 -i ${mod}_maxL.par -n 1000000 -C 10 -m -q -c 12
		#removed -0 argument (2/13/2021)
		# Fastsimcoal will generate a new folder called ${model}_maxL and write files in there
		
		# collect the lhood values (Note that >> appends to the file, whereas > would overwrite it)
		#sed -n '2,3p' ${mod}_maxL/${mod}_maxL.lhoods >> ${mod}.lhoods #sed version of extracting likelihood value
                awk 'NR == 2 {print $1}' ${mod}_maxL/${mod}_maxL.lhoods >> ${mod}.lhoods
		
		if [ $collect_sfs = "yes" ]
		then
		  cp ${mod}_maxL/${mod}_maxL_jointMAFpop1_0.txt $storage_path/$likelihood_results_dir/sfs_dir/${mod}_maxL_jointMAFpop1_0_iter${iter}.txt  
		fi
		
		# delete the folder with results
		rm -r ${mod}_maxL/
	done
	
	#move likelihood file to the final storage directory
	mv ${mod}.lhoods $storage_path/$likelihood_results_dir/likelihood_dir/

	if [ $collect_sfs = "yes" ]
        then
	  start_val=`awk 'NR == 1 {print $1}' $storage_path/$likelihood_results_dir/likelihood_dir/${mod}.lhoods`
	  top_iter=`awk -v first_row=$start_val 'BEGIN{val=first_row; row_num=1} {if ($1>val+0) {val=$1; row_num=FNR}} END{print row_num}' $storage_path/$likelihood_results_dir/likelihood_dir/${mod}.lhoods`

	  cd $storage_path/$likelihood_results_dir/sfs_dir/
	  
	  if [ $retain_all_sfs != "yes" ]
	  then
            find . -name "${mod}_maxL_jointMAFpop1_0_iter*.txt" -type f | grep -v "${mod}_maxL_jointMAFpop1_0_iter${top_iter}.txt" | xargs rm
	  else
	    mkdir ${mod}_sfs
	    mv ${mod}_maxL_jointMAFpop1_0_iter${top_iter}.txt ${mod}_maxL_jointMAFpop1_0_iter${top_iter}_toplikelihood.txt

	    #https://unix.stackexchange.com/questions/154818/how-to-integrate-mv-command-after-find-command
	    find . -name "${mod}_maxL_jointMAFpop1_0_iter*.txt" -type f -exec mv -t ${mod}_sfs/ {} +
	  fi

	fi
done
