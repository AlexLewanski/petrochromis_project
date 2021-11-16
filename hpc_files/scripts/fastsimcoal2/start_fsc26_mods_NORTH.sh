#!/bin/sh

### Job Name
#SBATCH --job-name start_fsc

### Declare an account for the job to run under (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_start_fsc26_mods
#SBATCH -e stderr_start_fsc26_mods

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0


### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=0-03:00:00

#loading modules
#module load gcc
#module load htslib



###############################
### SETTING UP THE ANALYSES ###
###############################

### Paths to directories and/or files ###
#what is the path to the script directory?
script_directory=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts

analysis_directory=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/two_pop_mods_maxmiss75_monomorph_7_7_2021_C10_popgrow

#sfs file with full file path
sfs_file=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/sfs_creation/easy_sfs_creation_maxmiss70_monomorph_7_7_2021/fastsimcoal2/north_allsites_polykaz_minq20_minDP5_maxDP75_miss0_jointMAFpop1_0.obs


#path to directory where tpl and est files are stored
tpl_est_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/template_inputs

#do you want to change the sample sizes in the tpl file?
#if YES, specify change_sample_sizes=YES and specify sample sizes as an array: sample_size_array=( pop0 pop1 pop2 pop3 ) 
change_sample_sizes=YES
sample_size_array=( 50 50 )


### Processing the array (containing the model names) used for launching the models ###
#if model names are not passed to script, then it will use the models named in internal_model_array

#mod release for monomorphic models (2/10/2021)
#internal_model_array=( no_geneflow geneflow_alt early_geneflow_alt recent_geneflow_alt change_geneflow_alt early_geneflow_symmetric recent_geneflow_symmetric geneflow_symmetric change_geneflow_symmetric early_geneflow_alt_unidirect0 early_geneflow_alt_unidirect1 geneflow_alt_unidirect0 geneflow_alt_unidirect1 recent_geneflow_alt_unidirect0 recent_geneflow_alt_unidirect1 change_geneflow_alt_unidirect0 change_geneflow_alt_unidirect1 )
internal_model_array=( continuous_geneflow_dualpopgrowth continuous_geneflow_symmetric_dualpopgrowth continuous_geneflow_unidirect0_dualpopgrowth continuous_geneflow_unidirect1_dualpopgrowth early_geneflow_dualpopgrowth early_geneflow_symmetric_dualpopgrowth early_geneflow_unidirect0_dualpopgrowth early_geneflow_unidirect1_dualpopgrowth no_geneflow_dualpopgrowth recent_geneflow_dualpopgrowth recent_geneflow_symmetric_dualpopgrowth recent_geneflow_unidirect0_dualpopgrowth recent_geneflow_unidirect1_dualpopgrowth )

if [ $# -eq 0 ]
then
    if [ ${#internal_model_array[@]} -ne 0 ]
    then
        final_model_array=( "${internal_model_array[@]}" )

    else
        echo "If you do not feed the model names to the script, you need to provide them within the script in the following array: internal_model_array"
        exit 1
    fi

else
    final_model_array=( "$@" )
fi



#######################
### RUNNING  MODELS ###
#######################

if [ ! -d $analysis_directory ]
then
        mkdir $analysis_directory
else 
        echo "The analysis directory already exists. Please specify a new analysis directory that doesn't exist."
        exit 1
fi

for mod in "${final_model_array[@]}"
do
    #copy and rename sfs file into the top-level analysis results directory (run_fsc26_mod.sh will move these into the specific modelling directory)
    cp ${sfs_file} ${analysis_directory}/${mod}_jointMAFpop1_0.obs

    #copy tpl and est files into the top-level analysis results directory
    cp ${tpl_est_path}/${mod}.est ${tpl_est_path}/${mod}.tpl  ${analysis_directory}/


    #optionally change sample sizes in tpl file
    if [ $change_sample_sizes = 'YES' ]
    then
        for i in "${!sample_size_array[@]}"
        do 
            line_num=$(( $i + 7 )) #first sample size is on line 7 in the tpl file (bash loops start indexing at 0)        
            sed -i.temp "${line_num}s/.*/${sample_size_array[$i]}/" ${analysis_directory}/${mod}.tpl && rm ${analysis_directory}/${mod}.tpl.temp
        done    
     fi


    #start the model runs
    sbatch ${script_directory}/run_fsc26_mod.sh ${mod} ${script_directory} ${analysis_directory}

    echo "started ${mod}"
    sleep 1
done



##############################################
### NOTES AND RESOURCES FOR WRITING SCRIPT ###
##############################################

#helpful sources for creating if statement above (mostly to correctly specify conditional statements):
#https://www.cyberciti.biz/faq/finding-bash-shell-array-length-elements/
#https://stackoverflow.com/questions/1378274/in-a-bash-script-how-can-i-exit-the-entire-script-if-a-certain-condition-occurs
#https://stackoverflow.com/questions/6482377/check-existence-of-input-argument-in-a-bash-shell-script
#https://opensource.com/article/18/5/you-dont-know-bash-intro-bash-arrays #info on working with arrays

### notes on creating portable in-place sed editing ###
#a more portable version of sed in-place editing that should work on sed versions available on macOS and Linux (GNU for many linux distributions)
#works by creating a temporary file (${mod}.tpl.temp) and then removes it after file is edited
#see: https://stackoverflow.com/questions/16745988/sed-command-with-i-option-in-place-editing-works-fine-on-ubuntu-but-not-mac
#info on && (from Stack Overflow): && lets you do something based on whether the previous command completed successfully (1/2)
#that's why you tend to see it chained as do_something && do_something_else_that_depended_on_something (2/2)
