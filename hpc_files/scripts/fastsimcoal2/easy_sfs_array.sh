#!/bin/sh

### Job Name
#SBATCH --job-name run_easySFS

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH --output=out_err_files/bs_sfs_%A_%a.out
#SBATCH --error=out_err_files/bs_sfs_%A_%a.err

### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --mem=130G

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=0-20:00:00

### Array specification
#SBATCH --array=1-100%25


#####################
### SCRIPT SET-UP ###
#####################

### DIRECTIONS ###
#To create each region's SFS, uncomment the lines below each region's header in SCRIPT SET-UP section and the RUNNING SCRIPT section

#load modules and conda environment
module load gcc
module load miniconda3
source activate easySFS

### OPTIONS ###
delete_bs_dataset="yes" #yes or no

### MID ###
#bs_dataset_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/mid_bs_creation_popgrow_altspec/bootstrap_vcf_dir
#prefix=recent_geneflow_dualpopgrowth_altspec
#sfs_output_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/mid_bs_creation_popgrow_altspec/sfs_creation
#pop_spec=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/sfs_creation/pop_spec_mid_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.txt
##16,16 --> projection

### NORTH ###
#bs_dataset_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/north_bs_creation_popgrow_altspec/bootstrap_vcf_dir
#prefix=recent_geneflow_dualpopgrowth_altspec
#sfs_output_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/north_bs_creation_popgrow_altspec/sfs_creation
#pop_spec=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/sfs_creation/pop_spec_north_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.txt
##50,50 --> projection



######################
### RUNNING SCRIPT ###
######################

### MID ###
#/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $bs_dataset_path/${prefix}.bs.$SLURM_ARRAY_TASK_ID.vcf.gz -p $pop_spec --proj 16,16 -o $sfs_output_path/bs_sfs$SLURM_ARRAY_TASK_ID -a -f

### NORTH ###
#/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $bs_dataset_path/${prefix}.bs.$SLURM_ARRAY_TASK_ID.vcf.gz -p $pop_spec --proj 50,50 -o $sfs_output_path/bs_sfs$SLURM_ARRAY_TASK_ID -a -f


#only remove the bootstrap vcf if the delete_bs_dataset variable is "yes" and the sfs was successful (as indicated by presence of datadict.txt)
if [ $delete_bs_dataset = "yes" ] && [ -f "$sfs_output_path/bs_sfs$SLURM_ARRAY_TASK_ID/datadict.txt" ]
then
  rm $bs_dataset_path/${prefix}.bs.$SLURM_ARRAY_TASK_ID.vcf.gz
fi
