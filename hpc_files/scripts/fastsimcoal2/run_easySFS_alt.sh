#!/bin/sh

### Job Name
#SBATCH --job-name run_easySFS

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_easySFS
#SBATCH -e stderr_easySFS

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
####SBATCH --mem=0
#SBATCH --mem=130G

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=3-12:00:00



#####################
### SCRIPT SET-UP ###
#####################

#loading libraries
module load gcc
module load miniconda3

source activate easySFS



#####################
### CREATING SFSs ###
#####################

### NORTH ###
echo 'started creating SFS for north region'

north_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/sfs_creation
north_path_output=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/north_region_models/sfs_creation/easy_sfs_creation_maxmiss70_monomorph_7_7_2021
north_vcf=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/allsites_datasets_6_16_2021/north_allsites_polykaz/north_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.vcf.gz
north_pop_spec=pop_spec_north_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.txt


if [ -d $north_path_output ]
then
  echo "The directory specified for the north sfs output already exists. Please specify a directory that doesn't already exist."
  exit 1
else
  mkdir $north_path_output
fi


#awk '{print $1, $2}' /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/final_variant_files/fsc_vcf_processing/subsampled_fsc_dataset_north.txt > $north_path/$north_pop_spec


/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $north_vcf -p $north_path/$north_pop_spec -a --proj 50,50 -o $north_path_output -f

echo 'finished creating SFS for north region'


### MID ###
echo 'started creating SFS for mid region'

mid_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/sfs_creation
mid_path_output=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/mid_region_models/sfs_creation/easy_sfs_creation_maxmiss70_monomorph_7_7_2021
mid_vcf=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/allsites_datasets_6_16_2021/mid_allsites_polykaz/mid_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.vcf.gz
mid_pop_spec=pop_spec_mid_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.txt


if [ -d $mid_path_output ]
then
  echo "The directory specified for the mid sfs output already exists. Please specify a directory that doesn't already exist."
  exit 1
else
  mkdir $mid_path_output
fi


#awk '{print $1, $2}' /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/final_variant_files/fsc_vcf_processing/subsampled_fsc_dataset_mid.txt > $mid_path/$mid_pop_spec

/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/sfs/easySFS/easySFS.py -i $mid_vcf -p $mid_path/$mid_pop_spec -a --proj 16,16 -o $mid_path_output -f

echo 'finished creating	SFS for	mid region'
