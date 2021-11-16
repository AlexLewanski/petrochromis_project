#!/bin/sh

### Job Name
#SBATCH --job-name fsc_boot

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_fsc_bootstrap
#SBATCH -e stderr_fsc_bootstrap

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
####SBATCH --ntasks-per-node=16
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=1-00:00:00

#LIBRARIES
#module load gcc
#module load htslib



##################################
### 1. SETTING UP THE ANALYSES ###
##################################

### VARIABLE DESCRIPTIONS ###
#vcf_path: path to original (gzipped) vcf file (don't include / at end)
#vcf_name: name of vcf file (don't include path in name because that is specified in vcf_path)
#working_dir: the directory where bootstrap dataset and (optional) sfs creation will be conducted (don't include / at end)
#bootstrap_vcf_dir: name of directory where bootstraps will be stored (don't include path in the name; don't include /	at end)
#sfs_dir_path: name of directory where sfs files will be stored (don't include path in the name; don't include /  at end)
#sfs_script_dir: path to location of sfs script (don't include / at end)
#number_blocks: number of blocks to break dataset into when forming the bootstrap datasets

### MID ###
#vcf_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/allsites_datasets_6_16_2021/mid_allsites_polykaz
#vcf_name=mid_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.vcf.gz
#working_dir=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/mid_bs_creation_popgrow_altspec
#bootstrap_vcf_dir=bootstrap_vcf_dir
#sfs_dir_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/mid_bs_creation_popgrow_altspec
#sfs_script_dir=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts
#number_blocks=100

### NORTH ####
#vcf_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/allsites_datasets_6_16_2021/north_allsites_polykaz
#vcf_name=north_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.vcf.gz
#working_dir=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/north_bs_creation_popgrow_altspec
#bootstrap_vcf_dir=bootstrap_vcf_dir
#sfs_dir_path=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/fsc_boostrap_creation/north_bs_creation_popgrow_altspec
#sfs_script_dir=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts
#number_blocks=100


PREFIX="recent_geneflow_dualpopgrowth_altspec"
bootstraps=100

#should the last site file be removed before creating the bootstrap vcf files? (options: "yes", "no")
#NOTE: this should only be yes if the last block has very few sites compared to the other files (the "leftover" sites)
remove_last_file="yes"
delete_building_blocks="yes"
create_sfs="no"
record_info="yes"



############################
### 2. RUNNING  ANALYSES ###
############################

### 1. FIRST STEPS TOWARDS CREATING BOOTSTRAP DATASETS ###
if [ ! -d "$working_dir" ]
then
  mkdir $working_dir/
  cp $vcf_path/$vcf_name $working_dir/
else
  echo "The working directory already exists. Please specify a new directory."
  exit 1
fi

cd $working_dir/


#calculate the block size
nsites=`zgrep -v '^#' $vcf_name | wc -l`
blsize=$(($nsites/$number_blocks))


### 2. CREATE BOOTSTRAP DATASETS AND (OPTIONALLY) LAUNCH SFS CREATION SCRIPT ###
# Get all lines with genomic data
zgrep -v "^#" $vcf_name | gzip -c > $PREFIX.allSites.gz

# Get the header
zgrep "^#" $vcf_name > header

# get 100 files with total_number_sites/100 sites each
#split -l $blsize $PREFIX.allSites $PREFIX.sites.
zcat $PREFIX.allSites.gz | split -l $blsize --filter='gzip > $FILE.gz' - $PREFIX.sites.

if [ $record_info = "yes" ]
then 
  trailing_file=`find . -maxdepth 1 -name "$PREFIX.sites.*" | sort | tail -1`
  nsites_finalblock=`zcat $trailing_file | wc -l`
fi


if [ $remove_last_file = "yes" ]
then
  #ls $PREFIX.sites.* | sort | tail -1 | xargs rm
  find . -maxdepth 1 -name "$PREFIX.sites.*" | sort | tail -1 | xargs rm
fi


# Make a new folder to store bs vcf files
mkdir $bootstrap_vcf_dir
cd $bootstrap_vcf_dir

if [ $create_sfs = "yes" ]
then
  if [ -d $sfs_dir_path/sfs_dir ]
  then
    echo "The directory where the sfs files are specified to be stored already exists. Please specify a new directory."
    exit 1
  fi

  mkdir $sfs_dir_path/sfs_dir
fi

for i in $(eval echo "{1..$bootstraps}")
do

  # Add the header to our new bootstrapped vcf file
  cat ../header > $PREFIX.bs.$i.vcf
  
  #Randomly add $number_blocks blocks to bootstrap dataset
  for r in $(eval echo "{1..$number_blocks}")
  do
    zcat `shuf -n1 -e ../$PREFIX.sites.*` >> ${PREFIX}.bs.$i.vcf
  done
  
  # Compress the fully created vcf file again
  gzip ${PREFIX}.bs.$i.vcf

  #LAUNCH SCRIPT TO CREATE SITE FREQUENCY SPECTRUM
  if [ $create_sfs = "yes" ]
  then
    sbatch $sfs_script_dir/sfs_creation.sh $working_dir/$bootstrap_vcf_dir/${PREFIX}.bs.$i.vcf.gz $pop_spec $sfs_dir_path/sfs_dir bs$i
    sleep 1
  fi

  echo "finished the following bootstrap dataset: $i"

done


### 3. FINAL STEPS: OPTIONAL CREATION OF INFORMATION FILE AND DELETION OF BLOCK BOOTSTRAP BUILDING BLOCKS ###
if [ $record_info = "yes" ]
then 
  echo "remove trailing block: $delete_building_blocks" > bootstrap_dataset_creation_notes.txt

  number_of_blocks=`find .. -maxdepth 1 -name "${PREFIX}.sites.*" -type f | wc -l`
  echo "number of blocks used in bs dataset creation: $number_of_blocks" >> bootstrap_dataset_creation_notes.txt

  echo "number of sites per block: $blsize" >> bootstrap_dataset_creation_notes.txt
  echo "number sites in trailing block: $nsites_finalblock" >> bootstrap_dataset_creation_notes.txt
fi

if [ $delete_building_blocks = "yes" ]
then
  cd ..
  find . -maxdepth 1 -name "${PREFIX}.sites.*" -delete
  rm header $PREFIX.allSites.gz
fi

