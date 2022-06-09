#!/bin/sh

### Job Name
#SBATCH --job-name start_pixy

### Declare an account for the job to run under
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_start_pixy
#SBATCH -e stderr_start_pixy

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
#SBATCH --time=0-00:10:00


module load gcc
#module load swset
module load r

[ ! -d full_species_kazumbe_subset_popfile_dir ] && mkdir full_species_kazumbe_subset_popfile_dir
[ ! -d pop_subset_kazumbe_subset_popfile_dir ] && mkdir pop_subset_kazumbe_subset_popfile_dir


Rscript create_subset_kazumbe_pop_files.R
