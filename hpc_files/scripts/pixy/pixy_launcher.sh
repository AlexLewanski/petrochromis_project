#!/bin/sh

### Job Name
#SBATCH --job-name pixy

### Declare an account for the job to run under
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_pixy
#SBATCH -e stderr_pixy

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#####SBATCH --cpus-per-task=4
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=5-00:00:00


###Load libraries
module load gcc
module load miniconda3

source activate pixy

pixy --stats pi dxy \
--vcf $1 \
--zarr_path $2 \
--populations $3 \
--window_size 500000 \
--bypass_filtration yes \
--outfile_prefix $4
