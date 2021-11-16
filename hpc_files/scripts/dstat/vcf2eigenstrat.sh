#!/bin/sh

### This is a general SLURM script. You'll need to make modifications for this to 
### work with the appropriate packages that you want. Remember that the .bashrc 
### file will get executed on each node upon login and any settings in this script
### will be in addition to, or will override, the .bashrc file settings. Users will
### find it advantageous to use only the specific modules they want or 
### specify a certain PATH environment variable, etc. If you have questions,
### please contact ARCC for help.

### Informational text is usually indicated by "###". Don't uncomment these lines.

### Lines beginning with "#SBATCH" are SLURM directives. They tell SLURM what to do.
### For example, #SBATCH --job-name my_job tells SLURM that the name of the job is "my_job".
### Don't remove the "#SBATCH".

### Job Name
#SBATCH --job-name vcf_2_eigen

### Declare an account for the job to run under
#SBATCH --account=WagnerLab

### By default, the standard output and error streams are sent to files in the current 
### working directory with names:
### job_name.osequence_number  (output stream)
### job_name.esequence_number  (error stream)
### where job_name is the name of the job and sequence_number is the job number assigned 
### when the job is submitted.
### Use the directives below to change the files to which the standard output and 
### error streams are sent.
#SBATCH -o stdout_file
#SBATCH -e stderr_file

### The directive below directs that the standard output and error streams are to be merged, 
### intermixed, as standard output. 
# # # This doesn't work yet....  SBATCH -j oe

### mailing options
#SBATCH --mail-type=END
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=1-00:00:00

### Change to working directory (location where qsub was done)
#  #  #  echo "Working Directory:  $PBS_O_WORKDIR"
#  #  # cd $PBS_O_WORKDIR

#module load swset/2018.05
module load gcc/7.3.0
module load python/2.7.15


### FILE PATHS ###
vcf=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/final_variant_files/dataset_processing_6_13_2021/biallelic_datasets/biallel_polykazoutgroups/biallel_polykazoutgroups_minq5_maxDP75_miss0.70.recode.vcf
dstat_typea="/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021/DStats_type_a"
dstat_typeb_kazumbe="/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021/DStats_type_b_kazumbe"
dstat_typeb_polyodon="/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021/DStats_type_b_polyodon"


### START THE JOB ###
#python vcf2eigenstrat.py -v /gscratch/rolson14/Tropheini/filtered_tropheini_3_nov_19.vcf -o out_root -i tropheini_nov_3_2019.ind
#python vcf2eigenstrat.py -v /project/ltcichlidgenomics/lt_cichlid_vcf/pundamilia_post_4.vcf -o out_root -i /project/ltcichlidgenomics/alewansk/Tropheine_Dstats/Dstats_12_9_19/tropheine_indfile_12_9_19.txt

python vcf2eigenstrat.py -v $vcf -o ${dstat_typea}/out_root -i ${dstat_typea}/indfile_TYPEA_colormorphs.txt

python vcf2eigenstrat.py -v $vcf -o ${dstat_typeb_kazumbe}/out_root -i ${dstat_typeb_kazumbe}/indfile_TYPEBKazumbe_colormorphs.txt

python vcf2eigenstrat.py -v $vcf -o ${dstat_typeb_polyodon}/out_root -i ${dstat_typeb_polyodon}/indfile_TYPEBPolyodon_colormorphs.txt
