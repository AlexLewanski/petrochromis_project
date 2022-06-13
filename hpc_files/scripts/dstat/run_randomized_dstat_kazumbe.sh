#!/bin/sh

####################
### SLURM SET-UP ###
###################
#SBATCH --job-name d_rand
#SBATCH --account=ltcichlidgenomics
#SBATCH --output=out_err_files/mod_%A_%a
#SBATCH --error=out_err_files/mod_%A_%a
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alewansk@uwyo.edu
#SBATCH --mem=0
#SBATCH --time=0-0:10:00
#SBATCH --array=1-1000%50



##########################
### SCRIPT SET-UP ###
##########################

### LIBRARIES ###
module load intel
module load intel-mkl
module load gsl
module load r


### FILE AND PATH INFO ###
dstat_dir_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/kazumbe_random
ind_file=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021/DStats_type_b_kazumbe/out_root.ind
pop_file=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021/DStats_type_b_kazumbe/popfile_TYPEBKazumbe_colormorphs.txt
geno=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021/DStats_type_b_kazumbe/out_root.geno
snp=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021/DStats_type_b_kazumbe/out_root.snp
par_dstat_file=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021/DStats_type_b_kazumbe/par_d_stat



##################################
### SETTING UP THE DSTAT CALCS ###
##################################

#[ ! -d $dstat_dir_path ] && mkdir $dstat_dir_path

#create directory
#if the job is requeued then the directory already exists (which causes problems) so if it already exists the directory is just deleted before being created
[ -d $dstat_dir_path/rand_dstat_kazumbe_$SLURM_ARRAY_TASK_ID ] && rm -rf $dstat_dir_path/rand_dstat_kazumbe_$SLURM_ARRAY_TASK_ID 
mkdir $dstat_dir_path/rand_dstat_kazumbe_$SLURM_ARRAY_TASK_ID

#create randomized ind file
Rscript randomize_ind.R $ind_file $dstat_dir_path/rand_dstat_kazumbe_$SLURM_ARRAY_TASK_ID/out_root.ind kazumbe

#copy par_dstat and popfile over to new directory
cp $pop_file $geno $snp $par_dstat_file $dstat_dir_path/rand_dstat_kazumbe_$SLURM_ARRAY_TASK_ID/

#enter into the new dstat directory
cd $dstat_dir_path/rand_dstat_kazumbe_$SLURM_ARRAY_TASK_ID

sed -i "s+DIR:.*+DIR:  $dstat_dir_path/rand_dstat_kazumbe_$SLURM_ARRAY_TASK_ID+" par_d_stat



##########################
### RUNNING DSTAT CALC ###
##########################

sh /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/run_dstat_alt.sh kazumbe

#get the file name
file_name=$(find . -name log_dstat* -exec basename {} \;) #file name with file extension (but with path removed)

file_strip_ext=$(echo "${file_name}" | cut -f 1 -d '.' | sed -e "s/^log_dstat_//") #file name without file extension and with "log_dstat_" removed
#create file to store processed results
#(1/2) NOTE: THIS TAXA ORDERING IS BASED ON THE UPDATE TO MY Dstat_input_files_function.R script where I order the taxa as: P1, P2, P3, outgroup
#(2/2) The previous version had the inverse ordering: outgroup, P3, P2, P1
echo -e P1"\t"P2"\t"P3"\t"outgroup"\t"D"\t"stderr"\t"z_score"\t"BABA"\t"ABBA"\t"n_snps > "dstat_results_${file_strip_ext}.txt"

#extract/process relevant rows from the log file and append to the results file
sed -n '/^result:/p' ${file_name} | sed 's/result://' >> "dstat_results_${file_strip_ext}.txt"

rm out_root.geno
rm out_root.snp
