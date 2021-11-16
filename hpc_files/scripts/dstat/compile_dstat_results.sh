#!/bin/sh

### Job Name
#SBATCH --job-name dstat_collect

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=ltcichlidgenomics

#SBATCH -o stdout_dstat_compile_results
#SBATCH -e stderr_dstat_compile_results

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=NONE
#SBATCH --mail-user=alewansk@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=0-01:00:00



### PATHS AND OTHER INFO ###
#array containing the directory names of the different D-stat calculations
dstat_array=( DStats_type_a DStats_type_b_kazumbe DStats_type_b_polyodon )

#path to the directories containing the D-stat results/working directories
dstat_top_dir_path="/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021"

#location where the output files will be located
output_dir_path="/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/DStats_9_19_2021"



### COLLECT AND PROCESS D-STAT RESULTS ###
for d_dir in ${dstat_array[*]}
do
    #enter into the d_stat result directory
    cd ${dstat_top_dir_path}/${d_dir}/   

    #get the file name
    file_name=$(find . -name log_dstat* -exec basename {} \;) #file name with file extension (but with path removed)

    file_strip_ext=$(echo "${file_name}" | cut -f 1 -d '.' | sed -e "s/^log_dstat_//") #file name without file extension and with "log_dstat_" removed

    #create file to store processed results
    #(1/2) NOTE: THIS TAXA ORDERING IS BASED ON THE UPDATE TO MY Dstat_input_files_function.R script where I order the taxa as: P1, P2, P3, outgroup
    #(2/2) The previous version had the inverse ordering: outgroup, P3, P2, P1
    echo -e P1"\t"P2"\t"P3"\t"outgroup"\t"D"\t"stderr"\t"z_score"\t"BABA"\t"ABBA"\t"n_snps > ${output_dir_path}/"dstat_results_${file_strip_ext}.txt"

    #extract/process relevant rows from the log file and append to the results file
    sed -n '/^result:/p' ${file_name} | sed 's/result://' >> ${output_dir_path}/"dstat_results_${file_strip_ext}.txt"
done



### MISC. COMMENTS/RESOURCES ###
#I was originally use to get the file name (stripped of of its file extension and "log_dstat_" because the name of the file was inexplicably starting with \n
#this problem was being caused by the file name being printed twice by the find function in the previous step, so tr -d '\n' is no longer needed
#file_strip_ext=$(echo "${file_name}" | cut -f 1 -d '.' | sed -e "s/^log_dstat_//" | tr -d '\n')

#Resources used for this:
#sed to delete characters: https://www.theunixschool.com/2014/08/sed-examples-remove-delete-chars-from-line-file.html
#sed to extract rows: https://www.theunixschool.com/2012/12/sed-10-examples-to-print-lines-from-file.html
#remove extensions from file name: https://stackoverflow.com/questions/12152626/how-can-i-remove-the-extension-of-a-filename-in-a-shell-script
#finding files with find: https://askubuntu.com/questions/1022172/save-the-result-of-find-as-a-variable-in-a-shell-script
#finding files with find: https://stackoverflow.com/questions/5456120/how-to-only-get-file-name-with-linux-find
