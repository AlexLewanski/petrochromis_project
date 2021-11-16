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
#SBATCH --job-name dstat

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
#SBATCH --mail-type=FAIL,END
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

module load intel
module load intel-mkl
module load gsl

### Start the job
### Command normally given on command line


# must set num_lines to number of lines in popfile
#num_lines=224


#use par_d_stat to extract number of lines in pop file
if [ -f par_d_stat ]
then
  pop_file_name=$(awk '$1 == "popfilename:" {print $2}' par_d_stat)

  if [ ! -f $pop_file_name ]
  then
    echo "the pop file in par_d_stat doesn't exist in the current directory"
    exit 1
  fi
  
  num_lines=$(wc -l $pop_file_name | awk '{print $1}')
else
  echo "par_d_stat not found"
  exit 1
fi


### TEST ###
echo "$pop_file_name"
echo "$num_lines"



batch_size=5

# # set logfile name for d or f states; ensure 
logfile="log_dstat_TYPE_A_colormorphs_9_19_21.txt"
# # ensure output file: 'logfile' does not exist as process appends to it
if [ -f "$logfile" ]
then
	echo "$logfile found"
	exit 1
else
	echo ""
fi
 
if [ $((num_lines%batch_size)) -ne 0 ]; then
	# loop based on lines in popfile / batch size
	for i in `seq 1 $(($num_lines/$batch_size+1))`
	do
		/project/ltcichlidgenomics/alewansk/lt_cichlid_fulldataset_exclude_sum2019/Admixtools/AdmixTools-master/bin/qpDstat -p par_d_stat -l $((($i-1)*batch_size+1)) -h $(($i*batch_size)) >> $logfile
	done
else
	# loop based on lines in popfile / batch size
	for i in `seq 1 $((num_lines/batch_size))`
	do
		/project/ltcichlidgenomics/alewansk/lt_cichlid_fulldataset_exclude_sum2019/Admixtools/AdmixTools-master/bin/qpDstat -p par_d_stat -l $((($i-1)*batch_size+1)) -h $(($i*batch_size)) >> $logfile
	done    
fi
