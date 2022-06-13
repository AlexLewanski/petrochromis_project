#!/bin/sh

species=$1
logfile="log_dstat_colormorphs_9_19_21_RANDOM_${species}.txt"
batch_size=5


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
#echo "$pop_file_name"
#echo "$num_lines"


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
