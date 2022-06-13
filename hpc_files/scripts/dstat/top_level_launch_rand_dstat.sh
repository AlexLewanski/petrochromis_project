#!/bin/bash

[ ! -d out_err_files ] && mkdir out_err_files

for SPECIES in kazumbe polyodon
do
  dstat_dir_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/DStats/DStats_Admixtools/${SPECIES}_random

  [ ! -d $dstat_dir_path ] && mkdir $dstat_dir_path
  sbatch run_randomized_dstat_${SPECIES}.sh
  
done
