#!/bin/bash

sfs_array=( mid_region north_region )

top_dir=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/stairway_plot/analyses


if [ -d $top_dir/stairway_9_3_2021 ]
then
  echo "analysis directory already exists"
  exit 1
fi

mkdir $top_dir/stairway_9_3_2021 #top level directory the analyses
mkdir $top_dir/stairway_9_3_2021/sfs_dir #directory storing sfs files
mkdir $top_dir/stairway_9_3_2021/results #directory for results

touch $top_dir/stairway_9_3_2021/sfs_file_name.txt


for i in ${sfs_array[@]}
do
 
  #polyodon
  sfs_full_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/${i}_models/sfs_creation/easy_sfs_creation_maxmiss70_monomorph_7_7_2021/fastsimcoal2/polyodon_MAFpop0.obs
  cp $sfs_full_path $top_dir/stairway_9_3_2021/sfs_dir/${i}_polyodon_MAFpop0.obs
  
  #sfs_name=$(basename $sfs_full_path)
  echo ${i}_polyodon_MAFpop0.obs >> $top_dir/stairway_9_3_2021/sfs_file_name.txt


  #kazumbe
  sfs_full_path=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/${i}_models/sfs_creation/easy_sfs_creation_maxmiss70_monomorph_7_7_2021/fastsimcoal2/kazumbe_MAFpop0.obs
  cp $sfs_full_path $top_dir/stairway_9_3_2021/sfs_dir/${i}_kazumbe_MAFpop0.obs

  #sfs_name=$(basename $sfs_full_path)
  echo ${i}_kazumbe_MAFpop0.obs >> $top_dir/stairway_9_3_2021/sfs_file_name.txt  

done

#launch the stairway plot launcher!
sbatch launch_stairway.sh $top_dir/stairway_9_3_2021/sfs_file_name.txt $top_dir/stairway_9_3_2021/sfs_dir $top_dir/stairway_9_3_2021/results
