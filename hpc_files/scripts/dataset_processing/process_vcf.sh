#!/bin/sh

##############################################################################################################
# AUTHOR: Alexander L. Lewanski
# DATE: June, 2021
# REQUIRED PACKAGES: bcftools, vcftools
#
# FUNCTIONALITY:
# From either vcf files or bgzipped vcf file vcf_processing_script.sh conducts the following processing steps:
# 1) subsets samples from vcf (optional)
# 2) filters vcf to particular region(s) (optional)
# 3) filters vcf by read depth, site-level missing data, map quality, allele count, and also removes indels
# 4) removes samples with too much missing data after filtering
#
# USAGE: 
# vcf_processing_script.sh -v vcf to be processed \
# -p directory where files will be stored \
# -s sample list (see bcftools for how include vs. exclude samples) \
# -q minimum quality score (default: 5) \
# -s sample list (see bcftools for how include vs. exclude samples) \
# -q minimum quality score (default: 5) \
# -i min alleles (default: 2) \
# -a max alleles (default: 2) \
# -f minor allele frequency (default: 0) \
# -d minimum mean read depth (default: 5) \
# -b maximum mean read depth (default: 100) \
# -o minimum read depth (default: 5) \
# -g maximum read depth (default: 100) \
# -h the regions to retain (default: all regions) \
# -x site-level missing data threshold (default: 0.50) \
# -l individual-level missing data threshold \
# -r inclusion r removes intermediate files \
# -n name of final vcf (default: final.vcf) \
# -z inclusion of z will bgzip the final file \
# -m include m flag if you want progress messages
##############################################################################################################


if [ "$#" -lt 1 ]
then
  echo "Performs filtering of vcf including (optional) removal and/or retention of samples,
  filtering (min and max alleles, minor allele frequency, min allele freq, min mean read
  depth), filter by region (e.g. chromosome), and filter by individual-level missing data.

  REQUIRED SOFTWARE
  bcftools, vcftools
  
  REQUIRED ARGUMENTS
  [-v] vcf to be processed (include path if vcf is not in the same directory as 
       this processing script)
  
  OPTIONAL ARGUMENTS
  [-p] directory where files will be stored
  [-s] sample list (see bcftools for how include vs. exclude samples)
  [-q] minimum quality score (default: 5)
  [-i] min alleles (default: 2)
  [-a] max alleles (default: 2)
  [-f] minor allele frequency (default: 0)
  [-d] minimum mean read depth (default: 5)
  [-b] maximum mean read depth (default: 100)
  [-o] minimum read depth (default: 5)
  [-g] maximum read depth (default: 100)
  [-h] the regions to retain (default: all regions)
  [-x] site-level missing data threshold (default: 0.50)
  [-l] individual-level missing data threshold, which occurs after the
       site-level missing data threshold step (default: 0.50)
  [-r] inclusion of r argument indicates that you want the intermediate vcf
       files to be removed
  [-n] name of final vcf (default: final.vcf)
  [-z] inclusion of -z argument indicates that you want the final vcf file
       to be bgzipped
  [-m] include m flag if you want progress messages"
  
else

  #####################
  ### SCRIPT SET-UP ###
  #####################
  
  while getopts p:v:s:q:i:a:f:d:b:o:g:h:x:l:rn:zm c
  do
    case $c in
      v) vcf=${OPTARG};;
      p) vcf_dir=${OPTARG};;
      s) samp_list=${OPTARG};;
      q) min_q=${OPTARG};;
      i) min_allele=${OPTARG};;
      a) max_allele=${OPTARG};;
      f) maf_val=${OPTARG};;
      d) min_mean_dp=${OPTARG};;
      b) max_mean_dp=${OPTARG};;
      o) min_dp=${OPTARG};;
      g) max_dp=${OPTARG};;
      h) chr_filt=${OPTARG};;
      x) site_missing=${OPTARG};;
      l) indiv_missing=${OPTARG};;
      r) remove_intermediate=remove;;
      n) output_name=${OPTARG};;
      z) zip_final=zip_file;;
      m) message=include_message;;
    esac
  done
  
  #set defaults
  vcf_dir=${vcf_dir:-processed_dir}
  samp_list=${samp_list:-none}
  min_q=${min_q:-5}
  min_allele=${min_allele:-2}
  max_allele=${max_allele:-2}
  maf_val=${maf_val:-0}
  min_mean_dp=${min_mean_dp:-5}
  max_mean_dp=${max_mean_dp:-100}
  min_dp=${min_dp:-5}
  max_dp=${max_dp:-100}
  chr_filt=${chr_filt:-no_chr_filt}
  site_missing=${site_missing:-0.50}
  indiv_missing=${indiv_missing:-0.50}
  output_name=${output_name:-final}
  
  
  
  ######################
  ### VCF PROCESSING ###
  ######################
  
  #PRE-STEP: check that the vcf exists and zipping the file if it is not already zipped
  
  #make directory to store the processed file (or exit if the directory already exists)
  if [ -d $vcf_dir ]
  then
    echo "please provide a different directory name to store the files because the one you specified already exists"
    exit 1
  else
    mkdir $vcf_dir/
  fi
  
  #check if vcf file exists and exit if it doesn't
  if [ ! -f $vcf ]
  then
    echo "The vcf file you specified does not exist. Check name and/or path of the file."
    exit 1
  fi
  

  #zipping the file if it is not already zipped
  extension=`echo $vcf | sed 's/^.*\.//'`
  #file $vcf | grep -q compressed; echo $?
  if [ $extension == "vcf" ]
  then
    if [ $message == "include_message" ]; then echo 'bgzipping vcf file'; fi
    vcf_name=`basename -s .vcf $vcf`
    bgzip -c $vcf > $vcf_dir/${vcf_name}_${output_name}.vcf.gz
    prestep_file=$vcf_dir/${vcf_name}_${output_name}.vcf.gz
    if [ $message == "include_message" ]; then echo 'finished bgzipping vcf file'; fi
  else
    prestep_file=$vcf
  fi

  if [ $message == "include_message" ]; then echo '##           (16.67%) finished initial set-up'; fi


  #STEP 1: filter sample step (optional)
  if [ "$samp_list" != "none" ]
  then
    samp_list_name=`echo $samp_list | sed 's/^\^//'`
    if [ -f "$samp_list_name" ]
    then
      bcftools view \
      -O z \
      -S $samp_list $prestep_file \
      -o $vcf_dir/sample_filter_${output_name}.vcf.gz

      filter_sample_name=$vcf_dir/sample_filter_${output_name}.vcf.gz
      
      if [ "$message" == "include_message" ]; then echo '####         (33.33%) finished filtering by samples'; fi
    else
      echo "If you want to filter the vcf by samples, you need to include a list of samples to exclude and/or include."
      exit 1
    fi
  else
    filter_sample_name=$prestep_file
  fi
  
  
  #STEP 2: filter by region (optional)
  if [ $chr_filt != "no_chr_filt" ]
  then
    #vcf file needs to be indexed prior to subsetting my region
    bcftools index -f $filter_sample_name
      
    bcftools view -O z \
    --regions $chr_filt $filter_sample_name \
    -o $vcf_dir/step2_${output_name}.vcf.gz
      
    filter_region_name=$vcf_dir/step2_${output_name}.vcf.gz
      
    if [ "$message" == "include_message" ]; then echo '######       (50%) finished quality filtering'; fi
    
    else
    filter_region_name=$filter_sample_name
  fi

  
  #STEP 3: quality filtering

  if [ "$message" == "include_message" ]
  then
    echo "Using the following command to filter vcf:
    vcftools --remove-indels 
    --minQ $min_q
    --min-alleles $min_allele 
    --max-alleles $max_allele 
    --maf $maf_val
    --max-missing $site_missing 
    --min-meanDP $min_mean_dp 
    --max-meanDP $max_mean_dp 
    --minDP $min_dp #currently excluded
    --maxDP $max_dp #currently excluded
    --gzvcf $filter_region_name 
    --recode --recode-INFO-all
    --stdout | gzip -c > $vcf_dir/step3_${output_name}.vcf.gz"
  fi
    
  vcftools --remove-indels \
  --minQ $min_q \
  --min-alleles $min_allele \
  --max-alleles $max_allele \
  --maf $maf_val \
  --max-missing $site_missing \
  --min-meanDP $min_mean_dp \
  --max-meanDP $max_mean_dp \
  --gzvcf $filter_region_name \
  --recode --recode-INFO-all \
  --stdout | gzip -c > $vcf_dir/step3_${output_name}.vcf.gz

#--min-meanDP $min_mean_dp \
#--max-meanDP $max_mean_dp \  
#  --minDP $min_dp \
#  --maxDP $max_dp \

  if [ "$message" == "include_message" ]; then echo '########     (66.67%) finished filter by region'; fi
    
  
  #STEP 4: remove samples with too much missing data
  vcftools --gzvcf $vcf_dir/step3_${output_name}.vcf.gz --missing-indv --out $vcf_dir/missing_data_info
  
  if [ ! -f $vcf_dir/missing_data_info.imiss ]
  then
    echo "stopping script because missing data file (missing_data_info.imiss) was not successfully created"
    exit 1
  fi

  awk -v thresh="$indiv_missing" 'NR > 1 && $5 > thresh {print $1}' $vcf_dir/missing_data_info.imiss > $vcf_dir/lowDP.indv

  if [ "$message" == "include_message" ]
  then
    echo "The following samples will be removed based on the missing data threshold (${indiv_missing}):
          $(cat $vcf_dir/lowDP.indv)"
  fi
  
  if [ "$zip_final" == "zip_file" ]
  then
    vcftools --gzvcf $vcf_dir/step3_${output_name}.vcf.gz \
    --remove $vcf_dir/lowDP.indv \
    --recode --recode-INFO-all \
    --stdout | gzip -c > $vcf_dir/${output_name}.vcf.gz
  else
    vcftools --gzvcf $vcf_dir/step3_${output_name}.vcf.gz \
    --remove $vcf_dir/lowDP.indv \
    --recode --recode-INFO-all \
    --out $vcf_dir/${output_name}
  fi
  
  if [ "$message" == "include_message" ]; then echo '##########   (83.33%) finished individual-level missing filter'; fi
  
    
  #STEP 5: REMOVE INTERMEDIATE FILES (OPTIONAL)
  if [ "$remove_intermediate" == "remove" ]
  then
    rm -f $vcf_dir/sample_filter_${output_name}.vcf.gz $vcf_dir/step2_${output_name}.vcf.gz $vcf_dir/step3_${output_name}.vcf.gz
    
    if [ $extension == "vcf" ]
    then
      rm $vcf_dir/${vcf_name}_${output_name}.vcf.gz
    fi
  fi
  
  if [ "$message" == "include_message" ]; then echo '############ (100%) processing complete'; fi

fi



#resources:
#https://stackoverflow.com/questions/13617843/unary-operator-expected-error-in-bash-if-condition
