#!/bin/sh

if [ "$#" -lt 1 ]
then
  echo "Creates the blueprint file for running Stairway Plot v2.1."

  REQUIRED ARGUMENTS
  [-q] number of sequences (twice the number of diploid individuals)
  [-l] total number of observed nucleic (polymorphic + monomorphic) sites
  [-s] the sfs
  [-r] number of random break points for each try (separated by white space)
  
  OPTIONAL ARGUMENTS
  [-n] name of output blueprint file and the popid (default: input)
  [-f] inclusion of f indicates that the sfs is folded
  [-m] the smallest size of SFS bin used for estimation (default: 1)
  [-b] the largest size of SFS bin used for estimation (default: nseq/2)
  [-p] percentage of sites used in training (default: 0.67)"
  
else

  #####################
  ### SCRIPT SET-UP ###
  #####################
  
  while getopts n:q:l:fs:m:b:p:r:i:k:t:e:u:g:h:x:y:z:a:o: c
  do
    case $c in
      n) name=${OPTARG};;
      q) nseq=${OPTARG};;
      l) length=${OPTARG};;
      f) folded=true;;
      s) sfs=${OPTARG};;
      m) smallest_bin_size=${OPTARG};;
      b) biggest_bin_size=${OPTARG};;
      p) pct_training=${OPTARG};;
      r) nrand=${OPTARG};;
      i) project_dir=${OPTARG};;
      k) stairway_plot_dir=${OPTARG};;
      t) ninput=${OPTARG};;
      e) random_seed=${OPTARG};;
      u) mu=${OPTARG};;
      g) year_per_generation=${OPTARG};;
      h) plot_title=${OPTARG};;
      x) xrange=${OPTARG};;
      y) yrange=${OPTARG};;
      z) xspacing=${OPTARG};;
      a) yspacing=${OPTARG};;
      o) fontsize=${OPTARG};;
    esac
  done

  #set defaults
  name=${name:-input}
  folded=${folded:-false}
  smallest_bin_size=${smallest_bin_size:-1}

  #beginning_biggest_bin_size=${biggest_bin_size:-"#"}
  [ -z "$biggest_bin_size" ] && beginning_biggest_bin_size="#"
  biggest_bin_size=${biggest_bin_size:-"default is nseq/2"}

  pct_training=${pct_training:-0.67}

  [ -z "$random_seed" ] && beginning_random_seed="#"
  random_seed=${random_seed:-"NONE"}

  stairway_plot_dir=${stairway_plot_dir:-stairway_plot_es}
  ninput=${ninput:-200}
  plot_title=${plot_title:-$name}
  xrange=${xrange:-"0,0"}
  yrange=${yrange:-"0,0"}
  xspacing=${xspacing:-1}
  yspacing=${yspacing:-1}
  fontsize=${fontsize:-12}


#############################
### OUTPUT BLUEPRINT FILE ###
#############################
      
  echo "###############################################
### BLUEPRINT FILE (for Stairway Plot v2.1) ###
###############################################
#input setting
popid: $name # id of the population (no white space)
nseq: $nseq # number of sequences
L: $length # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: $folded # whether the SFS is folded (true or false)
SFS:    $sfs # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)
smallest_size_of_SFS_bin_used_for_estimation: $smallest_bin_size # default is 1; to ignore singletons, uncomment this line and change this number to 2
${beginning_biggest_bin_size}largest_size_of_SFS_bin_used_for_estimation: $biggest_bin_size  # default is nseq/2 for folded SFS
pct_training: $pct_training # percentage of sites for training
nrand: $nrand # number of random break points for each try (separated by white space)
project_dir: $project_dir # project directory
stairway_plot_dir: $stairway_plot_dir # directory to the stairway plot files
ninput: $ninput # number of input files to be created for each estimation
${beginning_random_seed}random_seed: $random_seed
#output setting
mu: $mu # assumed mutation rate per site per generation (.e.g 1.2e-8)
year_per_generation: $year_per_generation # assumed generation time (in years)
#plot setting
plot_title: $plot_title # title of the plot
xrange: $xrange # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: $yrange # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: $xspacing # X axis spacing
yspacing: $yspacing # Y axis spacing
fontsize: $fontsize # Font size" > ${name}.blueprint

fi
