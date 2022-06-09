#!/bin/sh

module load bcftools

k2_processed_results=/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/entropy/entropy_6_13_2021/qmodel_fulldataset/processed_results/k2_processed_results.txt
vcf=/gscratch/alewansk/polyodon_kazumbe_bioinformatics/allsites_datasets_6_16_2021/allsites_polykaz/allsites_polykaz_minq20_minDP5_maxDP75_miss0.70.vcf.gz

vcf_basename=$(basename -s .vcf.gz $vcf)
bcftools query -l $vcf > SAMPLES_${vcf_basename}.txt

#species-level pop file
awk 'NR != 1 {print $1"\t"$4}' $k2_processed_results > INITIAL_full_species_7_1_2021.txt
awk 'FNR==NR {a[$1]; next}; $1 in a' SAMPLES_${vcf_basename}.txt INITIAL_full_species_7_1_2021.txt > full_species_7_1_2021.txt

if [ "$(wc -l < SAMPLES_${vcf_basename}.txt)" -ne "$(wc -l < full_species_7_1_2021.txt)" ]
then 
  echo "the species-level pop file does NOT include the same number of lines as the vcf"
  exit 1
else
  echo "the species-level pop file does include the same number of lines as the vcf"
fi


#location-specific pop file
awk 'NR != 1 {print $1"\t"$4"_"$5}' $k2_processed_results > INITIAL_population_subset_7_1_2021.txt
awk 'FNR==NR {a[$1]; next}; $1 in a' SAMPLES_${vcf_basename}.txt INITIAL_population_subset_7_1_2021.txt > population_subset_7_1_2021.txt

if [ "$(wc -l < SAMPLES_${vcf_basename}.txt)" -ne "$(wc -l < population_subset_7_1_2021.txt)" ]
then 
  echo "the pop-level pop file does NOT include the same number of lines as the vcf"
  exit 1
else
  echo "the pop-level pop file does include the same number of lines as the vcf"
fi


#resources:
#https://stackoverflow.com/questions/40785453/extract-rows-in-file-where-a-column-value-is-included-in-a-list
