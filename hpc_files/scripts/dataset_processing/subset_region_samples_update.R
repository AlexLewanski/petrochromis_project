#!/usr/bin/env Rscript


#####################
### SCRIPT SET-UP ###
#####################

#load libraries
library(optparse)

#parsing region input
arg_list <- list(
        make_option(c("-r", "--region"), type = "character", default = NULL),
        make_option(c("-m", "--metadata"), type = "character"),
        make_option(c("-f", "--file_name"), type = "character", default = NULL),
        make_option(c("-v", "--validate"), type = "logical", default = FALSE),
        make_option(c("-w", "--validate_file"), type = "character", default = NULL)
)


option <- parse_args(OptionParser(option_list=arg_list))
region <- option$region
metadata <- option$metadata
file_name <- option$file_name
validate <- option$validate
validate_file_name <- option$validate_file

#check inputs
if (is.null(validate_file_name) & !validate)
	stop('you need to provide a file name if you want to validate')


#region specific attributes (locations, output path)
region_list <- list(north = c('Nondwa', 'Kalala', 'Ska', 'Katongwe_S', 'Katongwe_N', 'Gombe_South'),
					mid = c('Hilltop', 'Jakob', 'Bangwe'))

#load metadata
cichlid_metadata <- read.csv(metadata, stringsAsFactors = FALSE)


###########################
### CREATE SPECIES LIST ###
###########################

#subset down to region
region_samples <- cichlid_metadata[ (cichlid_metadata$location %in% region_list[[region]]) & (cichlid_metadata$sciname2 %in% c('Petrochromis kazumbe', 'Petrochromis polyodon')),]$sample_code

#load samples found in the vcf
if (validate) {
	vcf_samples <- read.table(validate_file_name, stringsAsFactors = FALSE)
	processed_samples <- vcf_samples$V1[vcf_samples$V1 %in% region_samples]
} else {
	processed_samples <- region_samples
}


write.table(processed_samples, file = file_name,
        row.names = FALSE, col.names = FALSE, quote = FALSE)



					
#path_list <- list(north = '/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/entropy/entropy_ancestry_complement_may2021_north/',
#					mid = '/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/entropy/entropy_ancestry_complement_may2021_mid/')


#load samples found in the vcf
#vcf_samples <- read.table('SAMPLE_IDs_kazumbe_polyodon_pundamilia_8_19_2020_maf0.01_maxmissing0.70_minDP5_post_4.recode.txt', stringsAsFactors = FALSE)

#load metadata
#cichlid_metadata <- read.csv('/project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/entropy/monster_23jul20.csv', stringsAsFactors = FALSE)


#region_samples <- cichlid_metadata[ (cichlid_metadata$location %in% region_list[[region]]) & (cichlid_metadata$sciname2 %in% c('Petrochromis kazumbe', 'Petrochromis polyodon')),]$sample_code
#processed_samples <- vcf_samples$V1[vcf_samples$V1 %in% region_samples]
#write.table(processed_samples, file = paste0(path_list[[region]], file_name, '.txt'),
#	row.names = FALSE, col.names = FALSE, quote = FALSE)

