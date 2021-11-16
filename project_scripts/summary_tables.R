##############################################################################################
### SCRIPT NAME: summary_tables.R
### PURPOSE: create summary tables for genomic dataset attributes and sampling information
### PRODUCTS:
###     sampling_summary_table.txt: information about sampling including location name,
###                                 letter used for each location in the main text, latitude
###                                 and longitude, and the number of samples from each 
###                                 species. This table is included in the supplementary
###                                 material.
###     dataset_postprocessing_table.txt: attributes of each processed dataset including number
###                                       number of sites and read depth information. This 
###                                       table is included in the supplementary material.
###     dataset_specs_table.txt: filtering thresholds for each dataset including whether they
###                              contain only variants, min quality, minor allele frequency,
###                              minimum mean read depth, and maximum mean read depth. This 
###                              table is included in the supplementary material.
###     read_info: dataframe containing information about the details about sequencing
###                including the number of reads and amount of reads mapped to the reference
###                genome. This information is not written as a file but is reported in the
###                Results section of the main text.
##############################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(tidyverse)
library(kableExtra)
library(abind)
library(rhdf5)
library(here)


### Loading data ###

#metadata
cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata

### ENTROPY ###
entropy_samples_order <- read.delim(here('working_metadata_file', 'SAMPLES_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.txt'), header = FALSE) #order of samples in entropy results

#load K = 2 results 
mod_list <- list()
for (chain in 1:3) {
  mod_list[[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0('kazumbe_polyodon_qmod_6_18_2021_', 2, 'rep', chain, '.hdf5')), "q")
}

#metadata
all_taxa_samples <- read.delim(here('working_metadata_file', 'SAMPLES_biallel_polykazoutgroups_minq5_maxDP75_miss0.70.recode.txt'), header = FALSE)

#genomic dataset attributes (read information, number of sites, read depth info, etc...)
read_info <- read.table(here('sequencing_summary', 'dataset_summaries_9_20_2021', 'assembled_per_ind.txt'), header = TRUE)
vcf_summary <- read.table(here('sequencing_summary', 'dataset_summaries_9_20_2021', 'vcf_summary_file.txt'), header = TRUE)



##########################
### INITIAL PROCESSING ###
##########################

#combining entropy and metadata
k2_combined.q <- abind(mod_list$chain_1, mod_list$chain_2, mod_list$chain_3, along=1)
k2_q_estimate <- as.data.frame(t(apply(k2_combined.q, 2:3, mean)))

metadata_subset_entropy <- cichlid_metadata[cichlid_metadata$sample_code %in% entropy_samples_order$V1, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2')]
metadata_subset_entropy_reorder <- metadata_subset_entropy[match(as.character(entropy_samples_order$V1), metadata_subset_entropy$sample_code),]

metadata_subset_entropy_reorder$sciname_collapse <- gsub(" ", "_", metadata_subset_entropy_reorder$sciname2)
k2_q_estimate_with_metadata <- cbind(k2_q_estimate, metadata_subset_entropy_reorder)

#V2 greater than 0.5 --> majority kazumbe ancestry
#V2 less than 0.5 --> majority kazumbe ancestry
#there are no individuals with exactly 50/50 kazume/polyodon ancestry
k2_q_estimate_with_metadata$species_ID_entropy[k2_q_estimate_with_metadata$V2 > 0.5] <- 'P. kazumbe'
k2_q_estimate_with_metadata$species_ID_entropy[k2_q_estimate_with_metadata$V2 < 0.5] <- 'P. polyodon'



################################################
### SUMMARY TABLE #1: SAMPLING SUMMARY TABLE ###
################################################

#how many individuals per species?
table(k2_q_estimate_with_metadata$species_ID_entropy)

cichlid_metadata_processed <- k2_q_estimate_with_metadata %>%
  select(sample_code, location, species_ID_entropy, long, lat) %>% 
  #filter(sciname2 %in% c("Petrochromis kazumbe", "Petrochromis polyodon")) %>% 
  mutate("Location name" = case_when(location == 'Gombe_South' ~ 'Gombe S',
                                     location == 'Katongwe_N' ~ 'Katongwe N',
                                     location == 'Katongwe_S' ~ 'Katongwe S',
                                     location == 'Ska' ~ 'S Kagongo',
                                     location == 'Kalala' ~ 'Kalalangabo',
                                     location == 'Nondwa' ~ 'Nondwa',
                                     location == 'Hilltop' ~ 'Hilltop',
                                     location == 'Bangwe' ~ 'Bangwe',
                                     location == 'Jakob' ~ "Jakobsen's Beach",
                                     location == 'Ulombola' ~ 'Ulombola'),
         "Location letter" = factor(case_when(location == 'Gombe_South' ~ 'A',
                                              location == 'Katongwe_N' ~ 'B',
                                              location == 'Katongwe_S' ~ 'C',
                                              location == 'Ska' ~ 'D',
                                              location == 'Kalala' ~ 'E',
                                              location == 'Nondwa' ~ 'F',
                                              location == 'Hilltop' ~ 'G',
                                              location == 'Bangwe' ~ 'H',
                                              location == 'Jakob' ~ 'I',
                                              location == 'Ulombola' ~ 'J'),
                                    levels = (LETTERS[1:10] ))) %>% 
  mutate(Region = case_when(`Location letter` %in% c('A', 'B', 'C', 'D', 'E', 'F') ~ 'North',
                            `Location letter` %in% c('G', 'H', 'I') ~ 'Mid',
                            `Location letter` %in% 'J' ~ 'South'),
         Species = species_ID_entropy,
         Lon = round(long, 3),
         Lat = round(lat, 3) ) %>% 
  group_by(`Location letter`, Species) %>% 
  summarise(Region = first(Region),
            Species = first(Species),
            `Location name` = first(`Location name`),
            `Location letter` = first(`Location letter`),
            Lat = first(Lat),
            Lon = first(Lon),
            n = n()) %>% 
  pivot_wider(names_from = Species, values_from = n, values_fill = 0) %>%
  rename("P. sp. 'kazumbe'" = `P. kazumbe`,
         "P. cf. polyodon" = `P. polyodon`) %>%
  relocate(Region, .before = `Location letter`)

sampling_summary_table <- kbl(cichlid_metadata_processed, 'latex', booktabs = TRUE, align = "c") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  collapse_rows(1, latex_hline = "major")

cat(sampling_summary_table, file = here('tables', 'sampling_summary_table.txt'), append = FALSE)


################################################
### SUMMARY TABLE #2: SAMPLING SUMMARY TABLE ###
################################################

dataset_names <- c('focaltaxa_variant_missing70_maf1',
                   'alltaxa_variant_missing70_maf1',
                   'focaltaxa_allsites_missing70_maf0',
                   'mid_focaltaxa_allsites_missing70_maf0',
                   'north_focaltaxa_allsites_missing70_maf0')

dataset_name_conversion <- data.frame(
  new_names = dataset_names,
  old_names = c('biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode',
                'biallel_polykazoutgroups_minq5_maxDP75_miss0.70.recode',
                'allsites_polykaz_minq20_minDP5_maxDP75_miss0.70',
                'mid_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70',
                'north_allsites_polykaz_minq20_minDP5_maxDP75_miss0.70')
)

vcf_summary[, -1] <- round(vcf_summary[, -1], 2)

dataset_postprocessing_specs <- merge(vcf_summary, dataset_name_conversion, by.x = 'dataset', by.y = 'old_names') %>%
  select(!dataset & !mean_indiv_read_depth_sd & !mean_indiv_read_depth_mean) %>% 
  relocate(new_names, .before = n_sites) %>% 
  rename(Dataset = new_names,
         Sites = n_sites,
         `Average mean read depth` = mean_site_read_depth_mean,
         `Mean read depth SD` = mean_site_read_depth_sd)

dataset_postprocessing_specs_reorder <- dataset_postprocessing_specs[match(dataset_names, dataset_postprocessing_specs$Dataset),]
rownames(dataset_postprocessing_specs_reorder) <- NULL


dataset_postprocessing_table <- kbl(dataset_postprocessing_specs_reorder, 'latex', booktabs = TRUE, align = "c") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position") )

cat(dataset_postprocessing_table, file = here('tables', 'dataset_postprocessing_table.txt'), append = FALSE)



#################################################
### SUMMARY TABLE #3: DATASET FILTERING SPECS ###
#################################################

variant_only <- c('true', 'true', 'false', 'false', 'false')
min_quality <- rep(20, 5)
maf <- c(0.01, 0.01, 0, 0, 0)
min_dp <- rep(20, 5)
max_dp <- rep(75, 5)

dataset_specs <- data.frame("Dataset" = dataset_names,
                            "Variant-only" = variant_only,
                            "Min quality" = min_quality,
                            "MAF" = maf,
                            "Min mean DP" = min_dp,
                            "Max mean DP" = max_dp,
                            check.names = FALSE)

#adding information about the analyses that each dataset is used for
dataset_specs <- dataset_specs %>% 
  mutate(Analyses = case_when(Dataset == 'focaltaxa_variant_missing70_maf1' ~ "Fst, popvae, PCA, entropy",
                              Dataset == 'alltaxa_variant_missing70_maf1' ~ "D statistics",
                              Dataset == 'focaltaxa_allsites_missing70_maf0' ~ "nucleotide diversity, dxy",
                              Dataset == 'mid_focaltaxa_allsites_missing70_maf0' ~ "mid region demographic modeling",
                              Dataset == 'north_focaltaxa_allsites_missing70_maf0' ~ "north region demographic modeling"))

dataset_specs_table <- kbl(dataset_specs, 'latex', booktabs = TRUE, align = "c") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position") )

cat(dataset_specs_table, file = here('tables', 'dataset_specs_table.txt'), append = FALSE)



############################################################################
### SUMMARY INFORMATION: SEQUENCING AND REFERENCE GENOME MAPPING DETAILS ###
############################################################################

read_info_processed <- read_info[read_info$ind %in% all_taxa_samples$V1,]

(read_info <- data.frame(total_reads = sum(read_info_processed$raw),
                         assembled_reads = sum(read_info_processed$assembled),
                         prop_mapped_reads = sum(read_info_processed$assembled)/sum(read_info_processed$raw),
                         average_assembled_reads_per_ind = sum(read_info_processed$assembled)/nrow(read_info_processed),
                         sd_assembled_reads_per_ind = sd(read_info_processed$assembled)/nrow(read_info_processed)))

