##############################################################################################
### SCRIPT NAME: admixtools_input_files_creation.R
### PURPOSE: create input files for calculating D statistics with Admixtools. This includes
###          the pop file, which specifies the set of taxa involved in each D statistic 
###          calculation, and the ind file, which specifies each sample's taxon identity.
### PRODUCTS:
###     popfile_TYPEA_colormorphs.txt: pop file for D stats with the following setups:
###                                    1b: kazumbe_pop1, kazumbe_pop2, polyodon_pop2;
###                                    2b: polyodon_pop1, polyodon_pop2, kazumbe_pop2
###     indfile_TYPEA_colormorphs.txt: ind file for D stats with the following setups:
###                                    1b: kazumbe_pop1, kazumbe_pop2, polyodon_pop2;
###                                    2b: polyodon_pop1, polyodon_pop2, kazumbe_pop2
###     popfile_TYPEBKazumbe_colormorphs.txt: pop file for D stats with the following 
###                                           setup:
###                                           1a: kazumbe_pop1, kazumbe_pop2, polyodon_all
###     indfile_TYPEBKazumbe_colormorphs.txt: ind file for D stats with the following 
###                                           setup:
###                                           1a: kazumbe_pop1, kazumbe_pop2, polyodon_all
###     popfile_TYPEBPolyodon_colormorphs.txt: pop file for D stats with the following 
###                                            setup:
###                                            1a: polyodon_pop1, polyodon_pop2, kazumbe_all
###     indfile_TYPEBPolyodon_colormorphs.txt: ind file for D stats with the following 
###                                            setup:
###                                            1a: polyodon_pop1, polyodon_pop2, kazumbe_all
##############################################################################################


##########################
### SCRIPT PREPARATION ###
##########################

### Loading libraries ###
library(rhdf5)
library(abind)
library(tidyverse)
library(here)

#load dstat_inputfiles function
source(here('project_scripts', 'Dstat_input_files_function.R'))


### Load data ###
#metadata and sample information
cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE)
samples_in_vcf <- read.csv(here('working_metadata_file', 'SAMPLES_biallel_polykazoutgroups_minq5_maxDP75_miss0.70.recode.txt'), header = FALSE)
entropy_samples_order <- read.delim(here('working_metadata_file', 'SAMPLES_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.txt'), header = FALSE)

#entropy results (K = 2)
mod_list <- list()
for (chain in 1:3) {
  mod_list[[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0('kazumbe_polyodon_qmod_6_18_2021_', 2, 'rep', chain, '.hdf5')), "q")
}



#######################
### DATA PROCESSING ###
#######################

# process metadata --> limit metadata to samples in the vcf
cichlid_metadata_subset_to_vcf <- cichlid_metadata[cichlid_metadata$sample_code %in% samples_in_vcf$V1,]

### process outgroup samples ###
#outgroup taxa: Petrochromis green Simochromis diagramma
outgroup_info_final <- cichlid_metadata_subset_to_vcf %>% 
  filter((species2 %in% c('green', 'diagramma')) &  region == 'Kigoma') %>%
  mutate(species_final_id = species2) %>% 
  select("sample_code", "location", "species_final_id")

metadata_subset_entropy <- cichlid_metadata[cichlid_metadata$sample_code %in% entropy_samples_order$V1, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2')]
metadata_subset_entropy_reorder <- metadata_subset_entropy[match(as.character(entropy_samples_order$V1), metadata_subset_entropy$sample_code),]


### Processing entropy results and combining with metadata 
k2_combined.q <- abind(mod_list$chain_1, mod_list$chain_2, mod_list$chain_3, along = 1)
k2_q_estimate <- as.data.frame(t(apply(k2_combined.q, 2:3, mean)))
k2_q_estimate_with_metadata <- cbind(k2_q_estimate, metadata_subset_entropy_reorder)
#V2 greater than 0.5 --> majority kazumbe ancestry
#V2 less than 0.5 --> majority kazumbe ancestry
#there are no individuals with exactly 50/50 kazumbe/polyodon ancestry
k2_q_estimate_with_metadata$species_final_id[k2_q_estimate_with_metadata$V2 > 0.5] <- 'kazumbe'
k2_q_estimate_with_metadata$species_final_id[k2_q_estimate_with_metadata$V2 < 0.5] <- 'polyodon'


#combine metadata with entropy results
main_taxa_metadata_final <- cichlid_metadata_subset_to_vcf %>% 
  filter(sample_code %in% k2_q_estimate_with_metadata$sample_code) %>%
  left_join(., k2_q_estimate_with_metadata[, c('sample_code', 'species_final_id')], by = 'sample_code') %>% 
  select("sample_code", "location", "species_final_id")

#combine dataframes about focal taxa and outgroups
final_info <- rbind(main_taxa_metadata_final, outgroup_info_final)

#check that final info has the same samples as the initial list of samples in the vcf
if (!identical(sort(final_info$sample_code), sort(samples_in_vcf$V1)))
  stop("PROCESSED DATA HAS DIFFERENT SET OF SAMPLES COMPARED TO INITIAL LIST OF SAMPLES IN VCF")



########################################################
### GENERATING INPUT FILES FOR D-STATS IN ADMIXTOOLS ###
########################################################

### type a: all combinations of population-level taxa ###
dstat_type_a <- dstat_inputfiles(pop_info = final_info, 
                                 pop_column = 'location', 
                                 species_column = 'species_final_id', 
                                 sampleID_column = 'sample_code', 
                                 species1 = 'kazumbe', 
                                 species2 = 'polyodon', 
                                 outgroup_samples = final_info %>% 
                                   filter(species_final_id %in% c('green', 'diagramma')) %>% 
                                   pull(sample_code), 
                                 minimum_sample_size = 2,
                                 type = "a")

#export ind and pop files
write.table(dstat_type_a$popfile, 
            here('DStats', 'input_files', 'files_9_18_2021', 'type_a', 'popfile_TYPEA_colormorphs.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(dstat_type_a$indfile, 
            here('DStats', 'input_files', 'files_9_18_2021', 'type_a', 'indfile_TYPEA_colormorphs.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


### type b --kazumbe: each kazumbe population compared to the full polyodon dataset ###
dstat_type_b_kazumbe <- dstat_inputfiles(pop_info = final_info, 
                                         pop_column = 'location', 
                                         species_column = 'species_final_id', 
                                         sampleID_column = 'sample_code', 
                                         species1 = 'kazumbe', 
                                         species2 = 'polyodon', 
                                         outgroup_samples = final_info %>% 
                                           filter(species_final_id %in% c('green', 'diagramma')) %>% 
                                           pull(sample_code), 
                                         minimum_sample_size = 2,
                                         type = "b")

#export ind and pop files
write.table(dstat_type_b_kazumbe$popfile, 
            here('DStats', 'input_files', 'files_9_18_2021', 'type_b_kazumbe', 'popfile_TYPEBKazumbe_colormorphs.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(dstat_type_b_kazumbe$indfile, 
            here('DStats', 'input_files', 'files_9_18_2021', 'type_b_kazumbe', 'indfile_TYPEBKazumbe_colormorphs.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


### type b --polyodon: each polyodon population compared to the full kazumbe dataset ###
dstat_type_b_polyodon <- dstat_inputfiles(pop_info = final_info, 
                                          pop_column = 'location', 
                                          species_column = 'species_final_id', 
                                          sampleID_column = 'sample_code', 
                                          species1 = 'polyodon', 
                                          species2 = 'kazumbe', 
                                          outgroup_samples = final_info %>% 
                                            filter(species_final_id %in% c('green', 'diagramma')) %>% 
                                            pull(sample_code), 
                                          minimum_sample_size = 2,
                                          type = "b")

#export ind and pop files
write.table(dstat_type_b_polyodon$popfile, 
            here('DStats', 'input_files', 'files_9_18_2021', 'type_b_polyodon', 'popfile_TYPEBPolyodon_colormorphs.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(dstat_type_b_polyodon$indfile, 
            here('DStats', 'input_files', 'files_9_18_2021', 'type_b_polyodon', 'indfile_TYPEBPolyodon_colormorphs.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

