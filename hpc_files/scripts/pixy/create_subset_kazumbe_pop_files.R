#!/usr/bin/env Rscript

#####################
### SCRIPT SET-UP ###
#####################

### LOAD LIBRARIES ###
library(dplyr)

number_subset_datasets <- 5

### LOAD POP FILES ###
pop_subset <- read.table('population_subset_7_1_2021.txt', header = FALSE, stringsAsFactors = FALSE)
colnames(pop_subset) <- c('code', 'pop')

full_species <- read.table('full_species_7_1_2021.txt', header = FALSE, stringsAsFactors = FALSE)
colnames(full_species) <- c('code', 'species')



###################################################
### CREATE AND EXPORT LOCATION-SPECIFIC SUBSETS ###
###################################################

pop_subset <- pop_subset %>% 
    mutate(species = gsub("_[A-Za-z]*", "", pop), 
           location = gsub("^[A-Za-z]*_", "", pop))

for (i in seq_len(number_subset_datasets)) {
  
  freq_info <- pop_subset %>% 
    group_by(pop) %>% 
    summarize(freq = n(),
              location = first(location), 
              .groups = 'drop') %>% 
    group_by(location) %>% 
    filter(n() >= 2) %>%
    filter(freq == min(freq)) %>% 
    slice(n = 1) %>% 
    ungroup() %>% 
    select(location, freq)
  
  freq_info_list <- split(freq_info, freq_info$location)
  
  
  subset_df <- lapply(freq_info_list,  function(INFO, pop_subset) {
    
    pop_subset %>% 
      filter(location == INFO$location) %>% 
      group_by(species) %>% 
      sample_n(INFO$freq)
    
  }, pop_subset = pop_subset) %>% 
    bind_rows() %>% 
    as.data.frame() %>% 
    filter(species == 'kazumbe')
  
  pop_subset[, 'SUBSET'] <- 'ignore'
  pop_subset[, 'SUBSET'][pop_subset$code %in% subset_df$code] <- pop_subset$pop[pop_subset$code %in% subset_df$code]
  

  write.table(pop_subset[, c('code', 'SUBSET')], 
              file = paste0('pop_subset_kazumbe_subset_popfile_dir/population_subset_7_1_2021_SUBSAMP_', i, '.txt'), 
              row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
}



###############################################
### CREATE AND EXPORT SPECIES-LEVEL SUBSETS ###
###############################################

ulombola_codes <- pop_subset %>%
  filter(location == 'Ulombola') %>%
  pull(code)


for (i in seq_len(number_subset_datasets)) {
  
  freq_info <- full_species %>% 
    group_by(species) %>% 
    summarize(freq = n(), .groups = 'drop') %>% 
    filter(freq == min(freq)) %>% 
    slice(n = 1) %>% 
    ungroup()
  
  subset_df <- full_species %>% 
    filter(!(code %in% ulombola_codes)) %>% #remove ulombola kazumbe to make kazumbe subset directly comparable to polyodon (samples from same extent)
    filter(species == 'kazumbe') %>%
    sample_n(freq_info$freq)
  
  
  full_species[, 'SUBSET'] <- 'ignore'
  full_species[, 'SUBSET'][full_species$code %in% subset_df$code] <- full_species$species[full_species$code %in% subset_df$code]
  
  full_species$species[full_species$code %in% subset_df$code]
  
  write.table(full_species[, c('code', 'SUBSET')], 
              file = paste0('full_species_kazumbe_subset_popfile_dir/full_species_7_1_2021_SUBSAMP_', i, '.txt'), 
              row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
}
