#####################################################################################################
### SCRIPT NAME: popgen_stats.R
### PURPOSE: Calculate and process population genetic summary statistics including Fst, dxy, and
###          nucleotide diversity (pi).
### PRODUCTS:
###     fst_7_2021.csv: file of interspecific Fst estimates and associated 95% CIs for calculations
###                     based on the full dataset and at each sampling location with sampling for
###                     both species
###     intraspecific_fst_10_2021.csv: file of intraspecific Fst estimates and associated 95% CIs,
###                                    which were calculated between all pairs of kazumbe
###                                    populations and between all pairs of polyodon populations
###     intraspecific_fst.png: visualization of the intraspecific Fst estimates, which is included
###                            in the supplementary material
###     pi_7_2021.csv: file of speicies-specific nucleotide diversity (pi) estimates and 95% CIs 
###                    for kazumbe and polyodon, which were calculated based on all samples for 
###                    each species and for the samples from each species collected from each 
###                    sampling location (i.e. population)
###     dxy_7_2021.csv: file of interspecific dxy estimates and associated 95% CIs calculated
###                     based on the full dataset and at each sampling location with sampling for
###                     both species
###     popgen_stats_table.txt: latex table summarizing the results from the interspecific 
###                             population genetic summary statistic calculations. This table
###                             is included in the supplementary material.
###     intraspecific_fst_table.txt: latex table summarizing the results from the intraspecific 
###                                  Fst calculations. This table is included in the supplementary 
###                                  material.
###################################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### Notes ###
# location     long       lat Region number_label
# 1  Gombe_South 29.60158 -4.746605  North            A
# 2   Katongwe_N 29.60270 -4.803450  North            B
# 3   Katongwe_S 29.60633 -4.824055  North            C
# 4          Ska 29.60493 -4.825766  North            D
# 5       Kalala 29.60854 -4.843925  North            E
# 6       Nondwa 29.60833 -4.860178  North            F
# 7      Hilltop 29.61296 -4.886808    Mid            G
# 8       Bangwe 29.61221 -4.894394    Mid            H
# 9        Jakob 29.59809 -4.910006    Mid            I
# 10    Ulombola 29.76905 -5.011201  South            J


### Load libraries ###
library(adegenet)
library(tidyverse)
library(vcfR)
library(rhdf5)
library(abind)
library(SNPRelate)
library(dartR)
library(kableExtra)
library(here)


### Analytical decisions ###
number_of_bootstraps <- 500 #number of bootstraps to use for calculating confidence intervals

### Custom functions ###

# Load Fst function
source(here('project_scripts', 'reich_fst.R'))

diversity_calc_alt <- function(input_data,
                               calc_type = c('raw', 'per_site'),
                               no_sites_col,
                               count_diffs_col,
                               count_comparisons_col) {
  
  if (calc_type == 'raw') {
    return( sum(count_diffs_col, na.rm = TRUE)/sum(count_comparisons_col, na.rm = TRUE) )
  } else if (calc_type == 'per_site') {
    return( (sum(count_diffs_col, na.rm = TRUE)/sum(count_comparisons_col, na.rm = TRUE))/sum(no_sites_col, na.rm = TRUE) )
  }
}


diversity_calc <- function(input_data,
                           calc_type = c('raw', 'per_site'),
                           no_sites_col,
                           count_diffs_col,
                           count_comparisons_col) {

  if (calc_type == 'raw') {
    return( sum(input_data[, count_diffs_col], na.rm = TRUE)/sum(input_data[, count_comparisons_col], na.rm = TRUE) )
  } else if (calc_type == 'per_site') {
    return( (sum(input_data[, count_diffs_col], na.rm = TRUE)/sum(input_data[, count_comparisons_col], na.rm = TRUE))/sum(input_data[, no_sites_col], na.rm = TRUE) )
  }
}


conduct_bootstrap <- function(dataset, bootstrap_number = 10) {
  bootstrap_name_vec <- paste0('bs', 1:bootstrap_number)
  
  bs_vec <- sapply(split(bootstrap_name_vec, factor(bootstrap_name_vec, levels = bootstrap_name_vec) ), function(x, bs_dat) {
    diversity_calc(input_data = bs_dat[sample(1:nrow(bs_dat), nrow(bs_dat), replace = TRUE),],
                   calc_type = c('raw', 'per_site')[1],
                   no_sites_col = 'no_sites',
                   count_diffs_col = 'count_diffs',
                   count_comparisons_col = 'count_comparisons')
  },
  bs_dat = dataset)
  
  return(data.frame(bootstrap = names(bs_vec),
                    val = bs_vec))
}



#################
### LOAD DATA ###
#################

cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata

#vcf file (this file was already processed with vcftools)
petro_vcfR <- read.vcfR(here('vcf', 'biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.vcf'))
genlight_polyodon_kazumbe <- vcfR2genlight(petro_vcfR)

### ENTROPY ###
entropy_samples_order <- read.delim(here('working_metadata_file', 'SAMPLES_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.txt'), header = FALSE) #order of samples in entropy results

mod_list <- list()
for (chain in 1:3) {
  mod_list[[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0('kazumbe_polyodon_qmod_6_18_2021_', 2, 'rep', chain, '.hdf5')), "q")
}



#######################
### PROCESSING DATA ###
#######################

### Processing entropy data and metadata ###
k2_combined.q <- abind(mod_list$chain_1, mod_list$chain_2, mod_list$chain_3, along = 1)
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

polyodon_kazumbe_freq_dataframe <- as.data.frame(table(k2_q_estimate_with_metadata$location, k2_q_estimate_with_metadata$species_ID_entropy)) %>%
  spread(key = 'Var2', value = 'Freq')
colnames(polyodon_kazumbe_freq_dataframe) <- c('location', 'Petrochromis_kazumbe', 'Petrochromis_polyodon')


### Add population info to genlight object ###
meta_poly_kazumbe <- k2_q_estimate_with_metadata[k2_q_estimate_with_metadata$sample_code %in% genlight_polyodon_kazumbe@ind.names,]
meta_poly_kazumbe_reorder <- meta_poly_kazumbe[match(genlight_polyodon_kazumbe@ind.names, meta_poly_kazumbe$sample_code),]

#https://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html
#pop(genlight_polyodon_kazumbe) <- as.factor(gsub('. ', '_', meta_poly_kazumbe_reorder$species_ID_entropy) )
pop(genlight_polyodon_kazumbe) <- gsub('. ', '_', meta_poly_kazumbe_reorder$species_ID_entropy)

genlight_polyodon_kazumbe <- gl.compliance.check(genlight_polyodon_kazumbe)



#####################################
### CALCULATING INTERSPECIFIC FST ###
#####################################
sympatric_sites <- as.character(polyodon_kazumbe_freq_dataframe[which(polyodon_kazumbe_freq_dataframe$Petrochromis_kazumbe > 1 & polyodon_kazumbe_freq_dataframe$Petrochromis_polyodon > 1),]$location)

comparisons <- c('full_species', sympatric_sites)
Fst_dataframe <- data.frame(pop = comparisons, Fst = NA, lower_CI = NA, upper_CI = NA)

species_level_fst <- reich.fst(gl = genlight_polyodon_kazumbe, bootstrap = number_of_bootstraps, verbose = TRUE)


Fst_dataframe[which(Fst_dataframe$pop == 'full_species'), 'Fst'] <- species_level_fst$bootstraps[names(species_level_fst$bootstraps) == 'fst_estimate' ]
Fst_dataframe[which(Fst_dataframe$pop == 'full_species'), 'lower_CI'] <- species_level_fst$bootstraps[names(species_level_fst$bootstraps) == 'min_CI' ]
Fst_dataframe[which(Fst_dataframe$pop == 'full_species'), 'upper_CI'] <- species_level_fst$bootstraps[names(species_level_fst$bootstraps) == 'max_CI' ]

index <- 1
Fst_calc_progress_bar <- txtProgressBar(min = 0, max = length(Fst_dataframe$pop[Fst_dataframe$pop != 'full_species']), style = 3)

for (location in as.character(Fst_dataframe$pop[Fst_dataframe$pop != 'full_species'])) {
  location_samples <- meta_poly_kazumbe_reorder[which(meta_poly_kazumbe_reorder$location == location),]$sample_code
  
  genlight_object_subset <- genlight_polyodon_kazumbe[(genlight_polyodon_kazumbe@ind.names %in% location_samples)]
  
  pop_level_fst <- reich.fst(gl = genlight_object_subset, bootstrap = number_of_bootstraps, verbose = TRUE, subset_option = c("subset", "function")[1])
  
  Fst_dataframe[which(Fst_dataframe$pop == location), 'Fst'] <- pop_level_fst$bootstraps[names(pop_level_fst$bootstraps) == 'fst_estimate' ]
  Fst_dataframe[which(Fst_dataframe$pop == location), 'lower_CI'] <- pop_level_fst$bootstraps[names(pop_level_fst$bootstraps) == 'min_CI' ]
  Fst_dataframe[which(Fst_dataframe$pop == location), 'upper_CI'] <- pop_level_fst$bootstraps[names(pop_level_fst$bootstraps) == 'max_CI' ]
  
  setTxtProgressBar(Fst_calc_progress_bar, index)
  index <- index + 1
}

Fst_dataframe$Fst <- as.numeric(Fst_dataframe$Fst)


Fst_dataframe_updated <- Fst_dataframe %>%
  mutate(pop = replace(pop, pop == 'full_species', 'full'),
          population = factor(pop,
                             levels = rev(c('Gombe_South', 'Katongwe_N', 'Katongwe_S', 'Ska', 'Kalala', 'Nondwa', 'Hilltop', 'Bangwe', 'Jakob', 'Ulombola', 'full'))),
         Population = factor(case_when(population == 'Gombe_South' ~ 'A',
                                       population == 'Katongwe_N' ~ 'B',
                                       population == 'Katongwe_S' ~ 'C',
                                       population == 'Ska' ~ 'D',
                                       population == 'Kalala' ~ 'E',
                                       population == 'Nondwa' ~ 'F',
                                       population == 'Hilltop' ~ 'G',
                                       population == 'Bangwe' ~ 'H',
                                       population == 'Jakob' ~ 'I',
                                       population == 'Ulombola' ~ 'J',
                                       population == 'full' ~ 'Full'),
                             levels = (c('Full', LETTERS[1:10]) )))

write.csv(Fst_dataframe_updated, here('diversity_divergence_stats', 'processed_results', 'fst_7_2021.csv'), row.names = FALSE)



#####################################
### CALCULATING INTRASPECIFIC FST ###
#####################################

intraspecific_fst_list <- list()

for (SPECIES in c('P_kazumbe', 'P_polyodon')) {
  
  genlight_single_species <- genlight_polyodon_kazumbe[(genlight_polyodon_kazumbe@pop == SPECIES)]
  metadata_single_species_subset <- meta_poly_kazumbe_reorder[meta_poly_kazumbe_reorder$sample_code %in% genlight_single_species@ind.names,]
  pop(genlight_single_species) <- metadata_single_species_subset[match(genlight_single_species@ind.names, metadata_single_species_subset$sample_code),]$location
  genlight_single_species <- gl.compliance.check(genlight_single_species)
  
  species_col <- switch(SPECIES, P_kazumbe = 'Petrochromis_kazumbe', P_polyodon = 'Petrochromis_polyodon')
  fst_dataframe_single_species <- as.data.frame(t(combn(x = as.character(polyodon_kazumbe_freq_dataframe[polyodon_kazumbe_freq_dataframe[,species_col] > 1,]$location), 
                                                        m = 2)))
  
  fst_dataframe_single_species <- fst_dataframe_single_species %>% 
    rename(pop1 = V1,
           pop2 = V2) %>% 
    mutate(species = SPECIES, Fst = NA, lower_CI = NA, upper_CI = NA)
  
  index <- 1
  Fst_calc_progress_bar <- txtProgressBar(min = 0, max = nrow(fst_dataframe_single_species), style = 3)
  
  for (ROW in seq_len(nrow(fst_dataframe_single_species))) {
    species_name <- switch(SPECIES, P_kazumbe = 'P. kazumbe', P_polyodon = 'P. polyodon')
    pop1_samples <- meta_poly_kazumbe_reorder[meta_poly_kazumbe_reorder$location == fst_dataframe_single_species[ROW,]$pop1 & meta_poly_kazumbe_reorder$species_ID_entropy == species_name,]$sample_code
    pop2_samples <- meta_poly_kazumbe_reorder[meta_poly_kazumbe_reorder$location == fst_dataframe_single_species[ROW,]$pop2 & meta_poly_kazumbe_reorder$species_ID_entropy == species_name,]$sample_code
    
    genlight_object_subset <- genlight_single_species[(genlight_single_species@ind.names %in% c(pop1_samples, pop2_samples))]
    
    pop_compare_fst <- reich.fst(gl = genlight_object_subset, bootstrap = number_of_bootstraps, verbose = TRUE, subset_option = c("subset", "function")[1])
    genlight_object_subset$other$loc.metrics.flags$monomorphs
    genlight_single_species$other$loc.metrics.flags$monomorphs
    
    fst_dataframe_single_species[ROW, 'Fst'] <- pop_compare_fst$bootstraps[names(pop_compare_fst$bootstraps) == 'fst_estimate']
    fst_dataframe_single_species[ROW, 'lower_CI'] <- pop_compare_fst$bootstraps[names(pop_compare_fst$bootstraps) == 'min_CI']
    fst_dataframe_single_species[ROW, 'upper_CI'] <- pop_compare_fst$bootstraps[names(pop_compare_fst$bootstraps) == 'max_CI']
    
    setTxtProgressBar(Fst_calc_progress_bar, index)
    index <- index + 1
  }
  
  fst_dataframe_single_species$Fst <- as.numeric(fst_dataframe_single_species$Fst)
  
  intraspecific_fst_list[[SPECIES]] <- list(results = fst_dataframe_single_species,
                                            single_species_genlight = genlight_single_species)
  
  message('Finished ', SPECIES)
}

intraspecific_fst_df <- do.call(rbind, lapply(intraspecific_fst_list, function(x) x[['results']]))
rownames(intraspecific_fst_df) <- NULL

#save intraspecific results
write.csv(intraspecific_fst_df, here('diversity_divergence_stats', 'processed_results', 'intraspecific_fst_10_2021.csv'), row.names = FALSE)



##################################
### PLOTTING INTRASPECIFIC FST ###
##################################

if(!exists("intraspecific_fst_df")) intraspecific_fst_df <- read.csv(here('diversity_divergence_stats', 'processed_results', 'intraspecific_fst_10_2021.csv'))

intraspecific_fst_df %>%
  mutate(pop1_factor_letter = factor(case_when(pop1 == 'Gombe_South' ~ "A",
                                               pop1 == 'Katongwe_N' ~ "B",
                                               pop1 == 'Katongwe_S' ~ "C",
                                               pop1 == 'Ska' ~ "D",
                                               pop1 == 'Kalala' ~ "E",
                                               pop1 == 'Nondwa' ~ "F",
                                               pop1 == 'Hilltop' ~ "G",
                                               pop1 == 'Bangwe' ~ "H",
                                               pop1 == 'Jakob' ~ "I",
                                               pop1 == 'Ulombola' ~ "J"),
                                     levels= LETTERS[10:1]),
         pop2_factor_letter = factor(case_when(pop2 == 'Gombe_South' ~ "A",
                                               pop2 == 'Katongwe_N' ~ "B",
                                               pop2 == 'Katongwe_S' ~ "C",
                                               pop2 == 'Ska' ~ "D",
                                               pop2 == 'Kalala' ~ "E",
                                               pop2 == 'Nondwa' ~ "F",
                                               pop2 == 'Hilltop' ~ "G",
                                               pop2 == 'Bangwe' ~ "H",
                                               pop2 == 'Jakob' ~ "I",
                                               pop2 == 'Ulombola' ~ "J"),
                                     levels= LETTERS[10:1]),
         Species = case_when(species == 'P_kazumbe' ~ "P. sp. 'kazumbe'",
                             species == 'P_polyodon' ~ "P. cf. polyodon")) %>% 
  rowwise() %>% 
  mutate(pop1_factor_letter_ordered = if_else(which(LETTERS == pop1_factor_letter) < which(LETTERS == pop2_factor_letter), pop1_factor_letter, pop2_factor_letter),
         pop2_factor_letter_ordered = if_else(which(LETTERS == pop1_factor_letter) < which(LETTERS == pop2_factor_letter), pop2_factor_letter, pop1_factor_letter)) %>%
  ggplot(aes(x = pop1_factor_letter_ordered, y = pop2_factor_letter_ordered)) +
  facet_grid(cols = vars(Species)) +
  geom_raster(aes(fill = Fst)) +
  scale_fill_gradient(low = "#f3f3f6", high = "#c3352b",
                      breaks=c(0, 0.05, 0.1, 0.15),
                      limits=c(0, 0.15)) +
  xlab('Population') + ylab('Population') +
  theme_bw() +
  theme(strip.text.x = element_text(size = 15),
        axis.title = element_text(size = 15))

#ggsave(here('figures', 'intraspecific_fst.png'), 
#       width = 14, height = 7, bg = "white")
ggsave(here('figures', 'intraspecific_fst.pdf'), 
       width = 14, height = 7, bg = "white", device = 'pdf')


##EXPLORATORY VIOLIN PLOT FIGURE SHOWING JUMPS IN FST BETWEEN SAMPLING REGIONS
# intraspecific_fst_df_hist <- intraspecific_fst_df %>%
#   mutate(pop1_factor_letter = factor(case_when(pop1 == 'Gombe_South' ~ "A",
#                                                pop1 == 'Katongwe_N' ~ "B",
#                                                pop1 == 'Katongwe_S' ~ "C",
#                                                pop1 == 'Ska' ~ "D",
#                                                pop1 == 'Kalala' ~ "E",
#                                                pop1 == 'Nondwa' ~ "F",
#                                                pop1 == 'Hilltop' ~ "G",
#                                                pop1 == 'Bangwe' ~ "H",
#                                                pop1 == 'Jakob' ~ "I",
#                                                pop1 == 'Ulombola' ~ "J"),
#                                      levels= LETTERS[10:1]),
#          pop2_factor_letter = factor(case_when(pop2 == 'Gombe_South' ~ "A",
#                                                pop2 == 'Katongwe_N' ~ "B",
#                                                pop2 == 'Katongwe_S' ~ "C",
#                                                pop2 == 'Ska' ~ "D",
#                                                pop2 == 'Kalala' ~ "E",
#                                                pop2 == 'Nondwa' ~ "F",
#                                                pop2 == 'Hilltop' ~ "G",
#                                                pop2 == 'Bangwe' ~ "H",
#                                                pop2 == 'Jakob' ~ "I",
#                                                pop2 == 'Ulombola' ~ "J"),
#                                      levels= LETTERS[10:1]),
#          Species = case_when(species == 'P_kazumbe' ~ "P. sp. 'kazumbe'",
#                              species == 'P_polyodon' ~ "P. cf. polyodon")) %>% 
#   rowwise() %>% 
#   mutate(pop1_factor_letter_ordered = if_else(which(LETTERS == pop1_factor_letter) < which(LETTERS == pop2_factor_letter), pop1_factor_letter, pop2_factor_letter),
#          pop2_factor_letter_ordered = if_else(which(LETTERS == pop1_factor_letter) < which(LETTERS == pop2_factor_letter), pop2_factor_letter, pop1_factor_letter))
# 
# intraspecific_fst_df_hist %>% 
#   #rowwise() %>%
#   mutate(pop_compare = case_when(pop1_factor_letter_ordered %in% c('A', 'B', 'C', 'D', 'E', 'F') & pop2_factor_letter_ordered %in% c('A', 'B', 'C', 'D', 'E', 'F') ~ 'north_north',
#                    pop1_factor_letter_ordered %in% c('G', 'H', 'I') & pop2_factor_letter_ordered %in% c('G', 'H', 'I') ~ 'mid_mid',
#                    (pop1_factor_letter_ordered %in% c('A', 'B', 'C', 'D', 'E', 'F') | pop2_factor_letter_ordered %in% c('A', 'B', 'C', 'D', 'E', 'F')) & (pop1_factor_letter_ordered %in% c('G', 'H', 'I') | pop2_factor_letter_ordered %in% c('G', 'H', 'I')) ~ 'north_mid',
#                    (pop1_factor_letter_ordered %in% c('A', 'B', 'C', 'D', 'E', 'F') | pop2_factor_letter_ordered %in% c('A', 'B', 'C', 'D', 'E', 'F')) & (pop1_factor_letter_ordered %in% c('J') | pop2_factor_letter_ordered %in% c('J')) ~ 'north_south',
#                    (pop1_factor_letter_ordered %in% c('G', 'H', 'I') | pop2_factor_letter_ordered %in% c('G', 'H', 'I')) | (pop2_factor_letter_ordered %in% c('J')) ~ 'north_mid')) %>% 
#   ggplot() + 
#   geom_violin(aes(x = pop_compare, y = Fst))



############################################
### DIVERGENCE/DIVERSITY STATS FROM PIXY ###
############################################

#loading pixy data
pixy_file_list <- list.files(here('diversity_divergence_stats', 'pixy', 'results'))[grepl(pattern = ".*7_1_2021.*", list.files(here('diversity_divergence_stats', 'pixy', 'results')))]
results_file <- lapply(split(pixy_file_list, gsub(".txt", "", pixy_file_list)), function(x) {
 read.delim(paste0(here('diversity_divergence_stats', 'pixy', 'results'), '/', x))
})


### Initial processing of pixy data ###
results_file$output_full_species_7_1_2021_pi <- results_file$output_full_species_7_1_2021_pi %>% 
  mutate(population = 'full',
         species = pop)

results_file$output_population_subset_7_1_2021_pi <- results_file$output_population_subset_7_1_2021_pi %>% 
  mutate(population = str_remove(pop, '[A-Za-z]+_'),
         species = str_remove(pop, '_[A-Za-z]+.*'))


results_file$output_full_species_7_1_2021_dxy <- results_file$output_full_species_7_1_2021_dxy %>% 
  mutate(pop_pop1 = 'full',
         species_pop1 = pop1,
         pop_pop2 = 'full',
         species_pop2 = pop2,
         comparison_id = 'full')

results_file$output_population_subset_7_1_2021_dxy <- results_file$output_population_subset_7_1_2021_dxy %>%
  mutate(pop_pop1 = str_remove(pop1, '[A-Za-z]+_'),
         species_pop1 = str_remove(pop1, '_[A-Za-z]+.*'),
         pop_pop2 = str_remove(pop2, '[A-Za-z]+_'),
         species_pop2 = str_remove(pop2, '_[A-Za-z]+.*'),
         comparison_id = paste(pop1, pop2, sep = '_'))

output_population_subset_7_1_2021_dxy_subset <- results_file$output_population_subset_7_1_2021_dxy %>%
  filter(pop_pop1 == pop_pop2 & species_pop1 != species_pop2)

#save pi and dxy in a list
diversity_list <- list()
diversity_list[['pi']] <- rbind(results_file$output_full_species_7_1_2021_pi,
                                results_file$output_population_subset_7_1_2021_pi)

diversity_list[['dxy']] <- rbind(results_file$output_full_species_7_1_2021_dxy,
                                 output_population_subset_7_1_2021_dxy_subset)


### Processing and calculating bootstraps for pi ###
pi_bootstrap_list <- list()
pi_bootstrap_process_list <- list()
for (taxa in unique(diversity_list$pi$pop)) {
  pi_bootstrap_list[[taxa]] <- conduct_bootstrap(dataset = diversity_list$pi[which(diversity_list$pi$pop == taxa),], bootstrap_number = number_of_bootstraps)
  pi_bootstrap_list[[taxa]]$pop <- taxa
  
  pi_bootstrap_process_list[[taxa]] <- quantile(pi_bootstrap_list[[taxa]]$val, probs = c(0.025, 0.975))
  
  message('finished calculations for ', taxa)
}

pi_info <- diversity_list$pi %>% 
  group_by(pop) %>%
  summarize(mean = diversity_calc_alt(calc_type = c('raw', 'per_site')[1],
                                       no_sites_col = no_sites,
                                       count_diffs_col = count_diffs,
                                       count_comparisons_col = count_comparisons), .groups = 'drop') %>% 
  merge(., do.call(rbind, pi_bootstrap_process_list), by.x = 'pop', by.y = 0) %>% 
  rename(lower = "2.5%",
         upper = "97.5%") %>% 
  mutate(population = factor(case_when(str_detect(pop, '_') ~ str_remove(pop, '[A-Za-z]+_'),
                                str_detect(pop, '^[A-Za-z]+$') ~ 'full'),
                             levels = rev(c('Gombe_South', 'Katongwe_N', 'Katongwe_S', 'Ska', 'Kalala', 'Nondwa', 'Hilltop', 'Bangwe', 'Jakob', 'Ulombola', 'full'))),
         species = case_when(str_detect(pop, '_') ~ str_remove(pop, '_[A-Za-z]+.*'),
                             str_detect(pop, '^[A-Za-z]+$') ~ pop),
         Population = factor(case_when(population == 'Gombe_South' ~ 'A',
                                       population == 'Katongwe_N' ~ 'B',
                                       population == 'Katongwe_S' ~ 'C',
                                       population == 'Ska' ~ 'D',
                                       population == 'Kalala' ~ 'E',
                                       population == 'Nondwa' ~ 'F',
                                       population == 'Hilltop' ~ 'G',
                                       population == 'Bangwe' ~ 'H',
                                       population == 'Jakob' ~ 'I',
                                       population == 'Ulombola' ~ 'J',
                                       population == 'full' ~ 'Full'),
                             levels = (c('Full', LETTERS[10:1]) )))

#save processed pi results
write.csv(pi_info, here('diversity_divergence_stats', 'processed_results', 'pi_7_2021.csv'), row.names = FALSE)


### Processing and calculating bootstraps for dxy ###
dxy_bootstrap_list <- list()
dxy_bootstrap_process_list <- list()
for (pair in unique(diversity_list$dxy$comparison_id)) {
  dxy_bootstrap_list[[pair]] <- conduct_bootstrap(dataset = diversity_list$dxy[which(diversity_list$dxy$comparison_id == pair),], bootstrap_number = number_of_bootstraps)
  dxy_bootstrap_list[[pair]]$pop_pair <- pair
  
  dxy_bootstrap_process_list[[pair]] <- quantile(dxy_bootstrap_list[[pair]]$val, probs = c(0.025, 0.975))
  
  message('finished calculations for ', pair)
}

dxy_info <- diversity_list$dxy %>% 
  group_by(comparison_id) %>%
  summarize(mean =  diversity_calc_alt(calc_type = c('raw', 'per_site')[1],
                                       no_sites_col = no_sites,
                                       count_diffs_col = count_diffs,
                                       count_comparisons_col = count_comparisons), .groups = 'drop') %>% 
  merge(., do.call(rbind, dxy_bootstrap_process_list), by.x = 'comparison_id', by.y = 0) %>% 
  rename(lower = "2.5%",
         upper = "97.5%") %>% 
  #mutate(population = factor(str_remove(comparison_id, '[A-Za-z]*.*polyodon_'),
  mutate(population = factor(str_remove(comparison_id, '[A-Za-z]*.*(polyodon|kazumbe)_'),
                             levels = rev(c('Gombe_South', 'Katongwe_N', 'Katongwe_S', 'Ska', 'Kalala', 'Nondwa', 'Hilltop', 'Bangwe', 'Jakob', 'Ulombola', 'full'))),
         Population = factor(case_when(population == 'Gombe_South' ~ 'A',
                                       population == 'Katongwe_N' ~ 'B',
                                       population == 'Katongwe_S' ~ 'C',
                                       population == 'Ska' ~ 'D',
                                       population == 'Kalala' ~ 'E',
                                       population == 'Nondwa' ~ 'F',
                                       population == 'Hilltop' ~ 'G',
                                       population == 'Bangwe' ~ 'H',
                                       population == 'Jakob' ~ 'I',
                                       population == 'Ulombola' ~ 'J',
                                       population == 'full' ~ 'Full'),
                             levels = (c( 'Full', LETTERS[10:1]) )))

#save processed dxy results
write.csv(dxy_info, here('diversity_divergence_stats', 'processed_results', 'dxy_7_2021.csv'), row.names = FALSE)



###################################
### SUMMARY TABLE OF STATISTICS ###
###################################

### Loading data (if required) ###
if (!exists("cichlid_metadata")) cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata

if (!exists("Fst_dataframe_updated")) Fst_dataframe_updated <- read.csv(here('diversity_divergence_stats', 'processed_results', 'fst_7_2021.csv'))
if(!exists("pi_info")) pi_info <- read.csv(here('diversity_divergence_stats', 'processed_results', 'pi_7_2021.csv'))
if(!exists("dxy_info")) dxy_info <- read.csv(here('diversity_divergence_stats', 'processed_results', 'dxy_7_2021.csv'))
if(!exists("intraspecific_fst_df")) intraspecific_fst_df <- read.csv(here('diversity_divergence_stats', 'processed_results', 'intraspecific_fst_10_2021.csv'))


### Creating pop gen stats table (interspecific Fst, dxy, pi) ###
Fst_dataframe_updated_processed <- Fst_dataframe_updated %>% 
  mutate(Fst = paste0(round(Fst, 4), " (", round(lower_CI, 4), ", ", round(upper_CI, 4), ")")) %>% 
  select(Fst, Population)

pi_info_processed <- pi_info %>% 
  mutate(pi = paste0(round(mean, 4), " (", round(lower, 4), ", ", round(upper, 4), ")")) %>%
  rename(Species = species) %>%
  select(pi, Population, Species) %>% 
  pivot_wider(names_from = Species, values_from = pi, values_fill = "not calculated")

dxy_info_processed <- dxy_info %>% 
  mutate(dxy = paste0(round(mean, 4), " (", round(lower, 4), ", ", round(upper, 4), ")")) %>%
  select(dxy, Population)

#7/15/2022: order rows so that locations are N --> S (alphabetical order)
popgen_stats <- left_join(Fst_dataframe_updated_processed, dxy_info_processed, by = "Population") %>% 
  left_join(., pi_info_processed, by = "Population") %>% 
  relocate(Population, Fst, dxy, kazumbe, polyodon) %>% 
  mutate(Population = factor(Population, levels = (c(LETTERS[1:10], 'Full') )) ) %>% 
  arrange(Population) %>% 
  rename("P. sp. 'kazumbe'" = kazumbe,
         "P. cf. polyodon" = polyodon)

popgen_stats_table <- kbl(popgen_stats, "latex", booktabs = TRUE, align = "c", linesep = "") %>%
  add_header_above(c(" ", " ", " ", "$\\\\pi$" = 2), escape = FALSE) %>% 
  kable_styling(latex_options = c("scale_down"))

cat(popgen_stats_table, file = here('tables', 'popgen_stats_table.txt'), append = FALSE)


### Creating intraspecific Fst table ###
intraspecific_fst_df_processed <- intraspecific_fst_df %>%
  mutate(Pop1 = factor(case_when(pop1 == 'Gombe_South' ~ "A",
                                               pop1 == 'Katongwe_N' ~ "B",
                                               pop1 == 'Katongwe_S' ~ "C",
                                               pop1 == 'Ska' ~ "D",
                                               pop1 == 'Kalala' ~ "E",
                                               pop1 == 'Nondwa' ~ "F",
                                               pop1 == 'Hilltop' ~ "G",
                                               pop1 == 'Bangwe' ~ "H",
                                               pop1 == 'Jakob' ~ "I",
                                               pop1 == 'Ulombola' ~ "J"),
                                     levels= LETTERS[10:1]),
         Pop2 = factor(case_when(pop2 == 'Gombe_South' ~ "A",
                                               pop2 == 'Katongwe_N' ~ "B",
                                               pop2 == 'Katongwe_S' ~ "C",
                                               pop2 == 'Ska' ~ "D",
                                               pop2 == 'Kalala' ~ "E",
                                               pop2 == 'Nondwa' ~ "F",
                                               pop2 == 'Hilltop' ~ "G",
                                               pop2 == 'Bangwe' ~ "H",
                                               pop2 == 'Jakob' ~ "I",
                                               pop2 == 'Ulombola' ~ "J"),
                                     levels= LETTERS[10:1]),
         Species = case_when(species == 'P_kazumbe' ~ "P. sp. 'kazumbe'",
                             species == 'P_polyodon' ~ "P. cf. polyodon"),
         Fst = if_else(Fst < 0, 0, Fst),
         lower_CI = if_else(lower_CI < 0, 0, lower_CI),
         upper_CI = if_else(upper_CI < 0, 0, upper_CI),
         Fst = paste0(round(Fst, 4), " (", round(lower_CI, 4), ", ", round(upper_CI, 4), ")")) %>%
  select(Species, Pop1, Pop2, Fst) %>% 
  pivot_wider(values_from = Fst, names_from = Species) %>% 
  mutate(`P. cf. polyodon` = if_else(is.na(`P. cf. polyodon`), " ", `P. cf. polyodon`))

intraspecific_fst_table <- kbl(intraspecific_fst_df_processed, "latex", booktabs = TRUE, align = "c", linesep = "") %>%
  kable_styling(latex_options = c("hold_position")) %>% 
  add_header_above(c(" ", " ", "Fst" = 2))

cat(intraspecific_fst_table, file = here('tables', 'intraspecific_fst_table.txt'), append = FALSE)

#identifying max and min fst values in kazumbe and polyodon (reported in results)
intraspecific_fst_df %>% 
  group_by(species) %>% 
  summarize(max = max(Fst),
            min = min(Fst)) %>% 
  as.data.frame()



###################################
### SUMMARY TABLE OF STATISTICS ###
###################################

library(tidyverse)
library(here)

number_of_bootstraps <- 500

pi_list_kazumbe_subset <- list()

### load full dataset diversity estimates for subset kazumbe ###
full_species_kazumbe_subset_file_names <- list.files(here('diversity_divergence_stats', 'pixy', 'results', 'fullspecies_subset_kazumbe'), full.names = TRUE)
full_species_kazumbe_subset_file_names <- full_species_kazumbe_subset_file_names[grepl(pattern = 'pi.txt$', x = full_species_kazumbe_subset_file_names)]
pi_list_kazumbe_subset[['full_species']] <- lapply(split(full_species_kazumbe_subset_file_names, gsub(".txt", "", basename(full_species_kazumbe_subset_file_names) )), read.delim)

### load location-level dataset diversity estimates for subset kazumbe ###
location_kazumbe_subset_file_names <- list.files(here('diversity_divergence_stats', 'pixy', 'results', 'pop_subset_kazumbe'), full.names = TRUE)
pi_list_kazumbe_subset[['location']] <- lapply(split(location_kazumbe_subset_file_names, gsub(".txt", "", basename(location_kazumbe_subset_file_names) )), read.delim)


### initial processing of nuc. div. data (remove the calculations that are being ignored)
pi_list_kazumbe_subset_processed <- lapply(pi_list_kazumbe_subset, function(sub_list) {
  lapply(sub_list, function(x) x[x$pop != 'ignore',])
})



pi_bootstrap_full_list <- list()
for (calc_type in names(pi_list_kazumbe_subset_processed)) {
  for (replicate in names(pi_list_kazumbe_subset_processed[[calc_type]])) {
    
    replicate_label <- paste0('replicate_', gsub('^.*SUBSAMP_([0-9]+)_pi', '\\1', replicate ))
    
    pi_bootstrap_list <- list()
    pi_bootstrap_process_list <- list()
    for (taxa in unique(pi_list_kazumbe_subset_processed[[calc_type]][[replicate]]$pop)) {
      pi_bootstrap_list[[taxa]] <- conduct_bootstrap(dataset = pi_list_kazumbe_subset_processed[[calc_type]][[replicate]][which(pi_list_kazumbe_subset_processed[[calc_type]][[replicate]]$pop == taxa),], bootstrap_number = 15)
      pi_bootstrap_list[[taxa]]$pop <- taxa
      
      pi_bootstrap_process_list[[taxa]] <- quantile(pi_bootstrap_list[[taxa]]$val, probs = c(0.025, 0.975))
      
      message('finished calculations for ', taxa)
    }
    
    pi_bootstrap_full_list[[calc_type]][[replicate]] <- pi_list_kazumbe_subset_processed[[calc_type]][[replicate]] %>% 
      group_by(pop) %>%
      summarize(mean = diversity_calc_alt(calc_type = c('raw', 'per_site')[1],
                                          no_sites_col = no_sites,
                                          count_diffs_col = count_diffs,
                                          count_comparisons_col = count_comparisons), .groups = 'drop') %>% 
      merge(., do.call(rbind, pi_bootstrap_process_list), by.x = 'pop', by.y = 0) %>% 
      rename(lower = "2.5%",
             upper = "97.5%") %>% 
      mutate(population = factor(case_when(str_detect(pop, '_') ~ str_remove(pop, '[A-Za-z]+_'),
                                           str_detect(pop, '^[A-Za-z]+$') ~ 'full'),
                                 levels = rev(c('Gombe_South', 'Katongwe_N', 'Katongwe_S', 'Ska', 'Kalala', 'Nondwa', 'Hilltop', 'Bangwe', 'Jakob', 'Ulombola', 'full'))),
             species = case_when(str_detect(pop, '_') ~ str_remove(pop, '_[A-Za-z]+.*'),
                                 str_detect(pop, '^[A-Za-z]+$') ~ pop),
             Population = factor(case_when(population == 'Gombe_South' ~ 'A',
                                           population == 'Katongwe_N' ~ 'B',
                                           population == 'Katongwe_S' ~ 'C',
                                           population == 'Ska' ~ 'D',
                                           population == 'Kalala' ~ 'E',
                                           population == 'Nondwa' ~ 'F',
                                           population == 'Hilltop' ~ 'G',
                                           population == 'Bangwe' ~ 'H',
                                           population == 'Jakob' ~ 'I',
                                           population == 'Ulombola' ~ 'J',
                                           population == 'full' ~ 'Full'),
                                 levels = (c('Full', LETTERS[10:1]) )),
             replicate_val = replicate_label)
  }
  
  pi_bootstrap_full_list[[calc_type]] <- bind_rows(pi_bootstrap_full_list[[calc_type]])
  
}

pi_bootstrap_full_list <- bind_rows(pi_bootstrap_full_list)

write.csv(pi_bootstrap_full_list, 
          here('diversity_divergence_stats', 'processed_results', 'pi_kazumbe_subsample.csv'), row.names = FALSE)



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# library(tidyverse)
# library(RColorBrewer)
# library(egg)
# library(sf)
# library(sp)
# library(ggmap)
# library(ggsn)
# library(cowplot)
# library(ggrepel)
# library(rgeos)
# 
# output_dir <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/diversity_divergence_stats/processed_results/'
# table_dir <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/tables/'
# 
# 
# if (!exists("cichlid_metadata")) cichlid_metadata <- read.csv("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/working_metadata_file/monster_23jul20.csv", stringsAsFactors = FALSE)
# 
# if (!exists("Fst_dataframe_updated")) Fst_dataframe_updated <- read.csv(paste0(output_dir, 'fst_7_2021.csv'))
# if(!exists("pi_info")) pi_info <- read.csv(paste0(output_dir, 'pi_7_2021.csv'))
# if(!exists("dxy_info")) dxy_info <- read.csv(paste0(output_dir, 'dxy_7_2021.csv'))
# 
#
# 
# 
# ### COLOR PALETTES ###
# location_palette <- brewer.pal(11,"Set3")
# #names(location_palette) <- unique(metadata_subset_reorder$location_plotting)
# names(location_palette) <- c("Hilltop", "Jakobsen's S", "Nondwa", "Kalalangabo", "Ulombola", "Katongwe N", "Gombe S", "Katongwe S", "Bangwe", "S Kagango", "Harembe")
# location_palette[location_palette == '#FFFFB3'] <- '#F7DC6F' #switch the yellow out
# location_palette[location_palette == "#D9D9D9"] <- "#AAB7B8" #switch the gray out
# location_palette[names(location_palette) == "Hilltop"] <- "#996633"
# location_palette[names(location_palette) == "Ulombola"] <- "#0000e6"
# location_palette[names(location_palette) == "Nondwa"] <- "#4dd2ff"
# location_palette[names(location_palette) == "Gombe S"] <- "#66ff66"
# location_palette[names(location_palette) == "Harembe"] <- "#99ffeb"
# 
# #switch the colors of Katongwe N and Ulombola
# Ulombola_col_switch <- location_palette[names(location_palette) == "Ulombola"]
# Katongwe_N_col_switch <- location_palette[names(location_palette) == "Katongwe N"]
# 
# location_palette[names(location_palette) == "Ulombola"] <- Katongwe_N_col_switch
# location_palette[names(location_palette) == "Katongwe N"] <- Ulombola_col_switch
# 
# 
# #switch Jakobsen's S and Bangwe
# Jakob_col_switch <- location_palette[names(location_palette) == "Jakobsen's S"]
# Bangwe_col_switch <- location_palette[names(location_palette) == "Bangwe"]
# 
# location_palette[names(location_palette) == "Jakobsen's S"] <- Bangwe_col_switch
# location_palette[names(location_palette) == "Bangwe"] <- Jakob_col_switch
# 
# location_palette_letter <- location_palette
# 
# location_letter_df <- data.frame(location = c('Gombe S', 'Katongwe N', 'Katongwe S', 'S Kagango', 'Kalalangabo', 'Nondwa', 'Hilltop', 'Bangwe', "Jakobsen's S", 'Ulombola', 'Harembe'),
#                                  loc_letter = LETTERS[1:11])
# 
# for (i in seq_along(location_palette_letter)) {
#   focal_loc <- names(location_palette_letter)[i]
#   focal_let <- location_letter_df$loc_letter[location_letter_df$location == focal_loc]
#   
#   #replace name with associated letter
#   names(location_palette_letter)[i] <- focal_let
# }
# 
# 
# 
# 
# pi_plot <- ggplot(data = pi_info) +
#   #geom_pointrange(aes(x = population, y = mean, ymin = lower, ymax = upper, group = species, color = species), size = 1.2) +
#   geom_linerange(aes(x = population, ymin = lower, ymax = upper, group = species, color = species), size = 1.5) +
#   geom_point(aes(x = population, y = mean, group = species, color = species), size = 5) +
#   scale_color_manual(values = c("#DC7633", "#3498DB")) +
#   #ylim(0, max(pi_info$upper)) +
#   xlab("Population") + ylab("pi") +
#   theme_cowplot() +
#   theme(panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         legend.position = "none") +
#   coord_flip() +
#   ggtitle('pi')
# 
# 
# 
# ggplot(data = dxy_info) +
#   geom_linerange(aes(x = population, ymin = lower, ymax = upper), size = 1.5) +
#   geom_point(aes(x = population, y = mean), size = 4.5) +
#   theme_cowplot() +
#   coord_flip()
# 
# dxy_plot <- ggplot(data = dxy_info) +
#   geom_linerange(aes(x = population, ymin = lower, ymax = upper), size = 1.5) +
#   geom_point(aes(x = population, y = mean), size = 5) +
#   theme_cowplot() +
#   coord_flip() +
#   theme(panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank()) +
#   scale_x_discrete(drop = FALSE) +
#   ggtitle('dxy')
# 
# plot_grid(pi_plot, dxy_plot, ncol = 2, align = 'hv')
# 
# 
# #https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
# 
# #combined plot of pi and dxy
# divergence_diversity_plot <- ggplot(data = pi_info %>%
#          mutate(metric = case_when(species == 'kazumbe' ~ '\u03C0 (kazumbe)',
#                             species == 'polyodon' ~ '\u03C0 (polyodon)')) %>% 
#          select(-c(pop, species)) %>%
#          rbind(., dxy_info %>% 
#                  select(-comparison_id) %>% 
#                  mutate(metric = 'dxy')) %>% 
#            mutate(metric = factor(metric, levels = c('\u03C0 (polyodon)', '\u03C0 (kazumbe)', 'dxy')))) +
#   geom_linerange(aes(x = Population, ymin = lower, ymax = upper, group = metric, color = metric), size = 2) +
#   geom_point(aes(x = Population, y = mean, group = metric, color = metric), shape = 21, stroke = 2, size = 5.5) +
#   scale_color_manual(values = c("#3498DB", "#DC7633", "#6c757d"),
#                      labels = c(expression(paste('\u03C0 (', italic('P. polyodon'), ')')), expression(paste('\u03C0 (', italic('P. kazumbe'), ')')), 'dxy' ) ) +
#   #ylim(0, max(pi_info$upper)) +
#   xlab("Population") + ylab("Average sequence divergence") +
#   theme_cowplot() +
#   theme(plot.margin = margin(5.5, 0, 5.5, 5.5),
#         panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
#         legend.title = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_text(size = 20),
#         axis.title.x = element_text(size = 20),
#         axis.text.x = element_text(size = 15),
#         legend.position="bottom",
#         legend.spacing.x = unit(0.2, 'cm'),
#         legend.text = element_text(size = 15, margin = margin(r = 10, unit = "pt")),) +
#   #theme(panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#   #      legend.position = "none") +
#   coord_flip()
# 
# 
# 
# Fst_plot <- ggplot(data = Fst_dataframe_updated %>% 
#                      mutate(group = 'a',
#                             Population = factor(Population, levels = (c('Full', LETTERS[10:1]) )))) +
#   #geom_pointrange(aes(x = population, y = mean, ymin = lower, ymax = upper, group = species, color = species), size = 1.2) +
#   geom_linerange(aes(x = Population, ymin = lower_CI, ymax = upper_CI, color = group), size = 2) +
#   geom_point(aes(x = Population, y = Fst, color = group), shape = 21, stroke = 2, size = 5.5) +
#   #ylim(0, max(pi_info$upper)) +
#   xlab("Population") + ylab("Fst") +
#   theme_cowplot() +
#   theme(plot.margin = margin(5.5, 5.5, 0, 5.5),
#         panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.x = element_text(size = 20),
#         axis.text.x = element_text(size = 15),
#         legend.position = "bottom",
#         legend.text = element_text(color = "white"),
#         legend.title = element_text(color = "white"),
#         legend.key = element_rect(fill = "white")) +
#   scale_x_discrete(drop = FALSE) +
#   scale_color_manual(values = c('#5e60ce'),
#                      guide = guide_legend(override.aes = list(color = "white"))) +
#   #scale_color_discrete(guide = guide_legend(override.aes = list(color = "white"))) +
#   coord_flip()
# 
# diversity_multipanel <- ggarrange(divergence_diversity_plot, Fst_plot,
#           ncol = 2)
#   
# 
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/diversity_divergence_stats/pixy/dxy_pi_plot.png',
#         width = 7, height = 9)
# 
# 
# 
# shapefile_path <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/non_vcf_datafiles/shapefiles/'
# output_path <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/'
# africa_lake_shp <- read_sf(paste0(shapefile_path, 'africawaterbody/Africa_waterbody.shp'))
# 
# #africa_lake_shp_project <- spTransform(africa_lake_shp, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
# #africa_lake_shp_project <- st_transform(africa_lake_shp, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
# #africa_lake_shp_project <- st_transform(africa_lake_shp, "+proj=longlat +ellps=WGS84 +datum=WGS84")
# 
# 
# tz_river_shp <- read_sf(paste0(shapefile_path, 'AFRICOVER_TZ_RIVERS-shapefile/AFRICOVER_TZ_RIVERS.shp'))
# 
# #metadata
# cichlid_metadata <- read.csv("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/working_metadata_file/monster_23jul20.csv", stringsAsFactors = FALSE)
# 
# 
# sampling_locations <- cichlid_metadata %>% 
#   dplyr::filter(sciname2 %in% c('Petrochromis kazumbe', 'Petrochromis polyodon') & location != 'Harembe') %>% 
#   select(location, long, lat) %>%
#   distinct(location, .keep_all = TRUE) %>% 
#   mutate(Region = case_when(location == 'Ulombola' ~ 'South',
#                             location %in% c('Hilltop', 'Jakob', 'Bangwe') ~ 'Mid',
#                             location %in% c('Nondwa', 'Kalala', 'Ska', 'Katongwe_S', 'Katongwe_N', 'Gombe_South') ~ 'North')) %>% 
#   arrange(desc(lat)) %>% 
#   mutate(number_label = LETTERS[1:n()],
#          #number_label = 1:n(),
#          Region = factor(Region, levels = c('North', 'Mid', 'South')))
# 
# 
# 
# tz_river_points <- sf::st_point_on_surface(tz_river_shp)
# 
# # retrieve the coordinates
# tz_river_coords <- as.data.frame(sf::st_coordinates(tz_river_points))
# tz_river_coords$NAME <- tz_river_shp$AUTO_ID
# 
# 
# ggplot() +
#   geom_sf(data = africa_lake_shp, fill = '#ade8f4', color = '#ade8f4', size = 3) +
#   geom_sf(data = tz_river_shp, fill = '#ade8f4', color = '#ade8f4', size = 1.5) +
#   geom_text(data = tz_river_coords, aes(X, Y, label = NAME), colour = "black") +
#   xlim(29.55, 29.79) + ylim(-5.03, -4.72)
# 
# ggplot() +
#   geom_sf(data = africa_lake_shp, fill = '#ade8f4', color = '#ade8f4', size = 3) +
#   geom_sf(data = tz_river_shp[tz_river_shp$AUTO_ID %in% c('5153', '5070'),], fill = '#ade8f4', color = '#ade8f4', size = 1.5) +
#   xlim(29.55, 29.79) + ylim(-5.03, -4.72)
# 
# 
# #https://geodata.lib.berkeley.edu/catalog/AFRICOVER_TZ_RIVERS
# 
# ####
# #water_color <- '#ade8f4'
# water_color <- '#caf0f8'
# 
# 
# kigoma_region_map <- ggplot() +
#   geom_sf(data = africa_lake_shp, fill = water_color, color = water_color, size = 2) +
#   geom_sf(data = tz_river_shp[tz_river_shp$AUTO_ID %in% c('5153', '5070'),], fill = water_color, color = water_color, size = 1.5) +
#   ggrepel::geom_label_repel(data = sampling_locations,
#                             aes(x = long, y = lat, label = number_label),
#                             size = 7, alpha = 1,
#                             label.size = NA,
#                             fill = NA,
#                             #segment.color = "black", segment.size = 1,
#                             box.padding = 0.009, #0.01
#                             point.padding = 0,
#                             #min.segment.length = 0.01, nudge_x = 0.02,
#                             #min.segment.length = 1, nudge_x = 0.011,
#                             min.segment.length = 1, nudge_x = -0.0112,
#                             seed = 1003) +
#   geom_text(data = data.frame(label = c('Lake Tanganyika', 'Luiche River'),
#                               lon = c(29.61, 29.705),
#                               lat = c(-4.97, -4.847)),
#             aes(x = lon, y = lat, label = label, family = "Times", fontface = "plain"), size = 8, color = '#0077b6') +
#   geom_point(data = sampling_locations, aes(x = long, y = lat, fill = Region), size = 8.5, pch=21, color = '#e7e7e6') +
#   scale_fill_manual(values = c('#602437', '#b9375e', '#ff9ebb')) +
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid.major = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position = c(0.72, 1),  #c(0.95, 0.45) c(0.63, 0.173) --> lower left, next to scale bar
#         legend.justification = c("right", "top"),
#         legend.text = element_text(size = 18),
#         legend.title = element_text(size = 20),
#         legend.key=element_blank(),
#         legend.background=element_blank()) +
#   #coord_cartesian(clip = 'off') +
#   xlim(29.55, 29.79) + ylim(-5.03, -4.72) + 
#   scalebar(x.min = 29.55, x.max = 29.79, y.min = -5.03, y.max = -4.72, 
#            dist = 5, dist_unit = "km",
#            st.dist = .02,
#            st.size = 4.5,
#            height = .013,
#            border.size = 0.6,
#            transform = TRUE, model = "WGS84", st.bottom = FALSE, location = "bottomleft")
# 
# 
# #lake tanganyika label: 29.61, -4.97
# #Luiche River label: 29.72, -4.85
# 
# 
# kigoma_box <- st_as_sfc(rgeos::bbox2SP(n = -4.72, s = -5.03, w = 29.55, e = 29.79,
#                                        proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))
# 
# kigoma_inset_map <- ggplot() + 
#   geom_sf(data = africa_lake_shp[africa_lake_shp$NAME_OF_WA == "Tanganyika", ], fill = water_color, color = water_color, size = 0.1) +
#   geom_sf(data = kigoma_box, fill = NA, color = "red", size = 1) +
#   #xlim(28.6, 31.3) + 
#   ylim(-8.8, -3.3) +
#   scale_x_continuous(breaks = c(29, 30, 31)) + 
#   #  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + 
#   theme_minimal() +
#   theme(#panel.grid.major = element_blank(),
#     axis.text = element_text(size = 9),
#     panel.grid.major = element_line(color = '#e7e7e6'))
# 
# #water body color: #0077b6
# 
# tanganyika_sampling_map_final <- ggdraw() +
#   draw_plot(kigoma_region_map, width = 1, height = 1) +
#   draw_plot(kigoma_inset_map, x = 0.65, y = 0.64, width = 0.37, height = 0.37) +
#   theme(plot.margin = margin(5.5, -5, 5.5, 5.5))
#   
# 
# 
# 
# 
# 
# 
# multipanel_diversity_toprow <- plot_grid(tanganyika_sampling_map_final, ggdraw(diversity_multipanel),
#                                          ncol = 2, rel_widths = c(0.45,0.5), labels = c('(a)', '(b)'), label_size = 32, vjust = 0, hjust = c(-0.1, 0.1) ) +
#   theme(plot.margin = margin(32, 0, 8, 0))
# 
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/map_divergencestats_multipanel_5_12_2021.png',
#         width = 16, height = 9.5)
# 
# 
# 
# library(magick)
# #library(ggtext)
# 
# 
# fish_image_list <- list()
# 
# fish_image_list[['kazumbe']] <- image_read("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Fish_images/kazumbe_cropped.tif")
# fish_image_list[['polyodon']] <- image_read("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Fish_images/polyodon_pic_cropped.tif")
# 
# # ggdraw() +
# #   draw_image(
# #     fish_image_list$kazumbe, scale = 1, x = 1,
# #     hjust = 1, halign = 0.5, valign = 0
# #   ) +
# #   ggtitle('Petrochromis kazumbe')
# # draw_image(fish_image_list$kazumbe)
# 
# kazumbe_plot <- ggplot() +
#   ggtitle('Petrochromis kazumbe') +
#   draw_image(
#     fish_image_list$kazumbe, scale = 1, x = 1,
#     hjust = 1, halign = 0.5, valign = 0
#   ) +
#   theme_void() +
#   theme(#plot.margin=unit(c(0.75, 0.75, 0.4, 0),"cm"),
#         plot.margin=unit(c(0.4, 0, 0.2, 0),"cm"),
#         plot.title = element_text(size = 22, face = 'italic', hjust = 0.5)) +
#   coord_cartesian(clip = 'off')
# 
# polyodon_plot <- ggplot() +
#   ggtitle('Petrochromis polyodon') +
#   draw_image(
#     fish_image_list$polyodon, scale = 1, x = 1,
#     hjust = 1, halign = 0.5, valign = 0
#   ) +
#   theme_void() +
#   theme(#plot.margin=unit(c(0.65, 0.75, 0.5, 0),"cm"),
#         plot.margin=unit(c(0.4, 0, 0.2, 0),"cm"),
#         plot.title = element_text(size = 22, face = 'italic', hjust = 0.5)) +
#   coord_cartesian(clip = 'off')
# 
# fish_multipanel <- plot_grid(kazumbe_plot, polyodon_plot, nrow = 2) +
#   theme(plot.margin = unit(c(0.2, -0.2, 0.1, 0.5), "cm"))
# 
#   
# plot_grid(fish_multipanel, map_diversity_multipanel, nrow = 2, 
#           rel_heights = c(0.25, 0.75))
# 
# 
# 
# 
# plot_grid(tanganyika_sampling_map_final, fish_multipanel, ncol = 2)
# 
# 
# 
# 
# ### PCAs ###
# both_species_pca <- read.csv("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/kazumbe_polyodon_PCA_withmetadata.csv")
# polyodon_pca <- read.csv("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/polyodon_PCA_withmetadata.csv")
# kazumbe_pca <- read.csv("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/kazumbe_PCA_withmetadata.csv")
# 
# pca_importance_info <- read.csv("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/pca_importance_info.csv")
# 
# scaleFUN <- function(x) sprintf("%.2f", x)
# 
# #PCA containing both species
# species_PCA <- ggplot(data = both_species_pca, aes(x = PC1, y= PC2, color = species_ID_entropy)) +
#   #geom_point(size = 5) +
#   geom_point(size = 4.5, alpha = 0.5) +
#   geom_point(size = 4.5, shape = 21, alpha = 0.7, stroke = 1) +
#   scale_color_manual(name = "Species", values = c("#DC7633", "#3498DB")) +
#   ###scale_color_manual(name = "Species", values = c("#ffc966", "#ccccff")) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 13),
#         axis.title.y = element_text(margin = margin(0, 0, 0, 0)), 
#         #legend.title = element_text(size = 18),
#         legend.title = element_blank(),
#         legend.position = "bottom",
#         legend.margin = margin(-6, 0, 0, 0),
#         legend.text = element_text(size = 15, face = "italic"),
#         plot.margin = margin(5, 5, 5, 9)) +
#   scale_x_continuous(labels = scaleFUN) +
#   scale_y_continuous(labels = scaleFUN) +
#   xlab(paste0("PC1 (", round(pca_importance_info[pca_importance_info$dataset == "combined_kazumbe_polyodon" & pca_importance_info$pc == 'PC1',]$prop_variance*100, 2),"%)")) +
#   ylab(paste0("PC2 (", round(pca_importance_info[pca_importance_info$dataset == "combined_kazumbe_polyodon" & pca_importance_info$pc == 'PC2',]$prop_variance*100, 2),"%)"))
# 
# 
# kazumbe_pca_plus_locletters <- kazumbe_pca %>% 
#   mutate(location_factor_letter = factor(case_when(location == 'Gombe_South' ~ "A",
#                                                    location == 'Katongwe_N' ~ "B",
#                                                    location == 'Katongwe_S' ~ "C",
#                                                    location == 'Ska' ~ "D",
#                                                    location == 'Kalala' ~ "E",
#                                                    location == 'Nondwa' ~ "F",
#                                                    location == 'Hilltop' ~ "G",
#                                                    location == 'Bangwe' ~ "H",
#                                                    location == 'Jakob' ~ "I",
#                                                    location == 'Ulombola' ~ "J"),
#                                          levels= LETTERS[1:10]))
# 
# #PCA containing only kazumbe
# kazumbe_PCA_plot <- ggplot(data = kazumbe_pca_plus_locletters, 
#                            aes(x = PC1, y= PC2, color = location_factor_letter)) +
#   geom_point(size = 4.5, alpha = 0.5) +
#   geom_point(size = 4.5, shape = 1, alpha = 0.5) +
#   scale_color_manual(name = "Location", values = location_palette_letter[names(location_palette_letter) %in% unique(kazumbe_pca_plus_locletters$location_factor_letter) ]) +
#   ###scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(kazumbe_PCA_combinedinfo$location_plotting) ]) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 13),
#         axis.title.y = element_text(margin = margin(0, 0, 0, 0)), 
#         legend.position = "none",
#         plot.title = element_text(hjust = 0.5, size = 22, face="italic"),
#         plot.margin = margin(5, 5, 5, 9)) +
#   ggtitle("P. kazumbe") + 
#   scale_x_continuous(labels = scaleFUN) +
#   scale_y_continuous(labels = scaleFUN) +
#   xlab(paste0("PC1 (", round(pca_importance_info[pca_importance_info$dataset == "kazumbe" & pca_importance_info$pc == 'PC1',]$prop_variance*100, 2),"%)")) +
#   ylab(paste0("PC2 (", round(pca_importance_info[pca_importance_info$dataset == "kazumbe" & pca_importance_info$pc == 'PC2',]$prop_variance*100, 2),"%)"))
# 
# 
# 
# polyodon_pca_plus_locletters <- polyodon_pca %>% 
#   mutate(location_factor_letter = factor(case_when(location == 'Gombe_South' ~ "A",
#                                                    location == 'Katongwe_N' ~ "B",
#                                                    location == 'Katongwe_S' ~ "C",
#                                                    location == 'Ska' ~ "D",
#                                                    location == 'Kalala' ~ "E",
#                                                    location == 'Nondwa' ~ "F",
#                                                    location == 'Hilltop' ~ "G",
#                                                    location == 'Bangwe' ~ "H",
#                                                    location == 'Jakob' ~ "I",
#                                                    location == 'Ulombola' ~ "J"),
#                                          levels= LETTERS[1:10]))
# #PCA containing only polyodon
# polyodon_PCA_plot <- ggplot(data = polyodon_pca_plus_locletters, 
#                             aes(x = PC1, y= PC2, color = location_factor_letter)) +
#   geom_point(size = 4.5, alpha = 0.5) +
#   geom_point(size = 4.5, shape = 1, alpha = 0.5) +
#   scale_color_manual(name = "Location", values = location_palette_letter[names(location_palette_letter) %in% unique(polyodon_pca_plus_locletters$location_factor_letter) ]) +
#   ###scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(polyodon_PCA_combinedinfo$location_plotting) ]) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 20),
#         axis.text = element_text(size = 13),
#         axis.title.y = element_text(margin = margin(0, 0, 0, 0)), 
#         legend.position = "none",
#         plot.title = element_text(hjust = 0.5, size = 22, face= c("italic") ),
#         plot.margin = margin(5, 5, 5, 9)) +
#   ggtitle("P. polyodon") + 
#   scale_x_continuous(labels = scaleFUN) +
#   scale_y_continuous(labels = scaleFUN) +
#   #xlab(paste0("PC1 (", round(pca_importance_info[pca_importance_info$dataset == "polyodon" & pca_importance_info$pc == 'PC1',]$prop_variance*100, 2),"%)")) +
#   xlab(paste0("PC1 (", format(round(pca_importance_info[pca_importance_info$dataset == "polyodon" & pca_importance_info$pc == 'PC1',]$prop_variance*100, 2), nsmall = 2) ,"%)")) +
#   ylab(paste0("PC2 (", round(pca_importance_info[pca_importance_info$dataset == "polyodon" & pca_importance_info$pc == 'PC2',]$prop_variance*100, 2),"%)"))
# 
# 
# 
# both_species_pca_with_locletters <- both_species_pca %>% 
#   mutate(Location = factor(case_when(location == 'Gombe_South' ~ "A",
#                                      location == 'Katongwe_N' ~ "B",
#                                      location == 'Katongwe_S' ~ "C",
#                                      location == 'Ska' ~ "D",
#                                      location == 'Kalala' ~ "E",
#                                      location == 'Nondwa' ~ "F",
#                                      location == 'Hilltop' ~ "G",
#                                      location == 'Bangwe' ~ "H",
#                                      location == 'Jakob' ~ "I",
#                                      location == 'Ulombola' ~ "J"),
#                            levels= LETTERS[1:10]))
# 
# #extracting legend
# location_legend <- get_legend(ggplot(data = both_species_pca_with_locletters,  
#                                      aes(x = PC1, y= PC2, color = Location)) +
#                                 geom_point(size = 4.5, alpha = 0.5) +
#                                 geom_point(size = 4.5, shape = 1, alpha = 0.5) +
#                                 scale_color_manual(name = "Location", values = location_palette_letter[names(location_palette_letter) %in% unique(both_species_pca_with_locletters$Location) ]) +
#                                 ##scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(petro_PCA_combinedinfo$location_plotting) ]) +
#                                 theme_bw() +
#                                 theme(panel.grid.major = element_blank(),
#                                       panel.grid.minor = element_blank(),
#                                       axis.title = element_text(size = 18),
#                                       plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
#                                       legend.title = element_text(size = 18),
#                                       legend.text = element_text(size = 15),
#                                       #legend.text = element_text(size = 18, margin = margin(r = 22, unit = "pt")),
#                                       legend.position = "right",
#                                       #legend.spacing.x = unit(0.3, 'cm')
#                                       )
#                               )
# 
# 
# plot_grid(plot_grid(kazumbe_PCA_plot, polyodon_PCA_plot, ncol = 2),
#           location_legend, nrow = 2, rel_heights = c(0.9, 0.1))
# 
# 
# #species_pca_multipanel <- plot_grid(kazumbe_PCA_plot, polyodon_PCA_plot, location_legend, ncol = 3, rel_widths = c(0.4, 0.4, 0.2) )
# 
# pca_multipanel <- plot_grid(species_PCA, 
#                             kazumbe_PCA_plot, 
#                             polyodon_PCA_plot,
#                             ncol = 3, align = 'hv', axis = "bt")
# 
# 
# #ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/prac_pca_multipanel_5_12_2021.png',
# #        width = 18, height = 5)
# 
# empty_plot <- ggplot() + theme_void()
# 
# multipanel_diversity_bottomrow <- plot_grid(fish_multipanel, pca_multipanel, plot_grid(location_legend, empty_plot, nrow = 2, rel_heights = c(0.9, 0.107)), 
#                                             ncol = 3, rel_widths = c(0.35, 0.95, 0.125),
#                                             labels = c('(c)', '(d)',''), label_size = 32, vjust = 1, hjust = c(-0.1, -0.3, -0.5) ) +
#   theme(plot.margin = margin(15, 0, 0, 0))
# 
# plot_grid(multipanel_diversity_toprow,
#           multipanel_diversity_bottomrow, nrow = 2, rel_heights = c(0.7, 0.3))
# 
# 
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/prac_full_multipanel_5_12_2021.png',
#         width = 19, height = 16)
# 
# plot_grid(species_PCA, kazumbe_PCA_plot, ncol = 2, align = 'h', axis = "bt")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # diversity_calc <- function(input_data,
# #                            calc_type = c('raw', 'per_site'),
# #                            no_sites_col,
# #                            count_diffs_col,
# #                            count_comparisons_col) {
# # 
# #   if (calc_type == 'raw') {
# #     return( sum(input_data[, count_diffs_col], na.rm = TRUE)/sum(input_data[, count_comparisons_col], na.rm = TRUE) )
# #   } else if (calc_type == 'per_site') {
# #     return( (sum(input_data[, count_diffs_col], na.rm = TRUE)/sum(input_data[, count_comparisons_col], na.rm = TRUE))/sum(input_data[, no_sites_col], na.rm = TRUE) )
# #   }
# # }
# 
# #https://stackoverflow.com/questions/42438450/make-legend-invisible-but-keep-figure-dimensions-and-margins-the-same
# #https://stackoverflow.com/questions/47614314/how-to-put-plots-without-any-space-using-plot-grid
# #https://stackoverflow.com/questions/57168938/ggplot-italicize-part-of-legend-across-multiple-lines
# 
# 
# divergence_diversity_plot_horizontal <- ggplot(data = pi_info %>%
#                                       mutate(metric = case_when(species == 'kazumbe' ~ '\u03C0 (kazumbe)',
#                                                                 species == 'polyodon' ~ '\u03C0 (polyodon)')) %>% 
#                                       select(-c(pop, species)) %>%
#                                       rbind(., dxy_info %>% 
#                                               select(-comparison_id) %>% 
#                                               mutate(metric = 'dxy')) %>% 
#                                       mutate(metric = factor(metric, levels = c('\u03C0 (polyodon)', '\u03C0 (kazumbe)', 'dxy')))) +
#   geom_linerange(aes(x = Population, ymin = lower, ymax = upper, group = metric, color = metric), size = 2) +
#   geom_point(aes(x = Population, y = mean, group = metric, color = metric), size = 5.5, shape = 21, stroke = 1.5) +
#   scale_shape(solid = FALSE) +
#   scale_color_manual(values = c("#3498DB", "#DC7633", "#6c757d"), 
#                      labels = c( (expression(paste("\u03C0 (", italic("P. polyodon"), ")"))),
#                                  (expression(paste("\u03C0 (", italic("P. kazumbe"), ")"))),
#                                  (("dxy")) )) +
#   #ylim(0, max(pi_info$upper)) +
#   xlab("Population") + ylab("Avg. seq. div.") +
#   theme_cowplot() +
#   theme(plot.margin = margin(0, 10, 0, 13),
#         panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
#         legend.title = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 22),
#         axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 11, b = 0, l = 0)),
#         axis.text.y = element_text(size = 18),
#         legend.text = element_text(size = 18, margin = margin(r = 22, unit = "pt")),
#         legend.position = "bottom",
#         legend.spacing.x = unit(0.3, 'cm')) 
#   #theme(panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#   #      legend.position = "none") +
#   #coord_flip()
# 
# 
# Fst_plot_horizontal <- ggplot(data = Fst_dataframe_updated %>% 
#                      mutate(group = 'a')) +
#   #geom_pointrange(aes(x = population, y = mean, ymin = lower, ymax = upper, group = species, color = species), size = 1.2) +
#   geom_linerange(aes(x = Population, ymin = lower_CI, ymax = upper_CI, color = group), size = 2) +
#   geom_point(aes(x = Population, y = Fst, color = group), size = 5.5, shape = 21, stroke = 1.5) +
#   scale_shape(solid = FALSE) +
#   #ylim(0, max(pi_info$upper)) +
#   xlab("Population") + ylab("Fst") +
#   theme_cowplot() +
#   theme(plot.margin = margin(5.5, 10, 0, 13),
#         panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 18),
#         axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 11, b = 0, l = 0)),
#         legend.position = "none",
#         legend.text = element_text(color = "white"),
#         legend.title = element_text(color = "white"),
#         legend.key = element_rect(fill = "white")) +
#   scale_x_discrete(drop = FALSE) +
#   scale_color_manual(values = c('#5e60ce'),
#                      guide = guide_legend(override.aes = list(color = "white")))
#   #scale_color_discrete(guide = guide_legend(override.aes = list(color = "white"))) +
#   #coord_flip()
# 
# diversity_multipanel_horizontal <- ggarrange(Fst_plot_horizontal, divergence_diversity_plot_horizontal,
#                                   nrow = 2)
# 
# 
# 
# map_fish_multi <- plot_grid(tanganyika_sampling_map_final, fish_multipanel, ncol = 2, rel_widths = c(0.5, 0.5), labels = c('(a)', '(b)'), label_x = c(0, -0.03), label_size = 28)
# 
# plot_grid(map_fish_multi, ggdraw(diversity_multipanel_horizontal) + theme(plot.margin = unit(c(0, 0.5, 0.5, 0.8),"cm")), 
#           nrow = 2, rel_heights = c(0.65, 0.35), labels = c('', '(c)'), label_size = 28, align = "hv")
# 
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/map_divergencestats_fish_multipanel.png',
#         width = 18, height = 15)
# 
# 
# 
# map_fish_multi_alt <- plot_grid(fish_multipanel, tanganyika_sampling_map_final, ncol = 2, rel_widths = c(0.52, 0.48), labels = c('(a)', '(b)'), label_x = c(0, -0.02), label_size = 28)
# 
# plot_grid(map_fish_multi_alt, ggdraw(diversity_multipanel_horizontal) + theme(plot.margin = unit(c(0, 0.5, 0.5, 0.8),"cm")), 
#           nrow = 2, rel_heights = c(0.65, 0.35), labels = c('', '(c)'), label_size = 28)
# 
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/map_divergencestats_fish_multipanel_alt.png',
#         width = 18, height = 15)
# 


# reich.fst_alt <- function(gl, bootstrap=FALSE, plot=FALSE, verbose=TRUE) { 
#   if (!require("matrixStats",character.only=T, quietly=T)) {
#     install.packages("matrixStats")
#     library(matrixStats, character.only=T)
#   }
#   if (!require("dplyr",character.only=T, quietly=T)) {
#     install.packages("dplyr")
#     library(dplyr, character.only=T)
#   }
#   
#   nloc <- gl@n.loc
#   npop <- length(levels(gl@pop))
#   
#   fsts <- matrix(nrow=npop,
#                  ncol=npop,
#                  dimnames=list(levels(gl@pop),levels(gl@pop)))
#   
#   if (bootstrap != FALSE){
#     n.bs <- bootstrap
#     bs <- data.frame(matrix(nrow=nrow(combinat::combn2(levels(gl@pop))),
#                             ncol=n.bs+5))
#   }
#   
#   k <- 0
#   
#   for (p1 in levels(gl@pop)){
#     for (p2 in levels(gl@pop)){
#       if (which(levels(gl@pop) == p1) < which(levels(gl@pop) == p2)) {
#         k <- 1+k
#         
#         #pop1 <- gl.keep.pop(gl, p1, mono.rm=FALSE, v=0)
#         
#         pop1 <- popsub(
#           gid = gl,
#           sublist = p1,
#           exclude = NULL,
#           blacklist = NULL,
#           mat = NULL,
#           drop = FALSE
#         )
#         
#         
#         #pop2 <- gl.keep.pop(gl, p2, mono.rm=FALSE, v=0)
#         
#         pop2 <- popsub(
#           gid = gl,
#           sublist = p2,
#           exclude = NULL,
#           blacklist = NULL,
#           mat = NULL,
#           drop = FALSE
#         )
#         
#         a1 <- colSums2(as.matrix(pop1),na.rm=T)
#         a2 <- colSums2(as.matrix(pop2),na.rm=T)
#         n1 <- apply(as.matrix(pop1),2,function(x) 2*sum(!is.na(x)))
#         n2 <- apply(as.matrix(pop2),2,function(x) 2*sum(!is.na(x)))
#         
#         h1 <- (a1*(n1-a1))/(n1*(n1-1))
#         h2 <- (a2*(n2-a2))/(n2*(n2-1))
#         
#         N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
#         D <- N + h1 + h2
#         
#         F <- sum(N, na.rm=T)/sum(D, na.rm=T)
#         fsts[p2,p1] <- F
#         if (verbose == TRUE) {
#           print(paste("Pop1: ",p1,", Pop2: ",p2,", Reich FST: ",F,sep=""))
#         }
#         
#         if (bootstrap != FALSE) {
#           if (verbose == TRUE) {
#             print("beginning bootstrapping")
#           }
#           
#           bs[k,1:3] <- c(p2,p1,as.numeric(F))
#           
#           for (i in 1:n.bs){
#             loci <- sample((1:nloc), nloc, replace=TRUE)
#             
#             pop1.bs <- matrix(as.matrix(pop1)[,loci],
#                               ncol=length(loci))
#             pop2.bs <- matrix(as.matrix(pop2)[,loci],
#                               ncol=length(loci))
#             
#             a1 <- colSums2(as.matrix(pop1.bs),na.rm=T)
#             a2 <- colSums2(as.matrix(pop2.bs),na.rm=T)
#             n1 <- apply(as.matrix(pop1.bs),2,function(x) 2*sum(!is.na(x)))
#             n2 <- apply(as.matrix(pop2.bs),2,function(x) 2*sum(!is.na(x)))
#             
#             h1 <- (a1*(n1-a1))/(n1*(n1-1))
#             h2 <- (a2*(n2-a2))/(n2*(n2-1))
#             
#             N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
#             D <- N + h1 + h2
#             
#             F.bs <- sum(N, na.rm=T)/sum(D, na.rm=T)
#             bs[k,i+5] <- F.bs
#           }
#           if (verbose == TRUE){
#             print(paste("bootstrapping 95% CI: ",
#                         quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),"-",
#                         quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T)))
#           }
#           
#           bs[k,4:5] <- c(quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),
#                          quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T))
#         }
#         
#       }
#     }
#   }
#   
#   fsts[fsts < 0] <- 0
#   
#   if (bootstrap != FALSE){
#     colnames(bs)[1:5] <- c("pop1","pop2","fst_estimate","min_CI","max_CI")
#     fst.list <- list(fsts,bs)
#     names(fst.list) <- c("fsts","bootstraps")
#     
#     if (plot == TRUE){
#       print("drawing plot with bootstraps")
#       
#       if (!require("ggplot2",character.only=T, quietly=T)) {
#         install.packages("ggplot2")
#         library(ggplot2, character.only=T)
#       }
#       
#       plot.data <- bs[,1:5]
#       plot.data$fst_estimate <- as.numeric(plot.data$fst_estimate)
#       plot.data$min_CI <- as.numeric(plot.data$min_CI)
#       plot.data$max_CI <- as.numeric(plot.data$max_CI)
#       plot.data$pop_pair <- paste(plot.data$pop1,plot.data$pop2,sep="_")
#       plot.data$signif <- case_when(plot.data$min_CI > 0 ~ TRUE,
#                                     TRUE ~ FALSE)
#       
#       
#       bs.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate,col=signif)) + 
#         geom_point(size=2) + 
#         coord_flip() + 
#         geom_errorbar(aes(ymin=min_CI,ymax=max_CI),width=0.1,size=1) + 
#         geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
#         theme_minimal() + 
#         theme(legend.position="none")
#       
#       print(bs.plot)
#     }
#   } else {
#     fst.list <- list(fsts)
#     names(fst.list) <- "fsts"
#     
#     if (plot == TRUE){
#       print("drawing plot without bootstraps")
#       
#       if (!require("ggplot2",character.only=T, quietly=T)) {
#         install.packages("ggplot2")
#         library(ggplot2, character.only=T)
#       }
#       
#       plot.data <- data.frame(combinat::combn2(row.names(fsts)),
#                               fst_estimate=fsts[lower.tri(fsts)])
#       plot.data$pop_pair <- paste(plot.data$X1,plot.data$X2,sep="_")
#       
#       fst.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate)) + 
#         geom_point(size=2) + 
#         coord_flip() + 
#         geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
#         theme_minimal() + 
#         theme(legend.position="none")
#       
#       print(fst.plot)
#     }
#   }
#   
#   return(fst.list)
#   beepr::beep()
# }

