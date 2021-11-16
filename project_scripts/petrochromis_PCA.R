#################################################################################################
### SCRIPT NAME: petrochromis_PCA.R
### PURPOSE: Principal component analysis (PCA) of kazumbe and polydon. These results are then
###          visualized in the studysystem_overview_fig.R script.
### PRODUCTS:
###     kazumbe_polyodon_PCA_withmetadata_9_7_2021.csv: results of PCA that included all samples
###                                                     of both species (with metadata)
###     polyodon_PCA_withmetadata_9_7_2021.csv: results of PCA that included only samples of
###                                             polyodon (with metadata)
###     kazumbe_PCA_withmetadata_9_7_2021.csv: results of PCA that included only samples of
###                                             kazumbe (with metadata)
###     pca_importance_info.csv: information about the amount of variance explained by each PC
###                              for the three PCAs.
#################################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(tidyverse)
library(cowplot)
library(adegenet)
library(vcfR)
library(RColorBrewer)
library(rhdf5)
library(abind)
library(here)


### Load data ###

#metadata
cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata

#vcf file (this file was already processed with vcftools)
petro_vcfR <- read.vcfR(here('vcf', 'biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.vcf'))

#entropy
entropy_samples_order <- read.delim(here('working_metadata_file', 'SAMPLES_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.txt'), header = FALSE) #order of samples in entropy results

mod_list <- list()
for (chain in 1:3) {
  mod_list[[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0('kazumbe_polyodon_qmod_6_18_2021_', 2, 'rep', chain, '.hdf5')), "q")
}


### Custom PCA function ###
doPCA <- function(gmat, dims = 10){  ### inds in columns, loci in rows
  gmn <- apply(gmat, 1, mean, na.rm = T)  
  gmnmat <- matrix(gmn, nrow = nrow(gmat), ncol = ncol(gmat))
  gprime <- gmat - gmnmat ## remove mean
  gcovarmat <- matrix(NA, nrow = ncol(gmat), ncol = ncol(gmat))
  
  for(i in 1:ncol(gmat)){
    for(j in i:ncol(gmat)){
      if (i == j){
        gcovarmat[i, j] <- cov(gprime[, i], gprime[, j], use = "pairwise.complete.obs")
      }
      else{
        gcovarmat[i, j] <- cov(gprime[, i], gprime[, j], use = "pairwise.complete.obs")
        gcovarmat[j, i] <- gcovarmat[i, j]
      }
    }
  }
  
  rownames(gcovarmat) <- colnames(gprime)
  colnames(gcovarmat) <- colnames(gprime)
  
  ## pca on the genotype covariance matrix
  pcgcov <- prcomp(x = gcovarmat, center = TRUE, scale = FALSE)
  
  return(list(PCA = as.data.frame(round(pcgcov$x[, 1:dims], digits = 5)),
              imp = summary(pcgcov)))
}



#######################
### PROCESSING DATA ###
#######################

### Process genetic dataset ###
petro_gen <- vcfR2genlight(petro_vcfR) #convert vcf to genlight abject
petro_geno_matrix <- t(as.matrix(petro_gen)) #extract genotype matrix

### Process entropy results and metadata and combine together ###
metadata_subset_entropy <- cichlid_metadata[cichlid_metadata$sample_code %in% entropy_samples_order$V1, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2')]
metadata_subset_entropy_reorder <- metadata_subset_entropy[match(as.character(entropy_samples_order$V1), metadata_subset_entropy$sample_code),]

#entropy
k2_combined.q <- abind(mod_list$chain_1, mod_list$chain_2, mod_list$chain_3, along=1)
k2_q_estimate <- as.data.frame(t(apply(k2_combined.q, 2:3, mean)))

k2_q_estimate_with_metadata <- cbind(k2_q_estimate, metadata_subset_entropy_reorder)
#V2 greater than 0.5 --> majority kazumbe ancestry
#V2 less than 0.5 --> majority kazumbe ancestry
#there are no individuals with exactly 50/50 kazume/polyodon ancestry
k2_q_estimate_with_metadata$species_ID_entropy[k2_q_estimate_with_metadata$V2 > 0.5] <- 'P. kazumbe'
k2_q_estimate_with_metadata$species_ID_entropy[k2_q_estimate_with_metadata$V2 < 0.5] <- 'P. polyodon'

#reduce down to taxa found in the genotype matrix
subset_metadata <- k2_q_estimate_with_metadata[k2_q_estimate_with_metadata$sample_code %in% colnames(petro_geno_matrix), ]

#Adding a column of nice location names ordered from N --> S
subset_metadata$location_plotting[subset_metadata$location == 'Gombe_South'] <- "Gombe S"
subset_metadata$location_plotting[subset_metadata$location == 'Katongwe_N'] <- "Katongwe N"
subset_metadata$location_plotting[subset_metadata$location == 'Katongwe_S'] <- "Katongwe S"
subset_metadata$location_plotting[subset_metadata$location == 'Ska'] <- "S Kagango"
subset_metadata$location_plotting[subset_metadata$location == 'Kalala'] <- "Kalalangabo"
subset_metadata$location_plotting[subset_metadata$location == 'Nondwa'] <- "Nondwa"
subset_metadata$location_plotting[subset_metadata$location == 'Hilltop'] <- "Hilltop"
subset_metadata$location_plotting[subset_metadata$location == 'Bangwe'] <- "Bangwe"
subset_metadata$location_plotting[subset_metadata$location == 'Jakob'] <- "Jakobsen's S"
subset_metadata$location_plotting[subset_metadata$location == 'Ulombola'] <- "Ulombola"
#subset_metadata$location_plotting[subset_metadata$location == 'Harembe'] <- "Harembe"

subset_metadata$location_plotting <- factor(subset_metadata$location_plotting, 
                                                   levels= c("Gombe S", "Katongwe N", "Katongwe S", "S Kagango", "Kalalangabo", "Nondwa", "Hilltop", "Bangwe", "Jakobsen's S", "Ulombola"))


### Make sure that genotype matrix only includes those found in the entropy analyses ###
#only include P. kazumbe and P. polyodon
kazumbe_polyodon_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% k2_q_estimate_with_metadata[k2_q_estimate_with_metadata$species_ID_entropy %in% c('P. kazumbe', 'P. polyodon') , ]$sample_code]

#only include P. kazumbe
kazumbe_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% k2_q_estimate_with_metadata[which(k2_q_estimate_with_metadata$species_ID_entropy == 'P. kazumbe'), ]$sample_code]

#only include P. polyodon 
polyodon_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% k2_q_estimate_with_metadata[which(k2_q_estimate_with_metadata$species_ID_entropy == 'P. polyodon'), ]$sample_code]



####################
### PCA ANALYSES ###
####################

### PCAs ###
kazumbe_polyodon_PCA <- doPCA(g = kazumbe_polyodon_geno_matrix, dims = 10)
kazumbe_PCA <- doPCA(g = kazumbe_geno_matrix, dims = 10)
polyodon_PCA <- doPCA(g = polyodon_geno_matrix, dims = 10)


### Combine PCAs with metadata ###
kazumbe_polyodon_PCA_withmetadata <- merge(subset_metadata, kazumbe_polyodon_PCA$PCA, 
                                           by.x = 'sample_code', by.y = 'row.names', all.x = FALSE, all.y = TRUE)
kazumbe_PCA_withmetadata <- merge(subset_metadata, kazumbe_PCA$PCA, 
                                           by.x = 'sample_code', by.y = 'row.names', all.x = FALSE, all.y = TRUE)
polyodon_PCA_withmetadata <- merge(subset_metadata, polyodon_PCA$PCA, 
                                           by.x = 'sample_code', by.y = 'row.names', all.x = FALSE, all.y = TRUE)


### Variance explained by each PC ###
importance_info_kazumbe_polyodon <- as.data.frame(t(kazumbe_polyodon_PCA$imp$importance)) %>% 
  rownames_to_column() %>% 
  rename(pc = rowname,
         stdiv = 'Standard deviation',
         prop_variance = 'Proportion of Variance',
         cumulative_prop = 'Cumulative Proportion') %>% 
  mutate(dataset = 'combined_kazumbe_polyodon')

importance_info_polyodon <- as.data.frame(t(polyodon_PCA$imp$importance)) %>% 
  rownames_to_column() %>% 
  rename(pc = rowname,
         stdiv = 'Standard deviation',
         prop_variance = 'Proportion of Variance',
         cumulative_prop = 'Cumulative Proportion') %>% 
  mutate(dataset = 'polyodon')

importance_info_kazumbe <- as.data.frame(t(kazumbe_PCA$imp$importance)) %>% 
  rownames_to_column() %>% 
  rename(pc = rowname,
         stdiv = 'Standard deviation',
         prop_variance = 'Proportion of Variance',
         cumulative_prop = 'Cumulative Proportion') %>% 
  mutate(dataset = 'kazumbe')

combined_importance_info_df <- rbind(importance_info_kazumbe_polyodon, importance_info_polyodon, importance_info_kazumbe)


### Export the PCA results ###
write.csv(kazumbe_polyodon_PCA_withmetadata, here('PCAs', 'kazumbe_polyodon_PCA_withmetadata_9_7_2021.csv'), row.names = FALSE)
write.csv(polyodon_PCA_withmetadata, here('PCAs', 'polyodon_PCA_withmetadata_9_7_2021.csv'), row.names = FALSE)
write.csv(kazumbe_PCA_withmetadata, here('PCAs', 'kazumbe_PCA_withmetadata_9_7_2021.csv'), row.names = FALSE)

write.csv(combined_importance_info_df, here('PCAs', 'pca_importance_info.csv'), row.names = FALSE)



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# ### COLOR PALETTES ###
# #qualitative color palette
# location_palette <- brewer.pal(length(unique(kazumbe_polyodon_PCA_withmetadata$location_plotting)),"Set3")
# names(location_palette) <- unique(kazumbe_polyodon_PCA_withmetadata$location_plotting)
# location_palette[location_palette == '#FFFFB3'] <- '#F7DC6F' #switch the yellow out
# location_palette[location_palette == "#D9D9D9"] <- "#AAB7B8" #switch the gray out
# location_palette[names(location_palette) == "Hilltop"] <- "#996633"
# location_palette[names(location_palette) == "Ulombola"] <- "#0000e6"
# location_palette[names(location_palette) == "Nondwa"] <- "#4dd2ff"
# location_palette[names(location_palette) == "Gombe S"] <- "#66ff66"
# #location_palette[names(location_palette) == "Harembe"] <- "#99ffeb"
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
# 
# ### PLOTTING THE PCAS ###
# 
# #PCA containing both species
# species_PCA <- ggplot(data = kazumbe_polyodon_PCA_withmetadata, aes(x = PC1, y= PC2, color = species_ID_entropy)) +
#   #geom_point(size = 5) +
#   geom_point(size = 6, alpha = 0.5) +
#   geom_point(size = 6, shape = 21, alpha = 0.7, stroke = 1) +
#   scale_color_manual(name = "Species", values = c("#DC7633", "#3498DB")) +
#   #scale_color_manual(name = "Species", values = c("#ffc966", "#ccccff")) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 22),
#         axis.text = element_text(size = 13),
#         legend.title = element_text(size = 21),
#         legend.position="bottom",
#         legend.text=element_text(size = 18, face = "italic")) +
#   xlab(paste0("PC1 (", round(kazumbe_polyodon_PCA[[2]]$importance[2,1]*100, 2),"%)")) +
#   ylab(paste0("PC2 (", round(kazumbe_polyodon_PCA[[2]]$importance[2,2]*100, 2),"%)"))
# 
# 
# #PCA containing only kazumbe
# kazumbe_PCA_plot <- ggplot(data = kazumbe_PCA_withmetadata, aes(x = PC1, y= PC2, color = location_plotting)) +
#   geom_point(size = 4.5, alpha = 0.5) +
#   geom_point(size = 4.5, shape = 1, alpha = 0.5) +
#   scale_color_manual(name = "Location", values = location_palette[names(location_palette) %in% unique(kazumbe_PCA_withmetadata$location_plotting) ]) +
#   #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(kazumbe_PCA_combinedinfo$location_plotting) ]) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 22),
#         axis.text = element_text(size = 13),
#         legend.position = "none",
#         plot.title = element_text(hjust = 0.5, size = 27, face="italic")) +
#   ggtitle("P. kazumbe") +
#   xlab(paste0("PC1 (", round(kazumbe_PCA[[2]]$importance[2,1]*100, 2),"%)")) +
#   ylab(paste0("PC2 (", round(kazumbe_PCA[[2]]$importance[2,2]*100, 2),"%)"))
# 
# #PCA containing only polyodon
# polyodon_PCA_plot <- ggplot(data = polyodon_PCA_withmetadata, aes(x = PC1, y= PC2, color = location_plotting)) +
#   geom_point(size = 4.5, alpha = 0.5) +
#   geom_point(size = 4.5, shape = 1, alpha = 0.5) +
#   scale_color_manual(name = "Location", values = location_palette[names(location_palette) %in% unique(polyodon_PCA_withmetadata$location_plotting) ]) +
#   #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(polyodon_PCA_combinedinfo$location_plotting) ]) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 22),
#         axis.text = element_text(size = 13),
#         legend.position = "none",
#         plot.title = element_text(hjust = 0.5, size = 27, face= c("italic") )) +
#   ggtitle("P. polyodon") +
#   xlab(paste0("PC1 (", round(polyodon_PCA[[2]]$importance[2,1]*100, 2),"%)")) +
#   ylab(paste0("PC2 (", round(polyodon_PCA[[2]]$importance[2,2]*100, 2),"%)"))
# 
# #extracting legend
# location_legend <- get_legend(ggplot(data = kazumbe_polyodon_PCA_withmetadata, aes(x = PC1, y= PC2, color = location_plotting)) +
#                                 geom_point(size = 4.5, alpha = 0.5) +
#                                 geom_point(size = 4.5, shape = 1, alpha = 0.5) +
#                                 scale_color_manual(name = "Location", values = location_palette[names(location_palette) %in% unique(kazumbe_polyodon_PCA_withmetadata$location_plotting) ]) +
#                                 #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(petro_PCA_combinedinfo$location_plotting) ]) +
#                                 theme_bw() +
#                                 theme(panel.grid.major = element_blank(),
#                                       panel.grid.minor = element_blank(),
#                                       axis.title = element_text(size = 18),
#                                       plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
#                                       legend.title = element_text(size = 21),
#                                       legend.text=element_text(size = 18)))
# 
# 
# ### CREATE THE MULTI-PANEL PLOT ###
# multipanel_kazumbe_polyodon_PCAs <- plot_grid(kazumbe_PCA_plot, polyodon_PCA_plot, nrow = 2, labels = c("(b)", "(c)"), label_size = 26.5)
# multipanel_kazumbe_polyodon_PCAs_withlegend <- plot_grid(multipanel_kazumbe_polyodon_PCAs, 
#                                                          location_legend, ncol = 2, rel_widths = c(0.9, 0.25))
# 
# plot_grid(species_PCA, multipanel_kazumbe_polyodon_PCAs_withlegend, rel_widths = c(0.6, 0.4), labels = c('(a)', ''), label_size = 26.5)
# 
# ### Export the multi-panel plot ###
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/Figures/petro_multipanelPCA_8_19_2020.png", width = 22, height = 12)

