####################
### PCA analysis ###
####################


#########################
### PREPARATION STEPS ###
#########################

### LOAD RELEVANT LIBRARIES ###
library(tidyverse)
library(cowplot)
library(adegenet)
library(vcfR)
library(RColorBrewer)
library(rhdf5)
library(abind)



### LOAD DATA ###

#metadata
cichlid_metadata <- read.csv("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/working_metadata_file/monster_23jul20.csv", stringsAsFactors = FALSE)

#vcf file (this file was already processed with vcftools)
petro_vcfR <- read.vcfR("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/VCFs/kazumbe_polyodon_pundamilia_8_19_2020_maf0.01_maxmissing0.70_minDP5_post_4.recode.vcf")


### ENTROPY ###
entropy_samples_order <- read.delim("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_8_27_2020/NAMES_kazumbe_polyodon_pundamilia_8_19_2020_maf0.01_maxmissing0.70_minDP5_post_4.recode.txt", header = FALSE)
metadata_subset_entropy <- cichlid_metadata[cichlid_metadata$sample_code %in% entropy_samples_order$V1, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2')]
metadata_subset_entropy_reorder <- metadata_subset_entropy[match(as.character(entropy_samples_order$V1), metadata_subset_entropy$sample_code),]


hdf5_file_path <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_8_27_2020/hdf5_files/'
#K = 2
k2_rep1.q <- h5read(paste0(hdf5_file_path, 'kazumbe_polyodon_pundamilia_8_19_2020_2rep1.hdf5'), "q")
k2_rep2.q <- h5read(paste0(hdf5_file_path, 'kazumbe_polyodon_pundamilia_8_19_2020_2rep2.hdf5'), "q")
k2_rep3.q <- h5read(paste0(hdf5_file_path, 'kazumbe_polyodon_pundamilia_8_19_2020_2rep3.hdf5'), "q")
k2_combined.q <- abind(k2_rep1.q, k2_rep2.q, k2_rep3.q, along=1)
k2_q_estimate <- as.data.frame(t(apply(k2_combined.q, 2:3, mean)))


### FUNCTIONS ###

#PCA function
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

#convert vcf to genlight abject
petro_gen <- vcfR2genlight(petro_vcfR)

#extracting genotype matrix
petro_geno_matrix <- t(as.matrix(petro_gen))

### PROCESS METADATA ###
#cichlid_metadata_reduced <- cichlid_metadata[, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2')]
#cichlid_metadata_reduced$reduced_species_names <- cichlid_metadata_reduced$sciname2
#cichlid_metadata_reduced$reduced_species_names[cichlid_metadata_reduced$sciname2 == "Petrochromis kazumbe"] <- 'P. kazumbe'
#cichlid_metadata_reduced$reduced_species_names[cichlid_metadata_reduced$sciname2 == "Petrochromis polyodon"] <- 'P. polyodon'

k2_q_estimate_with_metadata <- cbind(k2_q_estimate, metadata_subset_entropy_reorder)
#V2 greater than 0.5 --> majority kazumbe ancestry
#V2 less than 0.5 --> majority kazumbe ancestry
#there are no individuals with exactly 50/50 kazume/polyodon ancestry
k2_q_estimate_with_metadata$species_ID_entropy[k2_q_estimate_with_metadata$V2 > 0.5] <- 'P. kazumbe'
k2_q_estimate_with_metadata$species_ID_entropy[k2_q_estimate_with_metadata$V2 < 0.5] <- 'P. polyodon'


#reduce down to taxa found in the genotype matrix
#subset_metadata <- cichlid_metadata_reduced[cichlid_metadata_reduced$sample_code %in% colnames(petro_geno_matrix), ]
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

#only include P. kazumbe and P. polyodon 
#kazumbe_polyodon_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% cichlid_metadata_reduced[cichlid_metadata_reduced$ %in% c('kazumbe', 'polyodon') , ]$sample_code]
kazumbe_polyodon_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% k2_q_estimate_with_metadata[k2_q_estimate_with_metadata$species_ID_entropy %in% c('P. kazumbe', 'P. polyodon') , ]$sample_code]


#write.table(as.data.frame(colnames(kazumbe_polyodon_geno_matrix)), '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/species_lists_and_metadata/NAMES_polyodon_kazumbe_vcfpost_4_8_19_2020.txt', 
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

#only include P. kazumbe
#kazumbe_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% cichlid_metadata_reduced[which(cichlid_metadata_reduced$species2 == 'kazumbe'), ]$sample_code]
kazumbe_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% k2_q_estimate_with_metadata[which(k2_q_estimate_with_metadata$species_ID_entropy == 'P. kazumbe'), ]$sample_code]

#only include P. polyodon 
#polyodon_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% cichlid_metadata_reduced[which(cichlid_metadata_reduced$species2 == 'polyodon'), ]$sample_code]
polyodon_geno_matrix <- petro_geno_matrix[, colnames(petro_geno_matrix) %in% k2_q_estimate_with_metadata[which(k2_q_estimate_with_metadata$species_ID_entropy == 'P. polyodon'), ]$sample_code]


################
### ANALYSES ###
################

### PCAs ###
#petro_PCA <- doPCA(g = petro_geno_matrix, dims = 10)
kazumbe_polyodon_PCA <- doPCA(g = kazumbe_polyodon_geno_matrix, dims = 10)
kazumbe_PCA <- doPCA(g = kazumbe_geno_matrix, dims = 10)
polyodon_PCA <- doPCA(g = polyodon_geno_matrix, dims = 10)

#petro_PCA_withmetadata <- merge(subset_metadata, petro_PCA$PCA, by.x = 'sample_code', by.y = 'row.names')
kazumbe_polyodon_PCA_withmetadata <- merge(subset_metadata, kazumbe_polyodon_PCA$PCA, 
                                           by.x = 'sample_code', by.y = 'row.names', all.x = FALSE, all.y = TRUE)
kazumbe_PCA_withmetadata <- merge(subset_metadata, kazumbe_PCA$PCA, 
                                           by.x = 'sample_code', by.y = 'row.names', all.x = FALSE, all.y = TRUE)
polyodon_PCA_withmetadata <- merge(subset_metadata, polyodon_PCA$PCA, 
                                           by.x = 'sample_code', by.y = 'row.names', all.x = FALSE, all.y = TRUE)



######################
### VISUALIZATIONS ###
######################

### COLOR PALETTES ###
#qualitative color palette
location_palette <- brewer.pal(length(unique(kazumbe_polyodon_PCA_withmetadata$location_plotting)),"Set3")
names(location_palette) <- unique(kazumbe_polyodon_PCA_withmetadata$location_plotting)
location_palette[location_palette == '#FFFFB3'] <- '#F7DC6F' #switch the yellow out
location_palette[location_palette == "#D9D9D9"] <- "#AAB7B8" #switch the gray out
location_palette[names(location_palette) == "Hilltop"] <- "#996633"
location_palette[names(location_palette) == "Ulombola"] <- "#0000e6"
location_palette[names(location_palette) == "Nondwa"] <- "#4dd2ff"
location_palette[names(location_palette) == "Gombe S"] <- "#66ff66"
#location_palette[names(location_palette) == "Harembe"] <- "#99ffeb"

#switch the colors of Katongwe N and Ulombola
Ulombola_col_switch <- location_palette[names(location_palette) == "Ulombola"]
Katongwe_N_col_switch <- location_palette[names(location_palette) == "Katongwe N"]

location_palette[names(location_palette) == "Ulombola"] <- Katongwe_N_col_switch
location_palette[names(location_palette) == "Katongwe N"] <- Ulombola_col_switch


#switch Jakobsen's S and Bangwe
Jakob_col_switch <- location_palette[names(location_palette) == "Jakobsen's S"]
Bangwe_col_switch <- location_palette[names(location_palette) == "Bangwe"]

location_palette[names(location_palette) == "Jakobsen's S"] <- Bangwe_col_switch
location_palette[names(location_palette) == "Bangwe"] <- Jakob_col_switch


### PLOTTING THE PCAS ###

#PCA containing both species
species_PCA <- ggplot(data = kazumbe_polyodon_PCA_withmetadata, aes(x = PC1, y= PC2, color = species_ID_entropy)) +
  #geom_point(size = 5) +
  geom_point(size = 6, alpha = 0.5) +
  geom_point(size = 6, shape = 21, alpha = 0.7, stroke = 1) +
  scale_color_manual(name = "Species", values = c("#DC7633", "#3498DB")) +
  #scale_color_manual(name = "Species", values = c("#ffc966", "#ccccff")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 21),
        legend.position="bottom",
        legend.text=element_text(size = 18, face = "italic")) +
  xlab(paste0("PC1 (", round(kazumbe_polyodon_PCA[[2]]$importance[2,1]*100, 2),"%)")) +
  ylab(paste0("PC2 (", round(kazumbe_polyodon_PCA[[2]]$importance[2,2]*100, 2),"%)"))


#PCA containing only kazumbe
kazumbe_PCA_plot <- ggplot(data = kazumbe_PCA_withmetadata, aes(x = PC1, y= PC2, color = location_plotting)) +
  geom_point(size = 4.5, alpha = 0.5) +
  geom_point(size = 4.5, shape = 1, alpha = 0.5) +
  scale_color_manual(name = "Location", values = location_palette[names(location_palette) %in% unique(kazumbe_PCA_withmetadata$location_plotting) ]) +
  #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(kazumbe_PCA_combinedinfo$location_plotting) ]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 13),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 27, face="italic")) +
  ggtitle("P. kazumbe") +
  xlab(paste0("PC1 (", round(kazumbe_PCA[[2]]$importance[2,1]*100, 2),"%)")) +
  ylab(paste0("PC2 (", round(kazumbe_PCA[[2]]$importance[2,2]*100, 2),"%)"))

#PCA containing only polyodon
polyodon_PCA_plot <- ggplot(data = polyodon_PCA_withmetadata, aes(x = PC1, y= PC2, color = location_plotting)) +
  geom_point(size = 4.5, alpha = 0.5) +
  geom_point(size = 4.5, shape = 1, alpha = 0.5) +
  scale_color_manual(name = "Location", values = location_palette[names(location_palette) %in% unique(polyodon_PCA_withmetadata$location_plotting) ]) +
  #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(polyodon_PCA_combinedinfo$location_plotting) ]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 13),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 27, face= c("italic") )) +
  ggtitle("P. polyodon") +
  xlab(paste0("PC1 (", round(polyodon_PCA[[2]]$importance[2,1]*100, 2),"%)")) +
  ylab(paste0("PC2 (", round(polyodon_PCA[[2]]$importance[2,2]*100, 2),"%)"))

#extracting legend
location_legend <- get_legend(ggplot(data = kazumbe_polyodon_PCA_withmetadata, aes(x = PC1, y= PC2, color = location_plotting)) +
                                geom_point(size = 4.5, alpha = 0.5) +
                                geom_point(size = 4.5, shape = 1, alpha = 0.5) +
                                scale_color_manual(name = "Location", values = location_palette[names(location_palette) %in% unique(kazumbe_polyodon_PCA_withmetadata$location_plotting) ]) +
                                #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(petro_PCA_combinedinfo$location_plotting) ]) +
                                theme_bw() +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.title = element_text(size = 18),
                                      plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
                                      legend.title = element_text(size = 21),
                                      legend.text=element_text(size = 18)))


### CREATE THE MULTI-PANEL PLOT ###
multipanel_kazumbe_polyodon_PCAs <- plot_grid(kazumbe_PCA_plot, polyodon_PCA_plot, nrow = 2, labels = c("(b)", "(c)"), label_size = 26.5)
multipanel_kazumbe_polyodon_PCAs_withlegend <- plot_grid(multipanel_kazumbe_polyodon_PCAs, 
                                                         location_legend, ncol = 2, rel_widths = c(0.9, 0.25))

plot_grid(species_PCA, multipanel_kazumbe_polyodon_PCAs_withlegend, rel_widths = c(0.6, 0.4), labels = c('(a)', ''), label_size = 26.5)

### Export the multi-panel plot ###
ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/Figures/petro_multipanelPCA_8_19_2020.png", width = 22, height = 12)


### Export the PCA dataframes ###
write.csv(kazumbe_polyodon_PCA_withmetadata, "/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/kazumbe_polyodon_PCA_withmetadata.csv", row.names = FALSE)
write.csv(polyodon_PCA_withmetadata, "/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/polyodon_PCA_withmetadata.csv", row.names = FALSE)
write.csv(kazumbe_PCA_withmetadata, "/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/PCAs/kazumbe_PCA_withmetadata.csv", row.names = FALSE)


