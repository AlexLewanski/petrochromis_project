##############################################################################################
### SCRIPT NAME: nucleotide_diversity_kazumbe_subset_plot.R
### PURPOSE: visualizing the nucleotide diversity estimates calculated on subset of kazumbe
### PRODUCTS:
###     pnucleotide_diversity_with_kazumbe_subset.pdf: figure of nucleotide estimates
##############################################################################################


#####################
### SCRIPT SET-UP ###
#####################
library(here)
library(tidyverse)
library(cowplot)

pi_info <- read.csv(here('diversity_divergence_stats', 'processed_results', 'pi_7_2021.csv'))
pi_info_kazumbe_subsample <- read.csv(here('diversity_divergence_stats', 'processed_results', 'pi_kazumbe_subsample.csv'))

pi_df_combined <- pi_info %>% 
  mutate(replicate_val = paste0('full_', species) ) %>% 
  rbind(., pi_info_kazumbe_subsample)

div_with_kazumbe_subset_plot <- ggplot(data = pi_df_combined %>%
         mutate(Population = factor(Population, levels = rev(c('Full', LETTERS[10:1]) )) ) %>% 
         mutate(metric = case_when(species == 'kazumbe' ~ '\u03C0 (kazumbe)',
                                   species == 'polyodon' ~ '\u03C0 (polyodon)'),
                type = if_else(replicate_val %in% c('full_kazumbe', 'full_polyodon'), 'full', 'subsamp')) %>% 
         filter(!(Population %in% c('H', 'J')))) +
  geom_linerange(aes(x = Population, ymin = lower, ymax = upper, group = replicate_val, color = species, alpha = type), 
                 size = 1.5, position = position_dodge(0.25)) +
  geom_point(aes(x = Population, y = mean, group = replicate_val, color = species, alpha = type), 
              size = 5, position = position_dodge(0.25)) +
  scale_color_manual(values = c("#DC7633", "#3498DB"),
                     labels = c(expression(paste(italic("P"), ". cf. ",  italic("polyodon"))),
                                expression(paste(italic("P"), ". sp. 'kazumbe'")))) +
  scale_alpha_manual(values = c(1, 0.3), guide = "none") +
  xlab("Population") + ylab("\u03C0") +
  theme_cowplot() +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5),
        panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
        panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
        legend.title = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 23),
        axis.text.x = element_text(size = 15),
        legend.position="bottom",
        legend.spacing.x = unit(0.2, 'cm'),
        legend.text = element_text(size = 13, margin = margin(r = 6, unit = "pt")))

#ggsave(filename = here('figures', 'nucleotide_diversity_with_kazumbe_subset.png'),
#       plot = div_with_kazumbe_subset_plot, 
#       width = 12, height = 7, bg = 'white')

#pdf does not like the pi symbol; it gives the following warning and prints out dots:
#18: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label),  ... :
#                            conversion failure on 'π' in 'mbcsToSbcs': dot substituted for <80>
#https://stackoverflow.com/questions/63985427/warning-in-grid-callc-textbounds-as-graphicsannotxlabel-xx-xy
#to fix, use cairo_pdf as the graphics device:
#https://github.com/wilkelab/cowplot/issues/174
#https://github.com/rstudio/rstudio/issues/3680
#https://stackoverflow.com/questions/61476002/cannot-save-ggplot-with-umlauts-in-expressions-with-cairo-pdf

ggsave2(filename = here('figures', 'nucleotide_diversity_with_kazumbe_subset.pdf'),
       plot = div_with_kazumbe_subset_plot, device = cairo_pdf,
       width = 12, height = 7, bg = 'white')
