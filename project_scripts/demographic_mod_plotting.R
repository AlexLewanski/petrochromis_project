###################################################################################################
### SCRIPT NAME: demographic_mod_plotting.R
### PURPOSE: Visualizing and summarizing fastsimcoal2 modeling results
### PRODUCTS:
###     fastsimcoal_mod_results_pluslegend_9_30_2021.png: main text figure showing visualizations
###                                                       of models (median bootstraps values) and
###                                                       relative fits of models (likelihood 
###                                                       distributions, AIC)
###     fsc_best_mods_params_table.txt: table of parameter estimates from the best fitting models
###                                     based on the parameter estimates of the bootstraps. This
###                                     table is included in the supplementary material.
###     fsc_mod_fit_table.txt: table of relative fit of all tested models including likelihood,
###                            AIC, and relative likelihood. This tabkle is included in the 
###                            supplementary material.
###     fastsimcoal_allmods.png: multipanel plot of visualizations of all tested models. This figure
###                              is included in the supplementary material.
###     fsc_param_info.txt: table of parameter search ranges used in modeling
###     fsc_north_region_topmod_2dsfs.png: visualization of how well the expected 2D SFS of the best 
###                                        fitting model in the north region fits compared to the
###                                        observed 2D SFS. This figure is included in the
###                                        supplementary material.
###     fsc_north_region_topmod_1dsfs.png: visualization of how well the expected marginal 1D SFSs of
###                                        of the best fitting model in the north region fits
###                                        compared to the observed marginal 1D SFSs. This figure is
###                                        included in the supplementary material.
###     fsc_mid_region_topmod_2dsfs.png: visualization of how well the expected 2D SFS of the best 
###                                      fitting model in the mid region fits compared to the observed
###                                      2D SFS. This figure is included in the supplementary material.
###     fsc_mid_region_topmod_1dsfs.png: visualization of how well the expected marginal 1D SFSs of
###                                      of the best fitting model in the mid region fits compared to
###                                      the observed marginal 1D SFSs. This figure is included in the
###                                      supplementary material.
###################################################################################################


####################################################
####################################################
### PART 1: VISUALIZATION OF FASTSIMCOAL2 MODELS ###
####################################################
####################################################

#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(cowplot)
library(tidyverse)
library(reshape2)
library(egg)
library(here)
library(scales)
library(kableExtra)


#load visualization functions
source(here('project_scripts', 'fsc_plotting_functions.R'))
source(here('project_scripts', 'sfs_compare_vis_funcs.R')) #needed for Part 1 for relative likelihood function


### Loading results (and some initial processing) ###
#par files
par_file_list <- lapply(setNames(nm = c('mid', 'north')), function(REGION) {
  par_file_names1 <- list.files(here('demographic_modelling', 'fastsimcoal', 'results_files', 'par', REGION, 'par_files'))
  par_file_names <- par_file_names1[grepl(pattern = 'continuous|no|updated', par_file_names1)]
  lapply(setNames(nm = par_file_names), function(FILE, REGION) {
    maxLpar <- readLines(paste0(here('demographic_modelling', 'fastsimcoal', 'results_files', 'par', REGION, 'par_files'), '/',FILE))
    return(parse_fsc_par(par = maxLpar, start = 2))
  }, REGION = REGION)
})

#AIC
aic_list <- lapply(setNames(nm = c('mid', 'north')), function(REGION) {
  #read.table(here('demographic_modelling', 'fastsimcoal', 'results_files', 'model_fit_results', REGION, 'allmodels_popgrow_altspec.AIC'), header = TRUE)
  read.table(here('demographic_modelling', 'fastsimcoal', 'results_files', 'model_fit_results', REGION, 'allmodels.AIC'), header = TRUE)
})


#likelihood values
lik_list <- lapply(setNames(nm = c('mid', 'north')), function(REGION) {
  lhood_vec1 <- list.files(here('demographic_modelling', 'fastsimcoal', 'results_files', 'model_fit_results', REGION, 'likelihood_dir', 'popgrow_alt_spec', 'likelihood_dir'))
  lhood_vec <- lhood_vec1[grepl('continuous|no|updated', lhood_vec1)]
  lhood_vec_named <- setNames(lhood_vec, nm = gsub(".lhoods", "", lhood_vec))
  
  lhood_load <-lapply(lhood_vec_named, function(name, REGION) {
    lhoods <- read.delim(here('demographic_modelling', 'fastsimcoal', 'results_files', 'model_fit_results', REGION, 'likelihood_dir', 'popgrow_alt_spec', 'likelihood_dir', name), header = FALSE)
    lhoods %>% 
      rename(likelihood = V1) %>% 
      mutate(model = gsub(".lhoods", "", name))
  }, REGION = REGION) %>% 
    bind_rows()
})

best_mod_info <- lapply(setNames(nm = c('north', 'mid')), function(REGION) {
  read.table(here('demographic_modelling', 'fastsimcoal', 'results_files', 'parameter_estimates', REGION, 'recent_geneflow_dualpopgrowth_altspec_updated.bestlhoods'), header = TRUE)
})

bootstrap_results <- lapply(setNames(nm = c('north', 'mid')), function(REGION) {
  read.table(here('demographic_modelling', 'fastsimcoal', 'results_files', 'parameter_estimates', REGION, 'bs_bestrun_param.txt'), header = TRUE)
})


#loading est files for reporting search ranges
all_files <- list.files(here('demographic_modelling', 'fastsimcoal', 'input_files', 'final_est'))
tpl_est_files_alt <- all_files[str_detect(all_files, pattern = "altspec.est$")]
tpl_est_file_list <- lapply(setNames(tpl_est_files_alt, nm = gsub(".est", "", tpl_est_files_alt)), function(FILE) {
  readLines(paste(here('demographic_modelling', 'fastsimcoal', 'input_files', 'final_est'), '/', FILE, sep = ""), warn = FALSE)
})


### Custom function ###
process_params <- function(raw_est) {
  
  start_pos <- which(str_detect(raw_est, '^// Priors and rules file')) + 5
  end_pos <- min(which(str_detect(raw_est, '^$'))[which(str_detect(raw_est, '^$')) > start_pos]) - 1
  
  estimated_params <- raw_est[start_pos:end_pos]
  
  param_df <- as.data.frame(do.call(rbind, str_split(gsub(" +", " ", gsub("\\t", " ", estimated_params)) , pattern = " ")))
  colnames(param_df) <- c("Value_Type", "Parameter", "Distribution", "Lower", "Upper", "include_in_output")
  
  return(
    param_df %>%
      mutate(Distribution = case_when(Distribution == 'logunif' ~ "Log-uniform",
                                      Distribution == 'unif' ~ "Uniform"),
             Value_Type = case_when(Value_Type == "1" ~ "Integer",
                                    Value_Type == "0" ~ "Float"),
             Bounded = "No") %>%
      select(!include_in_output) %>%
      rename(`Value Type` = Value_Type)
  )
}



####################################################
### Main fsc model results figure (in main text) ###
####################################################

#extract best fitting model for each region (the recent, asymmetric gene flow model)
best_fit_mods_par <- lapply(par_file_list, function(region) region[['recent_geneflow_dualpopgrowth_altspec_updated_maxL.par']])

#update the parameter values in the best fitting model par to reflect the median bootstrap values
#the complex parameters (CHANGM, growth rates) are recalculated based on the median parameter estimates
#of the simple paramters
best_fit_mods_par_bs_update <- lapply(setNames(nm = names(best_fit_mods_par)), function(REGION, PAR, BOOTSTRAP) {
  
  processed_bs_params <- as.data.frame(t(apply(BOOTSTRAP[[REGION]], 2, function(x) quantile(x, probs = 0.5 ))))
  
  PAR[[REGION]]$pop_size <- c(processed_bs_params$N_POP1, processed_bs_params$N_POP2)
  PAR[[REGION]]$growth_rate <- c(log(processed_bs_params$N_ANC_POP1/processed_bs_params$N_POP1)/processed_bs_params$TDIV, 
                                 log(processed_bs_params$N_ANC_POP2/processed_bs_params$N_POP2)/processed_bs_params$TDIV)
  PAR[[REGION]]$mig_mat_list[[1]][2, 1] <- processed_bs_params$MIG10
  PAR[[REGION]]$mig_mat_list[[1]][1, 2] <- processed_bs_params$MIG01
  
  #PAR[[REGION]]$hist_event_info$processed_hist_events$time[1] <- processed_bs_params$CHANGM
  PAR[[REGION]]$hist_event_info$processed_hist_events$time[1] <- processed_bs_params$TDIV*processed_bs_params$TPROP
  
  PAR[[REGION]]$hist_event_info$processed_hist_events$time[2] <- processed_bs_params$TDIV
  PAR[[REGION]]$hist_event_info$processed_hist_events$new_deme_size[2] <- processed_bs_params$RSANC
  
  return(PAR[[REGION]])
}, PAR = best_fit_mods_par, BOOTSTRAP = bootstrap_results)


# #visualize the best fitting models based on the median bootstrap estimates
# best_fit_mods_vis <- lapply(best_fit_mods_par_bs_update, function(par) {
#   plot_demographic_mod_popgrow(par = par,
#                                colors = TRUE,
#                                pop_colors = c("#DC7633", "#3498DB"),
#                                scale_pops = TRUE,
#                                pop_color = '#a3a3a3',
#                                generic_plotting = FALSE,
#                                label_size = 5.5,
#                                perspective = c("forward", "backward")[1])
# })


# annotated_demo_plot_list <- lapply(setNames(nm = c('north', 'mid')), function(REGION, plot_list) {
#   
#   region_upper_case <- switch(REGION, mid = {'Mid'}, north = {'North'})
#   
#   plot_list[[REGION]]$plot +
#     geom_segment(data = data.frame(x1 = plot_list[[REGION]]$right_side, x2 = plot_list[[REGION]]$right_side, 
#                                    y1 = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'TDIV', "Lower CI"]), 
#                                    y2 = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'TDIV', "Upper CI"])),
#                  aes(x = x1, y = y1, xend = x2, yend = y2),
#                  lineend = 'round', size = 2.3) +
#     geom_point(data = data.frame(x = plot_list[[REGION]]$right_side, y = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'TDIV', "Median"])), 
#                aes(x = x, y = y), size = 5) +
#     geom_segment(data = data.frame(x1 = plot_list[[REGION]]$right_side, x2 = plot_list[[REGION]]$right_side, 
#                                    y1 = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'CHANGM', "Lower CI"]), 
#                                    y2 = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'CHANGM', "Upper CI"])),
#                  aes(x = x1, y = y1, xend = x2, yend = y2),
#                  lineend = 'round', size = 2.3) +
#     geom_point(data = data.frame(x = plot_list[[REGION]]$right_side, y = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'CHANGM', "Median"])), 
#                aes(x = x, y = y), size = 5) +
#     ylim(0, 308000)
#   
# }, plot_list = best_fit_mods_vis)


aic_list_processed <- lapply(setNames(nm = names(aic_list)), function(REGION, aic_list) {
  aic_list[[REGION]] %>%
    mutate(region = REGION,
           first_word = gsub("_.*", '', model),
           deltaAIC = AIC - min(AIC),
           #deltaLik = MaxObsLhood - MaxEstLhood,
           relative_likelihood = relative_lik(AIC),
           relative_likelihood_updated = case_when(relative_likelihood == 0 ~ "0",
                                                   relative_likelihood < 1e-20 ~ "approx. 0",
                                                   relative_likelihood > 1e-20 ~ as.character(relative_likelihood)),
           mod_class = factor(case_when(first_word == 'no' ~ 'No gene flow',
                                        first_word == 'continuous' ~ 'Continuous gene flow',
                                        first_word == 'recent' ~ 'Recent gene flow',
                                        first_word == 'early' ~ 'Early gene flow'),
                              levels = c('Continuous gene flow', 'Early gene flow', 'Recent gene flow', 'No gene flow')),
           mod_name = factor(
             case_when(model == "no_geneflow_dualpopgrowth_altspec" ~ "no gf",
                       model == "continuous_geneflow_dualpopgrowth_altspec" ~ "contin. gf (asym.)",
                       model == "early_geneflow_dualpopgrowth_altspec_updated" ~ "early gf (asym.)",
                       model == "recent_geneflow_dualpopgrowth_altspec_updated" ~ "recent gf (asym.)",
                       model == "early_geneflow_symmetric_dualpopgrowth_altspec_updated" ~ "early gf (sym.)",
                       model == "recent_geneflow_symmetric_dualpopgrowth_altspec_updated" ~ "recent gf (sym.)",
                       model == "continuous_geneflow_symmetric_dualpopgrowth_altspec" ~ "contin. gf (sym.)",
                       model == "early_geneflow_unidirect0_dualpopgrowth_altspec_updated" ~ "early gf (Pk to Pp)",
                       model == "early_geneflow_unidirect1_dualpopgrowth_altspec_updated" ~ "early gf (Pp to Pk)",
                       model == "continuous_geneflow_unidirect0_dualpopgrowth_altspec" ~ "contin. gf (Pk to Pp)",
                       model == "continuous_geneflow_unidirect1_dualpopgrowth_altspec" ~ "contin. gf (Pp to Pk)",
                       model == "recent_geneflow_unidirect0_dualpopgrowth_altspec_updated" ~ "recent gf (Pk to Pp)",
                       model == "recent_geneflow_unidirect1_dualpopgrowth_altspec_updated" ~ "recent gf (Pp to Pk)"),
             levels = rev(c("contin. gf (sym.)", "contin. gf (asym.)", "contin. gf (Pk to Pp)", "contin. gf (Pp to Pk)", "early gf (sym.)", "early gf (asym.)", "early gf (Pk to Pp)", "early gf (Pp to Pk)", "recent gf (sym.)", "recent gf (asym.)", "recent gf (Pk to Pp)", "recent gf (Pp to Pk)", "no gf"))
             
           )) 
}, aic_list = aic_list)

aic_list_processed_df <- do.call(rbind, aic_list_processed)


lik_list_processed <- lapply(lik_list, function(REGION) {
  REGION %>%
    mutate(first_word = gsub("_.*", '', model)) %>% 
    filter(first_word != 'change') %>% 
    mutate(mod_class = factor(case_when(first_word == 'no' ~ 'No gene flow',
                                        first_word == 'continuous' ~ 'Continuous gene flow',
                                        first_word == 'recent' ~ 'Recent gene flow',
                                        first_word == 'early' ~ 'Early gene flow'),
                              levels = c('Continuous gene flow', 'Early gene flow', 'Recent gene flow', 'No gene flow')),
           mod_name = factor(
             case_when(model == "no_geneflow_dualpopgrowth_altspec" ~ "no gf",
                       model == "continuous_geneflow_dualpopgrowth_altspec" ~ "contin. gf (asym.)",
                       model == "early_geneflow_dualpopgrowth_altspec_updated" ~ "early gf (asym.)",
                       model == "recent_geneflow_dualpopgrowth_altspec_updated" ~ "recent gf (asym.)",
                       model == "early_geneflow_symmetric_dualpopgrowth_altspec_updated" ~ "early gf (sym.)",
                       model == "recent_geneflow_symmetric_dualpopgrowth_altspec_updated" ~ "recent gf (sym.)",
                       model == "continuous_geneflow_symmetric_dualpopgrowth_altspec" ~ "contin. gf (sym.)",
                       model == "early_geneflow_unidirect0_dualpopgrowth_altspec_updated" ~ "early gf (Pk to Pp)",
                       model == "early_geneflow_unidirect1_dualpopgrowth_altspec_updated" ~ "early gf (Pp to Pk)",
                       model == "continuous_geneflow_unidirect0_dualpopgrowth_altspec" ~ "contin. gf (Pk to Pp)",
                       model == "continuous_geneflow_unidirect1_dualpopgrowth_altspec" ~ "contin. gf (Pp to Pk)",
                       model == "recent_geneflow_unidirect0_dualpopgrowth_altspec_updated" ~ "recent gf (Pk to Pp)",
                       model == "recent_geneflow_unidirect1_dualpopgrowth_altspec_updated" ~ "recent gf (Pp to Pk)"),
             levels = rev(c("contin. gf (sym.)", "contin. gf (asym.)", "contin. gf (Pk to Pp)", "contin. gf (Pp to Pk)", "early gf (sym.)", "early gf (asym.)", "early gf (Pk to Pp)", "early gf (Pp to Pk)", "recent gf (sym.)", "recent gf (asym.)", "recent gf (Pk to Pp)", "recent gf (Pp to Pk)", "no gf"))
    ))
})


fsc_results_list <- list(par = best_fit_mods_par_bs_update,
                         aic = aic_list_processed,
                         likelihood = lik_list_processed)

mod_type_cols <- setNames(c('#cb5cbb', '#FFC60A', '#2A9D8F', '#E24D28'), 
                          nm = c('Continuous gene flow', 'Early gene flow', 'Recent gene flow', 'No gene flow'))

processed_param_list <- lapply(setNames(nm = names(best_mod_info)), function(REGION, param_list) {
  
  processed_df <- t(apply(param_list[[REGION]], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975) ))) %>%
    as.data.frame() %>%
    rename("Lower CI" = `2.5%`,
           Median = `50%`,
           "Upper CI" = `97.5%`) %>% 
    rownames_to_column(var = "Parameter") %>%
    filter(!(Parameter %in% c("MaxEstLhood", "MaxObsLhood"))) %>%
    mutate(Description = case_when(Parameter == "N_POP1" ~ "kazumbe pop size",
                                   Parameter == "N_POP2" ~ "polyodon pop size",
                                   Parameter == "N_ANC_POP1" ~ "initial kazumbe pop size",
                                   Parameter == "N_ANC_POP2" ~ "initial polyodon pop size",
                                   Parameter == "MIG01" ~ "mig. rate, kazumbe to polyodon",
                                   Parameter == "MIG10" ~ "mig. rate, polyodon to kazumbe",
                                   Parameter == "TPROP" ~ "prop of TDIV to mig. cessation",
                                   Parameter == "RSANC" ~ "anc. pop size (relative to sink deme)",
                                   Parameter == "TDIV" ~ "time to divergence",
                                   Parameter == "CHANGM" ~ "time to mig. cessation",
                                   Parameter == "GR_POP1" ~ "kazumbe growth rate",
                                   Parameter == "GR_POP2" ~ "polyodon growth rate"),
           Region = switch(REGION, north = "North", mid = "Mid") ) %>% 
    relocate(Region, .before = Parameter)
  
  #update the complex parameters based on the lower CIs, medians, and upper CIs of the simple parameters used to calculate them
  processed_df[processed_df$Parameter == 'GR_POP1', c("Lower CI", "Median", "Upper CI")] <- log(processed_df[processed_df$Parameter == 'N_ANC_POP1', c("Lower CI", "Median", "Upper CI")]/processed_df[processed_df$Parameter == 'N_POP1', c("Lower CI", "Median", "Upper CI")])/processed_df[processed_df$Parameter == 'TDIV', c("Lower CI", "Median", "Upper CI")]
  processed_df[processed_df$Parameter == 'GR_POP2', c("Lower CI", "Median", "Upper CI")] <- log(processed_df[processed_df$Parameter == 'N_ANC_POP2', c("Lower CI", "Median", "Upper CI")]/processed_df[processed_df$Parameter == 'N_POP2', c("Lower CI", "Median", "Upper CI")])/processed_df[processed_df$Parameter == 'TDIV', c("Lower CI", "Median", "Upper CI")]
  processed_df[processed_df$Parameter == 'CHANGM', c("Lower CI", "Median", "Upper CI")] <- processed_df[processed_df$Parameter == 'TDIV', c("Lower CI", "Median", "Upper CI")]*processed_df[processed_df$Parameter == 'TPROP', c("Lower CI", "Median", "Upper CI")]
  
  return(
    list(table_df = processed_df %>% 
           mutate(Median = formatC(Median, format = "e", digits = 2),
                  "Lower CI" = formatC(`Lower CI`, format = "e", digits = 2),
                  "Upper CI" = formatC(`Upper CI`, format = "e", digits = 2)),
         raw_values = processed_df%>% 
           mutate(Median = formatC(Median, format = "f", digits = 7),
                  "Lower CI" = formatC(`Lower CI`, format = "f", digits = 7),
                  "Upper CI" = formatC(`Upper CI`, format = "f", digits = 7)))
    
  )
}, param_list = bootstrap_results)


region_param_info <- lapply(processed_param_list, function(x) x$raw_values) %>% 
  bind_rows()


fsc_vis_components_list <- lapply(setNames(nm = c('north', 'mid')), function(REGION, results_list, region_param_info, col) {
  
  processed_results <- list()
  
  demog_plot <- plot_demographic_mod_popgrow(par = results_list[['par']][[REGION]],
                                             colors = TRUE,
                                             pop_colors = c("#DC7633", "#3498DB"),
                                             scale_pops = TRUE,
                                             pop_color = '#a3a3a3',
                                             generic_plotting = FALSE,
                                             label_size = 5.5,
                                             label_nudge = 0.05,
                                             perspective = c("forward", "backward")[1])
  
  text_shift <- 1.03
  
  processed_results[['par']] <- demog_plot$plot
  
  region_upper_case <- switch(REGION, mid = {'Mid'}, north = {'North'})
  alt_region <- switch(REGION, north = {'Mid'}, mid = {'North'})
  
  processed_results[['annotated_par']] <- processed_results[['par']] +
    geom_segment(data = data.frame(x1 = demog_plot$right_side, x2 = demog_plot$right_side, 
                                   y1 = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'TDIV', "Lower CI"]), 
                                   y2 = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'TDIV', "Upper CI"])),
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 lineend = 'round', size = 2.3) +
    geom_point(data = data.frame(x = demog_plot$right_side, y = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'TDIV', "Median"])), 
               aes(x = x, y = y), size = 5) +
    geom_segment(data = data.frame(x1 = demog_plot$right_side * text_shift, x2 = demog_plot$right_side * text_shift, 
                                   y1 = as.numeric(region_param_info[region_param_info$Region == alt_region & region_param_info$Parameter == 'TDIV', "Lower CI"]), 
                                   y2 = as.numeric(region_param_info[region_param_info$Region == alt_region & region_param_info$Parameter == 'TDIV', "Upper CI"])),
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 lineend = 'round', size = 1.8, color = "#bfbfbf", alpha=1) +
    geom_point(data = data.frame(x = demog_plot$right_side * text_shift, y = as.numeric(region_param_info[region_param_info$Region == alt_region & region_param_info$Parameter == 'TDIV', "Median"])), 
               aes(x = x, y = y), size = 3.5, color = "#bfbfbf", alpha=1) +
    geom_segment(data = data.frame(x1 = demog_plot$right_side, x2 = demog_plot$right_side, 
                                   y1 = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'CHANGM', "Lower CI"]), 
                                   y2 = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'CHANGM', "Upper CI"])),
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 lineend = 'round', size = 2.3) +
    geom_point(data = data.frame(x = demog_plot$right_side, y = as.numeric(region_param_info[region_param_info$Region == region_upper_case & region_param_info$Parameter == 'CHANGM', "Median"])), 
               aes(x = x, y = y), size = 5) +
    geom_segment(data = data.frame(x1 = demog_plot$right_side * text_shift, x2 = demog_plot$right_side * text_shift, 
                                   y1 = as.numeric(region_param_info[region_param_info$Region == alt_region & region_param_info$Parameter == 'CHANGM', "Lower CI"]), 
                                   y2 = as.numeric(region_param_info[region_param_info$Region == alt_region & region_param_info$Parameter == 'CHANGM', "Upper CI"])),
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 lineend = 'round', size = 1.8, color = "#bfbfbf", alpha=1) +
    geom_point(data = data.frame(x = demog_plot$right_side * text_shift, y = as.numeric(region_param_info[region_param_info$Region == alt_region & region_param_info$Parameter == 'CHANGM', "Median"])), 
               aes(x = x, y = y), size = 3.5, color = "#bfbfbf", alpha=1) +
    ylim(0, 256000)
  
  
  processed_results[['aic']] <- results_list[['aic']][[REGION]] %>%
    mutate(deltaAIC = deltaAIC*-1) %>% 
    ggplot(aes(x = deltaAIC, y = mod_name)) +
    geom_point(aes(color = mod_class), size = 5) +
    xlab('-\u0394AIC') +
    scale_color_manual(values = col) +
    theme_cowplot() +
    theme(plot.margin = margin(5.5, 0, 0, 10),
          panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
          panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 17, margin = margin(t = 2, r = 11, b = 0, l = 0)),
          axis.text.x = element_text(size = 13.5),
          legend.text = element_text(size = 18, margin = margin(r = 22, unit = "pt")),
          legend.position = "none",
          legend.spacing.x = unit(0.3, 'cm'))
  
  processed_results[['likelihood']] <- results_list[['likelihood']][[REGION]] %>% 
    ggplot(aes(x = likelihood, y = mod_name)) +
    geom_violin(aes(fill = mod_class, color = mod_class), size = 1) +
    #xlab('Likelihood') +
    xlab(expression(log[10]*"(Likelihood)")) +
    #scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    theme_cowplot() +
    theme(plot.margin = margin(5.5, 10, 0, 0),
          panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
          panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 13.5),
          axis.title.x = element_text(size = 17, margin = margin(t = 2, r = 11, b = 0, l = 0)),
          legend.position = "none",
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.key = element_rect(fill = "white"))
  
  return(processed_results)
  
}, results_list = fsc_results_list, region_param_info = region_param_info, col = mod_type_cols)


region_title <- lapply(setNames(c("North region", "Mid region"), nm = c('north', 'mid')), function(REGION) {
  multipanel_title <- ggdraw() + 
    draw_label(
      REGION,
      x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
      #x = 0, hjust = 0,
      fontface = "bold", size = 27
    ) #+ theme(plot.margin = margin(0, 0, 0, 7))
})


multipanel_list <- lapply(setNames(nm = c('north', 'mid')), function (REGION, vis_list, region_title_list) {
  
  multipanel <- cowplot::plot_grid(vis_list[[REGION]]$annotated_par + theme(plot.margin = margin(5, 55, 5, 5)), 
                     ggdraw(egg::ggarrange(vis_list[[REGION]]$aic, vis_list[[REGION]]$likelihood, ncol = 2)) + theme(plot.margin = margin(5, 5, 5, 35)),
                     labels = switch(REGION, north = {c('(a)', '(b)')}, mid = {c('(c)', '(d)')}) , 
                     rel_widths = c(0.38, 0.62), label_size = 22, hjust = c(0.3, -0.5))
  
  multipanel_plustitle <- plot_grid(region_title_list[[REGION]], multipanel, 
                                    nrow = 2, rel_heights = c(0.1, 0.9))
  
  return(multipanel_plustitle)
}, vis_list = fsc_vis_components_list, region_title_list = region_title)


### CREATING THE LEGEND ###
tree_df <- data.frame(x1 = c(2, 2, 2),
                      x2 = c(2, 2.5, 1.5),
                      y1 = c(1, 2, 2),
                      y2 = c(2, 4, 4))

find_xval <- function(yval, branch = c('left', 'right')) {
  switch(branch,
         left = {(yval - 10)/-4},
         right = {(yval + 6)/4})
}

arrow_yval <- function(ystart, yend, pos = c(0.25, 0.75)) {
  yrange <- abs(yend - ystart)
  bottom_val <- min(c(ystart, yend))
  return(bottom_val + (pos*yrange) )
}


#264653 --> dark blue
mod_type_cols <- setNames(c('#cb5cbb', '#FFC60A', '#2A9D8F', '#E24D28'), 
                          nm = c('Continuous gene flow', 'Early gene flow', 'Recent gene flow', 'No gene flow'))


#CONTINUOUS GENE FLOW
contin_gf_df <- data.frame(
  y = arrow_yval(ystart = 2, yend = 4, pos = c(0.4, 0.75)),
  x_start = find_xval(yval = arrow_yval(ystart = 2, yend = 4, pos = c(0.4, 0.75)), branch = 'left'),
  x_end = find_xval(yval = arrow_yval(ystart = 2, yend = 4, pos = c(0.4, 0.75)), branch = 'right')
)
contin_gf_df[2,] <- c(contin_gf_df[2,][1], contin_gf_df[2,][3], contin_gf_df[2,][2])


continuous_shaded_df <- data.frame(x = c(1.5, 2, 2.5),
                                   y = c(4, 2, 4))

contin_gf_legend_tree <- ggplot(data = tree_df) +
  geom_polygon(data = continuous_shaded_df, aes(x = x, y = y), fill = "#bfbfbf") +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               size = 4, lineend = "round", col = mod_type_cols[['Continuous gene flow']]) +
  geom_segment(data = contin_gf_df,
               aes(x = x_start, y = y, xend = x_end, yend = y),
               arrow = arrow(length = unit(0.23, "cm"), type = "closed"), 
               size = 1.5, col = "black", lineend = "round") +
  xlim(1, 3) +
  theme_void() +
  theme(plot.margin = margin(15, 0, 15, 0)) +
  coord_cartesian(clip = 'off')


#RECENT GENE FLOW
recent_gf_df <- data.frame(
  y = arrow_yval(ystart = 3, yend = 4, pos = c(0.25, 0.75)),
  x_start = find_xval(yval = arrow_yval(ystart = 3, yend = 4, pos = c(0.25, 0.75)), branch = 'left'),
  x_end = find_xval(yval = arrow_yval(ystart = 3, yend = 4, pos = c(0.25, 0.75)), branch = 'right')
)
recent_gf_df[2,] <- c(recent_gf_df[2,][1], recent_gf_df[2,][3], recent_gf_df[2,][2])


recent_shaded_df <- data.frame(x = c(1.5, find_xval(yval = 3, branch = 'left'), find_xval(yval = 3, branch = 'right'), 2.5),
                               y = c(4, 3, 3, 4))

recent_gf_legend_tree <- ggplot(data = tree_df) +
  geom_polygon(data = recent_shaded_df, aes(x = x, y = y), fill = "#bfbfbf") +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               size = 4, lineend = "round", col = mod_type_cols[['Recent gene flow']]) +
  geom_segment(data = recent_gf_df,
               aes(x = x_start, y = y, xend = x_end, yend = y),
               arrow = arrow(length = unit(0.23, "cm"), type = "closed"), 
               size = 1.5, col = "black", lineend = "round") +
  xlim(1, 3) +
  theme_void() +
  theme(plot.margin = margin(15, 0, 15, 0)) +
  coord_cartesian(clip = 'off')


#EARLY GENE FLOW
early_gf_df <- data.frame(
  y = arrow_yval(ystart = 2, yend = 3, pos = c(0.4, 0.75)),
  x_start = find_xval(yval = arrow_yval(ystart = 2, yend = 3, pos = c(0.5, 0.8)), branch = 'left'),
  x_end = find_xval(yval = arrow_yval(ystart = 2, yend = 3, pos = c(0.5, 0.8)), branch = 'right')
)
early_gf_df[2,] <- c(early_gf_df[2,][1], early_gf_df[2,][3], early_gf_df[2,][2])


early_shaded_df <- data.frame(x = c(find_xval(yval = 3, branch = 'left'), 2, find_xval(yval = 3, branch = 'right')),
                              y = c(3, 2, 3))

early_gf_legend_tree <- ggplot(data = tree_df) +
  geom_polygon(data = early_shaded_df, aes(x = x, y = y), fill = "#bfbfbf") +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               size = 4, lineend = "round", col = mod_type_cols[['Early gene flow']]) +
  geom_segment(data = early_gf_df,
               aes(x = x_start, y = y, xend = x_end, yend = y),
               arrow = arrow(length = unit(0.23, "cm"), type = "closed"), 
               size = 1.5, col = "black", lineend = "round") +
  xlim(1, 3) +
  theme_void() +
  theme(plot.margin = margin(15, 0, 15, 0)) +
  coord_cartesian(clip = 'off')


#NO GENE FLOW

no_gf_legend_tree <- ggplot(data = tree_df) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               size = 4.5, lineend = "round", col = mod_type_cols[['No gene flow']]) +
  xlim(1, 3) +
  theme_void() +
  theme(plot.margin = margin(15, 0, 15, 0)) +
  coord_cartesian(clip = 'off')

#Create full legend
tree_legend <- cowplot::plot_grid(contin_gf_legend_tree,
                                  early_gf_legend_tree,
                                  recent_gf_legend_tree,
                                  no_gf_legend_tree, 
                                  nrow = 4) +
  theme(plot.margin = margin(47, 0, 35, 0 ))


### Putting together the full plot ###
north_fsc_results <- cowplot::plot_grid(fsc_vis_components_list$north$par + theme(plot.margin = margin(5, 55, 5, 5)), 
                                        ggdraw(egg::ggarrange(fsc_vis_components_list$north$aic, fsc_vis_components_list$north$likelihood, ncol = 2)) + theme(plot.margin = margin(5, 5, 5, 35)),
                                        labels = c('(a)', '(b)'), rel_widths = c(0.35, 0.65), label_size = 18)

north_fsc_results_plustitle <- plot_grid(region_title$North, north_fsc_results, nrow = 2, rel_heights = c(0.1, 0.9))

multipanel1 <- cowplot::plot_grid(multipanel_list$north + theme(plot.margin = margin(0, 0, 6, 0)), 
                                  multipanel_list$mid + theme(plot.margin = margin(6, 0, 0, 0)), 
                                  nrow = 2)

cowplot::plot_grid(multipanel1, tree_legend, ncol = 2, rel_widths = c(0.92, 0.08)) + theme(plot.margin = margin(0, 0, 0, 12))

# ggsave(here('figures', 'fastsimcoal_mod_results_pluslegend_9_30_2021.png'), 
#        width = 2.25*18, height = 2.25*10.5, units = "cm", bg = "white")

ggsave(here('figures', 'fastsimcoal_mod_results_pluslegend_3_31_2022.png'), 
       width = 2.25*18, height = 2.25*10.5, units = "cm", bg = "white")



#################################
### TABLE OF MODEL PARAMETERS ###
#################################

#STRUCTURE OF MIGRATION MATRIX
# 0 MIG_01
# MIG_10 0

#The above-mentioned matrix states that, for each generation backward in 
#time, any gene from population 0 has probability MIG_01 to be sent to 
#population 1, and that a gene from population 1 has a probability MIG_10
#to move to population 0.

processed_param_list <- lapply(names(best_mod_info), function(REGION, param_list) {
  
  processed_df <- t(apply(param_list[[REGION]], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975) ))) %>%
    as.data.frame() %>%
    rename("Lower CI" = `2.5%`,
           Median = `50%`,
           "Upper CI" = `97.5%`) %>% 
    rownames_to_column(var = "Parameter") %>%
    filter(!(Parameter %in% c("MaxEstLhood", "MaxObsLhood"))) %>%
    mutate(Description = case_when(Parameter == "N_POP1" ~ "kazumbe pop size",
                                   Parameter == "N_POP2" ~ "polyodon pop size",
                                   Parameter == "N_ANC_POP1" ~ "initial kazumbe pop size",
                                   Parameter == "N_ANC_POP2" ~ "initial polyodon pop size",
                                   Parameter == "MIG01" ~ "mig. rate, kazumbe to polyodon",
                                   Parameter == "MIG10" ~ "mig. rate, polyodon to kazumbe",
                                   Parameter == "TPROP" ~ "prop of TDIV to mig. cessation",
                                   Parameter == "RSANC" ~ "anc. pop size (relative to sink deme)",
                                   Parameter == "TDIV" ~ "time to divergence",
                                   Parameter == "CHANGM" ~ "time to mig. cessation",
                                   Parameter == "GR_POP1" ~ "kazumbe growth rate",
                                   Parameter == "GR_POP2" ~ "polyodon growth rate"),
           Region = switch(REGION, north = "North", mid = "Mid") ) %>% 
    relocate(Region, .before = Parameter)
  
  #update the complex parameters based on the lower CIs, medians, and upper CIs of the simple parameters used to calculate them
  processed_df[processed_df$Parameter == 'GR_POP1', c("Lower CI", "Median", "Upper CI")] <- log(processed_df[processed_df$Parameter == 'N_ANC_POP1', c("Lower CI", "Median", "Upper CI")]/processed_df[processed_df$Parameter == 'N_POP1', c("Lower CI", "Median", "Upper CI")])/processed_df[processed_df$Parameter == 'TDIV', c("Lower CI", "Median", "Upper CI")]
  processed_df[processed_df$Parameter == 'GR_POP2', c("Lower CI", "Median", "Upper CI")] <- log(processed_df[processed_df$Parameter == 'N_ANC_POP2', c("Lower CI", "Median", "Upper CI")]/processed_df[processed_df$Parameter == 'N_POP2', c("Lower CI", "Median", "Upper CI")])/processed_df[processed_df$Parameter == 'TDIV', c("Lower CI", "Median", "Upper CI")]
  processed_df[processed_df$Parameter == 'CHANGM', c("Lower CI", "Median", "Upper CI")] <- processed_df[processed_df$Parameter == 'TDIV', c("Lower CI", "Median", "Upper CI")]*processed_df[processed_df$Parameter == 'TPROP', c("Lower CI", "Median", "Upper CI")]
  
  return(
    list(table_df = processed_df %>% 
           mutate(Median = formatC(Median, format = "e", digits = 2),
                  "Lower CI" = formatC(`Lower CI`, format = "e", digits = 2),
                  "Upper CI" = formatC(`Upper CI`, format = "e", digits = 2)),
         raw_values = processed_df%>% 
           mutate(Median = formatC(Median, format = "f", digits = 7),
                  "Lower CI" = formatC(`Lower CI`, format = "f", digits = 7),
                  "Upper CI" = formatC(`Upper CI`, format = "f", digits = 7)))
    
  )
}, param_list = bootstrap_results)

#used for more detailed reporting of the parameter estimates in the text of the Methods
nonscientific_param_df <- lapply(processed_param_list, function(x) x[['raw_values']]) %>% 
  bind_rows()

#used for reporting parameter estimates in supplementary table
processed_param_df <- lapply(processed_param_list, function(x) x[['table_df']]) %>% 
  bind_rows()


fsc_best_mods_params <- processed_param_df %>% 
  kbl('latex', booktabs = TRUE, align = "c", escape = TRUE) %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  collapse_rows(1, latex_hline = "major")

cat(fsc_best_mods_params, file = here('tables', 'fsc_best_mods_params_table.txt'), append = FALSE)



##################################
### TABLE OF MODEL FIT RESULTS ###
##################################

#the same order as the models in the main text figure
mod_order <- c('contin. gf (sym.)',
               'contin. gf (asym.)',
               'contin. gf (Pk to Pp)',
               'contin. gf (Pp to Pk)',
               'early gf (sym.)',
               'early gf (asym.)',
               'early gf (Pk to Pp)',
               'early gf (Pp to Pk)',
               'recent gf (sym.)',
               'recent gf (asym.)',
               'recent gf (Pk to Pp)',
               'recent gf (Pp to Pk)',
               'no gf')

fsc_mod_fit_processed_final <- aic_list_processed_df %>% 
  select(region, mod_class, mod_name, number_params, MaxEstLhood, deltaL, AIC, deltaAIC, relative_likelihood_updated) %>% 
  mutate(MaxEstLhood = round(MaxEstLhood, 2),
         deltaL = round(deltaL, 2),
         AIC = round(AIC, 2),
         deltaAIC = round(deltaAIC, 2),
         region = recode(region, mid = "Mid", north = "North")) %>% 
  rename(Region = region,
         "Model class" = mod_class,
         "Model name" = mod_name,
         "Number params." = number_params,
         "log10(Lhood)" = MaxEstLhood,
         "$\\Delta$Lhood" = deltaL,
         "$\\Delta$AIC" = deltaAIC,
         "Rel. Lhood" = relative_likelihood_updated) %>% 
  mutate(`Model name` = factor(`Model name`, levels = mod_order)) %>% 
  dplyr::arrange(desc(Region), `Model name`)

rownames(fsc_mod_fit_processed_final) <- NULL

fsc_mod_fit_table <- fsc_mod_fit_processed_final %>% 
  kbl('latex', booktabs = TRUE, align = "c", escape = FALSE) %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  collapse_rows(1, latex_hline = "major")

cat(fsc_mod_fit_table, file = here('tables', 'fsc_mod_fit_table.txt'), append = FALSE)



################################################
### MULTIPANEL PLOT OF ALL TESTED FSC MODELS ###
################################################

library(grid)

#source(here('project_scripts', 'fsc_plotting_functions.R'))

par_file_list_processed <- lapply(par_file_list, function(region) {
  region[!grepl(".*change_geneflow.*", names(region))]
})


### process names ###
#1. remove all underscores
names1 <- gsub(pattern = '_', replacement = ' ', names(par_file_list_processed$mid))

first_word <- unique(word(gsub(pattern = '_', replacement = ' ', names(par_file_list_processed$mid)), start = 1, end = 1))

par_multi_level_list <- lapply(setNames(nm = first_word), function(first_w, name_vec, par_list) {
  processed_names <- word(gsub(pattern = '_', replacement = ' ', name_vec), start = 1, end = 1)
  full_names <- name_vec[processed_names %in% first_w]
  
  par_list[names(par_list) %in% full_names]
}, name_vec = names(par_file_list_processed$mid), par_list = par_file_list_processed$mid)

new_names_vec <- data.frame(orig_name = names(par_multi_level_list)) %>% 
  mutate(new_names = case_when(orig_name == "early" ~ "Early gene flow",
                               orig_name == "continuous" ~ "Continuous gene flow",
                               orig_name == "no" ~ "No gene flow",
                               orig_name == "recent" ~ "Recent gene flow")) %>% 
  pull(new_names)

names(par_multi_level_list) <- new_names_vec

par_full_updated_names <- lapply(par_multi_level_list, function(mod_list) {
  
  names(mod_list)[!grepl(pattern = "unidirect0|symmetric|unidirect1", names(mod_list))] <- "Asymmetric"
  names(mod_list)[grepl(pattern = "symmetric_dualpopgrowth_altspec", names(mod_list))] <- "Symmetric"
  names(mod_list)[grepl(pattern = "unidirect0_dualpopgrowth_altspec", names(mod_list))] <- "Unidirect (Pk to Pp)"
  names(mod_list)[grepl(pattern = "unidirect1_dualpopgrowth_altspec", names(mod_list))] <- "Unidirect (Pp to Pk)"
  return(mod_list)
})

mod_type_cols <- setNames(c('#cb5cbb', '#FFC60A', '#2A9D8F', '#E24D28'), 
                          nm = c('Continuous gene flow', 'Early gene flow', 'Recent gene flow', 'No gene flow'))

mod_type_cols_background <- setNames(c('#efceea', '#fff3ce', '#d4ebe8', '#f6c9be'), 
                                     nm = c('Continuous gene flow', 'Early gene flow', 'Recent gene flow', 'No gene flow'))


plot_list_migmatcolor <- list()

for (CLASS in names(par_full_updated_names)) {
  for (MOD in names(par_full_updated_names[[CLASS]])) {
    
    plot_list_migmatcolor[[CLASS]][[MOD]] <- plot_demographic_mod(par = par_full_updated_names[[CLASS]][[MOD]],
                                                      colors = TRUE,
                                                      pop_colors = c("#DC7633", "#3498DB"),
                                                      scale_pops = FALSE,
                                                      pop_color = '#929292', #a3a3a3
                                                      mig_band_color = mod_type_cols[[CLASS]],
                                                      generic_plotting = TRUE,
                                                      generic_arrow_head_size = 0.09)
    
    #plot_list_migmatcolor[[CLASS]][[MOD]] <- plot_list_migmatcolor[[CLASS]][[MOD]] + 
      #theme(panel.background = element_rect(fill = mod_type_cols_background[[CLASS]], color = mod_type_cols_background[[CLASS]]))
    
    #if ( (CLASS != "No gene flow") & (MOD %in% c("Symmetric", "Unidirect (Pk to Pp)"))) {
    # if (CLASS  %in% c("Continuous gene flow", "Recent gene flow")) {
    #   plot_list[[CLASS]][[MOD]] <- plot_list[[CLASS]][[MOD]] + 
    #     theme(panel.background = element_rect(fill = mod_type_cols_background[[CLASS]], color = mod_type_cols_background[[CLASS]]))
    # }
    
  }
}



plot_list_migmatcolor_withpadding <- lapply(plot_list_migmatcolor, function(x) {
  lapply(x, function(y) y  + theme(plot.margin = unit(c(0.2, 0.7, 0.2, 0.7), "cm")) )
})


prac_multipanel <- plot_grid(plot_list_migmatcolor_withpadding$`Continuous gene flow`$Symmetric, plot_list_migmatcolor_withpadding$`Early gene flow`$Symmetric, plot_list_migmatcolor_withpadding$`Recent gene flow`$Symmetric,
          plot_list_migmatcolor_withpadding$`Continuous gene flow`$Asymmetric, plot_list_migmatcolor_withpadding$`Early gene flow`$Asymmetric, plot_list_migmatcolor_withpadding$`Recent gene flow`$Asymmetric,
          plot_list_migmatcolor_withpadding$`Continuous gene flow`$`Unidirect (Pk to Pp)`, plot_list_migmatcolor_withpadding$`Early gene flow`$`Unidirect (Pk to Pp)`, plot_list_migmatcolor_withpadding$`Recent gene flow`$`Unidirect (Pk to Pp)`,
          plot_list_migmatcolor_withpadding$`Continuous gene flow`$`Unidirect (Pp to Pk)`, plot_list_migmatcolor_withpadding$`Early gene flow`$`Unidirect (Pp to Pk)`, plot_list_migmatcolor_withpadding$`Recent gene flow`$`Unidirect (Pp to Pk)`,
          nrow = 4, labels = paste0('(', letters[1:12], ')') )

gf_timing_legend <- get_legend(
  data.frame(x = c(1, 2, 3), 
             "Gene flow timing" = factor(c('continuous', 'early', 'recent'), levels = c('continuous', 'early', 'recent')),
             check.names = FALSE) %>% 
    ggplot() +
    geom_histogram(aes(x = x, fill = `Gene flow timing`)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.margin = margin(-6, 0, 0, 0),
          legend.title = element_text(size = 16.5),
          legend.text = element_text(size = 13, margin = margin(l = 1, r = 8, unit = "pt"))) +
    scale_fill_manual(values = setNames(c('#cb5cbb', '#FFC60A', '#2A9D8F'), 
                                        nm = c('continuous', 'early', 'recent')))
  
)


arrow_legend <- data.frame(y = c(4, 3.1, 2.9, 2, 1)) %>% 
  ggplot() +
  geom_segment(
    aes(x = 2.5, xend = 3, y = y, yend = y), 
    arrow = arrow(length = unit(0.25, 'cm'), type = 'closed', ends = c("both", "first", "last", "last", "first") ),
    size = 0.7
  ) +
  geom_text(data = data.frame(y_text = 4:1,
                       text = c('symmetric', 'asymmetric', 'unidirect (Pk to Pp)', 'unidirect (Pp to Pk)')),
            aes(x = 3.5, y = y_text, label = text), size = 5) +
  theme_void() +
  ggtitle('Gene flow parameterization') +
  theme(plot.title = element_text(size = 20)) +
  xlim(0, 6)



y.grob_symmetric <- textGrob("Bidirectional\nsymmetric", 
                             gp=gpar(fontface="plain", col="black", fontsize=15), rot=90)

y.grob_asymmetric <- textGrob("Bidirectional\nasymmetric", 
                              gp=gpar(fontface="plain", col="black", fontsize=15), rot=90)

y.grob_uni01 <- textGrob("Unidirectional\n(Pk to Pp)", 
                              gp=gpar(fontface="plain", col="black", fontsize=15), rot=90)

y.grob_uni10 <- textGrob("Unidirectional\n(Pp to Pk)", 
                         gp=gpar(fontface="plain", col="black", fontsize=15), rot=90)


symmetric_plot <- grid.arrange(arrangeGrob(plot_grid(plot_list_migmatcolor_withpadding$`Continuous gene flow`$Symmetric,
                                                     plot_list_migmatcolor_withpadding$`Early gene flow`$Symmetric,
                                                     plot_list_migmatcolor_withpadding$`Recent gene flow`$Symmetric, 
                                                     nrow = 1,
                                                     labels = c("(a)", "(b)", "(c)"), label_size = 18), 
                                           left = y.grob_symmetric))

asymmetric_plot <- grid.arrange(arrangeGrob(plot_grid(plot_list_migmatcolor_withpadding$`Continuous gene flow`$Asymmetric,
                                                     plot_list_migmatcolor_withpadding$`Early gene flow`$Asymmetric,
                                                     plot_list_migmatcolor_withpadding$`Recent gene flow`$Asymmetric, 
                                                     nrow = 1,
                                                     labels = c("(d)", "(e)", "(f)"), label_size = 18), 
                                           left = y.grob_asymmetric))

unidirect_01 <- grid.arrange(arrangeGrob(plot_grid(plot_list_migmatcolor_withpadding$`Continuous gene flow`$`Unidirect (Pk to Pp)`,
                                                      plot_list_migmatcolor_withpadding$`Early gene flow`$`Unidirect (Pk to Pp)`,
                                                      plot_list_migmatcolor_withpadding$`Recent gene flow`$`Unidirect (Pk to Pp)`, 
                                                      nrow = 1,
                                                      labels = c("(g)", "(h)", "(i)"), label_size = 18), 
                                            left = y.grob_uni01))

unidirect_10 <- grid.arrange(arrangeGrob(plot_grid(plot_list_migmatcolor_withpadding$`Continuous gene flow`$`Unidirect (Pp to Pk)`,
                                                   plot_list_migmatcolor_withpadding$`Early gene flow`$`Unidirect (Pp to Pk)`,
                                                   plot_list_migmatcolor_withpadding$`Recent gene flow`$`Unidirect (Pp to Pk)`, 
                                                   nrow = 1,
                                                   labels = c("(j)", "(k)", "(l)"), label_size = 18), 
                                         left = y.grob_uni10))

gf_timing_legend_multipanel <- plot_grid(plot_grid(ggplot() + theme_void() + theme(panel.background = element_rect(fill = 'white', color = "white"))),
                      gf_timing_legend, 
                      ncol = 2, 
                      rel_widths = c(0.07, 0.93))

plot_grid(symmetric_plot,
          asymmetric_plot,
          unidirect_01,
          unidirect_10,
          gf_timing_legend_multipanel, nrow = 5, rel_heights = c(100, 100, 100, 100, 20)) +
  theme(plot.margin = unit(c(0, 0, 0, 0.23), "cm"))

ggsave(here('figures', 'fastsimcoal_allmods_testcolor.png'), 
       width = 1.2*20, height = 1.2*16.5, units = "cm", bg = "white")



####################################
### PARAMETER SEARCH RANGE TABLE ###
####################################

processed_param_df_full <- do.call(rbind, lapply(tpl_est_file_list, process_params))
rownames(processed_param_df_full) <- NULL

processed_param_info_final <- processed_param_df_full[!duplicated(processed_param_df_full$Parameter),] %>%
  filter(!(Parameter %in% c('MIG_01', 'MIG_10') )) %>% 
  mutate(Lower = recode(Lower, '0.000001' = '1e-06', '0.0001' = '1e-04'),
         Upper = recode(Upper, '250000' = '2.5e05'),
         Description = case_when(Parameter == "N_POP1" ~ "kazumbe pop size",
                                 Parameter == "N_POP2" ~ "polyodon pop size",
                                 Parameter == "N_ANC_POP1" ~ "initial kazumbe pop size",
                                 Parameter == "N_ANC_POP2" ~ "initial polyodon pop size",
                                 Parameter == "MIG01" ~ "mig. rate, kazumbe to polyodon",
                                 Parameter == "MIG10" ~ "mig. rate, polyodon to kazumbe",
                                 Parameter == "MIG" ~ "mig. rate (symmetric mig. models)",
                                 Parameter == "TPROP" ~ "prop of TDIV to mig. cessation",
                                 Parameter == "RSANC" ~ "anc. pop size (relative to sink deme)",
                                 Parameter == "TDIV" ~ "time to divergence",
                                 Parameter == "CHANGM" ~ "time to mig. cessation",
                                 Parameter == "GR_POP1" ~ "kazumbe growth rate",
                                 Parameter == "GR_POP2" ~ "polyodon growth rate"))

processed_param_info_final_table <- processed_param_info_final %>% 
  kbl('latex',  booktabs = TRUE, align = "c") %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" ", " ", " ", "Search Range" = 2, " ", " "))

cat(processed_param_info_final_table, file = here('tables', 'fsc_param_info.txt'), append = FALSE)



#####################################################
#####################################################
### PART 2: VISUALIZING AND SUMMARIZING MODEL FIT ###
#####################################################
#####################################################


#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(cowplot)
library(tidyverse)
library(reshape2)
library(egg)
library(here)
library(scales)

#load visualization functions
source(here('project_scripts', 'sfs_compare_vis_funcs.R'))

### NOTES ###
#vertical axis of SFS ends up being the horizontal axis in visualized sfs using plot_2d_sfs
#the rows output of plot_1d_sfs is the 

# Where the first pop specified is listed in the rows and the second pop
## specified is listed in the columns.


### Loading results (and some initial processing) ###
load_sfs_list <- lapply(setNames(nm = c('north', 'mid')), function(REGION) {
  exp_sfs_list_names <- list.files(here('demographic_modelling', 'fastsimcoal', 'results_files', 'model_fit_results', REGION, 'likelihood_dir', 'popgrow_alt_spec', 'recent_geneflow_dualpopgrowth_altspec_updated_sfs'))
  exp_sfs_list <- lapply(setNames(paste0(here('demographic_modelling', 'fastsimcoal', 'results_files', 'model_fit_results', REGION, 'likelihood_dir', 'popgrow_alt_spec', 'recent_geneflow_dualpopgrowth_altspec_updated_sfs'), '/',exp_sfs_list_names), nm = gsub(".txt", "", exp_sfs_list_names)), function(FILE) {
    as.matrix(read.table(FILE, header = T, row.names = 1))
  })
  
  obs_sfs <- as.matrix(read.table(here('demographic_modelling', 'fastsimcoal', 'sfs_files', paste0(REGION, '_allsites_polykaz_minq20_minDP5_maxDP75_miss0_jointMAFpop1_0.obs')), 
                                  skip = 1, header = TRUE, row.names = 1))
  best_exp <- exp_sfs_list[[which(grepl("toplikelihood$", names(exp_sfs_list)))]]
  
  return(list(obs = obs_sfs,
              exp_list = exp_sfs_list,
              best_exp = best_exp))
})


plot_sfs_fits <- lapply(setNames(nm = c('north', 'mid')), function(REGION, sfs_list) {
  
  plot_list <- list()
  
  ##############
  ### 2D VIZ ###
  ##############
  joint_sfs <- plot_2d_sfs(sfs_list = list(observed = sfs_list[[REGION]]$obs, expected = sfs_list[[REGION]]$best_exp),
                           plot_obs = TRUE,
                           plot_exp = TRUE,
                           plot_comparison = TRUE,
                           log10 = TRUE,
                           drop_monomorph = TRUE,
                           minimum_count = 10,
                           divergent_col_scheme = c('#4575b4', '#ffffbf', '#d73027'),
                           sequential_col_scheme = c('white', '#d73027'),
                           scale_resolution = c(obs_exp = 1, comparison = 0.5))
  
  
  plot_list[['twod_sfs']] <- plot_grid(joint_sfs$plots$observed +
                                         ggtitle('Observed SFS') +
                                         labs(fill='Log10(SFS count)') +
                                         theme(plot.title = element_text(size = 15, hjust = 0.5),
                                               legend.title = element_text(size = 8)) +
                                         ylab(expression(paste(italic("P"), ". sp. 'kazumbe'", sep = ""))) +
                                         xlab(expression(paste(italic("P"), ". cf. ",  italic("polyodon"), sep = ""))),
                                       joint_sfs$plots$expected +
                                         labs(fill='Log10(SFS count)') +
                                         ggtitle('Expected SFS') +
                                         theme(plot.title = element_text(size = 15, hjust = 0.5),
                                               legend.title = element_text(size = 8)) +
                                         ylab(expression(paste(italic("P"), ". sp. 'kazumbe'", sep = ""))) +
                                         xlab(expression(paste(italic("P"), ". cf. ",  italic("polyodon"), sep = ""))),
                                       joint_sfs$plots$comparison+
                                         labs(fill='Exp. SFS/Obs. SFS') +
                                         ggtitle('Exp. SFS/Obs. SFS') +
                                         theme(plot.title = element_text(size = 15, hjust = 0.5),
                                               legend.title = element_text(size = 8)) +
                                         ylab(expression(paste(italic("P"), ". sp. 'kazumbe'", sep = ""))) +
                                         xlab(expression(paste(italic("P"), ". cf. ",  italic("polyodon"), sep = ""))),
                                       ncol = 3
  )
  
  
  ###################
  ### MARG 1D VIZ ###
  ###################
  
  one_d_sfs <- plot_1d_sfs(sfs_list = list(observed = sfs_list[[REGION]]$obs, expected = sfs_list[[REGION]]$exp_list),
                           max_freq = ncol(sfs_list[[REGION]]$obs) - 1,
                           plot = TRUE,
                           log10 = FALSE)
  
  kazumbe_1d <- one_d_sfs$plots$column + 
    theme_minimal() +
    #scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values=c("#118ab2", "#ef476f")) +
    theme(panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 1),
          legend.position = "none",
          plot.title = element_text(size = 18, hjust = 0.5),
          axis.text = element_text(size = switch(REGION, north = 4, mid = 11))) +
    ylab("Number of sites") +
    xlab("Minor allele frequency") +
    ggtitle(expression(paste(italic("Petrochromis"), " sp. 'kazumbe'", sep = "")))
  
  polyodon_1d <- one_d_sfs$plots$row + 
    theme_minimal() +
    #scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values=c("#118ab2", "#ef476f")) +
    theme(panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 1),
          legend.position = "none",
          plot.title = element_text(size = 18, hjust = 0.5),
          axis.text = element_text(size = switch(REGION, north = 4, mid = 11))) +
    ylab("Number of sites") +
    xlab("Minor allele frequency") +
    ggtitle(expression(paste(italic("Petrochromis"), " cf. ",  italic("polyodon"), sep = "")))
  
  legend_1d_sfs <- get_legend(one_d_sfs$plots$row + 
                                theme_minimal() +
                                #scale_x_continuous(expand = c(0, 0)) +
                                scale_y_continuous(expand = c(0, 0)) +
                                scale_fill_manual(values=c("#118ab2", "#ef476f")) +
                                theme(panel.grid.minor = element_blank(),
                                      axis.line = element_line(color = "black", size = 1),
                                      legend.title = element_blank(),
                                      plot.title = element_text(size = 20, hjust = 0.5)) +
                                ylab("Number of sites") +
                                xlab("Minor allele frequency") +
                                ggtitle("polyodon"))
  
  plot_list[['oned_sfs']] <- plot_grid(kazumbe_1d, polyodon_1d, legend_1d_sfs, 
                                       ncol = 3, rel_widths = c(0.45, 0.45, 0.12))
  
  return(plot_list)
  
}, sfs_list = load_sfs_list)


plot_sfs_fits$north$twod_sfs
ggsave(here('figures', 'fsc_north_region_topmod_2dsfs.png'), 
       width = 1.5*18, height = 1.5*5.5, units = "cm", bg = "white")

plot_sfs_fits$north$oned_sfs
ggsave(here('figures', 'fsc_north_region_topmod_1dsfs.png'), 
       width = 1.5*16, height = 1.5*7, units = "cm", bg = "white")

plot_sfs_fits$mid$twod_sfs
ggsave(here('figures', 'fsc_mid_region_topmod_2dsfs.png'), 
       width = 1.5*18, height = 1.5*5.5, units = "cm", bg = "white")

plot_sfs_fits$mid$oned_sfs
ggsave(here('figures', 'fsc_mid_region_topmod_1dsfs.png'), 
       width = 1.5*16, height = 1.5*7, units = "cm", bg = "white")



#################################
### CODE CURRENTLY NOT IN USE ###
#################################

# generic_par <- function(par, generic_element = c('pop_size', 'mig_mat_list', 'hist_event_info')) {
#   
#   for (element in names(par)) {
#     
#     if (element == 'pop_size' && element %in% generic_element) {
#       par[[element]] <- rep(10, length(par[[element]]))
#     }
#     
#     if (element == 'mig_mat_list' && element %in% generic_element) {
#       par[[element]] <- lapply(par[[element]], function(y) {
#         y[y != 0] <- 10
#         return(y)
#       })
#     }
#     
#     if (element == 'hist_event_info' && element %in% generic_element) {
#       par[[element]]$processed_hist_events$time <- seq(10, length(par[[element]]$processed_hist_events$time)*10, by = 10)
#     }
#     
#   }
#   
#   return(par)
# }
# 
# 
# 
# time_interval <- function(mig_mat_list, historical_event, diverge_cushion) {
#   time_df_list <- list(data.frame(mig_mat = 0, start = 0, end = NA))
#   
#   new_migmat_pos <- which(!duplicated(historical_event$mig_mat_index) & historical_event$mig_mat_index != 0)
#   
#   if (length(new_migmat_pos) == 0) {
#     time_df_list[[1]][1, 'end']<- historical_event[historical_event$source != historical_event$sink,]$time - diverge_cushion
#   } else {
#     subset_historical_events <- historical_event[new_migmat_pos,]
#     for (i in seq_len(nrow(subset_historical_events))) {
#       if (subset_historical_events$source[i] == subset_historical_events$sink[i]) {
#         time_df_list[[i]]$end <- subset_historical_events$time[i]
#         time_df_list <- c(time_df_list, list(data.frame(mig_mat = i, start = subset_historical_events$time[i], end = NA)))
#       } else {
#         time_df_list[[i]]$end <- subset_historical_events$time[i] - diverge_cushion
#         time_df_list <- c(time_df_list, list(data.frame(mig_mat = i, start = subset_historical_events$time[i] - diverge_cushion, end = NA)))
#       }
#     }
#     
#   }
#   
#   return(time_df_list)
# }
# 
# migmat_arrow_info <- function(mig_mat,
#                               start_y, 
#                               end_y,
#                               offset = 0,
#                               arrow_replicates = 2,
#                               max_val = NULL,
#                               max_arrow_size = 2,
#                               min_arrow_size = 0.2,
#                               max_arrow_head_size = 0.15,
#                               min_arrow_head_size = 0.04,
#                               pop0_x,
#                               pop1_x,
#                               direction = c('forward', 'backward')) {
#   
#   if (all(mig_mat == 0) | any(mig_mat == 'NONE') ) {
#     return(data.frame(from = 0, to = 0, mig_est = 0, arrow_size = 0, arrow_head_size = 0, arrow_pos = 0, arrow_start = 0, arrow_end = 0))
#   }
#   
#   ### 1. argument checks ###
#   if (start_y >= end_y)
#     stop('start_y is more recent than end_y (migmat_arrow_info has a present --> past perspective) so start_y needs to be smaller than end_y')
#   direction <- match.arg(direction, several.ok = FALSE)
#   
#   
#   ### 2. process migration matrix ###
#   dimnames(mig_mat) <- replicate(2, 0:(ncol(mig_mat) - 1), simplify = FALSE) #add row and column names 
#   
#   #transform matrix to a "long" dataframe where each row represents an entry in the matrix (and remove zero entries)
#   melted_mig_mat <- melt(mig_mat) %>% 
#     filter(value > 0)
#   #Var1 --> from; Var2 --> to; value --> migration param
#   colnames(melted_mig_mat) <- c('from', 'to', 'mig_est')
#   
#   mig_params <- unique(melted_mig_mat$mig_est[melted_mig_mat$mig_est > 0]) #vector of unique non-zero migration parameters
#   number_unique_arrows <- length(mig_params) #how many unique migration parameters?
#   
#   
#   ### 3. arrow size specifications ###
#   if (is.null(max_val)) max_val <- max(mig_params) #set max value if one isn't provided
#   
#   #arrow size
#   melted_mig_mat$arrow_size <- (melted_mig_mat$mig_est/max_val)*max_arrow_size #standardize 
#   melted_mig_mat$arrow_size[melted_mig_mat$arrow_size < min_arrow_size] <- min_arrow_size
#   
#   #arrow head size
#   arrow_head_size_dif <- max_arrow_head_size - min_arrow_head_size
#   mig_mat_prop <- arrow_head_size_dif*(melted_mig_mat$mig_est/max_val)
#   melted_mig_mat$arrow_head_size <- min_arrow_head_size + mig_mat_prop
#   
#   
#   ### 4. arrow positions ###
#   start_end_dif <- abs(end_y - start_y)
#   updated_start <- start_y + (start_end_dif*offset)
#   updated_end <- end_y - (start_end_dif*offset)
#   updated_dif <- abs(updated_end - updated_start)
#   interval_vec <- seq(updated_start, updated_end, by = abs(updated_start - updated_end)/((arrow_replicates*nrow(melted_mig_mat)) + 1))
#   melted_mig_mat_updated <- do.call(rbind, replicate(arrow_replicates, melted_mig_mat, simplify = FALSE)) %>%
#     #mutate(arrow_pos = seq(updated_start, updated_end, length.out = nrow(melted_mig_mat)*arrow_replicates))
#     mutate(arrow_pos = interval_vec[-c(1, length(interval_vec))])
#   
#   
#   ### 5. arrow directionality -- forward or backward in time representation ###
#   if (direction == 'forward') {
#     melted_mig_mat_updated <- melted_mig_mat_updated %>% 
#       mutate(arrow_end = ifelse(from == 0, pop0_x, pop1_x),
#              arrow_start = ifelse(from == 0, pop1_x, pop0_x))
#   } else {
#     melted_mig_mat_updated <- melted_mig_mat_updated %>% 
#       mutate(arrow_start = ifelse(from == 0, pop0_x, pop1_x),
#              arrow_end = ifelse(from == 0, pop1_x, pop0_x))
#   }
#   
#   
#   ### 6. symmetric migration -- add arrow head on other side of arrow with symmetric migration ###
#   if ( (length(unique(melted_mig_mat_updated$mig_est)) == 1) && (nrow(melted_mig_mat) > 1) )  {
#     melted_mig_mat_opposite <- melted_mig_mat_updated %>% 
#       mutate(arrow_end1 = arrow_start,
#              arrow_start1 = arrow_end) %>% 
#       select(from, to, mig_est, arrow_size, arrow_head_size, arrow_pos, arrow_end = arrow_end1, arrow_start = arrow_start1)
#     
#     melted_mig_mat_updated <- rbind(melted_mig_mat_updated, melted_mig_mat_opposite)
#   }
#   
#   return(melted_mig_mat_updated)
# }
# 
# rescale <- function(x,min,max) {
#   (((x - min(data$pop_size)) * (max - min)) / diff(range(data$pop_size))) + min
# }
# 
# rescale_alt <- function(x, min, max, pop_size_vec) {
#   (((x - min(c(pop_size_vec, x) )) * (max - min)) / diff(range(c(pop_size_vec, x)))) + min
# }
# 
# theme_demovis <- function() {
#   theme(panel.background = element_rect(color = 'white', fill = 'white'),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.line.y = element_line(color = 'black', size = 1.25, lineend = 'square'),
#         axis.ticks.x = element_blank(),
#         axis.ticks.length.y = unit(2.5, "mm"),
#         axis.ticks.y = element_line(size = 1.25),
#         axis.text.y = element_text(color = 'black', size = 11),
#         legend.position = 'none')
# }
# 
# 
# 
# 
# plot_demographic_mod <- function(par,
#                                  colors = TRUE,
#                                  pop_colors = NULL,
#                                  scale_pops = TRUE,
#                                  pop_color = 'gray',
#                                  generic_plotting = FALSE,
#                                  max_arrow_head_size = 0.15,
#                                  min_arrow_head_size = 0.05,
#                                  generic_arrow_head_size = 0.15) {
#   
#   scale_pops <- scale_pops & !generic_plotting
#   
#   if (generic_plotting) {
#     data <- generic_par(par, generic_element = c('pop_size', 'mig_mat_list', 'hist_event_info')[c(1, 3)])
#   } else {
#     data <- par
#   }
#   
#   events <- data$hist_event_info$processed_hist_events
#   
#   if (!is.null(pop_colors) & colors & length(pop_colors) == data$number_samples) {
#     pop.cols <- pop_colors
#   } else {
#     pop.cols <- viridis::viridis(data$number_samples)
#   }
#   
#   ############################
#   # subsetting data into migration vs div events
#   ############################
#   mig_events <- events %>%
#     filter(migrants < 1) %>%
#     filter(source != sink)
#   div_events <- events %>%
#     filter(migrants == 1) %>%
#     filter(source != sink)
#   
#   # if (scale_pops) {
#   #   pops_rel_ne <- rescale(data$pop_size, 0.5, 1.2)
#   #   #root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[1], 0.5, 1.2)
#   #   root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1], 0.5, 1.2)
#   #   div_events$source_rel_ne <- pops_rel_ne[div_events$source + 1]
#   #   div_events$sink_rel_ne <- pops_rel_ne[div_events$sink + 1]
#   #   
#   #   mig_events$source_rel_ne <- pops_rel_ne[mig_events$source + 1]
#   #   mig_events$sink_rel_ne <- pops_rel_ne[mig_events$sink + 1]
#   # }
#   
#   if (scale_pops == TRUE) {
#     pops_rel_ne <- rescale_alt(data$pop_size, 0.5, 1.2, pop_size_vec = c(data$pop_size, div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1]) )
#     #root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[1], 0.5, 1.2)
#     root_rel_ne <- rescale_alt( div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1], 0.5, 1.2, pop_size_vec = c(data$pop_size, div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1]) )
#     div_events$source_rel_ne <- pops_rel_ne[div_events$source + 1]
#     div_events$sink_rel_ne <- pops_rel_ne[div_events$sink + 1]
#     
#     mig_events$source_rel_ne <- pops_rel_ne[mig_events$source + 1]
#     mig_events$sink_rel_ne <- pops_rel_ne[mig_events$sink + 1]
#   }
#   
#   
#   
#   ############################
#   # first, set up the skeleton using divergences
#   # it'd be good if we can make this more flexible (i.e. works with any number of div events)
#   ############################
#   tgap <- max(events$time)*0.10
#   troot <- max(events$time)*1.33
#   
#   p <- ggplot() +
#     xlim(0,length(data$pop_size)*1.5)+
#     ylim(0,troot)
#   
#   width <- 0.75
#   if (scale_pops) {
#     if (nrow(div_events) == 1){
#       p1 <- p + 
#         geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width*pops_rel_ne[1]/2,
#                               xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width*pops_rel_ne[2]/2,
#                               ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) + 
#         # root
#         geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width*root_rel_ne/2,
#                               xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width*root_rel_ne/2,
#                               ymin=div_events$time[1],
#                               ymax=troot), fill=pop_color) + 
#         # pop0
#         geom_rect(mapping=aes(xmin=div_events$source[1] + 1 - width*pops_rel_ne[1]/2,
#                               xmax=div_events$source[1] + 1 + width*pops_rel_ne[1]/2,
#                               ymin=0,
#                               ymax=div_events$time[1]), fill=pop_color) +
#         # pop1
#         geom_rect(mapping=aes(xmin=div_events$sink[1] + 1 - width*pops_rel_ne[2]/2,
#                               xmax=div_events$sink[1] + 1 + width*pops_rel_ne[2]/2,
#                               ymin=0,
#                               ymax=div_events$time[1]), fill=pop_color)
#       #print(p1)
#       
#     } else if (nrow(div_events) == 2) { 
#       p1 <- p + 
#         geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width*pops_rel_ne[1]/2,
#                               xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width*pops_rel_ne[2]/2,
#                               ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) +
#         # root
#         geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width*root_rel_ne/2,
#                               xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width*root_rel_ne/2,
#                               ymin=div_events$time[1],
#                               ymax=troot), fill=pop_color) + 
#         # pop0
#         geom_rect(mapping=aes(xmin=div_events$source[1] + 1 - width*pops_rel_ne[1]/2,
#                               xmax=div_events$source[1] + 1 + width*pops_rel_ne[1]/2,
#                               ymin=0,
#                               ymax=div_events$time[1]), fill=pop_color) +
#         # pop1
#         geom_rect(mapping=aes(xmin=div_events$sink[1] + 1 - width*pops_rel_ne[2]/2,
#                               xmax=div_events$sink[1] + 1 + width*pops_rel_ne[2]/2,
#                               ymin=0,
#                               ymax=div_events$time[1]), fill=pop_color) +
#         
#         geom_rect(mapping=aes(xmin=max(div_events$source[2],div_events$sink[2]) + 1 - width*pops_rel_ne[3]/2,
#                               xmax=max(div_events$source[2],div_events$sink[2]) + 1 + width*pops_rel_ne[3]/2,
#                               ymin=0,
#                               ymax=div_events$time[2]), fill=pop_color) +
#         
#         geom_rect(mapping=aes(xmin=min(div_events$source[2],div_events$sink[2]) + 1 - width*pops_rel_ne[2]/2,
#                               xmax=max(div_events$source[2],div_events$sink[2]) + 1 + width*pops_rel_ne[3]/2,
#                               ymin=div_events$time[2]-tgap,
#                               ymax=div_events$time[2]), fill=pop_color)
#       #print(p1)
#     }
#   } else {
#     if (nrow(div_events) == 1){
#       p1 <- p + 
#         geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width/2,
#                               xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width/2,
#                               ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) + 
#         # root
#         geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width/2,
#                               xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width/2,
#                               ymin=div_events$time[1],
#                               ymax=troot), fill=pop_color) + 
#         # pop0
#         geom_rect(mapping=aes(xmin=div_events$source[1] + 1 - width/2,
#                               xmax=div_events$source[1] + 1 + width/2,
#                               ymin=0,
#                               ymax=div_events$time[1]), fill=pop_color) +
#         # pop1
#         geom_rect(mapping=aes(xmin=div_events$sink[1] + 1 - width/2,
#                               xmax=div_events$sink[1] + 1 + width/2,
#                               ymin=0,
#                               ymax=div_events$time[1]), fill=pop_color)
#       #print(p1)
#       
#     } else if (nrow(div_events) == 2) { 
#       p1 <- p + 
#         geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width/2,
#                               xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width/2,
#                               ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) +
#         # root
#         geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width/2,
#                               xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width/2,
#                               ymin=div_events$time[1],
#                               ymax=troot), fill=pop_color) + 
#         # pop0
#         geom_rect(mapping=aes(xmin=div_events$source[1] + 1 - width/2,
#                               xmax=div_events$source[1] + 1 + width/2,
#                               ymin=0,
#                               ymax=div_events$time[1]), fill=pop_color) +
#         # pop1
#         geom_rect(mapping=aes(xmin=div_events$sink[1] + 1 - width/2,
#                               xmax=div_events$sink[1] + 1 + width/2,
#                               ymin=0,
#                               ymax=div_events$time[1]), fill=pop_color) +
#         
#         geom_rect(mapping=aes(xmin=max(div_events$source[2],div_events$sink[2]) + 1 - width/2,
#                               xmax=max(div_events$source[2],div_events$sink[2]) + 1 + width/2,
#                               ymin=0,
#                               ymax=div_events$time[2]), fill=pop_color) +
#         
#         geom_rect(mapping=aes(xmin=min(div_events$source[2],div_events$sink[2]) + 1 - width/2,
#                               xmax=max(div_events$source[2],div_events$sink[2]) + 1 + width/2,
#                               ymin=div_events$time[2]-tgap,
#                               ymax=div_events$time[2]), fill=pop_color)
#       #print(p1)
#     }
#   }
#   
#   
#   if (scale_pops) {
#     pop0_x_left <- div_events$source[1] + 1 - width*pops_rel_ne[1]/2
#     pop0_x_right <- div_events$source[1] + 1 + width*pops_rel_ne[1]/2
#     pop1_x_left <- div_events$sink[1] + 1 - width*pops_rel_ne[2]/2
#     pop1_x_right <-div_events$sink[1] + 1 + width*pops_rel_ne[2]/2
#     
#   } else {
#     pop0_x_left <- div_events$source[1] + 1 - width/2
#     pop0_x_right <- div_events$source[1] + 1 + width/2
#     pop1_x_left <- div_events$sink[1] + 1 - width/2
#     pop1_x_right <- div_events$sink[1] + 1 + width/2
#   }
#   
#   
#   ############################
#   # add colors for populations, if desired
#   ############################
#   if (colors){
#     if (scale_pops) {
#       p1.2 <- p1 +
#         geom_rect(data=div_events, mapping=aes(xmin=source + 1 - width*source_rel_ne/2,
#                                                xmax=source + 1 + width*source_rel_ne/2,
#                                                ymin=0,
#                                                ymax=tgap/2, 
#                                                fill=factor(source))) +
#         geom_rect(data=div_events, mapping=aes(xmin=sink + 1 - width*sink_rel_ne/2,
#                                                xmax=sink + 1 + width*sink_rel_ne/2,
#                                                ymin=0,
#                                                ymax=tgap/2,
#                                                fill=factor(sink))) +
#         scale_fill_manual(values=pop.cols)
#     } else {
#       p1.2 <- p1 +
#         geom_rect(data=div_events, mapping=aes(xmin=source + 1 - width/2,
#                                                xmax=source + 1 + width/2,
#                                                ymin=0,
#                                                ymax=tgap/2, 
#                                                fill=factor(source))) +
#         geom_rect(data=div_events, mapping=aes(xmin=sink + 1 - width/2,
#                                                xmax=sink + 1 + width/2,
#                                                ymin=0,
#                                                ymax=tgap/2,
#                                                fill=factor(sink))) +
#         scale_fill_manual(values=pop.cols)
#     }
#     
#   } else {
#     p1.2 <- p1
#   }
#   
#   
#   ############################  
#   # then, add migration events
#   # needs to be extended for > 2 pops
#   ############################
#   if (nrow(mig_events) > 0) {
#     if (length(pops) == 2) {
#       p2 <- p1.2 +
#         geom_segment(aes(x = if_else(mig_events$source < mig_events$sink,
#                                      (mig_events$source*2+2)/2+0.25,
#                                      (mig_events$source*2+2)/2-0.25), 
#                          xend = if_else(mig_events$source < mig_events$sink,
#                                         (mig_events$sink*2+2)/2-0.25,
#                                         (mig_events$sink*2+2)/2+0.25),
#                          y = mig_events$time,
#                          yend = mig_events$time),
#                      lineend = "round", # See available arrow types in example above
#                      linejoin = "round",
#                      size = 1.5, 
#                      arrow = arrow(length = unit(0.1, "inches")),
#                      colour = "black" # Also accepts "red", "blue' etc
#         ) 
#       print(p2)
#     } else if (length(pops) == 3) {
#       p2 <- p1.2 +
#         geom_segment(data=mig_events[mig_events$time > div_events$time[2],], aes(x = if_else(source < sink,
#                                                                                              (source*2+2)/2+0.25,
#                                                                                              (source*2+2)/2-0.25), 
#                                                                                  xend = if_else(source < sink,
#                                                                                                 (sink*2+2)/2-0.25,
#                                                                                                 (sink*2+2)/2+0.25),
#                                                                                  y = time,
#                                                                                  yend = time),
#                      lineend = "round", # See available arrow types in example above
#                      linejoin = "round",
#                      size = 1.5, 
#                      arrow = arrow(length = unit(0.1, "inches")),
#                      colour = "black" # Also accepts "red", "blue' etc
#         ) +
#         geom_segment(data=mig_events[mig_events$time < div_events$time[2],], 
#                      aes(x = if_else(source < sink,
#                                      (source*2+2)/2+0.24,
#                                      (source*2+2)/2-0.51), 
#                          xend = if_else(source < sink,
#                                         (sink*2+2)/2-0.51,
#                                         (sink*2+2)/2+0.24),
#                          y = time,
#                          yend = time),
#                      lineend = "round", # See available arrow types in example above
#                      linejoin = "round",
#                      size = 1.5, 
#                      arrow = arrow(length = unit(0.1, "inches")),
#                      colour = "black" # Also accepts "red", "blue' etc
#         ) 
#       #print(p2)
#     }
#   } else {
#     p2 <- p1.2
#   }
#   
#   events_sorted <- events %>%
#     arrange((time))
#   
#   #MIGRATION MATRICES
#   time_interval_list <- time_interval(mig_mat_list = data$mig_mat_list,
#                                       historical_event = events_sorted,
#                                       diverge_cushion = tgap)
#   
#   
#   migmat_arrow_info_list <- lapply(split(seq_len(length(data$mig_mat_list)), seq_len(length(data$mig_mat_list))), function(i, mig_mat, time_int_list) {
#     
#     if (all(!is.na(time_int_list[[i]]$end))) {
#       
#       migmat_arrow_info(mig_mat = mig_mat[[i]],
#                         start_y = time_int_list[[i]]$start,
#                         #end_y = div_events$time[1],
#                         end_y = time_int_list[[i]]$end,
#                         offset = 0.05,
#                         arrow_replicates = 1,
#                         max_val = NULL,
#                         max_arrow_size = 1,
#                         min_arrow_size = 0.2,
#                         max_arrow_head_size = max_arrow_head_size,
#                         min_arrow_head_size = min_arrow_head_size,
#                         pop0_x = pop0_x_right,
#                         pop1_x = pop1_x_left,
#                         direction = c('forward', 'backward')[2])
#     }
#     
#   }, time_int_list = time_interval_list, mig_mat = data$mig_mat_list)
#   
#   migmat_arrow_info_list_filtered <- migmat_arrow_info_list[sapply(migmat_arrow_info_list, function(x) !is.null(x) & !all(x == 0))]
#   
#   
#   if (generic_plotting) {
#     migmat_arrow_info_list_filtered <- lapply(migmat_arrow_info_list_filtered, function(x, arrow_head_size) {
#       x$arrow_size <- 1
#       x$arrow_head_size <- arrow_head_size
#       return(x)
#     }, arrow_head_size = generic_arrow_head_size)
#   }
#   
#   
#   mig_mat_eval <- sapply(data$mig_mat_list, function(x) !all(x == 0) & x != 'NONE')
#   if (any(mig_mat_eval)) {
#     #migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0) & x != 'NONE')][[1]]
#     migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0) & any(x != 'NONE'))][[1]]
#     
#     p3 <- p2 + 
#       geom_rect(data = migration_period,
#                 aes(xmin = pop0_x_right, 
#                     xmax = pop1_x_left, 
#                     ymin = start, ymax = end), 
#                 color = NA, fill = "#bfbfbf", alpha=0.5)
#     
#   } else {
#     p3 <- p2
#   }
#   
#   #migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0))][[1]]
#   
#   
#   # p3 <- p2 + 
#   #   geom_rect(data = migration_period,
#   #             aes(xmin = pop0_x_right, 
#   #                 xmax = pop1_x_left, 
#   #                 ymin = start, ymax = end), 
#   #             color = NA, fill = "#bfbfbf", alpha=0.5)
#   
#   
#   if (length(migmat_arrow_info_list_filtered) != 0) {
#     for (i in seq_len(length(migmat_arrow_info_list_filtered))) {
#       
#       p3 <- p3 +
#         geom_segment(aes(x = migmat_arrow_info_list_filtered[[i]]$arrow_start, 
#                          xend = migmat_arrow_info_list_filtered[[i]]$arrow_end, 
#                          y = migmat_arrow_info_list_filtered[[i]]$arrow_pos, # might not want to hard code this? or maybe it's fine?
#                          yend = migmat_arrow_info_list_filtered[[i]]$arrow_pos), # might not want to hard code this? or maybe it's fine?
#                      lineend = "round", # See available arrow types in example above
#                      linejoin = "round",
#                      size = migmat_arrow_info_list_filtered[[i]]$arrow_size, 
#                      arrow = arrow(length = unit(migmat_arrow_info_list_filtered[[i]]$arrow_head_size, "inches"), type = "closed"),
#                      lty=1,
#                      colour = "black" # Also accepts "red", "blue' etc
#         )
#     }
#   }
#   
#   
#   if (generic_plotting) {
#     return(p3 + theme_void() + theme(legend.position = 'none'))
#   } else {
#     
#     return(
#       p3 +
#         theme_demovis() +
#         ylab('Time (generations)') +
#         theme(axis.title.y = element_text(size = 16)) +
#         scale_x_continuous(expand = c(0.03, 0.03)) +
#         scale_y_continuous(labels = function(x) format(x, scientific = TRUE), expand = c(0, 0)) +
#         geom_segment(data=data$hist_event_info$processed_hist_events,
#                      aes(x = pop0_x_left, xend = pop1_x_right,
#                          y=time, yend=time),
#                      lty=2, size = 1.25, lineend = "round") +
#         geom_text(aes(x = pop1_x_right*1.02,
#                       y = data$hist_event_info$processed_hist_events$time, label = data$hist_event_info$processed_hist_events$time), size = 7, hjust = 0) +
#         coord_cartesian(clip = 'off')
#       #xlim(pop0_x_left, pop1_x_right*1.1)
#     )
#     
#   }
#   
# }



#mid_sfs_subsamp_exp <- as.matrix(read.table('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/demographic_modelling/fastsimcoal/mid_region_models/mods_maxmiss80_monomorph_subsample/recent_geneflow_alt_jointMAFpop1_0.txt', header=T, row.names = 1))


# mid_joint_sfs <- plot_2d_sfs(sfs_list = list(observed = mid_obs_sfs, expected = mid_bestexpsfs),
#                              plot_obs = TRUE,
#                              plot_exp = TRUE,
#                              plot_comparison = TRUE,
#                              log10 = TRUE,
#                              drop_monomorph = TRUE,
#                              minimum_count = 10,
#                              divergent_col_scheme = c('#4575b4', '#ffffbf', '#d73027'),
#                              sequential_col_scheme = c('white', '#d73027'),
#                              scale_resolution = c(obs_exp = 1, comparison = 0.5))
# 
# 
# mid_obs_exp_plots <- plot_grid(mid_joint_sfs$plots$observed +
#             ggtitle('Observed SFS') +
#               labs(fill='Log10(SFS count)') +
#             theme(plot.title = element_text(size = 15, hjust = 0.5),
#                   legend.title = element_text(size = 8)) +
#             ylab("kazumbe") +
#             xlab("polyodon"),
#           mid_joint_sfs$plots$expected +
#             labs(fill='Log10(SFS count)') +
#             ggtitle('Expected SFS') +
#             theme(plot.title = element_text(size = 15, hjust = 0.5),
#                   legend.title = element_text(size = 8)) +
#             ylab("kazumbe") +
#             xlab("polyodon"),
#           mid_joint_sfs$plots$comparison+
#             labs(fill='Exp. SFS/Obs. SFS') +
#             ggtitle('Exp. SFS/Obs. SFS') +
#             theme(plot.title = element_text(size = 15, hjust = 0.5),
#                   legend.title = element_text(size = 8)) +
#             ylab("kazumbe") +
#             xlab("polyodon"),
#           ncol = 3
#             )
# 
#
# mid_1d_sfs <- plot_1d_sfs(sfs_list =list(observed = mid_obs_sfs, expected = mid_exp_sfs_list_list),
#                           max_freq = 16,
#                           plot = TRUE,
#                           log10 = FALSE)
# 
# #073b4c
# #ef476f
# 
# kazumbe_1d <- mid_1d_sfs$plots$column + 
#   theme_minimal() +
#   #scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_manual(values=c("#118ab2", "#ef476f")) +
#   theme(panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black", size = 1),
#         legend.position = "none",
#         plot.title = element_text(size = 20, hjust = 0.5)) +
#   ylab("Number of sites") +
#   xlab("Minor allele frequency") +
#   ggtitle("kazumbe")
# 
# polyodon_1d <- mid_1d_sfs$plots$row + 
#   theme_minimal() +
#   #scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_manual(values=c("#118ab2", "#ef476f")) +
#   theme(panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black", size = 1),
#         legend.position = "none",
#         plot.title = element_text(size = 20, hjust = 0.5)) +
#   ylab("Number of sites") +
#   xlab("Minor allele frequency") +
#   ggtitle("polyodon")
# 
# legend_1d_sfs <- get_legend(mid_1d_sfs$plots$row + 
#   theme_minimal() +
#   #scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_manual(values=c("#118ab2", "#ef476f")) +
#   theme(panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black", size = 1),
#         legend.title = element_blank(),
#         plot.title = element_text(size = 20, hjust = 0.5)) +
#   ylab("Number of sites") +
#   xlab("Minor allele frequency") +
#   ggtitle("polyodon"))
# 
# plot_grid(kazumbe_1d, polyodon_1d, legend_1d_sfs, ncol = 3, rel_widths = c(0.45, 0.45, 0.12))
# 
# 
# mid_exp_sfs_list_path <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/demographic_modelling/fastsimcoal/results_files/model_fit_results/mid/likelihood_dir/popgrow_alt_spec/recent_geneflow_dualpopgrowth_altspec_sfs/'
# mid_exp_sfs_list_names <- list.files(mid_exp_sfs_list_path)
# mid_exp_sfs_list_list <- lapply(setNames(paste0(mid_exp_sfs_list_path, mid_exp_sfs_list_names), nm = gsub(".txt", "", mid_exp_sfs_list_names)), function(FILE) {
#   as.matrix(read.table(FILE, header = T, row.names = 1))
# })
# 
# mid_bestexpsfs <- mid_exp_sfs_list_list[[which(grepl("toplikelihood$", names(mid_exp_sfs_list_list)))]]


# bagley_sfs <- as.matrix(read.table('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/demographic_modelling/Demography_fastsimcoal2_Bagley_etal_2016/MainDataParameteterEstimation/ss14_01_2_mig_01_jointMAFpop1_0.obs',
#                                    skip = 1, header = TRUE, row.names = 1))
# 
# 
# bagley_sfs_vis <- plot_2d_sfs(sfs_list = list(observed = bagley_sfs),
#                                                plot_obs = TRUE,
#                                                plot_exp = FALSE,
#                                                plot_comparison = FALSE,
#                                                log10 = TRUE,
#                                                drop_monomorph = TRUE,
#                                                minimum_count = 10,
#                                                divergent_col_scheme = c('#4575b4', '#ffffbf', '#d73027'),
#                                                sequential_col_scheme = c('white', '#d73027'),
#                                                scale_resolution = c(obs_exp = 1, comparison = 0.5))




# north_fsc_results <- cowplot::plot_grid(fsc_vis_components_list$north$par + theme(plot.margin = margin(5, 55, 5, 5)), 
#                                         ggdraw(egg::ggarrange(fsc_vis_components_list$north$aic, fsc_vis_components_list$north$likelihood, ncol = 2)) + theme(plot.margin = margin(5, 5, 5, 35)),
#                                         labels = c('(a)', '(b)'), rel_widths = c(0.35, 0.65), label_size = 18)
# 
# 
# north_fsc_results_plustitle <- plot_grid(region_title$North, north_fsc_results, nrow = 2, rel_heights = c(0.1, 0.9))
# 
# mid_fsc_results <- cowplot::plot_grid(best_fit_mods_vis$north + theme(plot.margin = margin(5, 55, 5, 5)), 
#                                       ggdraw(egg::ggarrange(prac_plot1, lik_plot_north, ncol = 2)) + theme(plot.margin = margin(5, 5, 5, 35)),
#                                       labels = c('(a)', '(b)'), rel_widths = c(0.35, 0.65))
# 
# 
# multipanel1 <- cowplot::plot_grid(north_fsc_results_plustitle + theme(plot.margin = margin(0, 0, 5, 0)), 
#                    north_fsc_results_plustitle + theme(plot.margin = margin(5, 0, 0, 0)), nrow = 2)
# 
# 
# 
# 
# aic_plot_north <- aic_list_processed$north %>% 
#   ggplot(aes(x = deltaAIC, y = mod_name)) +
#   geom_point(aes(color = mod_class), size = 5) +
#   xlab('\u0394AIC') +
#   scale_color_manual(values = mod_type_cols) +
#   theme_cowplot() +
#   theme(plot.margin = margin(5.5, 0, 0, 10),
#         panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
#         legend.title = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_text(size = 14),
#         axis.title.x = element_text(size = 17, margin = margin(t = 0, r = 11, b = 0, l = 0)),
#         axis.text.x = element_text(size = 14),
#         legend.text = element_text(size = 18, margin = margin(r = 22, unit = "pt")),
#         legend.position = "none",
#         legend.spacing.x = unit(0.3, 'cm'))
# 
# lik_plot_north <- aic_list_processed$north %>% 
#   ggplot(aes(x = deltaAIC, y = mod_name)) +
#   geom_point(aes(color = mod_class), size = 5) +
#   xlab('Likelihood') +
#   scale_color_manual(values = mod_type_cols) +
#   theme_cowplot() +
#   theme(plot.margin = margin(5.5, 10, 0, 0),
#         panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(size = 14),
#         axis.title.x = element_text(size = 17, margin = margin(t = 0, r = 11, b = 0, l = 0)),
#         legend.position = "none",
#         legend.text = element_text(color = "white"),
#         legend.title = element_text(color = "white"),
#         legend.key = element_rect(fill = "white"))
# 
# north_fsc_results <- cowplot::plot_grid(best_fit_mods_vis$north + theme(plot.margin = margin(5, 55, 5, 5)), 
#                    ggdraw(egg::ggarrange(aic_plot_north, lik_plot_north, ncol = 2)) + theme(plot.margin = margin(5, 5, 5, 35)),
#                    labels = c('(a)', '(b)'), rel_widths = c(0.35, 0.65), label_size = 18)
# 
# 
# north_fsc_results_plustitle <- plot_grid(region_title$North, north_fsc_results, nrow = 2, rel_heights = c(0.1, 0.9))
# 
# mid_fsc_results <- cowplot::plot_grid(best_fit_mods_vis$north + theme(plot.margin = margin(5, 55, 5, 5)), 
#                                         ggdraw(egg::ggarrange(prac_plot1, lik_plot_north, ncol = 2)) + theme(plot.margin = margin(5, 5, 5, 35)),
#                                         labels = c('(a)', '(b)'), rel_widths = c(0.35, 0.65))
# 
# 
# cowplot::plot_grid(north_fsc_results_plustitle + theme(plot.margin = margin(0, 0, 5, 0)), 
#                    north_fsc_results_plustitle + theme(plot.margin = margin(5, 0, 0, 0)), nrow = 2)




# 
# 
# 
# 
# 
# pop.cols <- viridis::viridis(2)
# 
# 
# 
# best_fit_mods_par_bs_update$mid$hist_event_info$processed_hist_events
# events <- best_fit_mods_par_bs_update$mid$hist_event_info$processed_hist_events
# 
# 
# ############################
# # subsetting data into migration vs div events
# ############################
# mig_events <- events %>%
#   filter(migrants < 1) %>%
#   filter(source != sink)
# div_events <- events %>%
#   filter(migrants == 1) %>%
#   filter(source != sink)
# 
# pop_size_list <- list(current = best_fit_mods_par_bs_update$mid$pop_size,
#                       historical = calc_historical_popsize(current_pop = best_fit_mods_par_bs_update$mid$pop_size, 
#                                                            growth_rate = best_fit_mods_par_bs_update$mid$growth_rate, 
#                                                            time = div_events$time[1]))
# 
# # if (scale_pops) {
# #   pops_rel_ne <- rescale(data$pop_size, 0.5, 1.2)
# #   #root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[1], 0.5, 1.2)
# #   root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1], 0.5, 1.2)
# #   div_events$source_rel_ne <- pops_rel_ne[div_events$source + 1]
# #   div_events$sink_rel_ne <- pops_rel_ne[div_events$sink + 1]
# #   
# #   mig_events$source_rel_ne <- pops_rel_ne[mig_events$source + 1]
# #   mig_events$sink_rel_ne <- pops_rel_ne[mig_events$sink + 1]
# # }
# 
# if (scale_pops == TRUE) {
#   root_pop_size <- div_events$new_deme_size[div_events$time == max(div_events$time)]*pop_size_list$historical[div_events$sink + 1]
#   root_rel_ne <- rescale_alt(root_pop_size, 0.1, 1.2, c(unlist(pop_size_list), root_pop_size))
#   
#   pops_rel_ne <- rescale_alt(pop_size_list$current, 0.1, 1.2, pop_size_vec = c(unlist(pop_size_list), root_pop_size) )
#   hist_pops_rel_ne <- rescale_alt(pop_size_list$historical, 0.1, 1.2, pop_size_vec = c(unlist(pop_size_list), root_pop_size) )
#   
#   mig_events$source_rel_ne <- pops_rel_ne[mig_events$source + 1]
#   mig_events$sink_rel_ne <- pops_rel_ne[mig_events$sink + 1]
# }
# 
# width <- 0.75
# pop_poly_list <- lapply(1:2, function(X, pops_rel_ne, hist_pops_rel_ne, width, div_events) {
#   
#   data.frame(y = rep(c(0, div_events$time[1]), each  = 2),
#              # x = c(div_events[[switch(X, `1` = 'source', `2` = 'sink')]][1] + 1 - width*pops_rel_ne[X]/2,
#              #       div_events[[switch(X, `1` = 'source', `2` = 'sink')]][1] + 1 + width*pops_rel_ne[X]/2,
#              #       div_events[[switch(X, `1` = 'source', `2` = 'sink')]][1] + 1 + width*hist_pops_rel_ne[X]/2,
#              #       div_events[[switch(X, `1` = 'source', `2` = 'sink')]][1] + 1 - width*hist_pops_rel_ne[X]/2)
#              x = c((X - 1) + 1 - width*pops_rel_ne[X]/2,
#                    (X - 1) + 1 + width*pops_rel_ne[X]/2,
#                    (X - 1) + 1 + width*hist_pops_rel_ne[X]/2,
#                    (X - 1) + 1 - width*hist_pops_rel_ne[X]/2)
#   )
#   
# }, pops_rel_ne = pops_rel_ne, hist_pops_rel_ne = hist_pops_rel_ne, width = width, div_events = div_events)
# 
# 
# #return(pop_poly_list)
# #3/30/2022
# #pop_poly_list <- pop_poly_list[match(div_events[c('source', 'sink')], 0:1)]
# pop1_pos <- which(div_events[c('source', 'sink')] == 0)
# pop2_pos <- which(div_events[c('source', 'sink')] == 1)
# 
# 
# ############################
# # first, set up the skeleton using divergences
# # it'd be good if we can make this more flexible (i.e. works with any number of div events)
# ############################
# tgap <- max(events$time)*0.10
# troot <- max(events$time)*1.33
# 
# p <- ggplot() +
#   xlim(0,length(best_fit_mods_par_bs_update$mid$pop_size)*1.5)+
#   ylim(0,troot)
# 
# #width <- 0.75
# p1 <- p + 
#   geom_rect(mapping=aes(xmin=min(pop_poly_list[[1]][pop_poly_list[[1]]$y != 0,]$x),
#                        xmax=max(pop_poly_list[[2]][pop_poly_list[[1]]$y != 0,]$x),
#                        ymin=div_events$time[1],ymax=div_events$time[1]+tgap), fill=pop_color) +
#   # geom_rect(mapping=aes(xmin=min(pop_poly_list[[pop1_pos]][pop_poly_list[[pop1_pos]]$y != 0,]$x),
#   #                       xmax=max(pop_poly_list[[pop2_pos]][pop_poly_list[[pop2_pos]]$y != 0,]$x),
#   #                       ymin=div_events$time[1],ymax=div_events$time[1]+tgap), fill=pop_color) + 
#   # root
#   geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width*root_rel_ne/2,
#                         xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width*root_rel_ne/2,
#                         ymin=div_events$time[1],
#                         ymax=troot), fill=pop_color) + 
#   geom_polygon(data = pop_poly_list[[1]],
#                aes(x = x, y = y), fill=pop_color)  +
#   geom_polygon(data = pop_poly_list[[2]],
#                aes(x = x, y = y), fill=pop_color)
#   # geom_polygon(data = pop_poly_list[[pop1_pos]],
#   #              aes(x = x, y = y), fill=pop_color)  + 
#   # geom_polygon(data = pop_poly_list[[pop2_pos]],
#   #              aes(x = x, y = y), fill=pop_color) 



# ################################################
# ### MULTIPANEL PLOT OF ALL TESTED FSC MODELS ###
# ################################################
# 
# 
# par_file_list_processed <- lapply(par_file_list, function(region) {
#   region[!grepl(".*change_geneflow.*", names(region))]
# })
# 
# 
# ### process names ###
# #1. remove all underscores
# names1 <- gsub(pattern = '_', replacement = ' ', names(par_file_list_processed$mid))
# 
# first_word <- unique(word(gsub(pattern = '_', replacement = ' ', names(par_file_list_processed$mid)), start = 1, end = 1))
# 
# par_multi_level_list <- lapply(setNames(nm = first_word), function(first_w, name_vec, par_list) {
#   processed_names <- word(gsub(pattern = '_', replacement = ' ', name_vec), start = 1, end = 1)
#   full_names <- name_vec[processed_names %in% first_w]
#   
#   par_list[names(par_list) %in% full_names]
# }, name_vec = names(par_file_list_processed$mid), par_list = par_file_list_processed$mid)
# 
# new_names_vec <- data.frame(orig_name = names(par_multi_level_list)) %>% 
#   mutate(new_names = case_when(orig_name == "early" ~ "Early gene flow",
#                                orig_name == "continuous" ~ "Continuous gene flow",
#                                orig_name == "no" ~ "No gene flow",
#                                orig_name == "recent" ~ "Recent gene flow")) %>% 
#   pull(new_names)
# 
# names(par_multi_level_list) <- new_names_vec
# 
# par_full_updated_names <- lapply(par_multi_level_list, function(mod_list) {
#   
#   names(mod_list)[!grepl(pattern = "unidirect0|symmetric|unidirect1", names(mod_list))] <- "Asymmetric"
#   names(mod_list)[grepl(pattern = "symmetric_dualpopgrowth_altspec_maxL.par", names(mod_list))] <- "Symmetric"
#   names(mod_list)[grepl(pattern = "unidirect0_dualpopgrowth_altspec_maxL.par", names(mod_list))] <- "Unidirect (Pk to Pp)"
#   names(mod_list)[grepl(pattern = "unidirect1_dualpopgrowth_altspec_maxL.par", names(mod_list))] <- "Unidirect (Pp to Pk)"
#   return(mod_list)
# })
# 
# 
# 
# full_plot_list_processed <- lapply(setNames(nm = names(par_full_updated_names)), function(MULTIPANEL_NAME, par_list) {
#   multipanel_title <- ggdraw() + 
#     draw_label(
#       MULTIPANEL_NAME,
#       x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
#       fontface = "bold", size = 18
#     )
#   
#   plot_list1 <- lapply(par_list[[MULTIPANEL_NAME]], function(mod) {
#     plot_demographic_mod(par = mod,
#                          colors = TRUE,
#                          pop_colors = c("#DC7633", "#3498DB"),
#                          scale_pops = FALSE,
#                          pop_color = '#a3a3a3',
#                          generic_plotting = TRUE,
#                          generic_arrow_head_size = 0.07)
#   } )
#   
#   
#   processed_plot_list <- lapply(setNames(nm = names(plot_list1)), function(TITLE, plot_list, mod_class) {
#     if (mod_class == "No gene flow") {
#       plot_list[[TITLE]] + 
#         ggtitle(TITLE) +
#         theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white"))
#     } else {
#       plot_list[[TITLE]] + 
#         ggtitle(TITLE) +
#         theme(plot.title = element_text(hjust = 0.5, size = 17))
#     }
#     
#   }, plot_list = plot_list1, mod_class = MULTIPANEL_NAME)
#   
#   if (length(processed_plot_list) == 1) {
#     return(
#       plot_grid(multipanel_title, 
#                 plot_grid(processed_plot_list[[1]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")), 
#                           nrow = 4 ) +
#                   theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)),
#                 rel_heights = c(0.07, 1),
#                 nrow = 2)
#     )
#   } else {
#     return(
#       plot_grid(multipanel_title, 
#                 plot_grid(processed_plot_list[["Symmetric"]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           processed_plot_list[["Asymmetric"]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           processed_plot_list[["Unidirect (Pk to Pp)"]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           processed_plot_list[["Unidirect (Pp to Pk)"]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           nrow = 4) +
#                   theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)),
#                 rel_heights = c(0.07, 1),
#                 nrow = 2)
#     )
#   }
#   
#   return(processed_plot_list)
# }, par_list = par_full_updated_names)
# 
# 
# 
# plot_grid(
#   full_plot_list_processed$`Continuous gene flow`,
#   full_plot_list_processed$`Early gene flow`,
#   full_plot_list_processed$`Recent gene flow`,
#   full_plot_list_processed$`No gene flow`,
#   ncol = 4
# ) +
#   theme(plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm"))
# 
# ggsave(here('figures', 'fastsimcoal_allmods.png'), 
#        width = 1.4*20, height = 1.4*16.5, units = "cm", bg = "white")
# 
# 
# 
# ################################################
# ### MULTIPANEL PLOT OF ALL TESTED FSC MODELS ###
# ################################################
# 
# par_file_list_processed <- lapply(par_file_list, function(region) {
#   region[!grepl(".*change_geneflow.*", names(region))]
# })
# 
# 
# 
# ### process names ###
# #1. remove all underscores
# library(grid)
# names1 <- gsub(pattern = '_', replacement = ' ', names(par_file_list_processed$mid))
# 
# first_word <- unique(word(gsub(pattern = '_', replacement = ' ', names(par_file_list_processed$mid)), start = 1, end = 1))
# 
# par_multi_level_list <- lapply(setNames(nm = first_word), function(first_w, name_vec, par_list) {
#   processed_names <- word(gsub(pattern = '_', replacement = ' ', name_vec), start = 1, end = 1)
#   full_names <- name_vec[processed_names %in% first_w]
#   
#   par_list[names(par_list) %in% full_names]
# }, name_vec = names(par_file_list_processed$mid), par_list = par_file_list_processed$mid)
# 
# new_names_vec <- data.frame(orig_name = names(par_multi_level_list)) %>% 
#   mutate(new_names = case_when(orig_name == "early" ~ "Early gene flow",
#                                orig_name == "continuous" ~ "Continuous gene flow",
#                                orig_name == "no" ~ "No gene flow",
#                                orig_name == "recent" ~ "Recent gene flow")) %>% 
#   pull(new_names)
# 
# names(par_multi_level_list) <- new_names_vec
# 
# par_full_updated_names <- lapply(par_multi_level_list, function(mod_list) {
#   
#   names(mod_list)[!grepl(pattern = "unidirect0|symmetric|unidirect1", names(mod_list))] <- "Asymmetric"
#   names(mod_list)[grepl(pattern = "symmetric_dualpopgrowth_altspec_maxL.par", names(mod_list))] <- "Symmetric"
#   names(mod_list)[grepl(pattern = "unidirect0_dualpopgrowth_altspec_maxL.par", names(mod_list))] <- "Unidirect (Pk to Pp)"
#   names(mod_list)[grepl(pattern = "unidirect1_dualpopgrowth_altspec_maxL.par", names(mod_list))] <- "Unidirect (Pp to Pk)"
#   return(mod_list)
# })
# 
# par_full_updated_names$`Continuous gene flow`$Asymmetric$hist_event_info$processed_hist_events
# 
# 
# mod_type_cols <- setNames(c('#cb5cbb', '#FFC60A', '#2A9D8F', '#E24D28'), 
#                           nm = c('Continuous gene flow', 'Early gene flow', 'Recent gene flow', 'No gene flow'))
# 
# mod_type_cols_background <- setNames(c('#efceea', '#fff3ce', '#d4ebe8', '#f6c9be'), 
#                                      nm = c('Continuous gene flow', 'Early gene flow', 'Recent gene flow', 'No gene flow'))
# 
# 
# plot_list <- list()
# 
# for (CLASS in names(par_full_updated_names)) {
#   for (MOD in names(par_full_updated_names[[CLASS]])) {
#     
#     plot_list[[CLASS]][[MOD]] <- plot_demographic_mod(par = par_full_updated_names[[CLASS]][[MOD]],
#                                                       colors = TRUE,
#                                                       pop_colors = c("#DC7633", "#3498DB"),
#                                                       scale_pops = FALSE,
#                                                       pop_color = '#929292', #a3a3a3
#                                                       mig_band_color = "#bfbfbf",
#                                                       generic_plotting = TRUE,
#                                                       generic_arrow_head_size = 0.09)
#     
#     plot_list[[CLASS]][[MOD]] <- plot_list[[CLASS]][[MOD]] + 
#       theme(panel.background = element_rect(fill = mod_type_cols_background[[CLASS]], color = mod_type_cols_background[[CLASS]]))
#     
#     #if ( (CLASS != "No gene flow") & (MOD %in% c("Symmetric", "Unidirect (Pk to Pp)"))) {
#     # if (CLASS  %in% c("Continuous gene flow", "Recent gene flow")) {
#     #   plot_list[[CLASS]][[MOD]] <- plot_list[[CLASS]][[MOD]] + 
#     #     theme(panel.background = element_rect(fill = mod_type_cols_background[[CLASS]], color = mod_type_cols_background[[CLASS]]))
#     # }
#     
#   }
# }
# 
# 
# title_list <- lapply(setNames(nm = c("Symmetric", "Asymmetric", "Unidirect (Pk to Pp)", "Unidirect (Pp to Pk)")), function(title) {
#   textGrob(title, 
#            gp=gpar(fontface="plain", col="black", fontsize = 14), rot = 90)
# })
# 
# plot_multipanel_list <- lapply(setNames(nm = c("Symmetric", "Asymmetric", "Unidirect (Pk to Pp)", "Unidirect (Pp to Pk)")), function(x, title, plot_list) {
#   
#   multipanel1 <- grid.arrange(arrangeGrob(plot_grid(plot_list$`Continuous gene flow`[[x]],
#                                                     plot_list$`Early gene flow`[[x]],
#                                                     plot_list$`Recent gene flow`[[x]], 
#                                                     nrow = 1) +
#                                             theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)),
#                                           #theme(panel.border = element_rect(colour = "black", fill=NA, size = 1),
#                                           #       plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm")), 
#                                           left = title[[x]], padding = unit(1, "line")))
#   
#   return(multipanel1)
#   
#   # if (x == 'Asymmetric') {
#   #   multipanel2 <- ggdraw(grid.arrange(multipanel1, 
#   #                                      plot_list$`No gene flow`$Asymmetric, 
#   #                                      ncol = 2, widths = c(0.75, 0.25)))
#   # } else {
#   #   multipanel2 <- ggdraw(grid.arrange(multipanel1, 
#   #                                      ggplot() + theme_void(), 
#   #                                      ncol = 2, widths = c(0.75, 0.25)))
#   # }
#   # return(multipanel2)
# }, title = title_list, plot_list = plot_list)
# 
# multipanel_fsc_mods_wout_titles <- plot_grid(plot_multipanel_list$Symmetric,
#                                              plot_multipanel_list$Asymmetric,
#                                              plot_multipanel_list$`Unidirect (Pk to Pp)`,
#                                              plot_multipanel_list$`Unidirect (Pp to Pk)`, nrow = 4) +
#   theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
# 
# multipanel_fsc_mods_wout_titles
# 
# 
# ggsave(here('figures', 'fastsimcoal_allmods_maintext_wout_titles.png'),
#        plot = multipanel_fsc_mods_wout_titles,
#        width = 1.3*22, height = 1.3*13.5, units = "cm", bg = "white")
# 
# 
# continuous_gf_title <- ggdraw() + 
#   draw_label(
#     "Continuous gene flow",
#     hjust = 0.38,
#     size = 21
#   )
# 
# early_gf_title <- ggdraw() + 
#   draw_label(
#     "Early gene flow",
#     hjust = 0.415,
#     size = 21
#   )
# 
# recent_gf_title <- ggdraw() + 
#   draw_label(
#     "Recent gene flow",
#     hjust = 0.5,
#     size = 21
#   )
# 
# 
# 
# multipanel_fsc_mods_w_titles <- plot_grid(plot_grid(continuous_gf_title, early_gf_title, recent_gf_title, ncol = 3),
#                                           multipanel_fsc_mods_wout_titles, rel_heights = c(0.05, 0.95), nrow = 2)
# 
# 
# multipanel_fsc_results <- cowplot::plot_grid(multipanel1, tree_legend, 
#                                              ncol = 2, 
#                                              rel_widths = c(0.92, 0.08)) + theme(plot.margin = margin(0, 0, 0, 12))
# 
# cowplot::plot_grid(multipanel_fsc_results,
#                    multipanel_fsc_mods_w_titles,
#                    nrow = 2, rel_heights = c(0.5, 0.5))
# 
# ggsave(here('figures', 'fastsimcoal_allmods_maintext_w_titles.png'),
#        plot = multipanel_fsc_mods_w_titles,
#        width = 1.3*22, height = 1.3*13.5, units = "cm", bg = "white")
# 
# 
# grid.arrange(plot_grid(ggplot() +
#                          theme_void() +
#                          theme(panel.background = element_rect(fill = 'blue')),
#                        early_gf_title, 
#                        early_gf_title, 
#                        early_gf_title, 
#                        early_gf_title, ncol = 5, rel_widths = c(0.05, 0.5, 0.5, 0.5, 0.5)), 
#              prac_multipanel, heights = c(0.15, 0.85)) 
#
###############################################
#
# #add to plot
# 
# prac_arrange <- grid.arrange(arrangeGrob(plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#                                                    plot_list$`Early gene flow`$Symmetric,
#                                                    plot_list$`Recent gene flow`$Symmetric, nrow = 1) +
#                                            theme(panel.border = element_rect(colour = "black", fill=NA, size = 1),
#                                                  plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm")), 
#                                          left = y.grob))
# 
# 
# plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#           plot_list$`Early gene flow`$Symmetric,
#           plot_list$`Recent gene flow`$Symmetric, nrow = 1) +
#   theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))
# 
# no_gene_flow_plot <- grid.arrange(title_blank, 
#                                   plot_grid(plot_list$`No gene flow`$Asymmetric, nrow = 1), 
#                                   heights=c(0.15, 0.85))
# 
# 
# plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#           plot_list$`Early gene flow`$Symmetric,
#           plot_list$`Recent gene flow`$Symmetric, nrow = 1) +
#   theme(panel.border = element_rect(colour = "black", fill=NA, size = 1),
#         plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm"))
# 
# 
# 
# y.grob_symmetric <- textGrob("Bidirectional\nsymmetric", 
#                              gp=gpar(fontface="plain", col="black", fontsize=15), rot=90)
# 
# y.grob_asymmetric <- textGrob("Asymmetric gene flow", 
#                               gp=gpar(fontface="plain", col="black", fontsize=15), rot=90)
# 
# 
# #add to plot
# 
# prac_arrange <- grid.arrange(arrangeGrob(plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#                                                    plot_list$`Early gene flow`$Symmetric,
#                                                    plot_list$`Recent gene flow`$Symmetric, nrow = 1) +
#                                            theme(panel.border = element_rect(colour = "black", fill=NA, size = 1),
#                                                  plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm")), 
#                                          left = y.grob))
# 
# nogf_arrange <- grid.arrange(arrangeGrob(plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#                                                    plot_list$`Early gene flow`$Symmetric,
#                                                    plot_list$`Recent gene flow`$Symmetric, 
#                                                    nrow = 1) +
#                                            theme(panel.border = element_rect(colour = "black", fill=NA, size = 1),
#                                                  plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm")), 
#                                          left = y.grob))
# 
# 
# prac_multipanel <- ggdraw(grid.arrange(prac_arrange, 
#                                        plot_list$`No gene flow`$Asymmetric + theme(plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm")), 
#                                        ncol = 2, widths = c(0.75, 0.25)))
# 
# 
# prac_multipanel +
#   theme(axis.ticks = element_line(color = 'black'))
# 
# 
# grid.arrange(plot_grid(ggplot() +
#                          theme_void() +
#                          theme(panel.background = element_rect(fill = 'blue')),
#                        early_gf_title, 
#                        early_gf_title, 
#                        early_gf_title, 
#                        early_gf_title, ncol = 5, rel_widths = c(0.05, 0.5, 0.5, 0.5, 0.5)), 
#              prac_multipanel, heights = c(0.15, 0.85)) 
# 
# ggplot() +
#   theme_void() +
#   theme(panel.background = element_rect(fill = 'blue'))
# 
# early_gf_title <- ggdraw() + 
#   draw_label(
#     "Early gene flow",
#     hjust = 0.5
#   )
# 
# 
# grid.arrange(prac_arrange, 
#              ggplot() + theme_void(), 
#              ncol = 2, widths = c(0.75, 0.25))
# 
# plot_list$`Continuous gene flow`
# 
# 
# 
# 
# 
# 
# 
# library(grid)
# 
# title_list <- lapply(setNames(nm = c("Symmetric gene flow", "Asymmetric gene flow", "Unidirect (Pk to Pp)", "Unidirect (Pp to Pk)")), function(title) {
#   grobTree(rectGrob(gp=gpar(fill="black")),
#            textGrob(title, x=0.5, hjust=0.5,
#                     gp=gpar(col="white", cex=1.3)))
# })
# 
# title_blank <- grobTree(rectGrob(gp=gpar(fill = "white", col = "white")),
#                         textGrob('', x=0.5, hjust=0.5,
#                                  gp=gpar(col = "white", cex=1.3)))
# 
# 
# symmetric_multipanel <- grid.arrange(title_list$`Symmetric gene flow`, 
#                                      plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#                                                plot_list$`Early gene flow`$Symmetric,
#                                                plot_list$`Recent gene flow`$Symmetric, nrow = 1) +
#                                        theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)), 
#                                      heights=c(0.15, 0.85))
# 
# asymmetric_multipanel <- grid.arrange(title_list$`Asymmetric gene flow`, 
#                                       plot_grid(plot_list$`Continuous gene flow`$Asymmetric,
#                                                 plot_list$`Early gene flow`$Asymmetric,
#                                                 plot_list$`Recent gene flow`$Asymmetric, nrow = 1) +
#                                         theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)), 
#                                       heights=c(0.15, 0.85))
# 
# unidirect01_multipanel <- grid.arrange(title_list$`Unidirect (Pk to Pp)`, 
#                                        plot_grid(plot_list$`Continuous gene flow`$`Unidirect (Pk to Pp)`,
#                                                  plot_list$`Early gene flow`$`Unidirect (Pk to Pp)`,
#                                                  plot_list$`Recent gene flow`$`Unidirect (Pk to Pp)`, nrow = 1) +
#                                          theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)), 
#                                        heights=c(0.15, 0.85))
# 
# unidirect10_multipanel <- grid.arrange(title_list$`Unidirect (Pp to Pk)`, 
#                                        plot_grid(plot_list$`Continuous gene flow`$`Unidirect (Pp to Pk)`,
#                                                  plot_list$`Early gene flow`$`Unidirect (Pp to Pk)`,
#                                                  plot_list$`Recent gene flow`$`Unidirect (Pp to Pk)`, nrow = 1) +
#                                          theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)), 
#                                        heights=c(0.15, 0.85))
# 
# 
# no_gene_flow_plot <- grid.arrange(title_blank, 
#                                   plot_grid(plot_list$`No gene flow`$Asymmetric, nrow = 1), 
#                                   heights=c(0.15, 0.85))
# 
# 
# grid.arrange()
# no_gene_flow_plot
# 
# gf_multipanel <- grid.arrange(symmetric_multipanel,
#                               asymmetric_multipanel,
#                               unidirect01_multipanel,
#                               unidirect10_multipanel, nrow = 4)
# 
# nogf_multipanel <- grid.arrange(no_gene_flow_plot,
#                                 ggplot() + theme_void(),
#                                 ggplot() + theme_void(),
#                                 ggplot() + theme_void(), nrow = 4)
# 
# grid.arrange(gf_multipanel, nogf_multipanel, ncol = 2, widths = c(0.75, 0.25))
# 
# ggplot() + theme_void()
# 
# title_list
# 
# plot_grid(
#   plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#             plot_list$`Early gene flow`$Symmetric,
#             plot_list$`Recent gene flow`$Symmetric, nrow = 1) +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)),
#   plot_grid(plot_list$`Continuous gene flow`$Asymmetric,
#             plot_list$`Early gene flow`$Asymmetric,
#             plot_list$`Recent gene flow`$Asymmetric, nrow = 1) +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)),
#   plot_grid(plot_list$`Continuous gene flow`$`Unidirect (Pk to Pp)`,
#             plot_list$`Early gene flow`$`Unidirect (Pk to Pp)`,
#             plot_list$`Recent gene flow`$`Unidirect (Pk to Pp)`, nrow = 1) +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)),
#   plot_grid(plot_list$`Continuous gene flow`$`Unidirect (Pp to Pk)`,
#             plot_list$`Early gene flow`$`Unidirect (Pp to Pk)`,
#             plot_list$`Recent gene flow`$`Unidirect (Pp to Pk)`, nrow = 1) +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size = 1)),
#   nrow = 4
# ) +
#   theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))
# 
# 
# 
# 
# plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#           plot_list$`Early gene flow`$Symmetric,
#           plot_list$`Recent gene flow`$Symmetric, nrow = 1)
# 
# plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#           plot_list$`Early gene flow`$Symmetric,
#           plot_list$`Recent gene flow`$Symmetric, nrow = 1)
# 
# plot_grid(plot_list$`Continuous gene flow`$Symmetric,
#           plot_list$`Continuous gene flow`$Asymmetric,
#           plot_list$`Continuous gene flow`$`Unidirect (Pk to Pp)`,
#           plot_list$`Continuous gene flow`$`Unidirect (Pp to Pk)`, nrow = 4)
# 
# 
# 
# 
# plot_list <- lapply(par_full_updated_names, function(CLASS) {
#   
#   lapply(CLASS, function(MOD) {
#     return(
#       plot_demographic_mod(par = MOD,
#                            colors = TRUE,
#                            pop_colors = c("#DC7633", "#3498DB"),
#                            scale_pops = FALSE,
#                            pop_color = '#a3a3a3',
#                            mig_band_color = '#10FFCB',
#                            generic_plotting = TRUE,
#                            generic_arrow_head_size = 0.1) +
#         theme(panel.background = element_rect(fill = '#E0E0E0', color = '#E0E0E0'))
#     )
#   })
#   
# })
# 
# plot_list$`Continuous gene flow`$Symmetric
# plot_list$`Recent gene flow`$Asymmetric
# 
# 
# 
# full_plot_list_processed <- lapply(setNames(nm = names(par_full_updated_names)), function(MULTIPANEL_NAME, par_list) {
#   multipanel_title <- ggdraw() + 
#     draw_label(
#       MULTIPANEL_NAME,
#       x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
#       fontface = "bold", size = 18
#     )
#   
#   plot_list1 <- lapply(par_list[[MULTIPANEL_NAME]], function(mod) {
#     plot_demographic_mod(par = mod,
#                          colors = TRUE,
#                          pop_colors = c("#DC7633", "#3498DB"),
#                          scale_pops = FALSE,
#                          pop_color = '#a3a3a3',
#                          generic_plotting = TRUE,
#                          generic_arrow_head_size = 0.07)
#   } )
#   
#   
#   processed_plot_list <- lapply(setNames(nm = names(plot_list1)), function(TITLE, plot_list, mod_class) {
#     
#     plot_list[[TITLE]]
#     
#     # if (mod_class == "No gene flow") {
#     #   plot_list[[TITLE]] + 
#     #     ggtitle(TITLE) +
#     #     theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white"))
#     # } else {
#     #   plot_list[[TITLE]] + 
#     #     ggtitle(TITLE) +
#     #     theme(plot.title = element_text(hjust = 0.5, size = 17))
#     # }
#     
#   }, plot_list = plot_list1, mod_class = MULTIPANEL_NAME)
#   
#   if (length(processed_plot_list) == 1) {
#     return(
#       plot_grid(multipanel_title, 
#                 plot_grid(processed_plot_list[[1]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")), 
#                           nrow = 4 ) +
#                   theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)),
#                 rel_heights = c(0.07, 1),
#                 nrow = 2)
#     )
#   } else {
#     return(
#       plot_grid(multipanel_title, 
#                 plot_grid(processed_plot_list[["Symmetric"]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           processed_plot_list[["Asymmetric"]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           processed_plot_list[["Unidirect (Pk to Pp)"]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           processed_plot_list[["Unidirect (Pp to Pk)"]] + theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm")),
#                           nrow = 4) +
#                   theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)),
#                 rel_heights = c(0.07, 1),
#                 nrow = 2)
#     )
#   }
#   
#   return(processed_plot_list)
# }, par_list = par_full_updated_names)
# 
# full_plot_list_processed$`Continuous gene flow`
# 
# 
# plot_grid(
#   full_plot_list_processed$`Continuous gene flow`,
#   full_plot_list_processed$`Early gene flow`,
#   full_plot_list_processed$`Recent gene flow`,
#   full_plot_list_processed$`No gene flow`,
#   ncol = 4
# ) +
#   theme(plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm"))
# 
# ggsave(here('figures', 'fastsimcoal_allmods.png'), 
#        width = 1.4*20, height = 1.4*16.5, units = "cm", bg = "white")

