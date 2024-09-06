library(targets)
library(tarchetypes)
library(tidyr)
source("custom_functions.R")
source("arguments.R")

if(!exists("./figures")){
  dir.create("./figures")
}

if(!exists("./diagnostics")){
  dir.create("./diagnostics")
}

if(!exists("./tables")){
  dir.create("./tables")
}

# Set target-specific options such as packages.
tar_option_set(packages = c("Hmsc", "reshape2", "magrittr", "ggplot2", "abind", "dplyr", "dismo", "corrplot", "cooccur", "vegan", "tidyverse"))

values <- expand_grid(p_target = p, i_target = i)

mapped <- tar_map(values = values,
        unlist = FALSE,
        tar_target(sampled_community, get.sample(original_community, p_target, size= 500, i = i_target)),
        tar_target(same_sample_original_community, original.same.sample(sampled_community, original_community)),
        tar_target(hmsc_model, calibrate.hmsc(data=sampled_community, thin=thin, samples=samples,
                                              transient=transient, nChains=nchains, nParallel=nParallel)),
        tar_target(diagnostics_plots, plot.diagnostics(hmsc_model, p_target)),
        tar_target(explanatory_power, tuning.metrics(com_testing= sampled_community, com_training= sampled_community,
                                                      model= hmsc_model, sp_names= sp_names, p_target, i= i_target)),
        tar_target(predictive_power, tuning.predictive(com_testing= sampled_community, model= hmsc_model, sp_names= sp_names, p_target, i= i_target)),
        tar_target(same_sample_metrics, tuning.metrics(com_testing = same_sample_original_community, com_training= sampled_community,
                                                                  model= hmsc_model, sp_names= sp_names, p_target, i=i_target)),
        tar_target(individual_betas, BETAS(data = artificial_data, model = hmsc_model, p = p_target, i= i_target))
)

spcor_mapped <- tar_map(values = values,
                       unlist = FALSE,
                       tar_target(sampled_community_spcor, get.sample(original_community_spcor, p_target, size= 500, i= i_target)),
                       tar_target(same_sample_original_community_spcor, original.same.sample(sampled_community_spcor, original_community_spcor)),
                       tar_target(false_negatives, compare.matrices(sp_names, sampled_community_spcor, same_sample_original_community_spcor, i_target, p_target)),
                       tar_target(hmsc_model_spcor, calibrate.hmsc(data=sampled_community_spcor, thin=thin, samples=samples,
                                                         transient=transient, nChains=nchains, nParallel=nParallel)),
                       tar_target(to_plot_corr, to_plot(hmsc_model_spcor, p_target)),
                       tar_target(coo_to_ggplot, coo_pattern(sampled_community_spcor, sp_names, p_target, i_target, TRUE)),
                       tar_target(corr_to_ggplot, corr.extraction(hmsc_model_spcor)),
                       tar_target(corr_pattern,  coo_pattern(corr_to_ggplot, sample = p_target, i = i_target)),
                       tar_target(procrus, procrustes.analysis(i_target, p_target, corr_pattern, artificial_data_spcor)),
                       tar_target(spat_data, spatial.prediction(p = p_target, i = i_target, env_xy_data = original_community_spcor, model = hmsc_model_spcor, sp_names = sp_names))
                       )

list(
  
  ###############################PRIMER BLOQUE DEL TRABAJO#######################
  
  tar_target(artificial_data, generate.artificial.data(grid_size)),
  tar_target(original_community, cbind(artificial_data$data$COORDS, artificial_data$data$X, artificial_data$data$Y, artificial_data$data$partition)),
  tar_target(p_target, p),
  tar_target(i_target, i),
  mapped,
  tar_combine(
    combined_explanatory_power,
    mapped[["explanatory_power"]]),
  tar_combine(
    combined_predictive_power,
    mapped[["predictive_power"]]),
  tar_combine(
    combined_same_sample_metrics,
    mapped[["same_sample_metrics"]]),
  tar_target(performance_plot, comparative.plot(combined_explanatory_power, combined_predictive_power, combined_same_sample_metrics)),
  tar_combine(
    combined_betas, mapped[["individual_betas"]]),
  tar_target(combined_rmas, RMAS(combined_betas, p_target, i_target), pattern = cross(p_target, i_target)),
  tar_target(betas_figure, betas.plot(combined_rmas)),
  tar_target(supl_combined_rmas, supl.RMAS(combined_betas, p_target), pattern = map(p_target)),
  tar_target(supl_betas_figure, supl.betas.plot(combined_betas, supl_combined_rmas)),


  ########################SEGUNDO BLOQUE DEL TRABAJO#############################


  tar_target(artificial_data_spcor, generate.artificial.data(grid_size, spCor = speciesCor)),
  tar_target(original_community_spcor, cbind(artificial_data_spcor$data$COORDS, artificial_data_spcor$data$X, artificial_data_spcor$data$Y,
                                             artificial_data_spcor$data$partition)),
  tar_target(plotOrder, corrMatOrder(artificial_data_spcor$param$spCor, order="AOE")),
  spcor_mapped,
  tar_combine(combined_false_negatives, spcor_mapped[["false_negatives"]]),
  tar_target(false_negative_ratio, combined_false_negatives %>% group_by(p) %>% summarise(mean=mean(FN), sd=sd(FN)) %>%
               write.csv("./tables/false_negative_ratio.csv")),
  tar_combine(to_plot_corr_all, spcor_mapped[["to_plot_corr"]]),
  tar_combine(combined_procrus, spcor_mapped[["procrus"]]),
  tar_target(procrus_rmse, combined_procrus %>% group_by(p) %>% summarise(mean=mean(rmse), sd=sd(rmse), mean_sign= mean(rmse_sign), sd_sign= sd(rmse_sign)) %>%
               write.csv("./tables/procrus_rmse.csv")),
  tar_target(to_plot_corr_groups, to_plot_corr_all %>% group_by(p_target, Var1, Var2) %>% tar_group(), iteration = "group"),
  tar_target(mean_corr, to_plot_corr_groups %>% summarize(distinct(to_plot_corr_groups, Var1) ,
                                                          distinct(to_plot_corr_groups, Var2) ,
                                                          value = mean(value),
                                                          distinct(to_plot_corr_groups, p_target)),
             pattern = map(to_plot_corr_groups)),
  tar_target(array_corr, acast(mean_corr, Var1 ~ Var2 ~ p_target)),
  tar_target(residuals_plot, individual.corr.plot(p_target, array_corr, plotOrder), pattern = map(p_target)),


  tar_target(full_modelo_hmsc_spcor, calibrate.hmsc(data=original_community_spcor, thin=thin, samples=samples,
  transient=transient, nChains=nchains, nParallel=nParallel)),
  tar_target(full_residuals_plot, correlation.plot(full_modelo_hmsc_spcor, "full", plotOrder )),
  tar_target(original_residuals_plot, original.correlation.plot(artificial_data_spcor, plotOrder)),
  tar_target(original_pattern, coo_pattern(original_community_spcor, sp_names, 1, NA, TRUE)),
  tar_combine(combined_sampled_coo, spcor_mapped[["coo_to_ggplot"]]),
  tar_target(to_plot_coo_groups, combined_sampled_coo %>% group_by(`sp-sp`, sample) %>% tar_group(), iteration = "group"),
  tar_target(mean_coo, to_plot_coo_groups %>% summarize(distinct(to_plot_coo_groups, var1),
                                                        distinct(to_plot_coo_groups, var2),
                                                        prob = mean(prob),
                                                        distinct(to_plot_coo_groups, `sp-sp`),
                                                        distinct(to_plot_coo_groups, sample)),
             pattern = map(to_plot_coo_groups)),
 tar_target(combined_pattern, rbind(original_pattern[,pattern_colnames], mean_coo)),
 tar_target(coo_levels, original_pattern$`sp-sp` [order(original_pattern$prob)]),
 tar_target(plot_pattern, pattern.plot(combined_pattern, coo_levels)),
 tar_combine(combined_corr_pattern, spcor_mapped[["corr_pattern"]]),
 tar_target(corr_groups, combined_corr_pattern %>% group_by(`sp-sp`, sample) %>% tar_group(), iteration = "group"),
 tar_target(sampled_corr, corr_groups %>% summarize(distinct(corr_groups, var1),
                                                       distinct(corr_groups, var2),
                                                       prob = mean(prob),
                                                       distinct(corr_groups, `sp-sp`),
                                                       distinct(corr_groups, sample)),
            pattern = map(corr_groups)),
 tar_target(full_corr_to_ggplot, corr.extraction(full_modelo_hmsc_spcor)),
 tar_target(full_corr_pattern,  coo_pattern(full_corr_to_ggplot, sample = 1, i = NA)),
 tar_target(original_corr, coo_pattern(artificial_data_spcor$param$spCor, sample = "original", i= NA)),
 tar_target(corr_plot, rbind(sampled_corr, full_corr_pattern[,pattern_colnames], original_corr[,pattern_colnames])),
 tar_target(corr_levels, original_corr$`sp-sp` [order(original_corr$prob)]),
 tar_target(plot_pattern_corr, corr.pattern.plot(corr_plot, corr_levels, "figures/omegas.jpg")),
 # tar_target(plot_pattern_corr_sign, corr.pattern.plot(corr_plot, "figures/omegas-sign.pdf", TRUE)), 
 tar_combine(combined_spat_data, spcor_mapped[["spat_data"]]),
 tar_target(full_spat_data, spatial.prediction(p = 1, i = 1, env_xy_data = original_community_spcor, model = full_modelo_hmsc_spcor, sp_names = sp_names)),
 tar_target(to_plot_spat_data, combined_spat_data %>% bind_rows(full_spat_data) %>% group_by(sample, variable, x, y) %>% 
              summarise(value = mean(value))),
 tar_target(spatial_plot, spatial.clust.plot(to_plot_spat_data))
)
  
