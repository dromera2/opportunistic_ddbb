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
        tar_target(comunidad_muestreada, muestreo(comunidad_original, p_target, size= 500, i = i_target)),
        tar_target(comunidad_original_misma_muestra, original_same_sample(comunidad_muestreada, comunidad_original)),
        tar_target(comunidad_original_mismo_tamaño, original_same_size(comunidad_muestreada, comunidad_original)),
        tar_target(modelo_hmsc, calibrate_hmsc(muestreo=comunidad_muestreada, thin=thin, samples=samples,
                                                transient=transient, nChains=nchains, nParallel=nParallel)),
        tar_target(graficos_diagnosticos, plot.diagnostics(modelo_hmsc, p_target)),
        tar_target(poder_explicativo, tuning_metricas(com_testing= comunidad_muestreada, com_training= comunidad_muestreada,
                                                      modelo= modelo_hmsc, sp_names= sp_names, p_target, i= i_target)),
        tar_target(poder_predictivo, tuning_predictive(com_testing= comunidad_muestreada, modelo= modelo_hmsc, sp_names= sp_names, p_target, i= i_target)),
        tar_target(metricas_reales, tuning_metricas(com_testing= comunidad_original, com_training= comunidad_muestreada,
                                                    modelo= modelo_hmsc, sp_names= sp_names, p_target, i= i_target)),
        tar_target(metricas_reales_misma_muestra, tuning_metricas(com_testing = comunidad_original_misma_muestra, com_training= comunidad_muestreada,
                                                                  modelo= modelo_hmsc, sp_names= sp_names, p_target, i=i_target)),
        tar_target(metricas_reales_mismo_tamaño, tuning_metricas(com_testing = comunidad_original_mismo_tamaño, com_training = comunidad_muestreada,
                                                                 modelo= modelo_hmsc, sp_names= sp_names, p_target, i=i_target)),
        tar_target(individual_betas, BETAS(datos_artificiales = datos_artificiales, modelo = modelo_hmsc, p = p_target, i= i_target))
        
)

spcor_mapped <- tar_map(values = values,
                       unlist = FALSE,
                       tar_target(comunidad_muestreada_spcor, muestreo(comunidad_original_spcor, p_target, size= 500, i= i_target)),
                       tar_target(comunidad_original_misma_muestra_spcor, original_same_sample(comunidad_muestreada_spcor, comunidad_original_spcor)),
                       tar_target(falsos_negativos, compare.matrices(sp_names, comunidad_muestreada_spcor, comunidad_original_misma_muestra_spcor, i_target, p_target)),
                       tar_target(modelo_hmsc_spcor, calibrate_hmsc(muestreo=comunidad_muestreada_spcor, thin=thin, samples=samples,
                                                         transient=transient, nChains=nchains, nParallel=nParallel)),
                       tar_target(to_plot_corr, to_plot(modelo_hmsc_spcor, p_target)),
                       tar_target(coo_to_ggplot, coo_pattern(comunidad_muestreada_spcor, sp_names, p_target, i_target, TRUE)),
                       tar_target(corr_to_ggplot, corr_extraction(modelo_hmsc_spcor)),
                       tar_target(corr_pattern,  coo_pattern(corr_to_ggplot, sample = p_target, i = i_target)),
                       tar_target(procrus, procrustes.analysis(i_target, p_target, corr_pattern, datos_artificiales_spcor))
                       )

list(
  
  ###############################PRIMER BLOQUE DEL TRABAJO#######################
  
  tar_target(datos_artificiales, generate_artificial_data(GRID.SIZE)),
  tar_target(comunidad_original, cbind(datos_artificiales$data$COORDS, datos_artificiales$data$X, datos_artificiales$data$Y, datos_artificiales$data$partition)),
  tar_target(p_target, p),
  tar_target(i_target, i),
  mapped,
  tar_combine(
    combined_poder_explicativo,
    mapped[["poder_explicativo"]]),
  tar_combine(
    combined_poder_predictivo,
    mapped[["poder_predictivo"]]),
  tar_combine(
    combined_metricas_reales,
    mapped[["metricas_reales"]]),
  tar_combine(
    combined_metricas_reales_misma_muestra,
    mapped[["metricas_reales_misma_muestra"]]),
  tar_combine(
    combined_metricas_reales_mismo_tamaño,
    mapped[["metricas_reales_mismo_tamaño"]]),
  tar_target(performance_plot, comparative_plot(combined_poder_explicativo, combined_poder_predictivo, combined_metricas_reales, combined_metricas_reales_misma_muestra, combined_metricas_reales_mismo_tamaño)),
  tar_combine(
    combined_betas, mapped[["individual_betas"]]),
  tar_target(combined_rmas, RMAS(combined_betas, p_target, i_target), pattern = cross(p_target, i_target)),
  tar_target(betas_figure, betas_plot(combined_rmas)),
  tar_target(supl_combined_rmas, supl_RMAS(combined_betas, p_target), pattern = map(p_target)),
  tar_target(supl_betas_figure, supl_betas_plot(combined_betas, supl_combined_rmas)),
  
  
  ########################SEGUNDO BLOQUE DEL TRABAJO#############################
  
  
  tar_target(datos_artificiales_spcor, generate_artificial_data(GRID.SIZE, spCor = speciesCor)),
  tar_target(comunidad_original_spcor, cbind(datos_artificiales_spcor$data$COORDS, datos_artificiales_spcor$data$X, datos_artificiales_spcor$data$Y,
                                             datos_artificiales_spcor$data$partition)),
  tar_target(plotOrder, corrMatOrder(datos_artificiales_spcor$param$spCor, order="AOE")),
  spcor_mapped,
  tar_combine(combined_false_negatives, spcor_mapped[["falsos_negativos"]]),
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
  tar_target(residuals_plot, individual_corr_plot(p_target, array_corr, plotOrder), pattern = map(p_target)),


  tar_target(full_modelo_hmsc_spcor, calibrate_hmsc(muestreo=comunidad_original_spcor, thin=thin, samples=samples,
  transient=transient, nChains=nchains, nParallel=nParallel)),
  tar_target(full_residuals_plot, correlation_plot(full_modelo_hmsc_spcor, "full", plotOrder )),
  tar_target(original_residuals_plot, original_correlation_plot(datos_artificiales_spcor, plotOrder)), 
  tar_target(original_pattern, coo_pattern(comunidad_original_spcor, sp_names, 1, NA, TRUE)),
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
 tar_target(plot_pattern, pattern_plot(combined_pattern, coo_levels)), 
 tar_combine(combined_corr_pattern, spcor_mapped[["corr_pattern"]]),
 tar_target(corr_groups, combined_corr_pattern %>% group_by(`sp-sp`, sample) %>% tar_group(), iteration = "group"),
 tar_target(sampled_corr, corr_groups %>% summarize(distinct(corr_groups, var1),
                                                       distinct(corr_groups, var2),
                                                       prob = mean(prob),
                                                       distinct(corr_groups, `sp-sp`),
                                                       distinct(corr_groups, sample)),
            pattern = map(corr_groups)),
 tar_target(full_corr_to_ggplot, corr_extraction(full_modelo_hmsc_spcor)),
 tar_target(full_corr_pattern,  coo_pattern(full_corr_to_ggplot, sample = 1, i = NA)),
 tar_target(original_corr, coo_pattern(datos_artificiales_spcor$param$spCor, sample = "original", i= NA)),
 tar_target(corr_plot, rbind(sampled_corr, full_corr_pattern[,pattern_colnames], original_corr[,pattern_colnames])),
 tar_target(corr_levels, original_corr$`sp-sp` [order(original_corr$prob)]), 
 tar_target(plot_pattern_corr, corr_pattern_plot(corr_plot, corr_levels)),
 tar_target(plot_pattern_corr_sign, corr_pattern_plot(corr_plot, corr_levels, TRUE))
  )
  
