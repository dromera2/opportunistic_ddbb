library(targets)
source("packages.R")
source("functions.R")
source("arguments.R")

# Replace the target list below with your own:
list(
  tar_target(gbif_data, get.gbif.data(code)),
  tar_target(climate_data, get.climate.data(study_area), format="qs"),
  tar_target(soil_data, get.soil.data(pH, sand, temp), format="qs"),
  tar_target(model_input, prepare.data(gbif_data, climate_data, soil_data)),
  tar_target(hmsc_model, calibrate.hmsc(model_input, thin, samples, transient,
                                        nchains, nParallel, myspecies, var_names,
                                        coord_names, formula)),
  tar_target(diagnostics, plot.diagnostics(hmsc_model)),
  tar_target(residuals, omega.matrix(hmsc_model, myspecies)), 
  tar_target(prediction, spatial.prediction(pred_clim, pred_soil, myspecies, hmsc_model, rl = TRUE)),
  tar_target(hmsc_model_no_rl, calibrate.hmsc.no.rl(model_input, thin, samples, transient,
                                        nchains, nParallel, myspecies, var_names,
                                        coord_names, formula)),
  tar_target(prediction_no_rl, spatial.prediction(pred_clim, pred_soil, myspecies, hmsc_model_no_rl)),
  tar_target(combined_pred, combine.predictions(model_input, myspecies, var_names, hmsc_model,
                                                            hmsc_model_no_rl, prediction, prediction_no_rl)),
  tar_target(pred_map, prediction.map(combined_pred))
  )
