library(targets)
library(tarchetypes)
source("packages.R")
source("functions.R")
source("arguments.R")


values <- expand_grid(i_target = i, dataset_names = d_names)

mapped <- tar_map(values = values,
                  unlist = FALSE,
                  tar_target(input, prepare.data(dataset_names, gbif_data, complection_rast,
                                                 complection_data, 0.75, climate_data, soil_data,
                                                 TRUE, 500, i_target, partition)),
                  tar_target(model, calibrate.hmsc(input, thin, samples, transient,
                                                        nchains, nParallel, myspecies, var_names,
                                                        coord_names, formula)),
                  tar_target(diagnostics, plot.diagnostics(model, dataset_names, i_target)),
                  tar_target(betas, get.betas(model, myspecies, dataset_names, i_target)),
                  tar_target(omegas, get.omegas(model)),
                  tar_target(omegas_df, omegas_to_df(omegas, myspecies, dataset_names, i_target)),
                  tar_target(spatial_pred, spatial.prediction(climate_data, soil_data, myspecies, model,
                                                              rl = TRUE, i_target, dataset_names)),
                  tar_target(coords_pred, coords.prediction(community, myspecies, model, i_target, dataset_names)),
                  tar_target(coords_eval, eval.pred(coords_pred, community, myspecies, i_target, dataset_names))
)


# Replace the target list below with your own:
list(
  tar_target(gbif_data, get.gbif.data(code,
                                      countries = Forest_Europe_list,
                                      minlon = lon_min)),
  tar_target(indf_data, bdvis::format_bdvis(gbif_data, config = conf)),
  tar_target(mask, make.mask(indf_data, 0.041666)),
  tar_target(complection_data, completeness(indf_data, mask)),
  tar_target(complection_rast, rasterize_completeness(complection_data, "c", mask), format="qs"),
  tar_target(sprichness_rast, rasterize_completeness(complection_data, "Sobs", mask), format="qs"),
  tar_target(nrec_rast, rasterize(gbif_data[,c("decimalLongitude", "decimalLatitude")] %>% as.matrix, mask %>% unwrap, fun = "count") %>% wrap, format="qs"),
  tar_target(climate_data, get.climate.data(mask), format="qs"),
  tar_target(soil_data, get.soil.data(pH, sand, mask), format="qs"),
  tar_target(dataset_names, d_names),
  tar_target(partition, get.blocks(gbif_data)),
  tar_target(i_target, i),
  tar_target(community, prepare.data("full", gbif_data, complection_rast, complection_data, 0, climate_data, soil_data,
                                     cut = FALSE, 500, 1, partition)),
  mapped,
  tar_combine(
   betas_all, mapped[["betas"]]),
  tar_target(betas_plot, plot.betas(betas_all)),
  tar_combine(
  omegas_all, mapped[["omegas_df"]]),
  tar_target(omegas_plot, plot.omegas(omegas_all)),
  tar_combine(
    all_pred, mapped[["spatial_pred"]]),
  tar_combine(
    all_coords_eval, mapped[["coords_eval"]]),
  tar_target(coords_plot, plot.coords.eval(all_coords_eval)),
  tar_target(spatial_map, prediction.map(all_pred, all_coords_eval, myspecies, dataset_names)),
  tar_target(subset_map, subset.figure(spatial_map)),
  tar_target(composite_plot, composite(betas_plot, omegas_plot, spatial_map))
)
