

# remove records without coordinates
clean.coords <- function(myspecies_coords) {
  myspecies_coords <- myspecies_coords %>%
    filter(!is.na(decimalLongitude)) %>%
    filter(!is.na(decimalLatitude)) %>% 
    mutate(countryCode = replace(countryCode, countryCode == "", "SK"))
  
  
  
  #convert country code from ISO2c to ISO3c
  myspecies_coords$countryCode <-  countrycode(myspecies_coords$countryCode, 
                                               origin =  'iso2c',
                                               destination = 'iso3c',
                                               custom_match = c("XK" = "SRB"))
  
  flags <- clean_coordinates(x = myspecies_coords, 
                             lon = "decimalLongitude", 
                             lat = "decimalLatitude",
                             # countries = "countryCode",
                             species = "species",
                             tests = c("capitals", "centroids", "equal","gbif",
                                       "institutions",
                                       # "outliers", 
                                       "seas", 
                                       "zeros", "urban", "validity", 
                                       "duplicates")) # most test are on by default
  
  #Exclude problematic records
  dat_cl <- myspecies_coords[flags$.summary,]
  dat_cl
}



#Remove records with low coordinate precision
low.precision <- function(dat_cl) {
  dat_cl <- dat_cl %>%
    filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))
  dat_cl
}


individual.count <- function(dat_cl) {
  dat_cl <- dat_cl%>%
    filter(individualCount > 0 | is.na(individualCount))%>%
    filter(individualCount < 99 | is.na(individualCount))
  dat_cl
}


record.age <- function(dat_cl) {
  dat_cl <- dat_cl%>%
    filter(year > 1945) # remove records from before second world war
  dat_cl
}

long_lat.cut <- function(dat_cl) {
  dat_cl <- dat_cl%>%
    filter (decimalLatitude > 20 & decimalLatitude < 50)%>%
    filter (decimalLongitude > -25 & decimalLongitude < 25)
  dat_cl
}


get.gbif.data <- function(code) {
 
  occurrences <- occ_download_get(code) %>%
    occ_download_import()
  
  filtered_coordinates <- clean.coords(occurrences)
  filtered_precision <- low.precision(filtered_coordinates)
  filtered_counts <- individual.count(filtered_precision)
  filtered_age <- record.age(filtered_counts)
  filtered_longlat <- long_lat.cut(filtered_age)
  
  
  # write.csv(filtered_longlat, "./data/gbif/clean_occurrences.csv")
  
  filtered_longlat
  
}



get.climate.data <- function(study_area) {
  clim <- worldclim_global(var = "bio", res = 0.5, path = "./data/climate")
  clim_crop <- crop(clim, extent(study_area), snap = "out")
  
  # terra::writeRaster(clim_crop, "./data/climate/cropped/clim_data.tif", overwrite = TRUE)
  
  pca <- prcomp(as.data.frame(clim_crop, na.rm=TRUE), center = TRUE, scale. = TRUE)
  
  pca_prediction <- predict(clim_crop, pca)
  
  pca_3axis <- c(pca_prediction$PC1, pca_prediction$PC2, pca_prediction$PC3)
  
  # terra::writeRaster(pca_3axis, "./data/climate/pca3_data.tif", overwrite = TRUE)
  # terra::writeRaster(pca_3axis %>% subset("PC1"), "./data/climate/pca1_data.tif", overwrite = TRUE)
  
  pca_3axis %>% terra::wrap()
  
}


get.soil.data <- function(pH_list, sand_list, template) {
  
  pH <- rast(pH_list) %>% subset(c("phh2o_60-100cm_mean_1000",
                                   "phh2o_100-200cm_mean_1000"),
                                 negate = TRUE)
  
  sand <- rast(sand_list) %>% subset(c("sand_60-100cm_mean_1000",
                                       "sand_100-200cm_mean_1000"),
                                     negate = TRUE)
  
  pH <- pH %>% project(template, method = "bilinear")  %>% 
    # aggregate(fact = 5, fun = "mean", na.rm = TRUE) %>% 
    app(mean)
  
  names(pH) <- "pH"
  
  sand <- sand %>% project(template, method = "bilinear")  %>% 
    # aggregate(fact = 5, fun = "mean", na.rm = TRUE) %>%
    app(mean)
  
  names(sand) <- "sand"
  
  soil_data <- c(pH, sand)
  
  terra::writeRaster(soil_data, "./data/soil/soil_data.tif", overwrite = TRUE)
  
  soil_data %>% terra::wrap()
}

prepare.data <- function(occ_data, clim_data, soil_data) {
 points <- occ_data %>% dplyr::select(c(species, decimalLongitude, decimalLatitude))
  
  points <- points  %>%  na.omit() %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"))
  
  clim_data <- clim_data %>% terra::unwrap()
  
  template <- clim_data %>% subset(1)
  
  points <- points %>% mutate(pixelId = cellFromXY(template, points %>% st_coordinates()), occurrence = 1) %>% 
    dplyr::select(species, pixelId, occurrence) %>% st_drop_geometry() %>% 
    pivot_wider(names_from = "species", values_from = "occurrence", values_fn = max, values_fill = 0)
  
  
  coords <- template %>% xyFromCell(cell = points$pixelId) %>% as.data.frame()
  
  
  clim_data <- clim_data %>%  terra::extract(y = points$pixelId)
  
  soil_data <- soil_data %>% unwrap() %>% terra::extract(y = points$pixelId)
  
  model_data <- bind_cols(coords, clim_data, soil_data, points) %>% dplyr::select(-pixelId) %>% na.omit()
  
  # write.csv(model_data, "./data/model_data.csv")
  
  model_data <- model_data %>% slice_sample(n = 5000)
  model_data
}


calibrate.hmsc <- function(data, thin, samples, transient, nChains, nParallel, sp_names, var_names, coord_names, formula){
  
  set.model <- function(data, sp_names, var_names, coord_names, formula) {
    Y <- data[, sp_names] 
    XData <- data[, var_names]
    xy <- data[, coord_names]
    Y <- as.matrix(Y)
    rownames(Y) <- rownames(XData)
    XFormula <- formula
    studyDesign <- data.frame(units = rownames(Y))
    studyDesign$units <- factor(studyDesign$units, levels = rownames(Y))
    rl <- HmscRandomLevel(sData = xy, sMethod = 'NNGP', nNeighbours = 10, longlat = TRUE)
    rl = setPriors(rl, nfMin = 1, nfMax = 2)
    Hmsc(Y = Y,
         XData = XData,
         XFormula = XFormula,
         studyDesign = studyDesign,
         ranLevels = list("units"=rl),
         distr = "probit")
  }
  
  model.pa<- set.model(data, sp_names, var_names, coord_names, formula)
  
  # fitting the model
  sampleMcmc(model.pa, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = nParallel)
}

calibrate.hmsc.no.rl <- function(data, thin, samples, transient, nChains, nParallel, sp_names, var_names, coord_names, formula){
  
  set.model <- function(data, sp_names, var_names, coord_names, formula) {
    Y <- data[, sp_names] 
    XData <- data[, var_names]
    xy <- data[, coord_names]
    Y <- as.matrix(Y)
    rownames(Y) <- rownames(XData)
    XFormula <- formula
    # studyDesign <- data.frame(units = rownames(Y))
    # studyDesign$units <- factor(studyDesign$units, levels = rownames(Y))
    # rl <- HmscRandomLevel(sData = xy, sMethod = 'NNGP', nNeighbours = 10, longlat = TRUE)
    # rl = setPriors(rl, nfMin = 1, nfMax = 2)
    Hmsc(Y = Y,
         XData = XData,
         XFormula = XFormula,
         # studyDesign = studyDesign,
         # ranLevels = list("units"=rl),
         distr = "probit")
  }
  
  model.pa<- set.model(data, sp_names, var_names, coord_names, formula)
  
  # fitting the model
  sampleMcmc(model.pa, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = nParallel)
}


plot.diagnostics <- function(model) {
  mpost = convertToCodaObject(model)
  
  ess.beta = effectiveSize(mpost$Beta)
  
  ess.V = effectiveSize(mpost$V)
  
  psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  
  psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
  
  par(mfrow = c(2, 2))
  pdf("figures/diagnostics.pdf")
  hist(ess.beta, xlab = expression("Effective sample size" ~ beta ~ ""), main=NULL, col= "light green")
  hist(ess.V, xlab = expression("Effective sample size" ~ v ~ ""))
  hist(psrf.beta, xlab = expression("Potential scale reduction factor" ~ beta ~ ""), main=NULL, col= "light green")
  hist(psrf.V, xlab = expression("Potential scale reduction factor" ~ v ~ ""))
  dev.off()
  
}


omega.matrix <- function(model, plot_order) {
  OmegaCor = computeAssociations(model)
  supportLevel = 0.95
  
  OmegaPValue <- ((OmegaCor[[1]]$support>supportLevel) + (OmegaCor[[1]]$support<(1-supportLevel))>0)
  
  # plotOrder = corrMatOrder(OmegaCor[[1]]$mean,order="AOE")
  
  # plotOrder <- plot_order
  
  plotOrder <- c("Pinus pinaster", "Pinus nigra", "Pinus sylvestris",
                 "Pinus halepensis", "Pinus pinea", "Pinus mugo",
                 "Quercus petraea", "Quercus pubescens", "Quercus pyrenaica",
                 "Quercus suber", "Quercus faginea", "Quercus ilex", "Quercus robur")
  
  pdf("./figures/residuals.pdf", width = 7, height = 7)

  corrplot(OmegaCor[[1]]$mean[plotOrder,plotOrder],
           method = "color",
           p.mat = OmegaPValue[plotOrder,plotOrder],
           col = rev(COL2('RdBu', 200)),
           sig.level = 0.95,
           insig = "pch",
           pch = 0,
           pch.cex = 3.3)

  dev.off()
  
  ggcorrplot(OmegaCor[[1]]$mean[plotOrder,plotOrder], p.mat = (OmegaPValue[plotOrder,plotOrder]*1),
             colors = c("#009E73", "white", "#CC79A7"))
  
  ggsave("./figures/residuals_2.pdf", width = 7, height = 7)
  
  OmegaCor
}

spatial.prediction <- function(climate, soil, species, model, rl = NULL) {

  env_data <- c(climate, soil)
    
    XData.grid <- env_data %>% values(dataframe = TRUE)
    xy.grid <- env_data %>% xyFromCell(cell = 1:ncell(.)) %>% as.data.frame()
    
    complete.index <- which(complete.cases(XData.grid))
    Gradient <- prepareGradient(model, XDataNew = XData.grid[complete.index,], sDataNew = list(units=xy.grid[complete.index,]))
    
    if(!is.null(rl)){
      row.names(Gradient$rLNew$units$s) <- as.character(Gradient$rLNew$units$pi)
    }
    
    predY <- predict(model, Gradient = Gradient, expected = TRUE, predictEtaMean = TRUE, nParallel=2)
    
    EpredY <- matrix(NA, nrow(XData.grid), 13)
 
    library(abind)
    EpredY[complete.index,] <- apply(abind(predY,along=3),c(1,2),mean)
   
    mapData <- data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY)
    colnames(mapData) <- c("x", "y", species)
    library(reshape2)
    mapData <- melt(mapData, id.vars = c("x", "y"))
    
    mapData
  }

combine.predictions <- function(model_data, sp_names, var_names, model,
                                model_no_rl, pred, pred_no_rl) {
  thres <- NULL
  
  for(j in 1:length(sp_names)){
    print(j)
    prediction <- predict(model, XData=model_data %>% dplyr::select(var_names), expected = TRUE, nParallel=2)
    prediction <- (Reduce("+", prediction) / length(prediction)) %>% as.data.frame()
    
    eval <- dismo::evaluate(prediction %>% filter(get(sp_names[j], pos = model_data) == 1) %>% pull(sp_names[j]),
                            prediction %>% filter(get(sp_names[j], pos = model_data) == 0) %>% pull(sp_names[j]))
    if(j == 1){
      thres <- dismo::threshold(eval)
    } else{
      thres[j,] <- dismo::threshold(eval)
    }
  }
  
  thres$sp <- sp_names
  
  thres_no_rl <- NULL
  
  for(j in 1:length(sp_names)){
    print(j)
    prediction <- predict(model_no_rl, XData=model_data %>% dplyr::select(var_names), expected = TRUE, nParallel=2)
    prediction <- (Reduce("+", prediction) / length(prediction)) %>% as.data.frame()
    
    eval <- dismo::evaluate(prediction %>% filter(get(sp_names[j], pos = model_data) == 1) %>% pull(sp_names[j]),
                            prediction %>% filter(get(sp_names[j], pos = model_data) == 0) %>% pull(sp_names[j]))
    if(j == 1){
      thres_no_rl <- dismo::threshold(eval)
    } else{
      thres_no_rl[j,] <- dismo::threshold(eval)
    }
  }
  
  thres_no_rl$sp <- sp_names
  
  
  thr <- thres %>% pull(spec_sens)
  thr_no_rl <- thres_no_rl %>% pull(spec_sens)
  
  
  sp_list <- list()
  
  for(i in 1:length(sp_names)){
    
    sp_list[[i]] <- pred %>% filter(variable == sp_names[i]) %>% mutate(value = case_when(
      value < thr[i] ~ 0,
      value >= thr[i] ~ 1
    ))
  }
  
  sp_list_no_rl <- list()
  
  for(i in 1:length(sp_names)){
    sp_list_no_rl[[i]] <- pred_no_rl %>% filter(variable == sp_names[i]) %>% mutate(value = case_when(
      value < thr_no_rl[i] ~ 0,
      value >= thr_no_rl[i] ~ 1
    ))
  }
  
  recla_pred <- sp_list %>% bind_rows()
  recla_pred$value <- (recla_pred$value)*10
  
  recla_pred_no_rl <- sp_list_no_rl %>% bind_rows()
  
  combined_pred <- recla_pred %>% mutate(value = recla_pred$value + recla_pred_no_rl$value)
}

prediction.map <- function(comb_pred) {
  comb_pred %>% 
    ggplot(aes(x=x, y=y, fill=as.factor(value))) +
    geom_tile() +
    # scale_fill_manual(values = c("light green", "yellow", "orange", "red" ), na.value="royalblue") +
    scale_fill_manual(values = c("#D2D1CE", "#D8B404", "#3A93EF", "#2F521C" ), na.value="#FFFFFF") +
    theme(legend.position="bottom", aspect.ratio= 1) +
    facet_wrap(~variable) +
    theme_bw() +
    theme(legend.title=element_blank()) 
  
  ggsave("./figures/combined_spatial_pred.pdf", width = 10, height = 7)
}
