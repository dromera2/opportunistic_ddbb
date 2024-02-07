

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

long_lat.cut <- function(dat_cl, minlon) {
  dat_cl <- dat_cl%>%
    filter (decimalLongitude > minlon)
  dat_cl
}

filter.country <- function(dat_cl, countries){
  dat_cl <- dat_cl %>% 
    filter(countryCode %in% countries)
  dat_cl
}

get.gbif.data <- function(code, countries, minlon) {
  
  occurrences <- occ_download_get(code) %>%
    occ_download_import()
  
  filtered_coordinates <- clean.coords(occurrences)
  filtered_precision <- low.precision(filtered_coordinates)
  filtered_counts <- individual.count(filtered_precision)
  filtered_age <- record.age(filtered_counts)
  filtered_countries <- filter.country(filtered_age, countries)
  filtered_longlat <- long_lat.cut(filtered_countries, minlon)
  
  filtered_longlat
  
}

make.mask <- function(indf, resolution) {
  r <- terra::rast()
  ext(r) <- c(min(indf$Longitude) %>% floor(),
                   max(indf$Longitude) %>% ceiling(),
                   min(indf$Latitude) %>% floor(),
                   max(indf$Latitude) %>% ceiling())
  res(r) <- resolution
  crs(r) <- "epsg:4326"
  
  terra::wrap(r)
}
########blockcv

get.blocks <- function(data, n = 5, i = 100, size = 500000) {
  data <- data %>% select(decimalLongitude, decimalLatitude) 
  data$occ <- 1
  
  sf_data <- st_as_sf(x = data, 
                      coords = c("decimalLongitude", "decimalLatitude"),
                      crs = "epsg:4326")
  
  #Spatial blocking by specified range with random assignment
  set.seed(1)
  sb <- cv_spatial(x = sf_data,
                     column = "occ",
                     size = size,
                     k = n,
                     selection = "random",
                     iteration = i, # find evenly dispersed folds
  )
  blocks <- data.frame(partition = sb$folds_ids)
  
  blocks
}

##############

completeness <- function (indf, raster, recs = 10) {
  require(tidyverse)
  require(terra)  
  
  raster <- terra::unwrap(raster)
  indf <- as.data.frame(indf)
  
  coords <- indf %>% 
    select(Longitude, Latitude) 
  
  indf$cust_cell_id <- cellFromXY(raster, coords)
  
  dat1 <- indf %>% # Remove duplicated entries
    select(Scientific_name, Date_collected, cust_cell_id) %>% 
    distinct(Scientific_name, Date_collected, cust_cell_id)
  dat2 <- dat1 %>% group_by(cust_cell_id) %>% # Number of records in each gridcell
    summarize(cell_ct = n())
  dat3 <- dat2 %>% # Filter by minimum  number of records
    filter(cell_ct >= recs) 
  if (nrow(dat3) < 1) {
    stop("Too few data records to compute completeness")
  }
  retmat <- NULL
  for (i in 1:nrow(dat3)) {
    print(paste0(i, " / ", nrow(dat3)))
    cust_cell_id <- dat3$cust_cell_id[i] # The cell ID
    nrec <- dat3$cell_ct[i] # Number of records in the cell
    cset <- dat1[which(dat1$cust_cell_id == cust_cell_id), ] # Comm. matrix of the cell
    csum <- cset %>% group_by(Scientific_name) %>% summarize(sct = n()) # Number of records for each species in the cell
    Q1 <- csum %>% filter(sct == 1) %>% nrow() # Number of unique species
    Q2 <- csum %>% filter(sct == 2) %>% nrow() # Number of duplicate species
    m <- cset %>% distinct(Date_collected) %>% nrow() # Number of samples (times sampled the grid) 
    Sobs <- cset %>% distinct(Scientific_name) %>% nrow()
    if (Sobs > 0) {
      Sest <- Sobs + (((m - 1)/m) * ((Q1 * (Q1 - 1))/(2 * (Q2 + 1))))
      c <- Sobs/Sest
      retmat <- rbind(retmat, c(cust_cell_id, nrec, Q1, Q2, m, Sobs, Sest, c))
    }
  }
  retmat <- as.data.frame(retmat)
  names(retmat) <- c("Cell_id", "nrec", "Q1", "Q2", "m", "Sobs", "Sest", "c")
  
  retmat
}

rasterize_completeness <- function(comp_an, var, raster) {
  
  raster <- raster %>% terra::unwrap()
  
  new_values <- rep(NA, ncell(raster))
  new_values[comp_an$Cell_id] <- comp_an[, var]
  values(raster) <- new_values
  raster %>% terra::wrap()
}

sample.data <- function(complection_rast, gbif_data, complection_data,  threshold) {
  raster <- complection_rast %>% terra::unwrap()
  
  filtered_data <- gbif_data %>% as.data.frame()
  
  coords <- filtered_data %>% select(decimalLongitude, decimalLatitude)
  
  filtered_data$Cell_id <- cellFromXY(raster, coords)
  
    comp_data <- complection_data %>%
      filter(c >= threshold) %>% 
      select(Cell_id) 
  
  filtered_data <- filtered_data %>%
    filter(Cell_id %in% comp_data$Cell_id) 
  
  filtered_data
}

under.threshold <- function(raw_data, filtered_data){
  under_data <- raw_data %>% 
    filter(!(gbifID %in% filtered_data$gbifID))
  
  under_data
}

get.climate.data <- function(study_area) {
  
  study_area <- study_area %>% terra::unwrap()
  
  clim <- worldclim_global(var = "bio", res = 2.5, path = "./data/climate")
  clim_crop <- crop(clim, ext(study_area), snap = "out")
  
  pca <- prcomp(as.data.frame(clim_crop, na.rm=TRUE), center = TRUE, scale. = TRUE)
  
  pca_prediction <- predict(clim_crop, pca)
  
  pca_3axis <- c(pca_prediction$PC1, pca_prediction$PC2, pca_prediction$PC3)
  
  pca_3axis %>% terra::wrap()
  
}


get.soil.data <- function(pH_list, sand_list, template) {
  
  pH <- rast(pH_list) %>% subset(c("phh2o_60-100cm_mean_1000",
                                   "phh2o_100-200cm_mean_1000"),
                                 negate = TRUE)
  
  sand <- rast(sand_list) %>% subset(c("sand_60-100cm_mean_1000",
                                       "sand_100-200cm_mean_1000"),
                                     negate = TRUE)
  
  template <- template %>% 
    terra::unwrap() %>% 
    disagg(fact = 5, method = "bilinear")
  
  pH <- pH %>% project(template, method = "bilinear")  %>% 
    aggregate(fact = 5, fun = "mean", na.rm = TRUE) %>%
    app(mean)
  
  names(pH) <- "pH"
  
  sand <- sand %>% project(template, method = "bilinear")  %>% 
    aggregate(fact = 5, fun = "mean", na.rm = TRUE) %>%
    app(mean)
  
  names(sand) <- "sand"
  
  soil_data <- c(pH, sand)
  
  soil_data %>% terra::wrap()
}

filter.data <- function(type, data, comp_rast, comp_data, th){
  if(type %in% c("over", "under")){
    source_data <- data
    data <- sample.data(comp_rast, data, comp_data, th)
    if (type == "under"){
      data <- under.threshold(source_data, data)
    }
  }
  return(data)
}

prepare.data <- function(type, occ_data, comp_rast, comp_data, th, clim_data, soil_data, cut, sample = NULL, i, blocks) {
  
  occ_data$partition <- blocks$partition
  
  occ_data <- filter.data(type, occ_data, comp_rast, comp_data, th)

  points <- occ_data %>% dplyr::select(c(species, decimalLongitude, decimalLatitude, partition))
  
  points <- points  %>%  na.omit() %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"))
  
  clim_data <- clim_data %>% terra::unwrap()
  
  template <- clim_data %>% subset(1)
  
  points <- points %>% mutate(pixelId = cellFromXY(template, points %>% st_coordinates()), occurrence = 1) %>% 
    dplyr::select(species, pixelId, occurrence, partition) %>% st_drop_geometry() %>% 
    pivot_wider(names_from = "species", values_from = "occurrence", values_fn = max, values_fill = 0) %>% 
    distinct(pixelId, .keep_all = TRUE)

  coords <- template %>% xyFromCell(cell = points$pixelId) %>% as.data.frame()
  
  clim_data <- clim_data %>%  terra::extract(y = points$pixelId)
  
  soil_data <- soil_data %>% unwrap() %>% terra::extract(y = points$pixelId)
  
  model_data <- bind_cols(coords, clim_data, soil_data, points) %>% dplyr::select(-pixelId) %>% na.omit()
  
    if(cut == TRUE){
    set.seed(i)  
    model_data <- model_data %>% slice_sample(n = sample) 
  }
  
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



plot.diagnostics <- function(model, name, i) {
  mpost = convertToCodaObject(model)
  
  ess.beta = effectiveSize(mpost$Beta)
  
  ess.V = effectiveSize(mpost$V)
  
  psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  
  psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
  
  par(mfrow = c(2, 2))
  pdf(paste0("figures/diagnostics/diagnostics_", name, i, ".pdf"))
  hist(ess.beta, xlab = expression("Effective sample size" ~ beta ~ ""), main=NULL, col= "light green")
  hist(ess.V, xlab = expression("Effective sample size" ~ v ~ ""))
  hist(psrf.beta, xlab = expression("Potential scale reduction factor" ~ beta ~ ""), main=NULL, col= "light green")
  hist(psrf.V, xlab = expression("Potential scale reduction factor" ~ v ~ ""))
  dev.off()
  
}

spatial.prediction <- function(climate, soil, species, model, rl = NULL, i, dataset) {

  mask <- ext(c(-11, 10, 42.5, 54))
  climate <- climate %>% terra::unwrap() %>% crop(mask)
  soil <- soil %>% terra::unwrap() %>% crop(mask)
  
  env_data <- c(climate, soil)
  
  XData.grid <- env_data %>% terra::values(dataframe = TRUE)
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
  
  mapData$i <- i
  mapData$dataset <- dataset
  
  mapData
 
}

coords.prediction <- function(comm, species, model, i, dataset) {
  
  XData.grid <- comm %>% select(PC1, PC2, PC3, pH, sand)
  xy.grid <- comm %>% select(x, y)
  
  complete.index <- which(complete.cases(XData.grid))
  Gradient <- prepareGradient(model, XDataNew = XData.grid[complete.index,], sDataNew = list(units=xy.grid[complete.index,]))
  
  row.names(Gradient$rLNew$units$s) <- as.character(Gradient$rLNew$units$pi)
  
  predY <- predict(model, Gradient = Gradient, expected = TRUE, predictEtaMean = TRUE, nParallel=2)
  
  EpredY <- matrix(NA, nrow(XData.grid), 13)
  
  library(abind)
  EpredY[complete.index,] <- apply(abind(predY,along=3),c(1,2),mean)
  
  mapData <- data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY)
  colnames(mapData) <- c("x", "y", species)
  library(reshape2)
  mapData <- melt(mapData, id.vars = c("x", "y"))
  
  mapData$i <- i
  mapData$dataset <- dataset
  
  mapData
  
}


eval.pred <- function(pred, obse, sp_names, i, dataset) {

  pred <- pred %>% 
    pivot_wider(names_from = variable,
                values_from = value)
  
  data <- obse %>% left_join(pred, by = c("x", "y"))
  
  thres <- data.frame(kappa = NA,
                      spec_sens = NA,
                      no_omission = NA,
                      prevalence = NA,
                      equal_sens_spec = NA,
                      sensitivity = NA,
                      AUC = NA,
                      KAPPA = NA,
                      TSS = NA)
  
  for(j in 1:length(sp_names)){
    print(j)
    
    eval <- dismo::evaluate(data %>% filter(get(paste0(sp_names[j], ".x"), pos = data) == 1) %>% pull(paste0(sp_names[j], ".y")),
                            data %>% filter(get(paste0(sp_names[j], ".x"), pos = data) == 0) %>% pull(paste0(sp_names[j], ".y")))
    thres[j,] <- dismo::threshold(eval)
    thres[j, "AUC"] <- eval@auc
    thres[j, "KAPPA"] <- eval@kappa[first(which(eval@t >= thres$spec_sens))]
    thres[j, "TSS"] <- eval@TPR[first(which(eval@t >= thres$spec_sens))] + eval@TNR[first(which(eval@t >= thres$spec_sens))] - 1
  }
  
  thres$sp <- sp_names
  
  thres$i <- i
  thres$dataset <- dataset
  
  thres
  
}


prediction.map <- function(pred, th, sp_names, datasets) {
  
  th <- th %>%
    select(spec_sens, i, dataset, sp) %>% 
    group_by(dataset, sp) %>% 
    summarize(spec_sens = mean(spec_sens, na.rm = TRUE))
  
   pred <- pred %>% group_by(x, y, variable, dataset) %>% summarize(value = mean(value, na.rm = TRUE))
  
  results <- list()
  
  for(i in 1:length(sp_names)){
  
    temp_list <- list()
    
    for(j in 1:length(datasets)){
     
    th_value <- th %>%
      filter(sp == sp_names[i] &
      dataset == datasets[j]) %>% 
      pull(spec_sens)
    
    temp_list[[j]] <- pred %>% 
      filter(variable == sp_names[i] &
               dataset == datasets[j]) %>% 
      mutate(th_pred = as.numeric(value >= th_value))
    }
    
    names(temp_list) <- datasets
    
    results[[i]] <- temp_list
  }
  
names(results) <- sp_names  
  
recla_pred <- reshape2::melt(results, id.vars = c("x", "y", "variable", "dataset", "value", "th_pred"))

  recla_pred %>% 
    ggplot(aes(x=x, y=y, fill= th_pred, col = th_pred)) +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_colour_viridis_c() +
    theme(legend.position="bottom", aspect.ratio= 1) +
    facet_grid(variable ~ factor(dataset, levels = c("under", "full", "over"))) +
    theme_bw() +
    theme(legend.title=element_blank()) 
 
  
  ggsave("./figures/spatial_pred.jpg", width = 6, height = 18)
  
  recla_pred
}

get.betas <- function (model, species, type, i){
  postBeta <- getPostEstimate(model, parName="Beta")
  beta.mod <- data.frame(t(postBeta$mean))
  colnames(beta.mod) <- c("int", "PC1_1", "PC2_1", "PC3_1", "pH_1", "sand_1",
                          "PC1_2", "PC2_2", "PC3_2", "pH_2", "sand_2")
  
  beta.mod <- melt(beta.mod, measure.vars= c("int", "PC1_1", "PC2_1", "PC3_1", "pH_1", "sand_1",
                                             "PC1_2", "PC2_2", "PC3_2", "pH_2", "sand_2"))
  beta.mod$species <- species
  beta.mod$type <- type
  beta.mod$i <- i
  colnames(beta.mod) <- c("variable", "beta", "species", "type", "i")
  return(beta.mod)
}

get.omegas <- function(model){
  
  OmegaCor = computeAssociations(model)
  
  plotOrder <- c("Pinus pinaster", "Pinus nigra", "Pinus sylvestris",
                 "Pinus halepensis", "Pinus pinea", "Pinus mugo",
                 "Quercus petraea", "Quercus pubescens", "Quercus pyrenaica",
                 "Quercus suber", "Quercus faginea", "Quercus ilex", "Quercus robur")
  
  omega <- OmegaCor[[1]]$mean[plotOrder,plotOrder]
} 

omegas_to_df <- function(matrix, sp_names, type, i) {

  plotOrder <- c("Pinus pinaster", "Pinus nigra", "Pinus sylvestris",
                 "Pinus halepensis", "Pinus pinea", "Pinus mugo",
                 "Quercus petraea", "Quercus pubescens", "Quercus pyrenaica",
                 "Quercus suber", "Quercus faginea", "Quercus ilex", "Quercus robur")
  prob <- matrix[plotOrder, plotOrder]
  prob_matrix <- as.matrix(prob)
  prob_matrix[upper.tri(prob_matrix, diag = TRUE)] <- NA
  prob_matrix <- na.omit(melt(prob_matrix))
  prob_matrix <- cbind(prob_matrix, paste0(prob_matrix[,1], "-", prob_matrix[,2]))
  prob_matrix[,5] <- type
  prob_matrix[,"i"] <- i
  colnames(prob_matrix) <- c("var1", "var2", "prob", "sp-sp", "type", "i")
  prob_matrix
}

plot.omegas <- function(corr_plot) {
  
  corr_plot <- corr_plot %>% group_by(type, var1, var2) %>% summarize(prob = mean(prob))
  
  corr_plot$type <- corr_plot$type %>% str_to_sentence() %>% factor(levels = c("Under", "Full", "Over"))
  plotOrder <- c("P. pinaster", "P. nigra", "P. sylvestris",
                 "P. halepensis", "P. pinea", "P. mugo",
                 "Q. petraea", "Q. pubescens", "Q. pyrenaica",
                 "Q. suber", "Q. faginea", "Q. ilex", "Q. robur")
  
  corr_plot$var1 <- corr_plot$var1 %>% str_replace("Quercus", "Q.") %>% str_replace("Pinus", "P.") %>% factor(levels = plotOrder)
  corr_plot$var2 <- corr_plot$var2 %>% str_replace("Quercus", "Q.") %>% str_replace("Pinus", "P.") %>% factor(levels = rev(plotOrder))
  ggplot(corr_plot) +
    geom_raster(aes(x = var1, y = var2, fill = prob)) +
    scale_fill_gradientn(colors= c("#7f3b08", "#e08214", "#f5f5f5", "#8073ac", "#2d004b"), values = c(0, 0.4, 0.5, 0.6, 1), limits= c(-1, 1)) +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(legend.title= element_blank(), axis.text.x=element_text(angle = 45, hjust = 1)) +
    facet_wrap(.~type)
  filename <- "./figures/omegas.pdf"
  ggplot2::ggsave(filename, device = "pdf", width = 9.5, height = 3.5)
  
  corr_plot
}

plot.betas <- function(data){
  data <- data %>% 
    group_by(species, variable, type) %>% 
    summarize(beta = mean(beta)) %>% 
    filter(variable != "int")
  
  data %>% ggplot() +
    geom_tile(aes(x = factor(type, levels = c("under", "full", "over")), y = species, fill = beta)) +
    scale_fill_gradientn(colors= c("#8e0152", "#c51b7d", "#f5f5f5", "#4d9221", "#276419"), values = c(0, 0.4, 0.5, 0.6, 1), limits= c(-1, 1)) +
    facet_wrap(.~factor(variable, levels = c("PC1_1", "PC1_2", "PC2_1", "PC2_2", "PC3_1", "PC3_2",
                                             "sand_1", "sand_2", "pH_1", "pH_2"))) +
    xlab("") +
    theme_bw()
  
  ggsave("./figures/betas.pdf", width = 9, height = 6)
  
  data %>% 
    filter(variable %in% c("PC1_1", "PC1_2", "PC2_1", "PC2_2", "sand_1", "pH_1")) %>% 
    ggplot() +
    geom_tile(aes(x = factor(type, levels = c("under", "full", "over")), y = species, fill = beta)) +
    scale_fill_gradientn(colors= c("#8e0152", "#c51b7d", "#f5f5f5", "#4d9221", "#276419"), values = c(0, 0.4, 0.5, 0.6, 1), limits= c(-1, 1)) +
    facet_wrap(.~factor(variable, levels = c("PC1_1", "PC1_2", "PC2_1", "PC2_2", "sand_1", "pH_1"))) +
    xlab("") +
    theme_bw()
  
  ggsave("./figures/subset_betas.pdf", width = 9, height = 6)
  
  return(data)
}

##########

plot.coords.eval <- function(eval) {
  
  eval <- eval %>% select(AUC, KAPPA, TSS, sp, i, dataset) %>% 
    pivot_longer(cols = c(AUC, KAPPA, TSS),
                 names_to = "variable",
                 values_to = "value") 
  
  eval %>% 
    ggplot() +
    geom_boxplot(aes(x = factor(dataset, levels = c("under", "full", "over")), y = value, fill = dataset)) +
    xlab("") +
    theme_bw() +
    facet_wrap(~variable) +
    theme(legend.position = "none")
  
  ggsave("./figures/comm_eval.pdf", device = "pdf", width = 7, height = 2.5)
  
  eval
}
  
subset.figure <- function(map_data) {
  map_data %>% filter(variable == "Quercus suber" |
                        variable == "Pinus sylvestris" |
                        variable == "Pinus pinea") %>% 
    ggplot(aes(x = x, y = y, fill = th_pred, col = th_pred)) +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_colour_viridis_c() +
    theme(legend.position="bottom", aspect.ratio= 1) +
    facet_grid(factor(dataset, levels = c("under", "full", "over")) ~ variable) +
    theme_bw() +
    theme(legend.title=element_blank()) 
  
  
  ggsave("./figures/subset_spatial.jpg", width = 10, height = 7)
}


composite <- function(betas, omegas, spatial){
  
  betas$species <- betas$species %>%  str_replace("Quercus", "Q.") %>% str_replace("Pinus", "P.")
  
  betas_plot <- betas %>% 
    filter(variable %in% c("PC1_1", "PC1_2", "PC2_1", "PC2_2", "sand_1", "pH_1")) %>% 
    ggplot() +
    geom_tile(aes(x = factor(type, levels = c("under", "full", "over")), y = species, fill = beta)) +
    scale_fill_gradientn(colors= c("#8e0152", "#c51b7d", "#f5f5f5", "#4d9221", "#276419"), values = c(0, 0.4, 0.5, 0.6, 1), limits= c(-1, 1)) +
    facet_wrap(.~factor(variable, levels = c("PC1_1", "PC1_2", "PC2_1", "PC2_2", "sand_1", "pH_1"))) +
    xlab("") +
    ylab("") +
    theme_bw()
  
  omegas_plot <- omegas %>% 
    ggplot() +
    geom_raster(aes(x = var1, y = var2, fill = prob)) +
    scale_fill_gradientn(colors= c("#7f3b08", "#e08214", "#f5f5f5", "#8073ac", "#2d004b"), values = c(0, 0.35, 0.5, 0.65, 1), limits= c(-1, 1)) +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(legend.title= element_blank(), axis.text.x=element_text(angle = 45, hjust = 1)) +
    facet_wrap(.~type)
  
  spatial_plot <- spatial %>% filter(variable == "Quercus suber" |
                                        variable == "Pinus sylvestris" |
                                        variable == "Pinus pinea") %>% 
    ggplot(aes(x = x, y = y, fill = th_pred, col = th_pred)) +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_colour_viridis_c() +
    theme(legend.position="bottom", aspect.ratio= 1) +
    facet_grid(factor(dataset, levels = c("under", "full", "over")) ~ variable) +
    theme_bw() +
    theme(legend.title=element_blank()) 
  
  library(ggpubr)
  theme_set(theme_pubr())

  ggarrange(betas_plot, omegas_plot, spatial_plot, ncol = 1, nrow = 3, heights = c(5, 3.5, 6.5),  align = "hv", labels = "auto") %>%
    ggexport(filename = "./figures/composite_figure.pdf", width = 10, height = 15)
}

