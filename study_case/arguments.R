myspecies <- c("Pinus pinaster", "Pinus nigra", "Pinus sylvestris",
               "Pinus halepensis", "Pinus pinea", "Pinus mugo",
               "Quercus petraea", "Quercus ilex", "Quercus suber", 
               "Quercus robur","Quercus pyrenaica", "Quercus pubescens", 
               "Quercus faginea")

study_area <- read_sf("./data/mask/Mask.shp")

#it is not needed to run this code, the data code will be provided.
# install.packages("usethis")
# usethis::edit_r_environ()
# 
# GBIF_USER="example_user"
# GBIF_PWD="example_password"
# GBIF_EMAIL="example_mail"

# get.code <- function(myspecies) {
#   gbif_taxon_keys <- myspecies %>%
#     name_backbone_checklist() %>% # match to backbone
#     dplyr::pull(usageKey)
# 
# 
#   occurrence_data <- occ_download(
#     pred_in("taxonKey", gbif_taxon_keys),
#     #### pred("basisOfRecord", "PRESERVED SPECIMEN"),
#     pred("hasCoordinate", TRUE),
#     pred("hasGeospatialIssue", FALSE),
#     format = "SIMPLE_CSV"
#   )
# 
#   occurrence_data
# }
# 
# get.code(myspecies)

code <- '0033370-231002084531237'

pH <- list.files(path = "../Public/Data/SoilGrids/1km/phh2o/", pattern = ".tif", full.names = TRUE )
sand <- list.files(path = "../Public/Data/SoilGrids/1km/sand/", pattern = ".tif", full.names = TRUE )
temp <- rast("./data/climate/pca3_data.tif") %>% subset("PC1") 

samples = 100
nchains = 2
nParallel = 2
thin = 50
transient = 80*thin

var_names <- c(paste0("PC", 1:3), "pH", "sand")
coord_names <- c("x", "y")

formula <- ~poly(PC1, degree = 2, raw = TRUE) + 
  poly(PC2, degree = 2, raw = TRUE) +
  poly(PC3, degree = 2, raw = TRUE) +
  poly(pH, degree = 2, raw = TRUE) +
  poly(sand, degree = 2, raw = TRUE) 

#We load these coarse resolution data to make lighter the spatial predictions for example purposes.
pred_clim <- rast("./data/prediction/clim.tif")
pred_soil <- rast("./data/prediction/soil.tif")