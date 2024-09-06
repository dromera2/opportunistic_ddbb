myspecies <- c("Pinus pinaster", "Pinus nigra", "Pinus sylvestris",
               "Pinus halepensis", "Pinus pinea", "Pinus mugo",
               "Quercus petraea", "Quercus ilex", "Quercus suber", 
               "Quercus robur","Quercus pyrenaica", "Quercus pubescens", 
               "Quercus faginea")

Forest_Europe_list <- c("Albania", "Andorra", "Austria", "Belarus", "Belgium",
                        "Bosnia and Herzegovina", "Bulgaria", "Croatia",
                        "Cyprus", "Czech Republic", "Denmark", "Estonia",
                        "Finland", "France", "Georgia", "Germany", "Greece",
                        "Holy see", "Hungary", "Iceland", "Ireland", "Italy",
                        "Latvia", "Liechtenstein", "Lithuania", "Luxembourg",
                        "Northen Macedonia", "Malta", "Republic of Moldova",
                        "Monaco", "Montenegro", "Netherlands", "Norway",
                        "Poland", "Portugal", "Romania", "Serbia", "Slovakia",
                        "Slovenia", "Spain", "Sweden", "Switzerland", "Turkey",
                        "Ukraine", "United Kingdom") %>%
  countrycode(origin = "country.name", destination = "iso3c")

lon_min <- -12

# study_area <- read_sf("./data/mask/Mask.shp")

conf <- list(Latitude = "decimalLatitude",
             Longitude = "decimalLongitude",
             Date_collected = "eventDate",
             Scientific_name = "species")


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

# Download this data file from: https://www.gbif.org/es/occurrence/download/0033370-231002084531237
code <- '0033370-231002084531237'

# my_download_metadata <- occ_download_meta(code)
# gbif_citation(my_download_metadata)

# Download data from this source: https://files.isric.org/soilgrids/former/2017-03-10/aggregated/1km/
pH <- list.files(path = "../../Public/Data/SoilGrids/1km/phh2o/", pattern = ".tif", full.names = TRUE )
sand <- list.files(path = "../../Public/Data/SoilGrids/1km/sand/", pattern = ".tif", full.names = TRUE )

i <- 1:10
d_names <- c("full", "under", "over")

samples <- 100
nchains <- 4
nParallel <- 2
thin <- 50
transient <- 80*thin

var_names <- c(paste0("PC", 1:3), "pH", "sand")
coord_names <- c("x", "y")

formula <- ~poly(PC1, degree = 2, raw = TRUE) + 
  poly(PC2, degree = 2, raw = TRUE) +
  poly(PC3, degree = 2, raw = TRUE) +
  poly(pH, degree = 2, raw = TRUE) +
  poly(sand, degree = 2, raw = TRUE) 
