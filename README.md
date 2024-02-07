# opportunistic_ddbb
To execute the main project code, it is only needed to open _targets.R and run tar_make() from package "targets"
Dependencies = c("Hmsc", "reshape2", "magrittr", "ggplot2", "abind", "dplyr", "dismo", "corrplot", "cooccur", "vegan", "tidyverse", "countrycode", "CoordinateCleaner", "rgbif", "sf", "raster", "terra", "usdm", "geodata", "qs", "ggcorrplot", "bdvis", "sqldf", "tidyterra", "s2", "blockCV"))

To execute the study case code, it is needed to download the input data and change their route on arguments.R, after doing so, it will run completely openning _targets.R and runnning tar_make()

Occurrence data can be downloaded from: https://www.gbif.org/es/occurrence/download/0033370-231002084531237

Soil data need to be downloaded from: https://files.isric.org/soilgrids/former/2017-03-10/aggregated/1km/
