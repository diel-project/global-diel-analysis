library(sf)
sf::sf_use_s2(FALSE)
library(stars)
library(dplyr)

source("./scripts/reprojection_functions.R")

# I'm assuming you've downloaded all the necessary files from
# browseURL("https://wcshumanfootprint.org/data-access")

my_files <- list.files(
  "D:/GIS/GHF",
  recursive = TRUE,
  full.names = TRUE
)
# read in diel data
a_dat <- read.csv(
  "./data/analysis_units/diel_data.csv"
)
dat <- a_dat[,c("mean_lat", "mean_lon")]
dat <- dplyr::distinct(dat)
dat$spatial_ID <- 1:nrow(dat)

# join on the spatial_ID to join it back up the with original data
a_dat <- dplyr::inner_join(
  a_dat,
  dat,
  by = c("mean_lat", "mean_lon")
)

dat <- sf::st_as_sf(
  dat,
  coords = c("mean_lon", "mean_lat"),
  crs = 4326
  
)

# get the utms
dat$zone <- longlat_to_utm(
  dat
)

dat <- split(
  dat,
  factor(dat$zone)
)
for(i in 1:length(dat)){
  dat[[i]] <- sf::st_transform(
    dat[[i]],
    crs = unique(dat[[i]]$zone)
  )
  dat[[i]] <- sf::st_buffer(
    dat[[i]],
    dist = 1000 * 20
  )
  dat[[i]] <- sf::st_transform(
    dat[[i]],
    crs = 4326
  )
}
dat <- dplyr::bind_rows(dat)

ghf_dat <- matrix(
  NA,
  nrow = nrow(dat),
  ncol = length(my_files) + 1
)

ghf_dat <- data.frame(ghf_dat)
colnames(ghf_dat) <- c(
  "spatialID",
  paste0("ghf_", 2006:2020)
)
for(i in 1:length(my_files)){
  print(
    paste(i ,"of", length(my_files))
  )
  ghf <- stars::read_stars(
    my_files[i]
  )
  dat <- sf::st_transform(
    dat,
    sf::st_crs(ghf)
  )

  ghf_query <- stars::st_extract(
    ghf,
    dat,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  ghf_dat[,i+1] <- data.frame(ghf_query)[,2,drop = TRUE]
  
}

ghf_dat$spatialID <- dat$spatial_ID

# actually it'd be easier to do this in long format
ghf_list <- vector("list", length = length(my_files) + 2)
for(i in 1:length(ghf_list)){
  if(i <= length(my_files)){
  ghf_list[[i]] <- data.frame(
    spatialID = ghf_dat$spatialID,
    year = 2005 + i,
    ghf = ghf_dat[,i+1]
  )
  }else{
    ghf_list[[i]] <- data.frame(
      spatialID = ghf_dat$spatialID,
      year = 2005 + i,
      ghf = ghf_dat[,16]
    )
  }
}

ghf_list <- dplyr::bind_rows(
  ghf_list
)

a_dat$year <- lubridate::year(a_dat$min_date)


# make a data.frame based on spatialID and year

ghf_list <- dplyr::inner_join(
  ghf_list,
  a_dat[,c("spatial_ID", "year", "analysis_unit", "file_name")],
  by = c("spatialID" = "spatial_ID", "year")
)
ghf_list <- dplyr::distinct(ghf_list)
ghf_list <- ghf_list[order(ghf_list$analysis_unit, ghf_list$year),]

write.csv(
  ghf_list[,-which(colnames(ghf_list) == "spatialID")],
  "./data/analysis_units/global_human_footprint_2023_09_25.csv",
  row.names = FALSE
)
