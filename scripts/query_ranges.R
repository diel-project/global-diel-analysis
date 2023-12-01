###################################
#
# Pull range info
# Written by M. Fidino
#
#
####################################

### load packages
library(sf)
library(dplyr)
library(mapview)
library(smoothr)
sf::sf_use_s2(FALSE)

# NOTE: This script requires you to have downloaded
#         the terrestrial mammal data. As such,
#         the path to this dataset has been hardcoded.
#         If you are running this locally you will need
#         to modify any line referencing
#         "D:/TERRESTRIAL_MAMMALS"

#       There were also a couple of species where I had
#         to either look up or redownload their specific
#         range data as it was not in the original 
#         dataset. The one species that I had to
#         re-download was 'Cebus capucinus.' While the
#         script below assumes you have the file
#         D:/GIS/redlist_species_data" to get this number,
#         that part of the script can just be skipped past
#         as the actual value is not in the species range
#         dataset.

sql_IN <- function(x){
  to_return <- paste0("'",x,"'")
  to_return <- paste0(to_return, collapse = ", ")
  to_return <- paste0("(", to_return, ")")
  return(to_return)
}



one_row <- sf::read_sf(
  "D:/TERRESTRIAL_MAMMALS",
  query = "SELECT * FROM TERRESTRIAL_MAMMALS LIMIT 1"
)

iucn_species <- sf::read_sf(
  "D:/TERRESTRIAL_MAMMALS",
  query = "SELECT DISTINCT tm.binomial FROM TERRESTRIAL_MAMMALS tm"
)

# get the legend info to subset down 
leg <- sf::read_sf(
  "D:/TERRESTRIAL_MAMMALS",
  query = "SELECT DISTINCT tm.legend FROM TERRESTRIAL_MAMMALS tm"
)

# read in diel data
dat <- read.csv(
  "./data/analysis_units/diel_data_2023_01_10.csv"
)

diel_species <- unique(dat$scientificName)

sum(diel_species %in% iucn_species$tm.binomial)

# there are 6 that we need to find their range.

to_fix <- diel_species[
  which(!diel_species %in% iucn_species$tm.binomial)
]

sp_info <- dat[
  dat$scientificName %in% to_fix,
  c(
    "scientificName", "phylum", "class", "order", "family"
  )]
sp_info <- sp_info[!duplicated(sp_info),]

# I was able to pull these online from the IUCN website,
#  the last one I needed to download their shape file
sp_info$EOO <- c(
  7930000, # Euphractus sexcinctus,
  2500000, # Cabassous tatouay
  937000, # Tolypeutes tricinctus
  9750000, # Priodontes maximus
  780000, # Cabassous centralis
  NA # Cebus capucinus
)
# the last one we need to pull from a shape file
last_sp <- sf::read_sf(
  "D:/GIS/redlist_species_data"
)
# determine the extent stuff
unique(last_sp$LEGEND)

# it's all resident, so we just need to collapse and sum.
last_sp <- sf::st_intersection(last_sp)
new_area <- sum(sf::st_area(last_sp))
#units(new_area, "km")
units(new_area) <- "km^2"
new_area <- as.numeric(new_area)
sp_info$EOO[5] <- new_area

# we will add these after we grab the stuff for all the others



# figure out the types of polygons we want to query,
#  which is just the ones that start with extant
extant <- leg$tm.legend[
  grep("^Extant", leg$tm.legend)
]

sp_in_db <- diel_species[diel_species %in% iucn_species$tm.binomial]

species_area <- data.frame(
  scientificName = sp_in_db,
  EOO = NA
)

pb <- txtProgressBar(max = length(sp_in_db))
for(i in 1:length(sp_in_db)){
  setTxtProgressBar(pb, i)
  mammal_qry <- paste0(
    "SELECT tm.binomial AS binomial FROM TERRESTRIAL_MAMMALS tm\n",
    "WHERE tm.binomial = '", sp_in_db[i], "'\n",
    "AND tm.legend IN ", sql_IN(extant))
  
  mams <- sf::read_sf("D:/TERRESTRIAL_MAMMALS",
                      query = mammal_qry
  )
  
  parts <- sf::st_make_valid(
    mams
  )
  parts <- suppressMessages(
    sf::st_cast(st_union(parts),"POLYGON")
  )
  parts <- sf::st_make_valid(
    parts
  )
  
  new_area <- sum(
    sf::st_area(parts)
  )
  units(new_area) <- "km^2"
  species_area$EOO[i] <- as.numeric(new_area)
}


# add on the other species I had to grab

species_area <- dplyr::bind_rows(
  species_area, sp_info[,c("scientificName", "EOO")]
)

write.csv(
  species_area,
  "./data/analysis_units/species_distrib_range.csv",
  row.names = FALSE
)

