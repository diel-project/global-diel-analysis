library(nimble)
library(dplyr)
library(lubridate)
library(ape)
library(phytools)
library(parallel)


# >= 20 analysis units analysis

# analysis unit level
# 1. latitude
# 2. hours per day (drop if needed)
# 3. human footprint index?

# species trait level
# 1. body size (logged)
# 2. proportion diet that is vertebrate (drop if needed)
# 3. Check on size of distributional range of species.

# add latitude - add quadratic term to this.
# add hours of sunlight & a quadratic term
# possibly interaction?
#od <- read.csv("./data/Traditional.species.hyps.with.analysis.units_old.csv")
dat <- read.csv(
  "./data/Traditional.species.hyps.with.analysis.units.csv"
)

# read in analysis units to join other data
aunit <- read.csv(
  "./data/analysis_units/diel_data.csv"
)
# add some columns
dat <- dplyr::inner_join(
  dat,
  aunit[,c("analysis_unit", "mean_lat","min_date","max_date","scientificName","file_name")],
  by = c("analysis_unit", "species" = "scientificName")
)


# need to link the two based on project, which is messy


source("./nimble/hypothesis_covariates.R")



if(!"hours.per.day" %in% colnames(dat)){
  # read in the hours per day data.
  hpd <- read.csv(
    "./data/hours_per_day_lat.csv"
  )
  
  dat$min_date <- lubridate::mdy(
    dat$min_date
  )
  dat$max_date <- lubridate::mdy(
    dat$max_date
  )
  dat$mean_date <- as.Date(NA)
  
  for(i in 1:nrow(dat)){
    dat$mean_date[i] <- mean.Date(
      c(dat$min_date[i], dat$max_date[i])
    )
  }
  ly <- lubridate::leap_year(dat$mean_date)
  to_sub <- lubridate::yday(lubridate::ymd("2020-02-29"))
  
  dat$ord_date <- lubridate::yday(dat$mean_date)
  dat$ord_date[
    ly & dat$ord_date >= to_sub
  ] <- dat$ord_date[
    ly & dat$ord_date >= to_sub
  ] - 1
  
  dat$floor_lat <- round(
    dat$mean_lat,
    0
  )
  
  dat <- dplyr::inner_join(
    dat,
    hpd,
    by = c('ord_date' = 'Julien.Day', 'floor_lat' = 'latitude')
  )
}



# query the global human footprint data if it has not been
#  loaded. Unfortunately, because the ghf data is so large some
#  hard coding of file paths had to be done. Check ./scripts/ghf.R
#  for additional details.
if(!file.exists("./data/analysis_units/global_human_footprint_2023_09_25.csv")){
  source("./scripts/ghf.R")
}
if(!"ghf" %in% colnames(dat)){
  ghf <- read.csv(
    "./data/analysis_units/global_human_footprint_2023_09_25.csv"
  )
  dat$year <- lubridate::year(
    dat$min_date
  )
  
  dat <- dplyr::inner_join(
    dat,
    ghf,
    by = c("analysis_unit", "file_name", "year"),
    relationship = "many-to-many"
  )
  if(any(duplicated(dat))){
    dat <- dat[!duplicated(dat),]
  }
}

# start subsetting the data

dat <- dat[dat$p_hypothesis>0.8,]

# over 500 grams
dat <- dat[dat$mass > 500,]

dat <- dat[-which(dat$unit_type== "allday"),]

# three is arboreal
dat <- dat[-which(dat$foraging_strata == 3),]
sp_table <- table(dat$species)
loss_amount <- rep(NA, 30)
for(i in 1:length(loss_amount)){
  to_go <- names(sp_table)[sp_table < i]
  loss_amount[i] <- length(sp_table) - length(to_go)
  
}

to_go <- names(sp_table)[sp_table < 20]

dat <- dat[
  -which(dat$species %in% to_go),
]

# MEAN CENTERING

# do the same species mean scaling
dat <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(hpd_scale = hours.per.day - mean(hours.per.day),
                hpd_mean = mean(hours.per.day),
                ghf_scale = log(ghf+1) - mean(log(ghf+1)),
                ghf_mean = mean(log(ghf+1))) %>% 
  data.frame

dat$hpd_mean <- as.numeric(
  scale(
    dat$hpd_mean
  )
)

dat$ghf_mean <- as.numeric(
  scale(
    dat$ghf_mean
  )
)

# do species mean centering for latitude.
dat <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(lat_scale = abs(mean_lat) - mean(abs(mean_lat)),
                sp_mean_lat = mean(abs(mean_lat))) %>% 
  data.frame

# divide by a number so the covariate has an acceptable range
#  for default priors
dat$lat_scale <- dat$lat_scale / 20

dat$sp_mean_lat <- as.numeric(
  scale(
    dat$sp_mean_lat
  )
)


# make hypothesis a category now
dat$hyp_cat <- as.numeric(
  factor(
    dat$hypothesis,
    levels = c("Cathemeral", "Diurnal", "Nocturnal")
  )
)

# getting all the trait data together. Some of it is in dat
#  while others are in sp_traits. First, get the unique 
#  expected range size for a species.

sp_traits <- read.csv(
  "./data/phylo/Phylacine_Trait_data_sp_corrected.csv"
)

sp_traits <- sp_traits[,c("scientificName", "Mass.g")]

# remove duplicates 
sp_traits <- dplyr::distinct(
  sp_traits
)

# and get unique species range sizes
sp_area <- read.csv(
  "./data/analysis_units/species_distrib_range.csv"
)

# combine the two
sp_traits <- dplyr::inner_join(
  sp_traits,
  sp_area,
  by = "scientificName"
)

sp_traits <- sp_traits[which(sp_traits$scientificName %in% dat$species),]

sp_traits$Mass.g <- as.numeric(
  scale(
    log(sp_traits$Mass.g)
  )
)
sp_traits$EOO <- as.numeric(
  scale(
    log(sp_traits$EOO)
  )
)

dat <- dplyr::inner_join(
  dat,
  sp_traits,
  by = c("species" = "scientificName")
)

# Reorganize data to order by species, filename, and then
#  mean date to generate an auto-regressive term.

dat <- dat[order(
  dat$species,
  dat$file_name,
  dat$mean_date
),]

dat$ar1 <- 0

tmp_dat <- split(
  dat,
  factor(
    paste(
      dat$species,
      dat$file_name,
      sep = "-"
    )
  )
)

# generate term
for(i in 1:length(tmp_dat)){
  one_study <- tmp_dat[[i]]
  if(nrow(one_study)>1){
    day_diff <- diff(
      one_study$mean_date
    )
    day_diff <- as.numeric(day_diff)
    one_study$ar1[2:nrow(one_study)] <- day_diff 
  }
  tmp_dat[[i]] <- one_study
}
tmp_dat <- dplyr::bind_rows(
  tmp_dat
)
tmp_dat$ar1[tmp_dat$ar1>365] <- 0

# divide by 27, which is the median difference in days.
tmp_dat$ar1 <- tmp_dat$ar1 / 27

# now divide by reciprocal if > 0 so that longer time
#  lags are close to 0, the median time lag is 
#  1, and any shorter time lags are > 1.
tmp_dat$ar1[tmp_dat$ar1>0] <- 1 / tmp_dat$ar1[tmp_dat$ar1>0]

# overwrite dat
dat <- tmp_dat

unit_dm <- cbind(
  1,
  dat$lat_scale,
  dat$hpd_scale,
  dat$ghf_scale,
  dat$ar1
)

# species dm with inxs,
#  No need to add an intercept here
#  as this is just an extension to the
#  unit_dm.
trait_dm <- cbind(
  dat$Mass.g,
  dat$EOO,
  dat$sp_mean_lat,
  dat$hpd_mean,
  dat$ghf_mean,
  dat$lat_scale * dat$Mass.g,
  dat$lat_scale * dat$EOO ,
  dat$lat_scale * dat$sp_mean_lat,
  dat$hpd_scale * dat$Mass.g,
  dat$hpd_scale * dat$EOO,
  dat$hpd_scale * dat$hpd_mean,
  dat$ghf_scale * dat$Mass.g,
  dat$ghf_scale * dat$EOO,
  dat$ghf_scale * dat$ghf_mean
)

data_list <- list(
  y = dat$hyp_cat,
  trait_dm = trait_dm,
  unit_dm = unit_dm
)
constant_list <- list(
  species_vec = as.numeric(
    factor(
      dat$species
    )
  ),
  family_vec = as.numeric(
    factor(
      dat$family
    )
  ),
  nspecies = length(unique(dat$species)),
  ncov_trait = ncol(trait_dm),
  ncov_unit = ncol(unit_dm),
  ndata = nrow(dat)
)


# Create a function with all the needed code
run_MCMC_allcode <- function(
    seed,
    data,
    cons,
    nimble_model
) {
  library(nimble)
  my_inits <- function(){
    list(
      diur_unit_beta = matrix(
        rnorm(
          cons$nspecies *
          cons$ncov_unit
        ),
        nrow = cons$nspecies,
        ncol = cons$ncov_unit
      ),
      noct_unit_beta = matrix(
        rnorm(
          cons$nspecies *
            cons$ncov_unit
        ),
        nrow = cons$nspecies,
        ncol = cons$ncov_unit
      ),
      diur_trait_beta = rnorm(
        cons$ncov_trait
      ),
      noct_trait_beta = rnorm(
        cons$ncov_trait
      ),
      diur_beta_mu = rnorm(
        cons$ncov_unit
      ),
      noct_beta_mu = rnorm(
        cons$ncov_unit
      ),
      diur_sd = rgamma(
        cons$ncov_unit,
        1,
        1
      ),
      noct_sd = rgamma(
        cons$ncov_unit,
        1,
        1
      )
    )
  }
  set.seed(seed = seed)
  myModel <- nimble::nimbleModel(
    code = nimble_model,
    data = data,
    constants = cons,
    inits = my_inits()
  )
  
  CmyModel <- nimble::compileNimble(
    myModel
  )
  mconf <- nimble::configureMCMC(
    myModel,
    onlySlice = TRUE
  )
  # get parameters
  h <- my_inits()
  mconf$setMonitors(
    names(h)
  )
  
  myMCMC <- nimble::buildMCMC(
    mconf
  )
  CmyMCMC <- nimble::compileNimble(
    myMCMC,
    project = myModel
  )
  
  results <- nimble::runMCMC(
    CmyMCMC, 
    niter = 350000,
    nburnin = 300000,
    thin = 10,
    setSeed = seed
  )
  
  return(results)
}

this_cluster <- parallel::makeCluster(10)
my_start <- Sys.time()
chain_output <- parLapply(
  cl = this_cluster, 
  X = 1:10, 
  fun = run_MCMC_allcode, 
  data = data_list,
  cons = constant_list,
  nimble_model = my_model
)

# It's good practice to close the cluster when you're done with it.
stopCluster(this_cluster)

saveRDS(
  chain_output,
  "./nimble/unconstrained_results_20plus_suncalc_fix.RDS"
)
my_end <- Sys.time()
my_end - my_start
# clear everything and run again
rm(list = ls())
gc()
