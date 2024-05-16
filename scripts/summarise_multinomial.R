library(nimble)
library(dplyr)
library(bbplot)
source("./scripts/mcmc_utility.R")


#### read in data ####
{
# read in data and all that
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
  aunit[,c("analysis_unit", "scientificName",
           "mean_lat", "file_name",
           "min_date", "max_date")],
  by = c(c("species"= "scientificName","analysis_unit"))
)


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
    by = c("analysis_unit", "file_name", "year")
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

species_unscaled_mean <- data.frame(
  species = sort(unique(dat$species))
)
# do the same species mean scaling
dat <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(hpd_scale = hours.per.day - mean(hours.per.day),
                hpd_mean = mean(hours.per.day),
                ghf_scale = log(ghf+1) - mean(log(ghf+1)),
                ghf_mean = mean(log(ghf+1))) %>% 
  data.frame
species_unscaled_mean <- dplyr::inner_join(
  species_unscaled_mean,
  dplyr::distinct(dat[,c("species", "hpd_mean")]),
  by = "species"
)

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

species_unscaled_mean <- dplyr::inner_join(
  species_unscaled_mean,
  dplyr::distinct(dat[,c("species", "ghf_mean")]),
  by = "species"
)

# do species mean centering for latitude.
dat <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(lat_scale = abs(mean_lat) - mean(abs(mean_lat)),
                sp_mean_lat = mean(abs(mean_lat))) %>% 
  data.frame

species_unscaled_mean <- dplyr::inner_join(
  species_unscaled_mean,
  dplyr::distinct(dat[,c("species", "sp_mean_lat")]),
  by = "species"
)
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
# keep them unscaled as well for plotting
sp_traits_us <- sp_traits
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

unit_dm <- cbind(
  1,
  dat$lat_scale,
  dat$hpd_scale,
  dat$ghf_scale
)

unit_formula <- c("intercept", "lat", "hpd", "ghf")
# species dm with inxs
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

trait_formula <- c(
  "mass", "eoo", "mean_lat", "hpd_mean",
  "ghf_mean", "lat:mass", "lat:eoo", "lat:mean_lat",
  "hpd:mass", "hpd:eoo", "hpd:hpd_mean",
  "ghf:mass", "ghf:eoo", "ghf:ghf_mean"
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
  project_vec = as.numeric(
    factor(
      dat$file_name
    )
  ),
  nspecies = length(unique(dat$species)),
  ncov_trait = ncol(trait_dm),
  ncov_unit = ncol(unit_dm),
  ndata = nrow(dat)
)

#### read in mcmc ####
output <- readRDS(
  "./nimble/unconstrained_results_family_project_re.RDS"
)
mc <- do.call(
  "rbind",
  output
)
set.seed(232)
mc <- mc[sample(1:nrow(mc), 10000),]

# make split mcmc object
mc <- split_mcmc(
  mc
)

noct_unit <- t(apply(
  mc$noct_unit_beta[,,4],
  2,
  quantile,
  probs = c(0.05,0.5,0.95)
))

test <- which(abs(rowSums(sign(noct_unit))) == 3)

tmp <- noct_unit[which(abs(rowSums(sign(noct_unit))) == 3),]
row.names(tmp) <- sort(unique(dat$species))[test]
round(tmp,2)


diur_unit <- t(apply(
  mc$diur_unit_beta[,,4],
  2,
  quantile,
  probs = c(0.05,0.5,0.95)
))

test <- which(abs(rowSums(sign(diur_unit))) == 3)

tmp <- diur_unit[which(abs(rowSums(sign(diur_unit))) == 3),]
row.names(tmp) <- sort(unique(dat$species))[test]
round(tmp,2)
#### color scheme ####
my_cols <- c(
  "cathermal" =  "#3B3B3B",
  "diurnal" =  "#EFC000",
  "nocturnal" = "#0073C2"
)


#### community effects, unit level ####
diur_com <- t(
  apply(
    mc$diur_beta_mu,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
)

row.names(diur_com) <- unit_formula
diur_com <- round(diur_com,2)

noct_com <- t(
  apply(
    mc$noct_beta_mu,
    2,
    quantile,
    probs = c(0.025,  0.5, 0.975)
  )
)

row.names(noct_com) <- unit_formula
noct_com <- round(noct_com,2)


#### Nocturnality by trait ####
noct_traits <- apply(
  mc$noct_trait_beta,
  2,
  quantile, 
  probs = c(0.025,0.5,0.975)
)
noct_traits <- t(noct_traits)

row.names(noct_traits) <- trait_formula
round(noct_traits,2)
# Nocturnality by mass
range(dat$mass)

pvec <- seq(500,  4500000, length.out = 5000)
pvec <- c(500, 50000)

# for plotting purposes!
pvec2<- c(500,5000, 50000, 500000, 5000000)

# multiply by 1k because it's actually grams
svec <- log(pvec)
svec <- (svec - mean(log(dat$mass))) / (
  sd(log(dat$mass))
)

onepiece <- which(
  trait_formula %in% c("mass")
)
pred <- mc$noct_beta_mu[,1] + mc$noct_trait_beta[,onepiece] %*% rbind(svec)
pred2 <-mc$diur_beta_mu[,1] + mc$diur_trait_beta[,onepiece] %*% rbind(svec)
pred <- exp(pred)
pred2 <- exp(pred2)

denom <- 1 + pred + pred2
pred_prob <- pred / denom

mean_pred <- t(
  apply(
    pred_prob,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
)
rm(pred_prob, pred2, pred)
gc()

# do it for each species now as well. But here I think
#  we need all of the mean traits for that specific
#  species
mean_traits <- dat[,c("species", "mass", "Mass.g", "EOO",
                      "sp_mean_lat", "hpd_mean", "ghf_mean")]
mean_traits <- dplyr::distinct(
  mean_traits
)

mean_traits <- mean_traits[order(mean_traits$species),]
mean_trait_dm <- rbind(
  mean_traits$Mass.g,
  mean_traits$EOO,
  mean_traits$sp_mean_lat,
  mean_traits$hpd_mean,
  mean_traits$ghf_mean
)

# get family vector to index the correct family level intercept
tmp_sp <- dplyr::distinct(
  dat[,c("species", "family")]
)
tmp_sp <- tmp_sp[order(tmp_sp$species),]

if(!all(tmp_sp$species == mean_traits$species)){
  stop("Error in species names, fix.")
}
tmp_sp$family_vec <- as.numeric(
  factor(
    tmp_sp$family
  )
)

sp_pred_noct <- mc$noct_unit_beta[,,1] +
  mc$noct_trait_beta[,1:5] %*% mean_trait_dm
sp_pred_noct <- sweep(
  sp_pred_noct,
  1,
  mc$noct_beta_mu[,1],
  FUN = "+"
)
# add in the correct family-level variable
for(i in 1:ncol(sp_pred_noct)){
  sp_pred_noct[,i] <- sp_pred_noct[,i] + 
    mc$noct_family_beta[,tmp_sp$family_vec[i]]
}
sp_pred_diur <- mc$diur_unit_beta[,,1] +
  mc$diur_trait_beta[,1:5] %*% mean_trait_dm
sp_pred_diur <- sweep(
  sp_pred_diur,
  1,
  mc$diur_beta_mu[,1],
  FUN = "+"
)
for(i in 1:ncol(sp_pred_diur)){
  sp_pred_diur[,i] <- sp_pred_diur[,i] + 
    mc$diur_family_beta[,tmp_sp$family_vec[i]]
}
sp_pred_noct <- exp(sp_pred_noct)
sp_pred_diur <- exp(sp_pred_diur)
denom <- 1 + sp_pred_noct + sp_pred_diur

sp_pred <- sp_pred_noct / denom

sp_pred <- t(
  apply(
    sp_pred,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
)

windows(5,5)
par(mar = c(6,6,0.5,0.5), lend = 1)
bbplot::blank(
  xlim = c(6,15.5),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(
  1,
  at = log(pvec2),
  minor = FALSE
)
bbplot::axis_blank(
  2
)
bbplot::axis_text(
  side = 2,
  las = 1,
  line = 0.5
)
bbplot::axis_text(
  pvec2 / 1000,
  line = 0.5,
  side = 1,
  at = log(pvec2)
)
bbplot::axis_text(
  "Mass (kg)",
  side =1,
  line = 2.5,
  cex = 1.5
)
bbplot::axis_text(
  "Pr(Nocturnality)",
  side =2,
  line = 2.5,
  cex = 1.5
)


for(i in 1:nrow(sp_pred)){
  lines(
    x = rep(
      log(
        mean_traits$mass[i]
      ),
      2
    ),
    y = sp_pred[i,-2],
    col = "gray60"
  )
}

points(
  x = log(mean_traits$mass),
  y = sp_pred[,2],
  pch = 21,
  cex = 1.2,
  bg = "gray80"
)

dat$month <- lubridate::month(dat$min_date)


bbplot::ribbon(
  x = log(pvec),
  y = mean_pred[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = log(pvec),
  y = mean_pred[,2],
  col = my_cols[3],
  lwd = 5
)

# make pred for mean latitude

range(abs(dat$mean_lat))


pvec <- c(0, 20)

svec <- pvec
svec <- (svec - mean(abs(dat$mean_lat))) / (
  sd(abs(dat$mean_lat))
)

onepiece <- which(
  trait_formula %in% c("mean_lat")
)
pred <- mc$noct_beta_mu[,1] + mc$noct_trait_beta[,onepiece] %*% rbind(svec)
pred2 <-mc$diur_beta_mu[,1] + mc$diur_trait_beta[,onepiece] %*% rbind(svec)
pred <- exp(pred)
pred2 <- exp(pred2)

denom <- 1 + pred + pred2
pred_prob <- pred / denom

mean_pred <- t(
  apply(
    pred_prob,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
)
round(mean_pred,2)
median(pred_prob[,1] / pred_prob[,2])

#### summarise species level coefficients ####

diur_unit_sum <- array(
  NA,
  dim = c(3, 126, 3)
)

mean_trait_dm <- rbind(
  mean_traits$Mass.g,
  mean_traits$EOO,
  mean_traits$sp_mean_lat,
  mean_traits$hpd_mean,
  mean_traits$ghf_mean
)

coef_loc <-  list(
  lat = list(
    unit = c(2),
    trait = c(1:3, 6:8),
    trait_dm = rbind(
      mean_trait_dm[1:3,],
      mean_trait_dm[1:3,]
    )
  ),
  hpd = list(
    unit = c(3),
    trait = c(1:2,4, 9:11),
    trait_dm = rbind(
      mean_trait_dm[c(1:2,4),],
      mean_trait_dm[1,],
      mean_trait_dm[2,],
      mean_trait_dm[4,]
    )
  ),
  ghf = list(
    unit = c(4),
    trait = c(1:2,5, 12:14),
    trait_dm = rbind(
      mean_trait_dm[c(1:2,5),],
      mean_trait_dm[1,],
      mean_trait_dm[2,],
      mean_trait_dm[5,]
    )
  )
)

coef_loc <-  list(
  lat = list(
    unit = c(2),
    trait = c(1:3),
    trait_dm = rbind(
      mean_trait_dm[1:3,]
    )
  ),
  hpd = list(
    unit = c(3),
    trait = c(1:2,4),
    trait_dm = rbind(
      mean_trait_dm[c(1:2,4),]
    )
  ),
  ghf = list(
    unit = c(4),
    trait = c(1:2,5),
    trait_dm = rbind(
      mean_trait_dm[c(1:2,5),]
    )
  )
)


for(i in 1:3){
  tmp <- mc$diur_beta_mu[,
    coef_loc[[i]]$unit] + 
    mc$diur_unit_beta[,,
    coef_loc[[i]]$unit
  ] +
    mc$diur_trait_beta[,
      coef_loc[[i]]$trait
    ] %*% coef_loc[[i]]$trait_dm

  diur_unit_sum[,,i] <- apply(
    tmp,
    2,
    quantile,
    probs = c(0.05,0.5,0.95)
  )
}
diur_unit_sign <- matrix(
  NA,
  ncol = dim(diur_unit_sum)[3],
  nrow = dim(diur_unit_sum)[2]
)

for(i in 1:3){
  tmp <- sign(diur_unit_sum[,,i])
  diur_unit_sign[,i] <- apply(
    tmp,
    2,
    function(x) length(unique(x)) == 1
  )
}




diur_sign_plot <- diur_unit_sign[
  order(diur_unit_sum[2,,3]),
]

diur_sum_plot <- diur_unit_sum[,
                               order(diur_unit_sum[2,,3]),
]



}
tiff(
  "./plots/figure_s5_diurnality_slope_terms.tiff",
  height = 6,
  width = 7,
  units = "in",
  res = 600,
  compression = "lzw"
)


m <- matrix(1:3, ncol = 3)
layout(m)
my_pch <- rep(18,3)# c(18:20)
par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,0,0,0))
my_type <- c("A) Latitude",
             "B) Hours per day", 
             "C) Global human footprint"
)
my_range <- list(
  c(-10,10),
  c(-10, 10),
  c(-10,10)
)



for(i in 1:3){
  bbplot::blank(xlim = my_range[[i]], ylim = c(0,140),
                yaxs = "i")
  
  bbplot::axis_blank(1)
  bbplot::axis_text(side = 1, line = 0.75)
  lines(
    x = c(0,0),
    y = c(0, 126),
    lty = 2
  )
  
  for(j in 1:126){
    lines(
      x = diur_sum_plot[-2,j,i],
      y = rep(j,2),
      col = ifelse(
        diur_sign_plot[j,i],
        my_cols[1],
        scales::alpha( my_cols[3],0.5)
      ),
      lwd = 2
    )
  }
  
  points(
    x = diur_sum_plot[2,,i],
    y = 1:126,
    pch = my_pch[i],
    col = ifelse(
      diur_sign_plot[,i],
      my_cols[1],
      scales::alpha( my_cols[3],0.5)
    ),
    cex = 1.7
  )
  u <- par("usr")
  text(
    x = u[1], 
    y = 135,
    label = my_type[i],
    cex=1.4,
    pos = 4
  )
  if(i == 2){
    bbplot::axis_text(
      "Species-specific diurnality slope terms",
      outer = TRUE,
      line = 3,
      cex = 1.4,
      side = 1,
      at = 0.5
    )
  }
  
  
}


dev.off()

noct_unit_sum <- array(
  NA,
  dim = c(3, 126, 3)
)


for(i in 1:3){
  tmp <- 
    mc$noct_beta_mu[,coef_loc[[i]]$unit]+ 
    mc$noct_unit_beta[,,
                           coef_loc[[i]]$unit
  ] +
    mc$noct_trait_beta[,
                       coef_loc[[i]]$trait
    ] %*% coef_loc[[i]]$trait_dm
  noct_unit_sum[,,i] <- apply(
    tmp,
    2,
    quantile,
    probs = c(0.05,0.5,0.95)
  )
}
noct_unit_sign <- matrix(
  NA,
  ncol = dim(noct_unit_sum)[3],
  nrow = dim(noct_unit_sum)[2]
)

for(i in 1:3){
  tmp <- sign(noct_unit_sum[,,i])
  noct_unit_sign[,i] <- apply(
    tmp,
    2,
    function(x) length(unique(x)) == 1
  )
}





noct_sign_plot <- noct_unit_sign[
  order(noct_unit_sum[2,,3]),
]

noct_sum_plot <- noct_unit_sum[,
                               order(noct_unit_sum[2,,3]),
]


noct_sp <- data.frame(
  beta = round(noct_unit_sum[2,which(noct_unit_sign[,3]), 3],2),
  sp = sp_traits$scientificName[which(noct_unit_sign[,3])]
)
noct_sp <- noct_sp[order(noct_sp$beta, decreasing = TRUE),]

both_sides <- diur_unit_sign + noct_unit_sign

colSums(both_sides > 0) / nrow(both_sides)

tiff(
  "./plots/figure_s6_nocturnality_slope_terms.tiff",
  height = 6,
  width = 7,
  units = "in",
  res = 600,
  compression = "lzw"
)

m <- matrix(1:3, ncol = 3)
layout(m)
my_pch <- rep(18,3) 
par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,0,0,0))
my_type <- c("A) Latitude",
             "B) Hours per day", 
             "C) Global human footprint"
)
my_range <- list(
  c(-10,10),
  c(-10, 10),
  c(-10,10)
)



for(i in 1:3){
  bbplot::blank(xlim = my_range[[i]], ylim = c(0,140),
                yaxs = "i")
  
  bbplot::axis_blank(1)
  bbplot::axis_text(side = 1, line = 0.75)
  lines(
    x = c(0,0),
    y = c(0, 126),
    lty = 2
  )
  
  for(j in 1:126){
    lines(
      x = noct_sum_plot[-2,j,i],
      y = rep(j,2),
      col = ifelse(
        noct_sign_plot[j,i],
        my_cols[1],
        scales::alpha( my_cols[3],0.5)
      ),
      lwd = 2
    )
  }
  
  points(
    x = noct_sum_plot[2,,i],
    y = 1:126,
    pch = my_pch[i],
    col = ifelse(
      noct_sign_plot[,i],
      my_cols[1],
      scales::alpha( my_cols[3],0.5)
    ),
    cex = 1.7
  )
  u <- par("usr")
  text(
    x = u[1], 
    y = 135,
    label = my_type[i],
    cex=1.4,
    pos = 4
  )
  if(i == 2){
    bbplot::axis_text(
      "Species-specific nocturnality slope terms",
      outer = TRUE,
      line = 3,
      cex = 1.4,
      side = 1,
      at = 0.5
    )
  }
  
  
}

dev.off()


#### diurnality by trait ####

diur_traits <- apply(
  mc$diur_trait_beta,
  2,
  quantile, 
  probs = c(0.025,0.5,0.975)
)
diur_traits <- t(diur_traits)
row.names(diur_traits) <- trait_formula
round(diur_traits,2)
# Nocturnality by mass
sp_traits_us <- sp_traits_us[
  sp_traits_us$scientificName %in%
    sp_traits$scientificName,
]
range(sp_traits_us$EOO)
#7.240247e+07
#72500000
pvec <- seq(300,  72500000, length.out = 5000)
pvec <- c(90000, 450000)

# for plotting purposes!
pvec2<- c(500,5000, 50000, 500000, 5000000)

# multiply by 1k because it's actually grams
svec <- log(pvec)
svec <- (svec - mean(log(sp_traits_us$EOO))) / (
  sd(log(sp_traits_us$EOO))
)

onepiece <- which(
  trait_formula %in% c("eoo")
)
pred <-  mc$noct_trait_beta[,onepiece] %*% rbind(svec)
pred2 <- mc$diur_trait_beta[,onepiece] %*% rbind(svec)
pred <- exp(pred)
pred2 <- exp(pred2)

denom <- 1 + pred + pred2
pred_prob <-pred2 / denom

mean_pred <- t(
  apply(
    pred_prob,
    2,
    qw
  )
)
rm(pred_prob, pred2, pred)
gc()

# diurnality by latitude

# recalc mean lat

tmp_lat <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(
    mean_lat =  mean(abs(mean_lat))
  ) 



range(dat$sp_mean_lat)
#7.240247e+07
#72500000

dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(lat_scale = abs(mean_lat) - mean(abs(mean_lat)),
                sp_mean_lat = mean(abs(mean_lat))) %>% 
  data.frame

pvec <-c(0, 20)


# multiply by 1k because it's actually grams
svec <- pvec
svec <- (svec - mean(tmp_lat$mean_lat)) / (
  sd(tmp_lat$mean_lat)
)

onepiece <- which(
  trait_formula %in% c("mean_lat")
)
pred <-  mc$noct_beta_mu[,1] +  mc$noct_trait_beta[,onepiece] %*% rbind(svec)
pred2 <- mc$diur_beta_mu[,1] +  mc$diur_trait_beta[,onepiece] %*% rbind(svec)
pred <- exp(pred)
pred2 <- exp(pred2)

denom <- 1 + pred + pred2
pred_prob <- pred / denom

mean_pred <- t(
  apply(
    pred_prob,
    2,
    qw
  )
)
round(mean_pred,2)
rm(pred_prob, pred2, pred)
gc()


# do it for each species now as well. But here I think
#  we need all of the mean traits for that specific
#  species
mean_traits <- dat[,c("species", "mass", "Mass.g", "EOO",
                      "sp_mean_lat", "hpd_mean", "ghf_mean")]
mean_traits <- dplyr::distinct(
  mean_traits
)

mean_traits <- mean_traits[order(mean_traits$species),]
mean_trait_dm <- rbind(
  1,
  mean_traits$Mass.g,
  mean_traits$EOO,
  mean_traits$sp_mean_lat,
  mean_traits$hpd_mean,
  mean_traits$ghf_mean
)
sp_pred_noct <- mc$noct_unit_beta[,,1] +
  mc$noct_trait_beta[,1:6] %*% mean_trait_dm
sp_pred_noct <- sweep(
  sp_pred_noct,
  1,
  mc$noct_beta_mu[,1],
  FUN = "+"
)
sp_pred_diur <- mc$diur_unit_beta[,,1] +
  mc$diur_trait_beta[,1:6] %*% mean_trait_dm
sp_pred_diur <- sweep(
  sp_pred_diur,
  1,
  mc$diur_beta_mu[,1],
  FUN = "+"
)
sp_pred_noct <- exp(sp_pred_noct)
sp_pred_diur <- exp(sp_pred_diur)
denom <- 1 + sp_pred_noct + sp_pred_diur

sp_pred <- sp_pred_diur / denom

sp_pred <- t(
  apply(
    sp_pred,
    2,
    qw
  )
)

windows(5,5)
par(mar = c(6,6,0.5,0.5), lend = 1)
bbplot::blank(
  xlim = c(5.5,18.1),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(
  1,
  at = log(pvec2),
  minor = FALSE
)
bbplot::axis_blank(
  2
)
bbplot::axis_text(
  side = 2,
  las = 1,
  line = 0.5
)
bbplot::axis_text(
  pvec2/1000,
  line = 0.5,
  side = 1,
  at = log(pvec2)
)
bbplot::axis_text(
  "Distributional extent (1000 km sq.)",
  side =1,
  line = 2.5,
  cex = 1.4
)
bbplot::axis_text(
  "Pr(Cathermal)",
  side =2,
  line = 2.5,
  cex = 1.5
)


for(i in 1:nrow(sp_pred)){
  lines(
    x = rep(
      log(
        sp_traits_us$EOO[i]
      ),
      2
    ),
    y = sp_pred[i,-2],
    col = "gray60"
  )
}

points(
  x = log(sp_traits_us$EOO),
  y = sp_pred[,2],
  pch = 21,
  cex = 1.2,
  bg = "gray80" 
)



bbplot::ribbon(
  x = log(pvec),
  y = mean_pred[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = log(pvec),
  y = mean_pred[,2],
  col = my_cols[3],
  lwd = 5
)


# summarise species level coefficients

diur_unit_sum <- apply(
  mc$diur_unit_beta,
  c(2,3),
  quantile,
  probs = c(0.025,0.5,0.975)
)

diur_unit_sign <- sign(diur_unit_sum)

diur_summary <- apply(
  diur_unit_sign[-2,,],
  c(2,3),
  function(x) length(unique(x)) == 1
)

diur_mu <- apply(
  mc$diur_beta_mu,
  2,
  quantile,
  probs  =c(0.025,0.5,0.975)
)


windows(7,3)

m <- matrix(1:3, ncol = 3)
layout(m)
my_pch <- c(18:20)
par(mar = c(0.5,0.5,0.5,0.5), oma = c(5,0,4,0))
my_type <- c("A) Latitude",
             "B) Hours per day", 
             "C) Global human footprint"
)
for(i in 1:3){
  bbplot::blank(xlim = c(-3,3), ylim = c(0,140))
  
  bbplot::axis_blank(1)
  bbplot::axis_text(side = 1)
  bbplot::axis_text("Parameter estimate")
  rect(
    xleft = diur_mu[1,i+1],
    ybottom = 1,
    xright = diur_mu[3,i+1],
    ytop = 126,
    col = "gray70",
    border = NA
  )
  
  points(
    x = diur_unit_sum[2,,i+1],
    y = 1:126,
    pch = my_pch[i],
    col = scales::alpha("purple", 0.5),
    cex = 1.3
  )
  
  lines(
    x = rep(diur_mu[2,i+1],2),
    y = c(1,126),
    lwd = 2,
    lty = 2,
    col = "black"
  )
  text(
    x = -3, 
    y = 135,
    label = my_type[i],
    cex=1.2,
    pos = 4
  )
  
}

#### unit-level plots ####

# this will be a 3x3 plot where we
#  have the interacting variable 2 sd
#  above and below the mean value. 
#  Everything else will be kept
#  at it's mean value for now.

# interaction variable
invar <- c(-2,2)
jj <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(mu = mean(abs(mean_lat)))

# latitude
range(abs(dat$mean_lat))
# this is based on 'relative to a species mean'
#  location, so we are going above and below that.
pvec <- seq(-5,20, 0.25)
# this was how we scaled it in the model
svec <- pvec / 20

lat_mass <- gen_preds_unit_inxs(
  mcmc = mc,
  c(1,2),
  c(1,6),
  svec,
  c(-1.492494, 2) # not-2 to 2 because of bounded covariate
)
# provide x axis for plotting
lat_mass$pvec <- pvec



# lat EOO

lat_eoo <- gen_preds_unit_inxs(
  mcmc = mc,
  c(1,2),
  c(2,7),
  svec,
 invar
)
# provide x axis for plotting
lat_eoo$pvec <- pvec

# lat, mean lat
lat_ml<- gen_preds_unit_inxs(
  mcmc = mc,
  c(1,2),
  c(3,8),
  svec,
  invar
)
lat_ml$pvec <- pvec

# hours per day as a function of mass
range(dat$hpd_scale)
pvec <- seq(-6,6, length.out = 200)
# no need for scaled vector

hpd_mass <- gen_preds_unit_inxs(
  mc,
  c(1,3),
  c(1,9),
  pvec,
  invar
)

hpd_mass$pvec <- pvec
# hours per day eoo

hpd_eoo <- gen_preds_unit_inxs(
  mc,
  c(1,3),
  c(2,10),
  pvec,
  c(-1.246657,  0.177646)#invar
)
hpd_eoo$pvec <- pvec

# hours per day mean hpd
hpd_mh <- gen_preds_unit_inxs(
  mc,
  c(1,3),
  c(4,11),
  pvec,
  invar
)
hpd_mh$pvec <- pvec

# ghf

range(dat$ghf_scale)
pvec <- seq(-7,4, length.out = 200)
# no need for scaled vector
ghf_mass <- gen_preds_unit_inxs(
  mc,
  c(1,4),
  c(1,12),
  pvec,
  invar
)

ghf_mass$pvec <- pvec
# ghf eoo

ghf_eoo <- gen_preds_unit_inxs(
  mc,
  c(1,4),
  c(2,13),
  pvec,
  invar
)
ghf_eoo$pvec <- pvec

# ghf day mean ghf
ghf_mg <- gen_preds_unit_inxs(
  mc,
  c(1,4),
  c(5,14),
  pvec,
  invar
)
ghf_mg$pvec <- pvec

# plot out nocturnality



#### nocturnal plot start ####
#windows(12,12)
tiff(
  "./plots/analysis_unit_plasticity_nocturnality2sd.tiff",
  width = 12,
  height = 12,
  units = "in",
  res = 600,
  compression = "lzw"
)
m <- matrix(1:9, ncol =3, nrow = 3)
par(mar = c(1,1,1,1), oma = c(8,8,1,8), lend = 1)
layout(m) 
{
# lat plots
bbplot::blank(
  xlim = range(lat_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 2, las = 1, line = 0.5)
bbplot::ribbon(
  x = lat_mass$pvec,
  y = lat_mass$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = lat_mass$pvec,
  y = lat_mass$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = lat_mass$pvec,
  y = lat_mass$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = lat_mass$pvec,
  y = lat_mass$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")
par(xpd = NA)
text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "A) Mass",
  pos = 4,
  cex = 1.5
)

bbplot::blank(
  xlim = range(lat_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 2, las = 1, line = 0.5)
bbplot::ribbon(
  x = lat_eoo$pvec,
  y = lat_eoo$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = lat_eoo$pvec,
  y = lat_eoo$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = lat_eoo$pvec,
  y = lat_eoo$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = lat_eoo$pvec,
  y = lat_eoo$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "B) Distributional extent",
  pos = 4,
  cex = 1.5
)

bbplot::blank(
  xlim = range(lat_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 2, las = 1, line = 0.5)
bbplot::axis_text(side = 1, line = 0.85)
bbplot::ribbon(
  x = lat_ml$pvec,
  y = lat_ml$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = lat_ml$pvec,
  y = lat_ml$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = lat_ml$pvec,
  y = lat_ml$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = lat_ml$pvec,
  y = lat_ml$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "C) Mean distance from equator",
  pos = 4,
  cex = 1.5
)

bbplot::axis_text(
  "Distance from equator",
  side = 1,
  at = mean(lat_mass$pvec),
  line = 4
)
# HPD plot


# hpd plots
bbplot::blank(
  xlim = range(hpd_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::ribbon(
  x = hpd_mass$pvec,
  y = hpd_mass$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = hpd_mass$pvec,
  y = hpd_mass$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = hpd_mass$pvec,
  y = hpd_mass$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = hpd_mass$pvec,
  y = hpd_mass$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "D) Mass",
  pos = 4,
  cex = 1.5
)

bbplot::blank(
  xlim = range(hpd_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::ribbon(
  x = hpd_eoo$pvec,
  y = hpd_eoo$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = hpd_eoo$pvec,
  y = hpd_eoo$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = hpd_eoo$pvec,
  y = hpd_eoo$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = hpd_eoo$pvec,
  y = hpd_eoo$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "E) Distributional extent",
  pos = 4,
  cex = 1.5
)

bbplot::blank(
  xlim = range(hpd_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 1, line = 0.85)
bbplot::ribbon(
  x = hpd_mh$pvec,
  y = hpd_mh$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = hpd_mh$pvec,
  y = hpd_mh$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = hpd_mh$pvec,
  y = hpd_mh$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = hpd_mh$pvec,
  y = hpd_mh$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "F) Mean hours per day",
  pos = 4,
  cex = 1.5
)
bbplot::axis_text(
  "Hours per day",
  side = 1,
  at = mean(hpd_mass$pvec),
  line = 4
)

# GHF plot



# hpd plots
bbplot::blank(
  xlim = range(ghf_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::ribbon(
  x = ghf_mass$pvec,
  y = ghf_mass$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = ghf_mass$pvec,
  y = ghf_mass$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = ghf_mass$pvec,
  y = ghf_mass$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = ghf_mass$pvec,
  y = ghf_mass$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "G) Mass",
  pos = 4,
  cex = 1.5
)

bbplot::blank(
  xlim = range(ghf_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::ribbon(
  x = ghf_eoo$pvec,
  y = ghf_eoo$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = ghf_eoo$pvec,
  y = ghf_eoo$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = ghf_eoo$pvec,
  y = ghf_eoo$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = ghf_eoo$pvec,
  y = ghf_eoo$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "H) Distributional extent",
  pos = 4,
  cex = 1.5
)

bbplot::blank(
  xlim = range(ghf_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 1, line = 0.85)
bbplot::ribbon(
  x = ghf_mg$pvec,
  y = ghf_mg$nocturnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = ghf_mg$pvec,
  y = ghf_mg$nocturnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)
lines(
  x = ghf_mg$pvec,
  y = ghf_mg$nocturnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = ghf_mg$pvec,
  y = ghf_mg$nocturnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.05,
  label = "I) Global human footprint",
  pos = 4,
  cex = 1.5
)
bbplot::axis_text(
  "Global human footprint",
  side = 1,
  at = mean(ghf_mass$pvec),
  line = 4
)
bbplot::axis_text("Pr(Nocturnality)", side = 2, outer = TRUE, at = 0.5,
                  line = 4, cex = 1.5)
par(xpd = NA)
legend(x = 4.7, y = 1.75, c("-2 SD", "+2 SD"), col = my_cols[c(1,3)],
       lwd = 4, bty = "n", title = "Interacting\nvariable",
       cex = 1.5)

}
dev.off()

#windows(12,12)
tiff(
  "./plots/analysis_unit_plasticity_diurnality2sd.tiff",
  width = 12,
  height = 12,
  units = "in",
  res = 600,
  compression = "lzw"
)
m <- matrix(1:9, ncol =3, nrow = 3)
par(mar = c(1,1,1,1), oma = c(8,8,1,8), lend = 1)
layout(m) 

#### Diurnal plot start ####
{
  # lat plots
  bbplot::blank(
    xlim = range(lat_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, las = 1, line = 0.5)
  bbplot::ribbon(
    x = lat_mass$pvec,
    y = lat_mass$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_mass$pvec,
    y = lat_mass$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = lat_mass$pvec,
    y = lat_mass$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_mass$pvec,
    y = lat_mass$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  par(xpd = NA)
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "A) Mass",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(lat_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, las = 1, line = 0.5)
  bbplot::ribbon(
    x = lat_eoo$pvec,
    y = lat_eoo$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_eoo$pvec,
    y = lat_eoo$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = lat_eoo$pvec,
    y = lat_eoo$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_eoo$pvec,
    y = lat_eoo$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "B) Distributional extent",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(lat_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, las = 1, line = 0.5)
  bbplot::axis_text(side = 1, line = 0.85)
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = lat_ml$pvec,
    y = lat_ml$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_ml$pvec,
    y = lat_ml$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "C) Mean distance from equator",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::axis_text(
    "Distance from equator",
    side = 1,
    at = mean(lat_mass$pvec),
    line = 4
  )
  # HPD plot
  
  
  # hpd plots
  bbplot::blank(
    xlim = range(hpd_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::ribbon(
    x = hpd_mass$pvec,
    y = hpd_mass$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_mass$pvec,
    y = hpd_mass$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = hpd_mass$pvec,
    y = hpd_mass$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_mass$pvec,
    y = hpd_mass$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "D) Mass",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(hpd_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "E) Distributional extent",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(hpd_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.85)
  bbplot::ribbon(
    x = hpd_mh$pvec,
    y = hpd_mh$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_mh$pvec,
    y = hpd_mh$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = hpd_mh$pvec,
    y = hpd_mh$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_mh$pvec,
    y = hpd_mh$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "F) Mean hours per day",
    pos = 4,
    cex = 1.5
  )
  bbplot::axis_text(
    "Hours per day",
    side = 1,
    at = mean(hpd_mass$pvec),
    line = 4
  )
  
  # GHF plot
  
  
  
  # hpd plots
  bbplot::blank(
    xlim = range(ghf_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::ribbon(
    x = ghf_mass$pvec,
    y = ghf_mass$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_mass$pvec,
    y = ghf_mass$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = ghf_mass$pvec,
    y = ghf_mass$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_mass$pvec,
    y = ghf_mass$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "G) Mass",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(ghf_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::ribbon(
    x = ghf_eoo$pvec,
    y = ghf_eoo$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_eoo$pvec,
    y = ghf_eoo$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = ghf_eoo$pvec,
    y = ghf_eoo$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_eoo$pvec,
    y = ghf_eoo$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "H) Distributional extent",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(ghf_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.85)
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "F) Mean global human footprint",
    pos = 4,
    cex = 1.5
  )
  bbplot::axis_text(
    "Global human footprint",
    side = 1,
    at = mean(ghf_mass$pvec),
    line = 4
  )
  bbplot::axis_text("Pr(Diurnality)", side = 2, outer = TRUE, at = 0.5,
                    line = 4, cex = 1.5)
  par(xpd = NA)
  legend(x = 4.7, y = 1.75, c("-2 SD", "+2 SD"), col = my_cols[c(1,3)],
         lwd = 4, bty = "n", title = "Interacting\nvariable",
         cex = 1.5)
  
}
dev.off()


tiff(
  "./plots/analysis_unit_plasticity_cathemerality2sd.tiff",
  width = 12,
  height = 12,
  units = "in",
  res = 600,
  compression = "lzw"
)
m <- matrix(1:9, ncol =3, nrow = 3)
par(mar = c(1,1,1,1), oma = c(8,8,1,8), lend = 1)
layout(m) 

#### Cathemerality plot start ####
{
  # lat plots
  bbplot::blank(
    xlim = range(lat_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, las = 1, line = 0.5)
  bbplot::ribbon(
    x = lat_mass$pvec,
    y = lat_mass$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_mass$pvec,
    y = lat_mass$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = lat_mass$pvec,
    y = lat_mass$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_mass$pvec,
    y = lat_mass$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  par(xpd = NA)
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "A) Mass",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(lat_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, las = 1, line = 0.5)
  bbplot::ribbon(
    x = lat_eoo$pvec,
    y = lat_eoo$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_eoo$pvec,
    y = lat_eoo$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = lat_eoo$pvec,
    y = lat_eoo$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_eoo$pvec,
    y = lat_eoo$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "B) Distributional extent",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(lat_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, las = 1, line = 0.5)
  bbplot::axis_text(side = 1, line = 0.85)
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = lat_ml$pvec,
    y = lat_ml$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_ml$pvec,
    y = lat_ml$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "C) Mean distance from equator",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::axis_text(
    "Distance from equator",
    side = 1,
    at = mean(lat_mass$pvec),
    line = 5
  )
  # HPD plot
  
  
  # hpd plots
  bbplot::blank(
    xlim = range(hpd_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::ribbon(
    x = hpd_mass$pvec,
    y = hpd_mass$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_mass$pvec,
    y = hpd_mass$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = hpd_mass$pvec,
    y = hpd_mass$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_mass$pvec,
    y = hpd_mass$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "D) Mass",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(hpd_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "E) Distributional extent",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(hpd_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.85)
  bbplot::ribbon(
    x = hpd_mh$pvec,
    y = hpd_mh$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_mh$pvec,
    y = hpd_mh$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = hpd_mh$pvec,
    y = hpd_mh$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_mh$pvec,
    y = hpd_mh$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "F) Mean hours per day",
    pos = 4,
    cex = 1.5
  )
  bbplot::axis_text(
    "Hours per day",
    side = 1,
    at = mean(hpd_mass$pvec),
    line = 5
  )
  
  # GHF plot
  
  
  
  # hpd plots
  bbplot::blank(
    xlim = range(ghf_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::ribbon(
    x = ghf_mass$pvec,
    y = ghf_mass$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_mass$pvec,
    y = ghf_mass$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = ghf_mass$pvec,
    y = ghf_mass$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_mass$pvec,
    y = ghf_mass$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "G) Mass",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(ghf_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::ribbon(
    x = ghf_eoo$pvec,
    y = ghf_eoo$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_eoo$pvec,
    y = ghf_eoo$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = ghf_eoo$pvec,
    y = ghf_eoo$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_eoo$pvec,
    y = ghf_eoo$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "H) Distributional extent",
    pos = 4,
    cex = 1.5
  )
  
  bbplot::blank(
    xlim = range(ghf_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.85)
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.05,
    label = "F) Mean global human footprint",
    pos = 4,
    cex = 1.5
  )
  bbplot::axis_text(
    "Global human footprint",
    side = 1,
    at = mean(ghf_mass$pvec),
    line = 5
  )
  bbplot::axis_text("Pr(Cathemerality)", side = 2, outer = TRUE, at = 0.5,
                    line = 4, cex = 1.5)
  par(xpd = NA)
  legend(x = 4.7, y = 1.75, c("-2 SD", "+2 SD"), col = my_cols[c(1,3)],
         lwd = 4, bty = "n", title = "Interacting\nvariable",
         cex = 1.5)
  
}
dev.off()

#### mass_hpd  ####
tiff(
  "./plots/mass_hpd.tiff",
  width = 9, height = 3,
  res = 900,
  units = "in",
  compression = "lzw"
)
#windows(9, 3)
m <- matrix(1:3, ncol = 3)
par(mar = c(1,1,1,1), oma = c(6,6,1,8), lend = 1)
layout(m)
{
  
  
  
  bbplot::blank(
    xlim = range(hpd_mass$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  bbplot::ribbon(
    x = hpd_mass$pvec,
    y = hpd_mass$nocturnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_mass$pvec,
    y = hpd_mass$nocturnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = hpd_mass$pvec,
    y = hpd_mass$nocturnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_mass$pvec,
    y = hpd_mass$nocturnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  bbplot::axis_text(side = 2, las = 1,
                    line = 0.5,cex = 1.3)
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "A) Nocturnal",
    pos = 4,
    cex = 1.5
  )
bbplot::blank(
  xlim = range(hpd_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)

bbplot::axis_text(side = 1, line = 1, cex = 1.3)

bbplot::axis_text(
  "Difference from a species mean hours per day",
  side = 1,
  outer = TRUE,
  at = 0.5,
  cex = 1.3,
  line = 2.8
)
bbplot::axis_text(
  "Probability",
  side = 2,
  outer = TRUE,
  at = 0.5,
  cex = 1.3,
  line = 2.8
)
bbplot::ribbon(
  x = hpd_mass$pvec,
  y = hpd_mass$cathemeral$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = hpd_mass$pvec,
  y = hpd_mass$cathemeral$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)

lines(
  x = hpd_mass$pvec,
  y = hpd_mass$cathemeral$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = hpd_mass$pvec,
  y = hpd_mass$cathemeral$hi[,2],
  lwd = 4,
  col = my_cols[3]
)
u <- par("usr")
par(xpd = NA)

text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.09,
  label = "B) Cathemeral",
  pos = 4,
  cex = 1.5
)




bbplot::blank(
  xlim = range(hpd_mass$pvec),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 1, line = 1, cex = 1.3)
bbplot::ribbon(
  x = hpd_mass$pvec,
  y = hpd_mass$diurnal$lo[,-2],
  col = my_cols[1],
  alpha = 0.3
)
bbplot::ribbon(
  x = hpd_mass$pvec,
  y = hpd_mass$diurnal$hi[,-2],
  col = my_cols[3],
  alpha = 0.3
)

lines(
  x = hpd_mass$pvec,
  y = hpd_mass$diurnal$lo[,2],
  lwd = 4,
  col = my_cols[1]
)
lines(
  x = hpd_mass$pvec,
  y = hpd_mass$diurnal$hi[,2],
  lwd = 4,
  col = my_cols[3]
)


text(
  x = u[1] + ((u[2] - u[1]) * 0.025),
  y = 1.09,
  label = "C) Diurnal",
  pos = 4,
  cex = 1.5
)
}
legend(x = 9.4, y = 0.9, c("Large", "Small"), col = my_cols[c(3,1)],
       lwd = 4, bty = "n", title = "Species\nmass",
       cex = 1.5)
dev.off()




#### ghf_meanghf  ####
tiff(
  "./plots/ghf_meanghf.tiff",
  width = 9, height = 3,
  res = 900,
  units = "in",
  compression = "lzw"
)
#windows(9, 3)
m <- matrix(1:3, ncol = 3)
par(mar = c(1,1,1,1), oma = c(6,6,1,8), lend = 1)
layout(m)
{
  
  
  
  bbplot::blank(
    xlim = range(ghf_mg$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$nocturnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$nocturnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$nocturnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$nocturnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  bbplot::axis_text(side = 2, las = 1,
                    line = 0.5,cex = 1.3)
  u <- par("usr")
  par(xpd = NA)
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "A) Nocturnal",
    pos = 4,
    cex = 1.5
  )
  bbplot::blank(
    xlim = range(ghf_mg$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  
  bbplot::axis_text(
    "Global human footprint (scaled)",
    side = 1,
    outer = TRUE,
    at = 0.5,
    cex = 1.3,
    line = 2.8
  )
  bbplot::axis_text(
    "Probability",
    side = 2,
    outer = TRUE,
    at = 0.5,
    cex = 1.3,
    line = 2.8
  )
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  par(xpd = NA)
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "B) Cathemeral",
    pos = 4,
    cex = 1.5
  )
  
  
  
  
  bbplot::blank(
    xlim = range(ghf_mg$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = ghf_mg$pvec,
    y = ghf_mg$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = ghf_mg$pvec,
    y = ghf_mg$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "C) Diurnal",
    pos = 4,
    cex = 1.5
  )
}
legend(x = 5, y = 0.9, c("-2.3", "1.6"), col = my_cols[c(1,3)],
       lwd = 4, bty = "n", title = "Species\naverage\nglobal\nhuman\nfootprint",
       cex = 1.5)
dev.off()


#### lat_meanlat  ####
tiff(
  "./plots/lat_meanlat.tiff",
  width = 9, height = 3,
  res = 900,
  units = "in",
  compression = "lzw"
)
#windows(9, 3)
m <- matrix(1:3, ncol = 3)
par(mar = c(1,1,1,1), oma = c(6,6,1,8), lend = 1)
layout(m)
{
  
  
  
  bbplot::blank(
    xlim = range(lat_ml$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$nocturnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$nocturnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = lat_ml$pvec,
    y = lat_ml$nocturnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_ml$pvec,
    y = lat_ml$nocturnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  bbplot::axis_text(side = 2, las = 1,
                    line = 0.5,cex = 1.3)
  u <- par("usr")
  par(xpd = NA)
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "A) Nocturnal",
    pos = 4,
    cex = 1.5
  )
  bbplot::blank(
    xlim = range(lat_ml$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  
  bbplot::axis_text(
    "Global human footprint (scaled)",
    side = 1,
    outer = TRUE,
    at = 0.5,
    cex = 1.3,
    line = 2.8
  )
  bbplot::axis_text(
    "Probability",
    side = 2,
    outer = TRUE,
    at = 0.5,
    cex = 1.3,
    line = 2.8
  )
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = lat_ml$pvec,
    y = lat_ml$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_ml$pvec,
    y = lat_ml$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  u <- par("usr")
  par(xpd = NA)
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "B) Cathemeral",
    pos = 4,
    cex = 1.5
  )
  
  
  
  
  bbplot::blank(
    xlim = range(lat_ml$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = lat_ml$pvec,
    y = lat_ml$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = lat_ml$pvec,
    y = lat_ml$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = lat_ml$pvec,
    y = lat_ml$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "C) Diurnal",
    pos = 4,
    cex = 1.5
  )
}
legend(x = 19.9, y = 0.9, c("5", "20"), col = my_cols[c(1,3)],
       lwd = 4, bty = "n", title = "Species\naverage\ndistance\nfrom\nequator",
       cex = 1.5)
dev.off()


#### hpd_eoo  ####
tiff(
  "./plots/hpd_eoo.tiff",
  width = 9, height = 3,
  res = 900,
  units = "in",
  compression = "lzw"
)
#windows(9, 3)
m <- matrix(1:3, ncol = 3)
par(mar = c(1,1,1,1), oma = c(6,6,1,8), lend = 1)
layout(m)
{
  
  
  
  bbplot::blank(
    xlim = range(hpd_eoo$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$nocturnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$nocturnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$nocturnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$nocturnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  bbplot::axis_text(side = 2, las = 1,
                    line = 0.5,cex = 1.3)
  u <- par("usr")
  par(xpd = NA)
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "A) Nocturnal",
    pos = 4,
    cex = 1.5
  )
  bbplot::blank(
    xlim = range(hpd_eoo$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  
  bbplot::axis_text(
    "Difference from a species mean hours per day",
    side = 1,
    outer = TRUE,
    at = 0.5,
    cex = 1.3,
    line = 2.8
  )
  bbplot::axis_text(
    "Probability",
    side = 2,
    outer = TRUE,
    at = 0.5,
    cex = 1.3,
    line = 2.8
  )
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$cathemeral$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$cathemeral$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$cathemeral$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$cathemeral$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "B) Cathemeral",
    pos = 4,
    cex = 1.5
  )
  
  
  
  
  bbplot::blank(
    xlim = range(hpd_eoo$pvec),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 1, cex = 1.3)
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$diurnal$lo[,-2],
    col = my_cols[1],
    alpha = 0.3
  )
  bbplot::ribbon(
    x = hpd_eoo$pvec,
    y = hpd_eoo$diurnal$hi[,-2],
    col = my_cols[3],
    alpha = 0.3
  )
  
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$diurnal$lo[,2],
    lwd = 4,
    col = my_cols[1]
  )
  lines(
    x = hpd_eoo$pvec,
    y = hpd_eoo$diurnal$hi[,2],
    lwd = 4,
    col = my_cols[3]
  )
  
  
  text(
    x = u[1] + ((u[2] - u[1]) * 0.025),
    y = 1.09,
    label = "C) Diurnal",
    pos = 4,
    cex = 1.5
  )
}
legend(x = 6.4, y = 0.9, c("Large", "Small"), col = my_cols[c(3,1)],
       lwd = 4, bty = "n", title = "Distributional\nextent",
       cex = 1.5)
dev.off()


#### ghf with sp level traits ####

# get the logit-scale predictions


# unit level stuff
pu <- seq(-4,4, length.out = 300)

ghf_dm <- array(
  NA,
  dim = c(300,6,length(unique(dat$species)))
)
n <- 300
for(i in 1:dim(ghf_dm)[3]){
  ghf_dm[,,i] <- matrix(
    c(
      rep(mean_trait_dm[1,i], n),
      rep(mean_trait_dm[2,i], n),
      rep(mean_trait_dm[5,i], n),
      pu * mean_trait_dm[1,i],
      pu * mean_trait_dm[2,i],
      pu * mean_trait_dm[5,i]
    ),
    ncol = 6,
    nrow = n
  )
}

cloc <- list(
  unit = 4,
  trait = c(1,2,5, 12:14)
)

# get family vector to index the correct family level intercept
tmp_sp <- dplyr::distinct(
  dat[,c("species", "family")]
)
tmp_sp <- tmp_sp[order(tmp_sp$species),]

if(!all(tmp_sp$species == mean_traits$species)){
  stop("Error in species names, fix.")
}
tmp_sp$family_vec <- as.numeric(
  factor(
    tmp_sp$family
  )
)


my_diff <- rep(NA, 126)
y_ax <- rep(NA, 126)


pb <- txtProgressBar(max = 126)
for(i in 1:126){
  setTxtProgressBar(pb, i)
  noct_intercept <- mc$noct_beta_mu[,1] +
    mc$noct_unit_beta[,i,1] + 
    mc$noct_trait_beta[,cloc$trait[1:3]] %*% t(ghf_dm[,1:3,i]) +
    mc$noct_family_beta[,tmp_sp$family_vec[i]]
  noct_slope <- mc$noct_beta_mu[,4] %*%t(pu) +
    mc$noct_unit_beta[,i,4] %*% t(pu) +
    mc$noct_trait_beta[,cloc$trait[4:6]] %*% t(ghf_dm[,4:6,i])
  
  noct <- exp(noct_intercept + noct_slope)
  diur_intercept <- mc$diur_beta_mu[,1] +
    mc$diur_unit_beta[,i,1] + 
    mc$diur_trait_beta[,cloc$trait[1:3]] %*% t(ghf_dm[,1:3,i]) +
    mc$diur_family_beta[,tmp_sp$family_vec[i]]
  diur_slope <- mc$diur_beta_mu[,4] %*%t(pu) +
    mc$diur_unit_beta[,i,4] %*% t(pu) +
    mc$diur_trait_beta[,cloc$trait[4:6]] %*% t(ghf_dm[,4:6,i])
  diur <- exp(diur_intercept + diur_slope)
  denom <- 1 + noct + diur
  p_noct <- noct / denom
  p_diur <- diur / denom
  p_cath <- 1 / denom
  p_noct <- apply(
    p_noct,
    2,
    quantile,
    probs = c(0.05,0.5,0.955)
  )
  
  my_diff[i] <- p_noct[2,300] - p_noct[2,1]
  y_ax[i] <- p_noct[2,300]
}


tmp_dat <- data.frame(
  species = 
    sp_traits$scientificName[which(noct_unit_sign[,3])],
  diff = round(my_diff[which(noct_unit_sign[,3])],2),
  y = round(y_ax[which(noct_unit_sign[,3])],2),
  id = which(noct_unit_sign[,3])
)

tmp_dat <- tmp_dat[order(tmp_dat$diff, decreasing = TRUE),]
tmp_dat
# see which species have the largest positive change

my_sp <- c(67, 51, 119, 39, 17)
#my_sp <- c(37, 97,51, 24, 112 )
#my_sp <- c(97,17,56, 124,67)

top_species <- vector("list", length = length(my_sp))

# get the family vector for these species
top_sp <- tmp_sp$family_vec[my_sp]
for(i in 1:length(my_sp)){
  noct_intercept <- mc$noct_beta_mu[,1] +
    mc$noct_unit_beta[,my_sp[i],1] + 
    mc$noct_trait_beta[,cloc$trait[1:3]] %*% t(ghf_dm[,1:3,my_sp[i]]) +
    mc$noct_family_beta[,top_sp[i]]
  noct_slope <- mc$noct_beta_mu[,4] %*%t(pu) +
    mc$noct_unit_beta[,my_sp[i],4] %*% t(pu) +
    mc$noct_trait_beta[,cloc$trait[4:6]] %*% t(ghf_dm[,4:6,my_sp[i]])
  noct <- exp(noct_intercept + noct_slope)
  diur_intercept <- mc$diur_beta_mu[,1] +
    mc$diur_unit_beta[,my_sp[i],1] + 
    mc$diur_trait_beta[,cloc$trait[1:3]] %*% t(ghf_dm[,1:3,my_sp[i]]) +
    mc$diur_family_beta[,top_sp[i]]
  diur_slope <- mc$diur_beta_mu[,4] %*%t(pu) +
    mc$diur_unit_beta[,my_sp[i],4] %*% t(pu) +
    mc$diur_trait_beta[,cloc$trait[4:6]] %*% t(ghf_dm[,4:6,my_sp[i]])
  diur <- exp(diur_intercept + diur_slope)
  denom <- 1 + noct + diur
  p_noct <- noct / denom
  top_species[[i]] <- apply(
    p_noct,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
}

sp_traits$scientificName[my_sp]
sp_names <- c(
  "M. mephitis",
  "L. americanus",
  "U. cinereoargenteus",
  "E. dorsatum",
  "C. thous"
)

#sp_names <- c("P. lotor", "C. thous", "L. rufus", "U. arctos", "M. meles")

ci_ord <- rep(NA, 5)
for(i in 1:length(my_sp)){
  ci_ord[i] <- mean(top_species[[i]][3,] - top_species[[i]][1,])
}

obs_range <- dat %>% 
  dplyr::filter(species %in% sp_traits$scientificName[my_sp]) %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(min = min(ghf_scale),
                   max = max(ghf_scale))

plot_ord <- order(ci_ord)

dcol <- c(my_cols, "#F55536", "#1B998B")
tiff(
  "./plots/ghf_nocturnality.tiff",
  height = 4,
  width = 5.25,
  units = "in",
  res = 600,
  compression = "lzw"
)

par(mar = c(5,5,1,8), lend = 1)

bbplot::blank(
  xlim = c(-4,4),
  ylim = c(0,1),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 1, line = 0.7)
bbplot::axis_text(side = 2, line = 0.7, las = 1)
bbplot::axis_text(
  "Global human footprint (scaled)",
  side = 1,
  line = 2.5,
  cex = 1.3
)
bbplot::axis_text(
  "Pr(Nocturnality)",
  side = 2,
  line = 2.8,
  cex = 1.3
)


#for(i in 1:5){
#  j <- plot_ord[i]
#  bbplot::ribbon(
#    x = pu,
#   y =t(top_species[[j]][-2,]),
#    col = dcol[j],
#    alpha = 0.25
#  )
#}

for(i in 1:5){
  j <- plot_ord[i]
  lines(
    x = pu,
    y =top_species[[j]][2,],
    col = dcol[j],
    lwd = 2
  )

}
for(i in 1:5){
  j <- plot_ord[i]
  
  to_keep <- which(
    pu >=obs_range$min[j] &
      pu <= obs_range$max[j]
  )
  lines(
    x = pu[to_keep],
    y =top_species[[j]][2,to_keep],
    col = dcol[j],
    lwd = 8
  )
}
par(xpd = NA)
for(i in 1:5){
  text(
    x = 3.9,
    y = top_species[[i]][2,n] - 0.005 + ifelse(i == 3, 0.025,0 ),
    pos = 4,
    labels = bquote(paste(italic(.(sp_names[i]))))
  )
}

dev.off()





#### quantify plasticity ####


# get range of each covariate

unit_df  <- data.frame(
  mean_lat = seq(-5,20, length.out = 200) / 20,
  hpd = seq(-6,6, length.out = 200),
  ghf =seq(-7,4, length.out = 200)
)
  
# this was how we scaled it in the model


mean_trait_dm <- rbind(
  mean_traits$Mass.g,
  mean_traits$EOO,
  mean_traits$sp_mean_lat,
  mean_traits$hpd_mean,
  mean_traits$ghf_mean
)

# get family vector to index the correct family level intercept
tmp_sp <- dplyr::distinct(
  dat[,c("species", "family")]
)
tmp_sp <- tmp_sp[order(tmp_sp$species),]

if(!all(tmp_sp$species == mean_traits$species)){
  stop("Error in species names, fix.")
}
tmp_sp$family_vec <- as.numeric(
  factor(
    tmp_sp$family
  )
)


coef_loc <-  list(
  lat = list(
    unit = c(1,2),
    trait = c(1:3, 6:8),
    trait_dm = mean_trait_dm[1:3,],
    unit_dm = rbind(
      1,
      unit_df$mean_lat
    )
  ),
  hpd = list(
    unit = c(1,3),
    trait = c(1:2,4, 9:11),
    trait_dm = mean_trait_dm[c(1:2,4),],
    unit_dm = rbind(
      1,
      unit_df$hpd
    )
  ),
  ghf = list(
    unit = c(1,4),
    trait = c(1:2,5, 12:14),
    trait_dm = mean_trait_dm[c(1:2,5),],
    unit_dm = rbind(
      1,
      unit_df$ghf
    )
  )
)

sp_sd <- array(
  NA,
  dim = c(3,200,3, 129, 3)
)

sd_posterior <- vector(
  "list",
  length = 126
)




pb <- txtProgressBar(max = 126)
for(j in 1:126){
  setTxtProgressBar(pb, j)
  # mcmc by (diur, noct, cath) by (lat, hpd, ghf)
  # 
  sd_posterior[[j]] <- array(
    NA,
    dim = c(10000,3,3)
  )
  for(i in 1:3){
    tmp_trait_dm <- rbind(
      matrix(
        coef_loc[[i]]$trait_dm[,j],
        ncol = length(coef_loc[[i]]$unit_dm[2,]),
        nrow = 3
      ),
      sweep(
        matrix(
          coef_loc[[i]]$unit_dm[2,],
          ncol = length(coef_loc[[i]]$unit_dm[2,]),
          nrow = 3,
          byrow = TRUE
        ),
        1,
        coef_loc[[i]]$trait_dm[,j],
        FUN = "*"
      )
    )
  
  tmp_diur <-  (
    mc$diur_beta_mu[, coef_loc[[i]]$unit] +
    mc$diur_unit_beta[,j,coef_loc[[i]]$unit]
    ) %*% coef_loc[[i]]$unit_dm  +
    mc$diur_trait_beta[,coef_loc[[i]]$trait] %*%
    tmp_trait_dm
  tmp_diur <- sweep(
    tmp_diur,
     1,
    mc$diur_family_beta[,tmp_sp$family_vec[j]],
    "+"
  )
    
  tmp_noct <-  (
    mc$noct_beta_mu[, coef_loc[[i]]$unit] +
      mc$noct_unit_beta[,j,coef_loc[[i]]$unit]
  ) %*% coef_loc[[i]]$unit_dm  +
    mc$noct_trait_beta[,coef_loc[[i]]$trait] %*%
    tmp_trait_dm + 
    mc$noct_family_beta[,tmp_sp$family_vec[j]]
  tmp_noct <- sweep(
    tmp_noct,
    1,
    mc$noct_family_beta[,tmp_sp$family_vec[j]],
    "+"
  )
  
  tmp_diur <- exp(tmp_diur)
  tmp_noct <- exp(tmp_noct)
  denom <- 1 + tmp_diur + tmp_noct
  prob_diur <- tmp_diur / denom
  prob_noct <- tmp_noct / denom  
  prob_cath <- 1 / denom
  sd_diur <- apply(
    prob_diur,
    1, 
    sd
  )
  sd_noct <- apply(
    prob_noct,
    1,
    sd
  )
  sd_cath <- apply(
    prob_cath,
    1, 
    sd
  )
  sd_posterior[[j]][,1,i] <- sd_diur
  sd_posterior[[j]][,2,i] <- sd_noct
  sd_posterior[[j]][,3,i] <- sd_cath
  prob_diur <- apply(
    prob_diur,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  prob_noct <- apply(
    prob_noct,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  prob_cath <- apply(
    prob_cath,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  sp_sd[,,1,j,i] <- prob_diur
  sp_sd[,,2,j,i] <- prob_noct
  sp_sd[,,3,j,i] <- prob_cath
  
  }
}

# DO SUBSETTING HERE MASON
# get observed range of scaled covariates
obs_range <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(
    min_lat = min(lat_scale),
    max_lat = max(lat_scale),
    min_hpd = min(hpd_scale),
    max_hpd = max(hpd_scale),
    min_ghf = min(ghf_scale),
    max_ghf = max(ghf_scale)
  )


sd_sum <- matrix(
  NA,
  ncol = 3,
  nrow = 126
)
sd_sum_subset <- matrix(
  NA,
  ncol = 3,
  nrow = 126
)
for(j in 1:3){
  for(k in 1:126){
    tmp <- sp_sd[2,,,k,j]
    tmp_sd <- apply(
      tmp,
      2,
      sd
    )
    #get max sd to represent most change
    sd_sum[k,j] <- max(tmp_sd)
  }
}





my_order <- order(
  sd_sum[,1],
  decreasing = TRUE
)


#species_switch <- sp_sd[2,,3,my_order[i],1] - sp_sd[2,,1,my_order[1],1]
#species_switch <- which(species_switch < 0)[1]
#unit_df$mean_lat[species_switch] * 20

tiff(
  "./plots/figure_s12_most_plastic_lat.tiff",
  height = 10,
  width = 10,
  units = "in",
  res = 600,
  compression = "lzw"
)

m <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
layout(m)
par(mar = c(1,1,2,1), oma = c(7,7,0,9))



for(i in 1:9){
  bbplot::blank(
    xlim = c(-5, 20),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(side = 1)
  bbplot::axis_blank(side = 2)
  u <-par("usr")
  par(xpd = NA)
  text(
    paste0(
      LETTERS[i],") ", sp_traits$scientificName[my_order[i]]
    ),
    x = u[1],
    y = u[4] + 0.02,
    pos = 4,
    cex = 1.5
  )
  for(j in 1:3){
    bbplot::ribbon(
      x = unit_df$mean_lat * 20,
      y = t(sp_sd[-2,,j,my_order[i],1]),
      col = my_cols[c(2,3,1)][j],
      alpha = 0.3
    )
  }
  for(j in 1:3){
    lines(
      x = unit_df$mean_lat * 20,
      y = sp_sd[2,,j,my_order[i],1],
      col = my_cols[c(2,3,1)][j],
      lwd = 4
    )
  }
  if(i %in% c(7:9)){
    bbplot::axis_text(side = 1, line = 1)
  }
  if(i %in% c(1,4,7)){
    bbplot::axis_text(side = 2, line = 1, las = 1)
  }
  if(i == 8){
    bbplot::axis_text(
      "Distance to equator (decimal degrees)",
      side = 1,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 4){
    bbplot::axis_text(
      "Probability",
      side = 2,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 6){
    legend(
      x = 19.8,
      y = 0.67,
      c("Cathemeral", 
        "Diurnal",
        "Nocturnal"),
      fill = my_cols,
      bty = "n",
      cex = 1.6
    )
  }
}

dev.off()




my_order <- order(
  sd_sum[,2],
  decreasing = TRUE
)


tiff(
  "./plots/figure_s13_most_plastic_hpd.tiff",
  height = 10,
  width = 10,
  units = "in",
  res = 600,
  compression = "lzw"
)

m <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
layout(m)
par(mar = c(1,1,2,1), oma = c(7,7,0,9))



for(i in 1:9){
  bbplot::blank(
    xlim = c(-6, 6),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(side = 1)
  bbplot::axis_blank(side = 2)
  u <-par("usr")
  par(xpd = NA)
  text(
    paste0(
      LETTERS[i],") ", sp_traits$scientificName[my_order[i]]
    ),
    x = u[1],
    y = u[4] + 0.02,
    pos = 4,
    cex = 1.5
  )
  for(j in 1:3){
    bbplot::ribbon(
      x = unit_df$hpd,
      y = t(sp_sd[-2,,j,my_order[i],2]),
      col = my_cols[c(2,3,1)][j],
      alpha = 0.3
    )
  }
  for(j in 1:3){
    lines(
      x = unit_df$hpd,
      y = sp_sd[2,,j,my_order[i],2],
      col = my_cols[c(2,3,1)][j],
      lwd = 4
    )
  }
  if(i %in% c(7:9)){
    bbplot::axis_text(side = 1, line = 1)
  }
  if(i %in% c(1,4,7)){
    bbplot::axis_text(side = 2, line = 1, las = 1)
  }
  if(i == 8){
    bbplot::axis_text(
      "Hours per day",
      side = 1,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 4){
    bbplot::axis_text(
      "Probability",
      side = 2,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 6){
    legend(
      x = 5.8,
      y = 0.67,
      c("Cathemeral", 
        "Diurnal",
        "Nocturnal"),
      fill = my_cols,
      bty = "n",
      cex = 1.6
    )
  }
}

dev.off()




my_order <- order(
  sd_sum[,3],
  decreasing = TRUE
)


tiff(
  "./plots/figure_s14_most_plastic_ghf.tiff",
  height = 10,
  width = 10,
  units = "in",
  res = 600,
  compression = "lzw"
)

m <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
layout(m)
par(mar = c(1,1,2,1), oma = c(7,7,0,9))

for(i in 1:9){
  bbplot::blank(
    xlim = c(-7,4),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(side = 1)
  bbplot::axis_blank(side = 2)
  u <-par("usr")
  par(xpd = NA)
  text(
    paste0(
      LETTERS[i],") ", sp_traits$scientificName[my_order[i]]
    ),
    x = u[1],
    y = u[4] + 0.02,
    pos = 4,
    cex = 1.5
  )
  for(j in 1:3){
    bbplot::ribbon(
      x = unit_df$ghf,
      y = t(sp_sd[-2,,j,my_order[i],3]),
      col = my_cols[c(2,3,1)][j],
      alpha = 0.3
    )
  }
  for(j in 1:3){
    lines(
      x = unit_df$ghf,
      y = sp_sd[2,,j,my_order[i],3],
      col = my_cols[c(2,3,1)][j],
      lwd = 4
    )
  }
  if(i %in% c(7:9)){
    bbplot::axis_text(side = 1, line = 1)
  }
  if(i %in% c(1,4,7)){
    bbplot::axis_text(side = 2, line = 1, las = 1)
  }
  if(i == 8){
    bbplot::axis_text(
      "Global human footprint",
      side = 1,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 4){
    bbplot::axis_text(
      "Probability",
      side = 2,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 6){
    legend(
      x = 3.8,
      y = 0.67,
      c("Cathemeral", 
        "Diurnal",
        "Nocturnal"),
      fill = my_cols,
      bty = "n",
      cex = 1.6
    )
  }
}

dev.off()



my_order <- order(
  sd_sum[,1],
  decreasing = FALSE
)


tiff(
  "./plots/figure_s15_least_plastic_lat.tiff",
  height = 10,
  width = 10,
  units = "in",
  res = 600,
  compression = "lzw"
)

m <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
layout(m)
par(mar = c(1,1,2,1), oma = c(7,7,0,9))



for(i in 1:9){
  bbplot::blank(
    xlim = c(-5, 20),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(side = 1)
  bbplot::axis_blank(side = 2)
  u <-par("usr")
  par(xpd = NA)
  text(
    paste0(
      LETTERS[i],") ", sp_traits$scientificName[my_order[i]]
    ),
    x = u[1],
    y = u[4] + 0.02,
    pos = 4,
    cex = 1.5
  )
  for(j in 1:3){
    bbplot::ribbon(
      x = unit_df$mean_lat * 20,
      y = t(sp_sd[-2,,j,my_order[i],1]),
      col = my_cols[c(2,3,1)][j],
      alpha = 0.3
    )
  }
  for(j in 1:3){
    lines(
      x = unit_df$mean_lat * 20,
      y = sp_sd[2,,j,my_order[i],1],
      col = my_cols[c(2,3,1)][j],
      lwd = 4
    )
  }
  if(i %in% c(7:9)){
    bbplot::axis_text(side = 1, line = 1)
  }
  if(i %in% c(1,4,7)){
    bbplot::axis_text(side = 2, line = 1, las = 1)
  }
  if(i == 8){
    bbplot::axis_text(
      "Distance to equator (decimal degrees)",
      side = 1,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 4){
    bbplot::axis_text(
      "Probability",
      side = 2,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 6){
    legend(
      x = 19.8,
      y = 0.67,
      c("Cathemeral", 
        "Diurnal",
        "Nocturnal"),
      fill = my_cols,
      bty = "n",
      cex = 1.6
    )
  }
}

dev.off()




my_order <- order(
  sd_sum[,2],
  decreasing = FALSE
)


tiff(
  "./plots/figure_s16_least_plastic_hpd.tiff",
  height = 10,
  width = 10,
  units = "in",
  res = 600,
  compression = "lzw"
)

m <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
layout(m)
par(mar = c(1,1,2,1), oma = c(7,7,0,9))



for(i in 1:9){
  bbplot::blank(
    xlim = c(-6, 6),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(side = 1)
  bbplot::axis_blank(side = 2)
  u <-par("usr")
  par(xpd = NA)
  text(
    paste0(
      LETTERS[i],") ", sp_traits$scientificName[my_order[i]]
    ),
    x = u[1],
    y = u[4] + 0.02,
    pos = 4,
    cex = 1.5
  )
  for(j in 1:3){
    bbplot::ribbon(
      x = unit_df$hpd,
      y = t(sp_sd[-2,,j,my_order[i],2]),
      col = my_cols[c(2,3,1)][j],
      alpha = 0.3
    )
  }
  for(j in 1:3){
    lines(
      x = unit_df$hpd,
      y = sp_sd[2,,j,my_order[i],2],
      col = my_cols[c(2,3,1)][j],
      lwd = 4
    )
  }
  if(i %in% c(7:9)){
    bbplot::axis_text(side = 1, line = 1)
  }
  if(i %in% c(1,4,7)){
    bbplot::axis_text(side = 2, line = 1, las = 1)
  }
  if(i == 8){
    bbplot::axis_text(
      "Hours per day",
      side = 1,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 4){
    bbplot::axis_text(
      "Probability",
      side = 2,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 6){
    legend(
      x = 5.8,
      y = 0.67,
      c("Cathemeral", 
        "Diurnal",
        "Nocturnal"),
      fill = my_cols,
      bty = "n",
      cex = 1.6
    )
  }
}

dev.off()




my_order <- order(
  sd_sum[,3],
  decreasing = FALSE
)


tiff(
  "./plots/figure_s17_least_plastic_ghf.tiff",
  height = 10,
  width = 10,
  units = "in",
  res = 600,
  compression = "lzw"
)

m <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
layout(m)
par(mar = c(1,1,2,1), oma = c(7,7,0,9))

for(i in 1:9){
  bbplot::blank(
    xlim = c(-7,4),
    ylim = c(0,1),
    bty = "l"
  )
  bbplot::axis_blank(side = 1)
  bbplot::axis_blank(side = 2)
  u <-par("usr")
  par(xpd = NA)
  text(
    paste0(
      LETTERS[i],") ", sp_traits$scientificName[my_order[i]]
    ),
    x = u[1],
    y = u[4] + 0.02,
    pos = 4,
    cex = 1.5
  )
  for(j in 1:3){
    bbplot::ribbon(
      x = unit_df$ghf,
      y = t(sp_sd[-2,,j,my_order[i],3]),
      col = my_cols[c(2,3,1)][j],
      alpha = 0.3
    )
  }
  for(j in 1:3){
    lines(
      x = unit_df$ghf,
      y = sp_sd[2,,j,my_order[i],3],
      col = my_cols[c(2,3,1)][j],
      lwd = 4
    )
  }
  if(i %in% c(7:9)){
    bbplot::axis_text(side = 1, line = 1)
  }
  if(i %in% c(1,4,7)){
    bbplot::axis_text(side = 2, line = 1, las = 1)
  }
  if(i == 8){
    bbplot::axis_text(
      "Global human footprint",
      side = 1,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 4){
    bbplot::axis_text(
      "Probability",
      side = 2,
      outer = TRUE,
      line = 3,
      cex = 1.4,
      at = 0.5
    )
  }
  if(i == 6){
    legend(
      x = 3.8,
      y = 0.67,
      c("Cathemeral", 
        "Diurnal",
        "Nocturnal"),
      fill = my_cols,
      bty = "n",
      cex = 1.6
    )
  }
}

dev.off()

#### quantify total variation ####

tmp <- mc$diur_sd / rowSums(mc$diur_sd)

round(
  apply(
    tmp, 
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  ),
  2
)

round(apply(
  mc$diur_sd,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
),2)

tmp2 <- mc$noct_sd / rowSums(mc$noct_sd)

round(
  apply(
    tmp2, 
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  ),
  2
)
round(apply(
  mc$noct_sd,
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
),2)



#### Figure S5 ####
# a 5 x 3 figure with estimates across each
#  phenotype / species specific estimates.

covariate <- c(
  "mass", "eoo", "mean_lat",
  "hpd_mean", "ghf_mean"
)

# eoo stuff
sp_area <- sp_area[sp_area$scientificName %in% sp_traits$scientificName,]



pvec_list <- list(
  mass = seq(500,  4500000, length.out = 5000),
  eoo = seq(300, 72500000, length.out = 5000),
  mean_lat = seq(0, 60, length.out = 5000),
  hpd_mean = seq(9, 14, length.out = 5000),
  ghf_mean = seq(-4.5, 2, length.out = 5000)
)


sp_traits_us <- sp_traits_us[sp_traits_us$scientificName %in% sp_traits$scientificName,]
species_unscaled_mean$mass <- sp_traits_us$Mass.g
species_unscaled_mean$eoo <- sp_traits_us$EOO
species_unscaled_mean <- species_unscaled_mean[
  c("mass", "eoo", "sp_mean_lat", "hpd_mean", "ghf_mean")
]
svec_list <- pvec_list
svec_list$mass <- 
  (log(pvec_list$mass) - mean(log(sp_traits_us$Mass.g))) /
  sd(log(sp_traits_us$Mass.g))

svec_list$eoo <-
  (log(pvec_list$eoo) - mean(log(sp_traits_us$EOO))) /
  sd(log(sp_traits_us$EOO))

tmp <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(sp_mean_lat = mean(abs(mean_lat))) %>% 
  data.frame %>% 
  dplyr::ungroup() %>% 
  dplyr::summarise(
    mu = mean(sp_mean_lat),
    sd = sd(sp_mean_lat)
  )

svec_list$mean_lat <-
  (svec_list$mean_lat - tmp$mu) / tmp$sd

tmp <- dat %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(
                hpd_mean = mean(hours.per.day)
  ) %>% 
  data.frame %>% 
  dplyr::ungroup() %>% 
  dplyr::summarise(
    mu = mean(hpd_mean),
    sd = sd(hpd_mean)
  )

svec_list$hpd_mean <- 
  (svec_list$hpd_mean - tmp$mu) / tmp$sd

lapply(svec_list, range)

results_list <- vector(
  "list",
  length = length(covariate)
)
for(i in 1:length(covariate)){
  print(i)
  
  onepiece <- which(
    trait_formula %in% c(covariate[i])
  )
  pred <- mc$noct_beta_mu[,1] +
    mc$noct_trait_beta[,onepiece] %*% rbind(svec_list[[i]])
  pred2 <-mc$diur_beta_mu[,1] +
    mc$diur_trait_beta[,onepiece] %*% rbind(svec_list[[i]])
  pred <- exp(pred)
  pred2 <- exp(pred2)
  
  denom <- 1 + pred + pred2
  noct_pred <- pred / denom
  diur_pred <- pred2 / denom
  cath_pred <- 1 / denom 
  
  noct_pred <- t(
    apply(
      noct_pred,
      2,
      quantile,
      probs = c(0.025,0.5,0.975)
    )
  )
  diur_pred <- t(
    apply(
      diur_pred,
      2,
      quantile,
      probs = c(0.025,0.5,0.975)
    )
  )
  cath_pred <- t(
    apply(
      cath_pred,
      2,
      quantile,
      probs = c(0.025,0.5,0.975)
    )
  )
  
  rm(pred2, pred)
  gc()
  
  # do it for each species now as well. But here I think
  #  we need all of the mean traits for that specific
  #  species
  mean_traits <- dat[,c("species", "mass", "Mass.g", "EOO",
                        "sp_mean_lat", "hpd_mean", "ghf_mean")]
  mean_traits <- dplyr::distinct(
    mean_traits
  )
  
  mean_traits <- mean_traits[order(mean_traits$species),]
  mean_trait_dm <- rbind(
    mean_traits$Mass.g,
    mean_traits$EOO,
    mean_traits$sp_mean_lat,
    mean_traits$hpd_mean,
    mean_traits$ghf_mean
  )
  # get family vector to index the correct family level intercept
  tmp_sp <- dplyr::distinct(
    dat[,c("species", "family")]
  )
  tmp_sp <- tmp_sp[order(tmp_sp$species),]
  
  if(!all(tmp_sp$species == mean_traits$species)){
    stop("Error in species names, fix.")
  }
  tmp_sp$family_vec <- as.numeric(
    factor(
      tmp_sp$family
    )
  )
  
  sp_pred_noct <- mc$noct_unit_beta[,,1] +
    mc$noct_trait_beta[,1:5] %*% mean_trait_dm
  sp_pred_noct <- sweep(
    sp_pred_noct,
    1,
    mc$noct_beta_mu[,1],
    FUN = "+"
  )
  # add on family-level random effect
  for(j in 1:126){
    sp_pred_noct[,j] <- sp_pred_noct[,j] + 
      mc$noct_family_beta[,tmp_sp$family_vec[j]]
  }
  sp_pred_diur <- mc$diur_unit_beta[,,1] +
    mc$diur_trait_beta[,1:5] %*% mean_trait_dm
  sp_pred_diur <- sweep(
    sp_pred_diur,
    1,
    mc$diur_beta_mu[,1],
    FUN = "+"
  )
  # add on family-level random effect
  for(j in 1:126){
    sp_pred_diur[,j] <- sp_pred_diur[,j] + 
      mc$diur_family_beta[,tmp_sp$family_vec[j]]
  }
  sp_pred_noct <- exp(sp_pred_noct)
  sp_pred_diur <- exp(sp_pred_diur)
  denom <- 1 + sp_pred_noct + sp_pred_diur
  
  sp_pred_noct <- sp_pred_noct / denom
  sp_pred_diur <- sp_pred_diur / denom
  sp_pred_cath <- 1 / denom
  
  sp_pred_noct <- t(
    apply(
      sp_pred_noct,
      2,
      quantile,
      probs = c(0.025,0.5,0.975)
    )
  )
  sp_pred_diur <- t(
    apply(
      sp_pred_diur,
      2,
      quantile,
      probs = c(0.025,0.5,0.975)
    )
  )
  sp_pred_cath <- t(
    apply(
      sp_pred_cath,
      2,
      quantile,
      probs = c(0.025,0.5,0.975)
    )
  )
  results_list[[i]] <- list(
    mean = list(
      noct = noct_pred,
      diur = diur_pred,
      cath = cath_pred
    ),
    species = list(
      noct = sp_pred_noct,
      diur = sp_pred_diur,
      cath = sp_pred_cath
    )
  )
}



tiff(
  "./plots/figure_5_mean_pred_with_species_points.tiff",
  height = 12,
  width = 9,
  units = "in",
  res = 900,
  compression = "lzw"
)

m <- matrix(
  1:15,
  ncol = 3,
  nrow = 5,
  byrow = TRUE
)
layout(m)
par(mar = c(6,2,1,0.2), oma = c(0,6,2,1))
for(i in 1:5){
  for(j in 1:3){
    my_xlim <- range(pvec_list[[i]])
    if(i %in% c(1,2)){
      if(i == 1){
        pvec2 <- c(500,5000,50000, 500000, 5000000)
      }
      if(i == 2){
        pvec2 <- c(500, 7250, 725000, 72500000)
      }
      my_xlim <- my_xlim / 1000
    
    bbplot::blank(
      xlim = log(my_xlim),
      ylim = c(0,1),
      bty = "l"
    )
    if(i == 1){
      bbplot::axis_text(
        c("Nocturnal", "Diurnal", "Cathemeral")[j],
        side = 3,
        line = 0.9,
        cex = 1.5
      )
    }
      points(
        x = log(species_unscaled_mean[,i] / 1000),
        y = results_list[[i]]$species[[j]][,2],
        pch = 21,
        bg = "gray"
      )
      bbplot::axis_blank(
        1,
        at = log(pvec2/1000),
        minor = FALSE
      )
      bbplot::axis_text(
        text = pvec2/1000,
        at = log(pvec2/1000),
        line = 0.7,
        side = 1
      )
      if(j == 2){
      bbplot::axis_text(
        text = c("Mass (kg)", "Distributional extent (1000 km sq.)")[i],
        at = log(pvec2/1000)[3],
        line = 3.6,
        cex = 1.5,
        side = 1
      )
      }
      
      bbplot::ribbon(
        x = log(pvec_list[[i]] / 1000),
        y = results_list[[i]]$mean[[j]][,-2],
        col = my_cols[4 - j],
        alpha= 0.5
      )
      lines(
        x = log(pvec_list[[i]] / 1000),
        y = results_list[[i]]$mean[[j]][,2],
        col = my_cols[4-j],
        lwd = 3
      )
    }else{
      
      bbplot::blank(
        xlim = my_xlim,
        ylim = c(0,1),
        bty = "l"
      )
      points(
        x = species_unscaled_mean[,i],
        y = results_list[[i]]$species[[j]][,2],
        pch = 21,
        bg = "gray"
      )
      bbplot::axis_blank(
        1
      )
      bbplot::axis_text(
        line = 0.7,
        side = 1
      )
      if(j == 2){
        bbplot::axis_text(
          text = c(NA, NA, 
                   "Mean distance from equator (decimal degrees)",
                   "Mean hours per day",
                   "Mean global human footprint")[i],
          at = mean(pvec_list[[i]]),
          line = 3.6,
          cex = 1.5,
          side = 1
        )
      }
      
      bbplot::ribbon(
        x = pvec_list[[i]],
        y = results_list[[i]]$mean[[j]][,-2],
        col = my_cols[4 - j],
        alpha= 0.5
      )
      lines(
        x = pvec_list[[i]],
        y = results_list[[i]]$mean[[j]][,2],
        col = my_cols[4-j],
        lwd = 3
      )
    }
    bbplot::axis_blank(side=2)
    if(j == 1){
      bbplot::axis_text(side = 2,las = 1, line = 0.7)
      bbplot::axis_text(
        "Probability",
        side = 2,
        at = 0.5,
        line = 4,
        cex = 1.5
        
      )
    }
  }
}
dev.off()
