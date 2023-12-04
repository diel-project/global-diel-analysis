#### load packages ####
library(dplyr)
library(bbplot)
library(caper)
library(geiger)
library(phytools)
library(picante)
library(RColorBrewer)


tmp <- read.csv(
  "./data/analysis_units/diel_data.csv"
)

length(unique(tmp$family))
length(unique(tmp$analysis_unit))

#### Uniform reference ####

# now that we have everything set up accordingly let's
#  bring in some of the csvs.
trad <- read.csv(
  "./data/diel.niche_results/uniform/Traditionalspecies.hyp.ref.csv"
)


# read in sample size
ss <- read.csv(
  "./data/diel.niche_results/total.sample.size.by.analysis.unit.csv"
)

all(trad$species == ss$dat.scientificName)

trad$samp_size <- ss$ss

# trad$ref_prob <- ifelse(
#   trad$Hypothesis == trad$Ref.Activity,
#   trad$P.Hypothesis.,
#   1 - trad$P.Hypothesis.
# )

trad <- trad %>% 
  dplyr::group_by(species, hyp.ref) %>%
  dplyr::mutate(
    samp_prob = samp_size / sum(samp_size) * prob 
  ) %>% 
  data.frame()

trad <- trad %>% 
  dplyr::group_by(species, hyp.ref) %>% 
  dplyr::summarise(
    samp_prob = sum(samp_prob)
  )


trad_mean <- trad %>% 
  dplyr::group_by(hyp.ref) %>% 
  dplyr::summarise(
    n = length(samp_prob),
    mean_prob = mean(samp_prob),
    sd_prob = sd(samp_prob)
  ) %>% 
  data.frame()

# plot this out
windows(4.7,3)
tiff(
  "./plots/uniform_prior_overall_accuracy.tiff",
  width = 4.7,
  height = 3,
  res = 900,
  units = "in",
  compression = "lzw"
  
)
par(
  mar = c(4,6.5,1,1),
  lend = 2
)

{
bbplot::blank(
  xlim = c(0,1),
  ylim = c(0.5,4.5),
  bty = "l",
  yaxs = "i"
)
bbplot::axis_blank(1)
bbplot::axis_blank(
  2,
  at = 1:4,
  minor = FALSE
)
bbplot::axis_text(side = 1, line = 0.5)
bbplot::axis_text(
  "Probability of reference diel phenotype",
  side = 1,
  line = 2.25,
  cex = 1.15
)

bbplot::axis_text(
  rev(trad_mean$hyp.ref),
  at = 1:4,
  las = 1,
  side = 2,
  cex = 1.15,
  line = .5
)
}
# add number to test
my_cols <- c("#3B3B3BFF", "#A73030FF", "#EFC000FF", "#0073C2FF")
trad$y <- as.numeric(
  factor(
    trad$hyp.ref,
    levels = rev(trad_mean$hyp.ref)
  )
)
trad$col <- my_cols[5 - trad$y]
set.seed(3)
trad$y <- jitter(trad$y, 1.75)


points(
  x = trad$samp_prob,
  y = trad$y,
  pch = 19,
  col =  scales::alpha(trad$col, 0.5)
)


for(i in 1:4){
  j <- 5 - i
  my_lo <- trad_mean$mean_prob[i] - trad_mean$sd_prob[i]
  my_hi <- trad_mean$mean_prob[i] + trad_mean$sd_prob[i]
  if(my_lo<0){
    my_lo <- 0
  }
  if(my_hi>1){
    my_hi <- 1
  }
  lines(
    y = c(j,j),
    x = c(my_lo, my_hi),
    lwd = 8,
    col = "white"
  )
  lines(
    y = c(j,j),
    x = c(my_lo, my_hi),
    lwd = 4,
    col = my_cols[i]
  )
  points(
    x = trad_mean$mean_prob[i],
    y = j,
    pch = 19,
    col = "white",
    cex = 2.5
  )
  points(
    x = trad_mean$mean_prob[i],
    y = j,
    pch = 19,
    col = my_cols[i],
    cex = 1.5
  )
}
dev.off()
#### Informed reference ####

trad <-read.csv(
  "./data/diel.niche_results/informed/Traditionalspecies.hyp.ref.PRIOR90.csv"
)

# read in sample size
ss <- read.csv(
  "./data/diel.niche_results/total.sample.size.by.analysis.unit.csv"
)

trad$samp_size <- ss$ss

# trad$ref_prob <- ifelse(
#   trad$Hypothesis == trad$Ref.Activity,
#   trad$P.Hypothesis.,
#   1 - trad$P.Hypothesis.
# )

trad <- trad %>% 
  dplyr::group_by(species, hyp.ref) %>%
  dplyr::mutate(
    samp_prob = samp_size / sum(samp_size) * prob 
  ) %>% 
  data.frame()

trad <- trad %>% 
  dplyr::group_by(species, hyp.ref) %>% 
  dplyr::summarise(
    samp_prob = sum(samp_prob)
  )


trad_mean <- trad %>% 
  dplyr::group_by(hyp.ref) %>% 
  dplyr::summarise(
    n = length(samp_prob),
    mean_prob = mean(samp_prob),
    sd_prob = sd(samp_prob)
  ) %>% 
  data.frame()

# plot this out
#windows(4.7,3)
tiff(
  "./plots/informed_prior_overall_accuracy.tiff",
  width = 4.7,
  height = 3,
  units = "in", 
  res = 900,
  compression = "lzw"
)
par(
  mar = c(4,6.5,1,1),
  lend = 2
)

{
  bbplot::blank(
    xlim = c(0,1),
    ylim = c(0.5,4.5),
    bty = "l",
    yaxs = "i"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(
    2,
    at = 1:4,
    minor = FALSE
  )
  bbplot::axis_text(side = 1, line = 0.5)
  bbplot::axis_text(
    "Probability of reference diel phenotype",
    side = 1,
    line = 2.25,
    cex = 1.15
  )
  
  bbplot::axis_text(
    rev(trad_mean$hyp.ref),
    at = 1:4,
    las = 1,
    side = 2,
    cex = 1.15,
    line = .5
  )
}

# add number to test
my_cols <- c("#3B3B3BFF", "#A73030FF", "#EFC000FF", "#0073C2FF")
trad$y <- as.numeric(
  factor(
    trad$hyp.ref,
    levels = rev(trad_mean$hyp.ref)
  )
)
trad$col <- my_cols[5 - trad$y]
set.seed(3)
trad$y <- jitter(trad$y, 1.75)


points(
  x = trad$samp_prob,
  y = trad$y,
  pch = 19,
  col =  scales::alpha(trad$col, 0.5)
)


for(i in 1:4){
  j <- 5 - i
  my_lo <- trad_mean$mean_prob[i] - trad_mean$sd_prob[i]
  my_hi <- trad_mean$mean_prob[i] + trad_mean$sd_prob[i]
  if(my_lo<0){
    my_lo <- 0
  }
  if(my_hi>1){
    my_hi <- 1
  }
  lines(
    y = c(j,j),
    x = c(my_lo, my_hi),
    lwd = 8,
    col = "white"
  )
  lines(
    y = c(j,j),
    x = c(my_lo, my_hi),
    lwd = 4,
    col = my_cols[i]
  )
  points(
    x = trad_mean$mean_prob[i],
    y = j,
    pch = 19,
    col = "white",
    cex = 2.5
  )
  points(
    x = trad_mean$mean_prob[i],
    y = j,
    pch = 19,
    col = my_cols[i],
    cex = 1.5
  )
}
dev.off()



#### prop by family ####



# read in all analysis units
t2 <- read.csv(
  "./data/diel.niche_results/uniform/Traditionalspecies.hyp.ref.csv"
)
# t2 <- read.csv(
#   "D:/diel_data/inform_prior/Informative Prior on Ref Activity/All Analysis Units/Traditionalspecies.hyp.ref.PRIOR90.csv"
# )
ss <- read.csv(
  "./data/diel.niche_results/total.sample.size.by.analysis.unit.csv"
)
t2$ss <- ss$ss


# get proportion agreement
trad_sum <- t2 %>% 
  dplyr::group_by(species, hyp.ref) %>% 
  dplyr::summarise(
    prop_agree = sum((ss / sum(ss)) * prob),
    n = length(hyp.ref)
  ) %>% 
  data.frame()

# read in family info
fin <- read.csv(
  "./data/Traditional.species.hyps.with.analysis.units.csv"
)
trad_sum <- dplyr::inner_join(
  trad_sum,
  dplyr::distinct(fin[,c("species", "family")]),
  by = "species"
)

trad_mean <- trad_sum %>% 
  dplyr::group_by(Family) %>% 
  dplyr::summarise(
    mu = mean(prop_agree),
    sd = sd(prop_agree),
    n = length(unique(species))
  )
trad_mean <- trad_mean[trad_mean$n > 4,]
trad_mean <- trad_mean[order(trad_mean$mu, decreasing = TRUE),]

#windows(4.5,8)
tiff(
  "./plots/uniform_prior_family_accuracy.tiff",
  height = 8,
  width = 4.5,
  units = "in",
  compression = "lzw",
  res = 600
)
par(mar = c(4,8,1,3))
{
  bbplot::blank(
    xlim = c(0,1),
    ylim = c(0.5,26.5),
    bty = "l",
    yaxs = "i"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(
    2,
    at = 1:26,
    minor = FALSE
  )
  bbplot::axis_text(side = 1, line = 0.5)
  bbplot::axis_text(
    "Probability of\nreference diel phenotype",
    side = 1,
    line = 2.75,
    cex = 1.15
  )
  
  bbplot::axis_text(
    rev(trad_mean$Family),
    at = 1:26,
    las = 1,
    side = 2,
    cex = 1.15,
    line = .5
  )
  bbplot::axis_text(
    rev(
      paste0("n = ", trad_mean$n)
    ),
    at = 1:26,
    las = 1,
    side = 4,
    cex = 1,
    line = 0.5
  )
}
# add number to test
my_cols <- c("#3B3B3BFF", "#A73030FF")
#c("#3B3B3BFF", "#A73030FF", "#EFC000FF", "#0073C2FF")
trad_sum <- trad_sum[trad_sum$Family %in% trad_mean$Family,]
trad_sum$y <- as.numeric(
  factor(
    trad_sum$Family,
    levels = rev(trad_mean$Family)
  )
)

trad_sum$col <- my_cols[(trad_sum$y %%2)+1]
set.seed(3)
trad_sum$y <- jitter(trad_sum$y, 1.75)


points(
  x = trad_sum$prop_agree,
  y = trad_sum$y,
  pch = 19,
  col =  scales::alpha(trad_sum$col, 0.5)
)

trad_mean$col <- rep(my_cols, 13)
for(i in 1:26){
  j <- 27 - i
  my_lo <- trad_mean$mu[i] - trad_mean$sd[i]
  my_hi <- trad_mean$mu[i] + trad_mean$sd[i]
  if(my_lo<0){
    my_lo <- 0
  }
  if(my_hi>1){
    my_hi <- 1
  }
  lines(
    y = c(j,j),
    x = c(my_lo, my_hi),
    lwd = 8,
    col = "white"
  )
  lines(
    y = c(j,j),
    x = c(my_lo, my_hi),
    lwd = 4,
    col = trad_mean$col[i]
  )
  points(
    x = trad_mean$mu[i],
    y = j,
    pch = 19,
    col = "white",
    cex = 2.5
  )
  points(
    x = trad_mean$mu[i],
    y = j,
    pch = 19,
    col = trad_mean$col[i],
    cex = 1.5
  )
}
dev.off()


# same thing but informed prior

t2 <- read.csv(
  "./data/diel.niche_results/informed/Traditionalspecies.hyp.ref.PRIOR90.csv"
)
ss <- read.csv(
  "./data/diel.niche_results/total.sample.size.by.analysis.unit.csv"
)
t2$ss <- ss$ss


# get proportion agreement
trad_sum <- t2 %>% 
  dplyr::group_by(species, hyp.ref) %>% 
  dplyr::summarise(
    prop_agree = sum((ss / sum(ss)) * prob),
    n = length(hyp.ref)
  ) %>% 
  data.frame()

# read in family info
fin <- read.csv(
  "./data/diel.niche_results/Traditional.species.hyps.with.analysis.units.csv"
)
trad_sum <- dplyr::inner_join(
  trad_sum,
  dplyr::distinct(fin[,c("species", "family")]),
  by = "species"
)

trad_mean <- trad_sum %>% 
  dplyr::group_by(Family) %>% 
  dplyr::summarise(
    mu = mean(prop_agree),
    sd = sd(prop_agree),
    n = length(unique(species))
  )
trad_mean <- trad_mean[trad_mean$n > 4,]
trad_mean <- trad_mean[order(trad_mean$mu, decreasing = TRUE),]

#windows(4,8)
tiff(
  "./plots/informed_prior_family_accuracy.tiff",
  height = 8,
  width = 4.5,
  units = "in",
  compression = "lzw",
  res = 600
)
par(mar = c(4,8,1,3))
{
  bbplot::blank(
    xlim = c(0,1),
    ylim = c(0.5,26.5),
    bty = "l",
    yaxs = "i"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(
    2,
    at = 1:26,
    minor = FALSE
  )
  bbplot::axis_text(side = 1, line = 0.5)
  bbplot::axis_text(
    "Probability of\nreference diel phenotype",
    side = 1,
    line = 2.75,
    cex = 1.15
  )
  
  bbplot::axis_text(
    rev(trad_mean$Family),
    at = 1:26,
    las = 1,
    side = 2,
    cex = 1.15,
    line = .5
  )
  bbplot::axis_text(
    rev(
      paste0("n = ", trad_mean$n)
    ),
    at = 1:26,
    las = 1,
    side = 4,
    cex = 1,
    line = 0.5
  )
}
# add number to test
my_cols <- c("#3B3B3BFF", "#A73030FF")
#c("#3B3B3BFF", "#A73030FF", "#EFC000FF", "#0073C2FF")
trad_sum <- trad_sum[trad_sum$Family %in% trad_mean$Family,]
trad_sum$y <- as.numeric(
  factor(
    trad_sum$Family,
    levels = rev(trad_mean$Family)
  )
)

trad_sum$col <- my_cols[(trad_sum$y %%2)+1]
set.seed(3)
trad_sum$y <- jitter(trad_sum$y, 1.75)


points(
  x = trad_sum$prop_agree,
  y = trad_sum$y,
  pch = 19,
  col =  scales::alpha(trad_sum$col, 0.5)
)

trad_mean$col <- rep(my_cols, 13)
for(i in 1:26){
  j <- 27 - i
  my_lo <- trad_mean$mu[i] - trad_mean$sd[i]
  my_hi <- trad_mean$mu[i] + trad_mean$sd[i]
  if(my_lo<0){
    my_lo <- 0
  }
  if(my_hi>1){
    my_hi <- 1
  }
  lines(
    y = c(j,j),
    x = c(my_lo, my_hi),
    lwd = 8,
    col = "white"
  )
  lines(
    y = c(j,j),
    x = c(my_lo, my_hi),
    lwd = 4,
    col = trad_mean$col[i]
  )
  points(
    x = trad_mean$mu[i],
    y = j,
    pch = 19,
    col = "white",
    cex = 2.5
  )
  points(
    x = trad_mean$mu[i],
    y = j,
    pch = 19,
    col = trad_mean$col[i],
    cex = 1.5
  )
}
dev.off()


