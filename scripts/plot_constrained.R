#### load packages ####
library(dplyr)
library(bbplot)
library(caper)
library(geiger)
library(phytools)
library(picante)
library(RColorBrewer)

# first off I need to recompile the file
#  structure from what I got off google

# my_files <- vector("list", length = 5)
# dir.create("Species")
# for(i in 1:5){
#   my_files[[i]] <- list.files(
#     paste0(
#       "D:/diel_data/Species-20230201T164740Z-00",
#       i
#     ),
#     recursive = TRUE,
#     full.names = TRUE
#   )
# }
# 
# # create new paths
# for(i in 2:5){
#   tmp <- my_files[[i]]
#   tmp <- strsplit(
#     tmp,
#     "/"
#   )
#   tmp <- lapply(
#     tmp,
#     function(x){x[-3]}
#   )
#   # get directories to create
#   to_create <- lapply(
#     tmp,
#     function(x){x[4:5]}
#   )
#   to_create <- sapply(
#     to_create,
#     function(x){paste0(
#       x, collapse = "/")}
#   )
#   to_create <- paste0(
#     "D:/diel_data/Species/",
#     to_create
#   )
#   for(j in 1:length(to_create)){
#     dir.create(
#       to_create[j],
#       recursive = TRUE
#     )
#   }
#   tmp <- sapply(
#     tmp,
#     function(x) {paste0(x, collapse = "/")}
#   )
#   file.copy(
#     my_files[[i]],
#     tmp
#   )
# }

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


#### Pr(reference) with phylogeny ####

#Load phylogeny data

trad <-read.csv(
  "./data/Traditional.species.hyps.df.with.analysis.units.csv"
)
# remove some stuff from the plotting
trad <- trad[-which(trad$Foraging_strata == 3),]

trad <- trad[trad$Body_mass_IM>= 500,]

trad <- trad[-which(
  trad$unit_type == "allday"
),]

# read in all analysis units
t2 <- read.csv(
  "./data/diel.niche_results/uniform/Traditionalspecies.hyp.ref.csv"
)
ss <- read.csv(
  "./data/diel.niche_results/total.sample.size.by.analysis.unit.csv"
)
t2$ss <- ss$ss

# subset down to the species in trad
t2 <- t2[t2$species %in% trad$Species,]

# read in phylogentic tree
faurby_tree<-read.nexus(
  "D:/GIS/phylacine/DataS1/Complete_phylogeny.nex"
)
# get one of the trees
faurby_tree <- faurby_tree[[1]]

# get proportion agreement
trad_sum <- t2 %>% 
  dplyr::group_by(species, hyp.ref) %>% 
  dplyr::summarise(
    prop_agree = sum((ss / sum(ss)) * prob),
    n = length(hyp.ref)
  ) %>% 
  data.frame()

names <- data.frame(Species=trad_sum$species)
names <- dplyr::distinct(names)
names$Species <- gsub(
  " ",
  "_",
  names$Species
)

rownames(names) <- names$Species

# subset down to the species of interest
overlap_species<-geiger::name.check(
  faurby_tree,
  names
)

mammal_tree <- ape::drop.tip(
  faurby_tree,
  overlap_species$tree_not_data
)


# Format trait data for creating phylogeny
tmat <- matrix(
  trad_sum$prop_agree,
  ncol = 1,
  nrow = length(mammal_tree$tip.label)
)
row.names(tmat) <-names$Species


## Using Pagel's Lambda for continuous traits ###
signal_there <- phytools::phylosig(
  mammal_tree,
  tmat,
  method="lambda",
  test=TRUE
)
signal_there

### Now use contmap to plot continuous traits ###
traitnewvector <- tmat[,1]
names(traitnewvector) <- row.names(tmat)
### Crappy Tree ###
obj <- phytools::contMap(
  mammal_tree,
  traitnewvector,
  fsize=c(0.2,1),
  lwd=1,
  plot = FALSE
)
obj<-phytools::setMap(
  obj,
  colors=c(
    "firebrick4", "brown2", "sienna2", "lightgoldenrod1",
    "seagreen3", "deepskyblue3", "darkslateblue"
  )
)

pdf("phy_signal_uniform.pdf", height = 15, width = 15)
svg(
  "phy_signal_uniform.svg",
  height = 15,
  width = 15
)
tiff(
  "phy_signal_uniform.tiff",
  height = 15,
  width = 15,
  res = 200,
  units = "in",
  compression = "lzw"
)

plot(obj, type="fan", lwd=4, outline=T, fsize=c(1, 1), res = 600)
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
  dplyr::distinct(fin[,c("Species", "Family")]),
  by = c("species" = "Species")
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
  dplyr::distinct(fin[,c("Species", "Family")]),
  by = c("species" = "Species")
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

#### family phylogeny ####


trad <-read.csv(
  "D:/diel_data/Traditionalspecies.hyps.df.csv"
)

# remove some stuff from the plotting
trad <- trad[-which(trad$Foraging_strata == "Arboreal"),]

trad <- trad[trad$Body_mass..g.>= 500,]

trad <- trad[-which(
  trad$Unit.Type == "allday"
),]

library(bdc)
unq <- unique(trad$Species)

tmp <- bdc::bdc_query_names_taxadb(
  sci_name = unq
)
tmp <- tmp[complete.cases(tmp),]
tmp <- tmp[,c("kingdom", "phylum", "class", "order", "family")]
tmp <- dplyr::distinct(tmp)
tmp$order <- factor(tmp$order)
tmp$family <- factor(tmp$family)

tr <- as.phylo(
  ~order/family, data = tmp, collapse=FALSE)
tr$edge.length <- rep(0.25, nrow(tr$edge))


#read in all analysis units
t2 <- read.csv(
  "D:/diel_data/uniform_prior/Uniform Prior on Hyp Set/All Analysis Units/Traditionalspecies.hyp.ref.csv"
)
ss <- read.csv(
  "D:/diel_data/total.sample.size.by.analysis.unit.csv"
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
  "D:/diel_data/Traditionalspecies.hyps.df.csv"
)
trad_sum <- dplyr::inner_join(
  trad_sum,
  dplyr::distinct(fin[,c("Species", "Family")]),
  by = c("species" = "Species")
)

trad_mean <- trad_sum %>% 
  dplyr::group_by(Family) %>% 
  dplyr::summarise(
    mu = mean(prop_agree),
    sd = sd(prop_agree),
    n = length(unique(species))
  )
trad_mean <- trad_mean[trad_mean$Family %in% tmp$family,]



names <- data.frame(family = trad_mean$Family)
names <- dplyr::distinct(names)


rownames(names) <- names$family

# subset down to the species of interest
mammal_tree <- tr


# Format trait data for creating phylogeny
tmat <- matrix(
  trad_mean$mu,
  ncol = 1,
  nrow = length(mammal_tree$tip.label)
)
row.names(tmat) <-names$family


traitnewvector <- tmat[,1]
names(traitnewvector) <- row.names(tmat)
### Crappy Tree ###
obj <- phytools::contMap(
  mammal_tree,
  traitnewvector,
  fsize=c(0.2,1),
  lwd=1,
  plot = FALSE
)
obj<-phytools::setMap(
  obj,
  colors=c(
    "firebrick4", "brown2", "sienna2", "lightgoldenrod1",
    "seagreen3", "deepskyblue3", "darkslateblue"
  )
)

pdf("phy_signal_uniform.pdf", height = 15, width = 15)
svg(
  "phy_signal_uniform.svg",
  height = 15,
  width = 15
)
tiff(
  "family_signal_uniform.tiff",
  height = 20,
  width = 6,
  res = 200,
  units = "in",
  compression = "lzw"
)

plot(obj, 
     type="phylogram", lwd=8, 
     outline=T, fsize=c(2.5, 1), res = 600,
     use.edge.length = FALSE)
dev.off()





#### look at plasticity ####

trad <- read.csv(
  "./data/diel.niche_results/Traditional.species.hyps.with.analysis.units.csv"
)

hm <- trad[trad$P.Hypothesis.>0.8,] %>% 
  dplyr::group_by(Species, Ref.Activity) %>% 
  dplyr::summarise(
    nsupport = length(unique(Hypothesis)),
    n_units = length(Species)
  )

# tack backon species details
hm <- dplyr::inner_join(
  hm,
  dplyr::distinct(trad[,c("Species", "Body_mass..g.")]),
  by = "Species"
)

# read in distrib range
dr <- read.csv(
  "./data/analysis_units/species_distrib_range.csv"
)

hm <- dplyr::inner_join(
  hm,
  dr,
  by = c("Species" = "scientificName")
)
hm <- data.frame(hm)
hm <- hm[hm$n_units>1,]

hm$log_mass <- log(hm$Body_mass..g.)
hm$log_eoo <- log(hm$EOO)

yo <- hm
yo$nsupport <- factor(yo$nsupport)

yo$log_eoo <- as.numeric(scale(yo$log_eoo))

yo$log_mass <- as.numeric(scale(yo$log_mass))
yo$su <- as.numeric(scale(yo$n_units))
longshot <- polr(
  nsupport ~ log_mass + log_eoo + su,
  data = yo
)

cli <- confint(longshot)

# make predictions for body mass

range(hm$log_mass)
mass_dat <- data.frame(
  log_mass = seq(2,15,length.out = 500),
  log_eoo = 0,
  su = 0
)
mds <- mass_dat
mds$log_mass <- (mds$log_mass - mean(hm$log_mass)) / sd(hm$log_mass)


mass_pred <- predict(
  longshot,
  newdata = mds,
  type = "probs",
  se.fit = TRUE
)

nboot <- 1000
boot_array <- array(
  NA, dim = c(nrow(mass_dat), 3, nboot)
)
for(i in 1:nboot){
  tmp_dat <- yo[
    sample(1:nrow(yo), replace = TRUE),
  ]
  l2 <- polr(
    nsupport ~ log_mass + log_eoo + su,
    data = tmp_dat
  )
  mass_pred <- predict(
    l2,
    newdata = mds,
    type = "probs",
    se.fit = TRUE
  )
  boot_array[,,i] <- mass_pred
}
boot_est <- apply(
  boot_array,
  c(1,2),
  quantile,
  probs = c(0.025,0.5,0.975)
)

# make some pretty x axis labels
my_xvals <- c(10, 100, 1000, 10000, 500000)
my_xtitles <- paste0(
  my_xvals/1000
)


windows(6, 4)
tiff("body_mass_support.tiff",
     height = 4, width = 6,
     units = "in", res = 200,
     compression = "lzw")
par(mar = c(4,4,1,1))
bbplot::blank(
  xlim = c(2,15),
  ylim = c(0,1),
  bty = "l",
  xaxs = "i",
  yaxs = "i"
  )
bbplot::axis_blank(1, at = log(my_xvals), minor = FALSE)
bbplot::axis_blank(2)
bbplot::axis_text(my_xtitles,side = 1, line = 0.5,at = log(my_xvals), )
bbplot::axis_text(side = 2, line = 0.5, las = 1)
bbplot::axis_text(
  "Body mass (kg)",
  side = 1,
  line = 2.5,
  cex = 1.25
)
bbplot::axis_text(
  "Probability",
  side = 2,
  line = 2.5,
  cex = 1.25
)
my_cols <-  c("#3B3B3B", "#A73030", "#0073C2")

for(i in 1:3){
  bbplot::ribbon(
    x = mass_dat$log_mass,
    y = t(boot_est[-2,,i]),
    col = my_cols[i],
    alpha = 0.25
  )
}
for(i in 1:3){
  
  lines(
    y = boot_est[2,,i],
    x = mass_dat$log_mass,
    col = my_cols[i],
    lty = i,
    lwd = 5
  )
}
u <- par("usr")
lines(x = c(u[1:2]), y = c(0,0))
lines(x = c(u[1],u[1]), y = c(0,1))
legend(
  "topright",
  legend = c("one", "two", "three"),
  title = "Hypotheses supported",
  bty = "n",
  lty = 1:3,
  col = my_cols[1:3],
  seg.len = 4.5,
  lwd = 5
)
dev.off()

# and now distributional area

range(hm$log_eoo)
eoo_dat <- data.frame(
  log_mass = 0,
  log_eoo = seq(5,18,length.out = 500),
  su = 0
)
mds <- eoo_dat
mds$log_eoo <- (mds$log_eoo - mean(hm$log_eoo)) / sd(hm$log_eoo)


eoo_pred <- predict(
  longshot,
  newdata = mds,
  type = "probs",
  se.fit = TRUE
)
nboot <- 1000
boot_array <- array(
  NA, dim = c(nrow(eoo_dat), 3, nboot)
)
pb <- txtProgressBar(max = nboot)
for(i in 1:nboot){
  setTxtProgressBar(pb,i)
  tmp_dat <- yo[
    sample(1:nrow(yo), replace = TRUE),
  ]
  l2 <- polr(
    nsupport ~ log_mass + log_eoo + su,
    data = tmp_dat
  )
  eoo_pred <- predict(
    l2,
    newdata = mds,
    type = "probs"
  )
  boot_array[,,i] <- eoo_pred
}
boot_est <- apply(
  boot_array,
  c(1,2),
  quantile,
  probs = c(0.025,0.5,0.975)
)


# make some pretty x axis labels
my_xvals <- c(1000, 10000, 500000,50000000)
my_xtitles <- paste0(
  my_xvals/1000
)


tiff("eoo_support.tiff",
     height = 4, width = 6,
     units = "in", res = 200,
     compression = "lzw")
par(mar = c(4,4,1,1.25))
bbplot::blank(
  xlim = c(5,18),
  ylim = c(0,1),
  bty = "l",
  xaxs = "i",
  yaxs = "i"
)
bbplot::axis_blank(1, at = log(my_xvals), minor = FALSE)
bbplot::axis_blank(2)
bbplot::axis_text(my_xtitles,side = 1, line = 0.5,at = log(my_xvals), )
bbplot::axis_text(side = 2, line = 0.5, las = 1)
mtext(
  expression(
    paste("Extent of species distribution (km"^"2", ")")
  ),
  side = 1,
  line= 2.5,
  cex = 1.25,
  at = mean(c(5,18))
)


bbplot::axis_text(
  "Probability",
  side = 2,
  line = 2.5,
  cex = 1.25
)
my_cols <-  c("#3B3B3B", "#A73030", "#0073C2")

for(i in 1:3){
  bbplot::ribbon(
    x = eoo_dat$log_eoo,
    y = t(boot_est[-2,,i]),
    col = my_cols[i],
    alpha = 0.25
  )
}
for(i in 1:3){
  
  lines(
    y = boot_est[2,,i],
    x = eoo_dat$log_eoo,
    col = my_cols[i],
    lty = i,
    lwd = 5
  )
}
u <- par("usr")
lines(x = c(u[1:2]), y = c(0,0))
lines(x = c(u[1],u[1]), y = c(0,1))
legend(
  "topright",
  legend = c("one", "two", "three"),
  title = "Hypotheses supported",
  bty = "n",
  lty = 1:3,
  col = my_cols[1:3],
  seg.len = 4.5,
  lwd = 5
)
dev.off()

#### all types multi-plot ####

the_files <- c(
  "D:/diel_data/uniform_prior/Uniform Prior on Hyp Set/Aggregated/Traditionalspecies.hyp.ref.AGGREGATED.csv",
  "D:/diel_data/uniform_prior/Uniform Prior on Hyp Set/Aggregated/Generalspecies.hyp.ref.AGGREGATED.csv",
  "D:/diel_data/uniform_prior/Uniform Prior on Hyp Set/Aggregated/Maximizingspecies.hyp.ref.AGGREGATED.csv",
  "D:/diel_data/uniform_prior/Uniform Prior on Hyp Set/All Analysis Units/Traditionalspecies.hyp.ref.csv",
  "D:/diel_data/uniform_prior/Uniform Prior on Hyp Set/All Analysis Units/Generalspecies.hyp.ref.csv",
  "D:/diel_data/uniform_prior/Uniform Prior on Hyp Set/All Analysis Units/Maximizingspecies.hyp.ref.csv"
)


jj <- read.csv("")

trad <- lapply(
  the_files,
  read.csv
)
names(trad) <- c(
  "aggregated_traditional",
  "aggregated_general",
  "aggregated_maximize",
  "all_units_traditional",
  "all_units_general",
  "all_units_maximize"
)
trad <-read.csv(
  "D:/diel_data/uniform_prior/Uniform Prior on Hyp Set/All Analysis Units/Traditionalspecies.hyp.ref.csv"
)

# read in sample size
ss <- read.csv(
  "D:/diel_data/total.sample.size.by.analysis.unit.csv"
)
# add on sample size non-aggregated
for(i in 4:6){
  trad[[i]]$samp_size <- ss$ss
}

# trad$ref_prob <- ifelse(
#   trad$Hypothesis == trad$Ref.Activity,
#   trad$P.Hypothesis.,
#   1 - trad$P.Hypothesis.
# )

# get average for the non-aggregated
for(i in 4:6){
  tmp <- trad[[i]]
  tmp <- tmp[complete.cases(tmp),]
  tmp <- tmp %>% 
    dplyr::group_by(species, hyp.ref) %>%
    dplyr::mutate(
      samp_prob = samp_size / sum(samp_size) * prob 
    ) %>% 
    data.frame()
  
  tmp <- tmp %>% 
    dplyr::group_by(species, hyp.ref) %>% 
    dplyr::summarise(
      samp_prob = sum(samp_prob)
    )
  trad[[i]] <- tmp
  
}



# get the average 

trad_mean <- vector(
  "list",
  length = length(trad)
)

for(i in 1:6){
  tmp <- trad[[i]]
  tmp <- tmp[complete.cases(tmp),]
  if(i %in% 1:3){
    trad_mean[[i]] <- tmp %>% 
      dplyr::group_by(hyp.ref) %>% 
      dplyr::summarise(
        n = length(prob),
        mean_prob = mean(prob),
        sd_prob = sd(prob)
      ) %>% 
      data.frame()
  }else{
      trad_mean[[i]] <- tmp %>% 
        dplyr::group_by(hyp.ref) %>% 
        dplyr::summarise(
          n = length(samp_prob),
          mean_prob = mean(samp_prob),
          sd_prob = sd(samp_prob)
        ) %>% 
        data.frame()
  }

}

windows(9,9)

m <- matrix(
  c(0,0,0,
    rep(1:3,3),
    rep(4:6,3)
    ),
  ncol =7,
  nrow = 3
)
layout(m)
par(
  mar = c(4,1,2,1),
  lend = 2
)
for(k in 1:6){

  {
    bbplot::blank(
      xlim = c(0,1),
      ylim = c(0.5,4.5),
      bty = "l",
      yaxs = "i"
    )
    bbplot::axis_text(
      names(trad)[k],
      line = 1
    )
    bbplot::axis_blank(1)
    bbplot::axis_blank(
      2,
      at = 1:4,
      minor = FALSE
    )
    bbplot::axis_text(side = 1, line = 1)
    if(k %in% c(3,6)){
    bbplot::axis_text(
      "Pr(Reference diel modality)",
      side = 1,
      line = 2.25,
      cex = 1.15
    )
    }
    if(k %in% 1:3){
     bbplot::axis_text(
       rev(trad_mean[[k]]$hyp.ref),
       at = 1:length(trad_mean[[k]]$hyp.ref),
       las = 1,
       side = 2,
       cex = 1.15,
       line = .5
     )
    }
  }
  # add number to test
  my_cols <- c("#3B3B3BFF", "#A73030FF", "#EFC000FF", "#0073C2FF")
  my_levs <- c("Nocturnal", "Diurnal", "Crepuscular",
               "Cathemeral")
  trad[[k]]$y <- as.numeric(
    factor(
      trad[[k]]$hyp.ref,
      levels = my_levs
    )
  )
  trad[[k]]$col <- my_cols[5 - trad[[k]]$y]
  set.seed(3)
  trad[[k]]$y <- jitter(trad[[k]]$y, 1.75)
  
  if("prob" %in% colnames(trad[[k]])){
    trad[[k]]$samp_prob <- trad[[k]]$prob
  }
  points(
    x = trad[[k]]$samp_prob,
    y = trad[[k]]$y,
    pch = 19,
    col =  scales::alpha(trad[[k]]$col, 0.5)
  )
  
  
  for(i in 1:4){
    if(k %in% c(3,6)){
      if(i == 1){
        next
      }
    }
    j <- 5 - i
    lev_loc <- which(
      trad_mean[[k]]$hyp.ref == rev(my_levs)[i]
    )
    my_lo <- trad_mean[[k]]$mean_prob[lev_loc] - trad_mean[[k]]$sd_prob[lev_loc]
    my_hi <- trad_mean[[k]]$mean_prob[lev_loc] + trad_mean[[k]]$sd_prob[lev_loc]
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
      x = trad_mean[[k]]$mean_prob[i],
      y = j,
      pch = 19,
      col = "white",
      cex = 2.5
    )
    points(
      x = trad_mean[[k]]$mean_prob[i],
      y = j,
      pch = 19,
      col = my_cols[i],
      cex = 1.5
    )
  }
}

#### plasticity 2 #####



trad <- read.csv(
  "D:/diel_data/Traditionalspecies.hyps.df.csv"
)

hm <- trad[trad$P.Hypothesis.>0.8,] %>% 
  dplyr::group_by(Species, Ref.Activity) %>% 
  dplyr::summarise(
    Diurnal = sum(Hypothesis == "Diurnal"),
    Cathemeral = sum(Hypothesis == "Cathemeral"),
    Crepuscular = sum(Hypothesis == "Crepuscular"),
    Nocturnal = sum(Hypothesis == "Nocturnal"),
    n_units = length(Species)
  ) %>% 
  data.frame()
hm$Hypothesis <- factor(hm$Hypothesis)

# tack backon species details
hm <- dplyr::inner_join(
  hm,
  dplyr::distinct(trad[,c("Species", "Body_mass..g.")]),
  by = "Species"
)

# read in distrib range
dr <- read.csv(
  "./data/analysis_units/species_distrib_range.csv"
)

hm <- dplyr::inner_join(
  hm,
  dr,
  by = c("Species" = "scientificName")
)
hm <- data.frame(hm)

hm$log_mass <- log(hm$Body_mass..g.)
hm$log_eoo <- log(hm$EOO)

yo <- hm


yo$log_eoo <- as.numeric(scale(yo$log_eoo))

yo$log_mass <- as.numeric(scale(yo$log_mass))

my_response <- as.matrix(
  yo[,c("Cathemeral", "Diurnal", "Nocturnal")]
)
yo$Hypothesis <- factor(
  yo$Hypothesis
)
longshot <- multinom(
  Hypothesis ~ log_mass + log_eoo,
  data = yo
)


ack <- data.frame(
  log_mass = 
    seq(
      -3.83,
      2.94,
      length.out = 100
    ),
  log_eoo = 0
)

my_x <- ack$log_mass
batt <- predict(
  longshot, newdata = ack, "probs"
)
plot(batt[,1] ~ my_x, type = "l", ylim = c(0,1),
     xlab = "scaled log body mass", lwd = 2,
     bty = "l", las = 1, ylab = "Probability of diel modality")
lines(batt[,2] ~ my_x, lty = 2,lwd = 2)
lines(batt[,3] ~ my_x, lty = 3, lwd = 2)

legend("topleft", legend = levels(yo$Hypothesis),
       lty = 1:3, bty = "n", lwd = 2)

yo$su <- as.numeric(scale(yo$n_units))
longshot <- polr(
  nsupport ~ log_mass + log_eoo + su,
  data = yo
)

cli <- confint(longshot)

# make predictions for body mass

range(hm$log_mass)
mass_dat <- data.frame(
  log_mass = seq(2,15,length.out = 500),
  log_eoo = 0,
  su = 0
)
mds <- mass_dat
mds$log_mass <- (mds$log_mass - mean(hm$log_mass)) / sd(hm$log_mass)


mass_pred <- predict(
  longshot,
  newdata = mds,
  type = "probs",
  se.fit = TRUE
)

nboot <- 1000
boot_array <- array(
  NA, dim = c(nrow(mass_dat), 3, nboot)
)
for(i in 1:nboot){
  tmp_dat <- yo[
    sample(1:nrow(yo), replace = TRUE),
  ]
  l2 <- polr(
    nsupport ~ log_mass + log_eoo + su,
    data = tmp_dat
  )
  mass_pred <- predict(
    l2,
    newdata = mds,
    type = "probs",
    se.fit = TRUE
  )
  boot_array[,,i] <- mass_pred
}
boot_est <- apply(
  boot_array,
  c(1,2),
  quantile,
  probs = c(0.025,0.5,0.975)
)

# make some pretty x axis labels
my_xvals <- c(10, 100, 1000, 10000, 500000)
my_xtitles <- paste0(
  "ln(", my_xvals/1000,")"
)


windows(6, 4)
tiff("body_mass_support.tiff",
     height = 4, width = 6,
     units = "in", res = 200,
     compression = "lzw")
par(mar = c(4,4,1,1))
bbplot::blank(
  xlim = c(2,15),
  ylim = c(0,1),
  bty = "l",
  xaxs = "i",
  yaxs = "i"
)
bbplot::axis_blank(1, at = log(my_xvals), minor = FALSE)
bbplot::axis_blank(2)
bbplot::axis_text(my_xtitles,side = 1, line = 0.5,at = log(my_xvals), )
bbplot::axis_text(side = 2, line = 0.5, las = 1)
bbplot::axis_text(
  "Body mass (kg)",
  side = 1,
  line = 2.5,
  cex = 1.25
)
bbplot::axis_text(
  "Probability",
  side = 2,
  line = 2.5,
  cex = 1.25
)
my_cols <-  c("#3B3B3B", "#A73030", "#0073C2")

for(i in 1:3){
  bbplot::ribbon(
    x = mass_dat$log_mass,
    y = t(boot_est[-2,,i]),
    col = my_cols[i],
    alpha = 0.25
  )
}
for(i in 1:3){
  
  lines(
    y = boot_est[2,,i],
    x = mass_dat$log_mass,
    col = my_cols[i],
    lty = i,
    lwd = 5
  )
}
u <- par("usr")
lines(x = c(u[1:2]), y = c(0,0))
lines(x = c(u[1],u[1]), y = c(0,1))
legend(
  "topright",
  legend = c("one", "two", "three"),
  title = "Hypotheses supported",
  bty = "n",
  lty = 1:3,
  col = my_cols[1:3],
  seg.len = 4.5,
  lwd = 5
)
dev.off()

# and now distributional area

range(hm$log_eoo)
eoo_dat <- data.frame(
  log_mass = 0,
  log_eoo = seq(5,18,length.out = 500),
  su = 0
)
mds <- eoo_dat
mds$log_eoo <- (mds$log_eoo - mean(hm$log_eoo)) / sd(hm$log_eoo)


eoo_pred <- predict(
  longshot,
  newdata = mds,
  type = "probs",
  se.fit = TRUE
)
nboot <- 1000
boot_array <- array(
  NA, dim = c(nrow(eoo_dat), 3, nboot)
)
pb <- txtProgressBar(max = nboot)
for(i in 1:nboot){
  setTxtProgressBar(pb,i)
  tmp_dat <- yo[
    sample(1:nrow(yo), replace = TRUE),
  ]
  l2 <- polr(
    nsupport ~ log_mass + log_eoo + su,
    data = tmp_dat
  )
  eoo_pred <- predict(
    l2,
    newdata = mds,
    type = "probs"
  )
  boot_array[,,i] <- eoo_pred
}
boot_est <- apply(
  boot_array,
  c(1,2),
  quantile,
  probs = c(0.025,0.5,0.975)
)


# make some pretty x axis labels
my_xvals <- c(1000, 10000, 500000,50000000)
my_xtitles <- paste0(
  "ln(", my_xvals/1000,")"
)


tiff("eoo_support.tiff",
     height = 4, width = 6,
     units = "in", res = 200,
     compression = "lzw")
par(mar = c(4,4,1,1.25))
bbplot::blank(
  xlim = c(5,18),
  ylim = c(0,1),
  bty = "l",
  xaxs = "i",
  yaxs = "i"
)
bbplot::axis_blank(1, at = log(my_xvals), minor = FALSE)
bbplot::axis_blank(2)
bbplot::axis_text(my_xtitles,side = 1, line = 0.5,at = log(my_xvals), )
bbplot::axis_text(side = 2, line = 0.5, las = 1)
mtext(
  expression(
    paste("Extent of species distribution (km"^"2", ")")
  ),
  side = 1,
  line= 2.5,
  cex = 1.25,
  at = mean(c(5,18))
)


bbplot::axis_text(
  "Probability",
  side = 2,
  line = 2.5,
  cex = 1.25
)
my_cols <-  c("#3B3B3B", "#A73030", "#0073C2")

for(i in 1:3){
  bbplot::ribbon(
    x = eoo_dat$log_eoo,
    y = t(boot_est[-2,,i]),
    col = my_cols[i],
    alpha = 0.25
  )
}
for(i in 1:3){
  
  lines(
    y = boot_est[2,,i],
    x = eoo_dat$log_eoo,
    col = my_cols[i],
    lty = i,
    lwd = 5
  )
}
u <- par("usr")
lines(x = c(u[1:2]), y = c(0,0))
lines(x = c(u[1],u[1]), y = c(0,1))
legend(
  "topright",
  legend = c("one", "two", "three"),
  title = "Hypotheses supported",
  bty = "n",
  lty = 1:3,
  col = my_cols[1:3],
  seg.len = 4.5,
  lwd = 5
)
dev.off()





# read in phylogentic tree
faurby_tree<-read.nexus(
  "D:/GIS/phylacine/DataS1/Complete_phylogeny.nex"
)
# get one of the trees
faurby_tree <- faurby_tree[[1]]



names <- data.frame(Species=tp$species)
names <- dplyr::distinct(names)
names$Species <- gsub(
  " ",
  "_",
  names$Species
)

rownames(names) <- names$Species

# subset down to the species of interest
overlap_species<-geiger::name.check(
  faurby_tree,
  names
)

mammal_tree <- ape::drop.tip(
  faurby_tree,
  overlap_species$tree_not_data
)


# Format trait data for creating phylogeny
tmat <- matrix(
  tp$var,
  ncol = 1,
  nrow = length(mammal_tree$tip.label)
)
row.names(tmat) <-names$Species


## Using Pagel's Lambda for continuous traits ###
signal_there <- phytools::phylosig(
  mammal_tree,
  tmat,
  method="lambda",
  test=TRUE
)
signal_there

### Now use contmap to plot continuous traits ###
traitnewvector <- tmat[,1]
names(traitnewvector) <- row.names(tmat)
### Crappy Tree ###
obj <- phytools::contMap(
  mammal_tree,
  traitnewvector,
  fsize=c(0.2,1),
  lwd=1,
  plot = FALSE
)
obj<-phytools::setMap(
  obj,
  colors=c(
    "firebrick4", "brown2", "sienna2", "lightgoldenrod1",
    "seagreen3", "deepskyblue3", "darkslateblue"
  )
)
windows(15,15)
pdf("phy_signal_uniform.pdf", height = 15, width = 15)
svg(
  "phy_signal_uniform.svg",
  height = 15,
  width = 15
)
tiff(
  "phy_signal_uniform.tiff",
  height = 15,
  width = 15,
  res = 200,
  units = "in",
  compression = "lzw"
)

plot(obj, type="fan", lwd=4, outline=T, fsize=c(1, 1), res = 600)
dev.off()
