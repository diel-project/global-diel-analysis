# Modeling and plotting of the analysis unit level model probabilities for each 
# species' reference/literature diel categorization

# We want to estimate the weighted (by sample size) species-level P(ref hypothesis)
# from each analysis unit while controlling for family and project. Specifically, we want species 
# nested within family nested in project. 

# First, We want to do this separately for each species reference categorization (diurnal, etc.) to look
# at species level variation within categorization.

# Second, We want to do this jointly for all reference categorization to look at family-wise agreement

########################################################
# Before modeling the data, here is an example of using a linear model function to do a weighted regression
# Specifically, where the weights reflect the sample size and what its weighting is the 
# observations within the likelihood. 

#Example probabilities and sample sizes
  diurnal = c(0.2,0.2,0.8,0.8,0.8)
  n = c(2,2,10,10,10)

# When not weighting by sample size the mean is perhaps not realistic, given the sample
# size differences and probability differences by sample size
  mean(diurnal)

#Get weighted mean via weighted.mean function (R Stats package) and then use lm (linear model)
  weighted.mean(diurnal, w=n)
  coef(lm(diurnal~1, weights=n))
#We get the same result

#Do this again with very different sample sizes
  n=c(2,2,100,100,100)

#These approaches provide the same sensible result
  weighted.mean(diurnal, w=n)
  coef(lm(diurnal~1, weights=n))
########################################################
# Get data and fit models
  
# Load libraries
  library(dplyr)
  library(brms)
  library(ordbetareg)
  #Install package if needed
  #devtools::install_github("dapperstats/bbplot")  
  library(bbplot)

# We need to fit this code for two datasets- where the model probabilities for each
# analysis unit was provided a uniform prior probabilities on each hypothesis within 
# the hypothesis set or whether the reference diel categorization was provided more weight. 
    
  
##################################  
# Go grab data needed for our analysis. Make this TRUE or FALSE to switch datasets
  model.probs.uniform=TRUE
  
  if(model.probs.uniform){
  #uniform data  
    dat <- read.csv(
      "./data/diel.niche_results/uniform/Traditionalspecies.hyp.ref.csv"
    )
  }else{
  # informed model probabilities
    dat <- read.csv(
      "./data/diel.niche_results/informed/Traditionalspecies.hyp.ref.PRIOR90.csv"
    )
  }
  
  head(dat)
  dim(dat)
  
##################################  
# Get additional information to append

  info <- read.csv(
    "./data/Traditional.species.hyps.with.analysis.units.csv"
  )

  dim(info)
  head(info)

# Check if info and dat are in the same order-TEST
  all(dat$species==info$species)
  which(dat$species=="Tupaia belangeri")
  which(info$species=="Tupaia belangeri")

# Combine the data
  dat2=data.frame(dat,info)
  head(dat2)

#Bring in the diel_data set that has additional information we need
  diel.data <- read.csv(
    "./data/analysis_units/diel_data.csv"
  )
  dim(diel.data)

#  Test if its in the same order - no it is not
  all(diel.data$scientificName==dat2$species)
  
#Join the two datasets by species and analysis unit number  
  dat3 <- dplyr::full_join(
    dat2,
    diel.data,
    by = c("analysis_unit", "species" = "scientificName")
  )
  
  head(dat3)  

# Create sample size column
  dat3$sample.size=dat3$twilight+dat3$day+dat3$night

# How many do we have?
  length(unique(dat3$species))
  length(unique(dat3$family.x))
  length(unique(dat3$file_name))


#Split data into diel reference categories- need to fit a model to each
  dat3.diurnal=dat3[which(dat3$hyp.ref=="Diurnal"),]
  dat3.noct=dat3[which(dat3$hyp.ref=="Nocturnal"),]
  dat3.crep=dat3[which(dat3$hyp.ref=="Crepuscular"),]
  dat3.cath=dat3[which(dat3$hyp.ref=="Cathemeral"),]

# For our modeling choices, we would like to 1) use a probability density function that 
# appropriately models the [0,1] data, 2) that allows nested random effects, and 
# 3) allows likelihood weighting. We can do all three using Ordered beta regression (Kubinec et al. 2023).
# This model is implemented in the R package ordbetareg, which provides wrapper functions to fit the relevanat
# model in brms
  
#  Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting model for continuous data with lower and upper bounds. 
#  Political Analysis, 31(4), 519-536.
  
# Fit oredered beta regression model, weighted by sample size with a nested random effects structure.
# Note that, an intercept only model leads to a coding error. However, we can get around this issue
# by specificing a zeroed out covariate. See the Github Issue and response by the package developer here,
# https://github.com/saudiwin/ordbetareg_pack/issues/21

# Create dummy covaraite  
  dat3.diurnal$zero_cov=0
  dat3.noct$zero_cov=0
  dat3.cath$zero_cov=0
  dat3$zero_cov=0
  
#Diurnal
    brm.diurnal.ordbeta <- ordbetareg(prob|weights(sample.size) ~ 1 + zero_cov  + (1|family.x/species/file_name), 
                                    data=dat3.diurnal,
                                    cores=4,
                                    chains=4,
                                    iter=10000
                                    ) 
    if(model.probs.uniform){filename="brm.diurnal.ordbeta"}else{filename="brm.diurnal.informed.ordbeta"}
    save(brm.diurnal.ordbeta,file=filename)

  
#Nocturnal  
  brm.noct.ordbeta <- ordbetareg(prob|weights(sample.size) ~ 1 + zero_cov + (1|family.x/species/file_name), 
                                 data=dat3.noct,
                                 cores=4,
                                 chains=4,
                                 iter=10000
                                 ) 
  if(model.probs.uniform){filename="brm.noct.ordbeta"}else{filename="brm.noct.informed.ordbeta"}
  save(brm.noct.ordbeta,file=filename)

#Cathemeral  
  brm.cath.ordbeta <- ordbetareg(prob|weights(sample.size) ~ 1 + zero_cov + (1|family.x/species/file_name), 
                                 data=dat3.cath,
                                 cores=4,
                                 chains=4,
                                 iter=10000) 
  if(model.probs.uniform){filename="brm.cath.ordbeta"}else{filename="brm.cath.informed.ordbeta"}
  save(brm.cath.ordbeta,file=filename)
  
  
#Crepuscular
  #Fitting these data are not informative as 99% are 0
  table(round(dat3.crep$prob,digits=2))
  #The species-level predictions will all be zero
  
    
#Examine traceplots    
  plot(brm.diurnal.ordbeta)
  plot(brm.noct.ordbeta)
  plot(brm.cath.ordbeta)

###################################################  
###################################################
# Predict species specific probabilities for each reference categorization
# We will do this by marginalizing over the predictions for each species 
# and across projects. 

#Get full samples of predictions
  pred.diurnal=predict(brm.diurnal.ordbeta,type="response",summary=FALSE)
  pred.nocturnal=predict(brm.noct.ordbeta,type="response",summary=FALSE)
  pred.cath=predict(brm.cath.ordbeta,type="response",summary=FALSE)

  dim(pred.diurnal)#mcmc by sample size

#Define species as columns
  colnames(pred.diurnal)=dat3.diurnal$species
  colnames(pred.nocturnal)=dat3.noct$species
  colnames(pred.cath)=dat3.cath$species

#Storage space for each species-prediction (averaged across projects)
  species.diurnal.preds=matrix(NA, 
                                   nrow=nrow(pred.diurnal),
                                   ncol=length(unique(dat3.diurnal$species))
                               )
  species.nocturnal.preds=matrix(NA, 
                               nrow=nrow(pred.nocturnal),
                               ncol=length(unique(dat3.noct$species))
  )
  species.cath.preds=matrix(NA, 
                               nrow=nrow(pred.cath),
                               ncol=length(unique(dat3.cath$species))
  )
  
#Loop across species and get marginalized predictions
#Diurnal  
  for(i in 1:length(unique(dat3.diurnal$species))){
    index=which(colnames(pred.diurnal)==unique(dat3.diurnal$species)[i])
    if(length(index)>1){
      species.level.mean=apply(pred.diurnal[,index],1,median)
    }else{
      species.level.mean=pred.diurnal[,index]
    }
    species.diurnal.preds[,i]=species.level.mean
  }
#Nocturnal  
  for(i in 1:length(unique(dat3.noct$species))){
    index=which(colnames(pred.nocturnal)==unique(dat3.noct$species)[i])
    if(length(index)>1){
      species.level.mean=apply(pred.nocturnal[,index],1,median)
    }else{
      species.level.mean=pred.nocturnal[,index]
    }
    species.nocturnal.preds[,i]=species.level.mean
  }
#Cathemeral
  for(i in 1:length(unique(dat3.cath$species))){
    index=which(colnames(pred.cath)==unique(dat3.cath$species)[i])
    if(length(index)>1){
      species.level.mean=apply(pred.cath[,index],1,median)
    }else{
      species.level.mean=pred.cath[,index]
    }
    species.cath.preds[,i]=species.level.mean
  }
  

  colnames(species.diurnal.preds)=unique(dat3.diurnal$species)
  colnames(species.nocturnal.preds)=unique(dat3.noct$species)
  colnames(species.cath.preds)=unique(dat3.cath$species)

#What is the posterior mean for each species
  post.mean.species.diurnal = data.frame(apply(species.diurnal.preds,2,mean))
  colnames(post.mean.species.diurnal)="Post.Mean"
  post.mean.species.diurnal = round(post.mean.species.diurnal,digits=4)
  
  post.mean.species.nocturnal = data.frame(apply(species.nocturnal.preds,2,mean))
  colnames(post.mean.species.nocturnal)="Post.Mean"
  post.mean.species.nocturnal = round(post.mean.species.nocturnal,digits=4)
  
  post.mean.species.cath = data.frame(apply(species.cath.preds,2,mean))
  colnames(post.mean.species.cath)="Post.Mean"
  post.mean.species.cath = round(post.mean.species.cath,digits=4)
  
  
  
#What is the mean across species - this is the marginal mean from the sample

  dirunal.across.species.average=apply(species.diurnal.preds,1,mean)
  dirunal.across.species.sd=apply(species.diurnal.preds,1,sd)
  mean(dirunal.across.species.average)
  
  #SE of the mean
  sd(dirunal.across.species.average)
  
  mean(dirunal.across.species.sd)
  
    
  noct.across.species.average=apply(species.nocturnal.preds,1,mean)
  noct.across.species.sd=apply(species.nocturnal.preds,1,sd)
  mean(noct.across.species.average)
  mean(noct.across.species.sd)
  
  cath.across.species.average=apply(species.cath.preds,1,mean)
  cath.across.species.sd=apply(species.cath.preds,1,sd)
  mean(cath.across.species.average)
  mean(cath.across.species.sd)
  
#This is the conditional posterior mean- across all species, families, and projects
# We won't be using this type of mean, prefering the marginal mean 
  diurnal.samples=as_draws(brm.diurnal.ordbeta)
  hist(diurnal.samples[[1]]$b_Intercept)
  mean(diurnal.samples[[1]]$b_Intercept)
  quantile(diurnal.samples[[1]]$b_Intercept,probs = c(0.025,0.975))

  noct.samples=as_draws(brm.noct.ordbeta)
  hist(noct.samples[[1]]$b_Intercept)
  mean(noct.samples[[1]]$b_Intercept)
  quantile(noct.samples[[1]]$b_Intercept,probs = c(0.025,0.975))
  
  
  cath.samples=as_draws(brm.cath.ordbeta)
  hist(cath.samples[[1]]$b_Intercept)
  mean(cath.samples[[1]]$b_Intercept)
  quantile(cath.samples[[1]]$b_Intercept,probs = c(0.025,0.975))
  
#########################################################  
# Organize information together

  diurnal.output=list(species=post.mean.species.diurnal,
                      mean.sp=mean(dirunal.across.species.average),
                      SE.mean.sp=sd(dirunal.across.species.average),
                      SD.across.sp=mean(dirunal.across.species.sd)
                      )
  nocturnal.output=list(species=post.mean.species.nocturnal,
                      mean.sp=mean(noct.across.species.average),
                      SE.mean.sp=sd(noct.across.species.average),
                      SD.across.sp=mean(noct.across.species.sd)
  )
  
  cathemeral.output=list(species=post.mean.species.cath,
                        mean.sp=mean(cath.across.species.average),
                        SE.mean.sp=sd(cath.across.species.average),
                        SD.across.sp=mean(cath.across.species.sd)
  )
  
# Additional setup for plotting in bbplot
  trad = data.frame(species = c(rownames(cathemeral.output$species),
                                unique(dat3.crep$species),
                                rownames(diurnal.output$species),
                                rownames(nocturnal.output$species)
                                ),
                    hyp.ref = c(rep("Cathemeral",length(unique(dat3.cath$species))),
                                rep("Crepuscular",length(unique(dat3.crep$species))),
                                rep("Diurnal",length(unique(dat3.diurnal$species))),
                                rep("Nocturnal",length(unique(dat3.noct$species)))
                              ),
                    samp_prob = c(cathemeral.output$species$Post.Mean,
                                  rep(0,length(unique(dat3.crep$species))),
                                  diurnal.output$species$Post.Mean,
                                  nocturnal.output$species$Post.Mean
                                  )
                  )
  
  
  
  trad_mean = data.frame(hyp.ref=c("Cathemeral","Crepuscular","Diurnal","Nocturnal"),
                        n = c(length(unique(dat3.cath$species)),
                              length(unique(dat3.crep$species)),
                              length(unique(dat3.diurnal$species)),
                              length(unique(dat3.noct$species))
                              ),
                        mean_prob = c(
                                      mean(cath.across.species.average),
                                      0,
                                      mean(dirunal.across.species.average),
                                      mean(noct.across.species.average)
                                      ),
                        sd_prob = c(
                                    mean(cath.across.species.sd),
                                    0.001,
                                    mean(dirunal.across.species.sd),
                                    mean(cath.across.species.sd)
                          
                        )
                        )

###############################################
# Plotting
  if(model.probs.uniform){filename="./plots/uniform_prior_overall_accuracy.model.based2.tiff"}else{filename="./plots/informed_prior_overall_accuracy.model.based2.tiff"}
  #windows(4.7,3)
  tiff(
    filename =filename,
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
  
  
################################################################
################################################################
################################################################
################################################################
  
# Second, we fit a model for all reference categories to investigate species
# within family
  brm.family.ordbeta <- brm(formula=prob|weights(sample.size) ~ 1 + zero_cov+ (1|family.x/species/file_name),
                            data=dat3, cores=4,chains=4,iter=10000)
  if(model.probs.uniform){filename="brm.family.ordbeta"}else{filename="brm.family.informed.ordbeta"}
  save(brm.family.ordbeta,file=filename)


#Predict posteriors for each analysis unit
  pred.family=predict(brm.family.ordbeta,type="response",summary=FALSE)


# We need to marginalize species across projects and then across species within family
# Loop over all species and get average across projects
  
  species.preds =   matrix(NA, 
                           nrow=nrow(pred.family),
                           ncol=length(unique(dat3$species))
  )
  colnames(species.preds)=unique(dat3$species)
  
  colnames(pred.family)=dat3$species
  
  for(i in 1:length(unique(dat3$species))){
    index = which(colnames(pred.family)==unique(dat3$species)[i])
    if(length(index)>1){
      species.level.mean=apply(pred.family[,index],1,median)
    }else{
      species.level.mean=pred.family[,index]
    }
    species.preds[,i]=species.level.mean
  }

# Now, we need to loop over all species to get family-level means
  colnames(species.preds)
  hist(species.preds[,1])
  
#Get family name for each species
   these=match(unique(dat3$species),dat3$species)
   familiy.by.species=dat3$family.x[these]
  
  
# posterior means by species
  post.sp.mean=data.frame(species=colnames(species.preds),
                          prob.mean=apply(species.preds,2,mean),
                          family=familiy.by.species
                          )


#Next, we need to get the family posterior mean by marginalizing over the species within each family
  
  family.preds=data.frame(matrix(NA, 
                                 nrow=length(unique(dat3$family.x)),
                                 ncol=4
                                 )
                          )
                      
  colnames(family.preds)=c("family_mu","family_sd","n_species","family")
                 
  family.preds.temp=species.preds
  colnames(family.preds.temp)=familiy.by.species
                       
    
  for(i in 1:length(unique(colnames(family.preds.temp)))){
             index=which(colnames(family.preds.temp)==unique(colnames(family.preds.temp))[i])
             family.preds[i,3]=length(index)
             family.preds[i,4]=unique(colnames(family.preds.temp))[i]
           if(length(index)>1){
             family.preds.mean=mean(apply(family.preds.temp[,index],2,mean))
             family.sd=sd(apply(family.preds.temp[,index],2,mean))
           }else{
                 family.preds.mean=mean(family.preds.temp[,index])
                 family.sd=NA
                }
            family.preds[i,1]=family.preds.mean
            family.preds[i,2]=family.sd           
    }

    
# We need to drop families that have <5 species
     index.keep=which(family.preds$n_species>=5)
     
     family.preds.limited=family.preds[index.keep,]
     colnames(family.preds.limited)
     
#Which families were dropped     
     unique(familiy.by.species)[-index.keep]

#for the species preds - only keep the families that have enough species     
     head(post.sp.mean)
     keep.these.species.rows=which(!is.na(match(post.sp.mean$family,family.preds.limited$family)))
     post.sp.mean.limited=post.sp.mean[keep.these.species.rows,]               
     
     #should be 26
     length(unique(post.sp.mean.limited$family))
   
     
###############################
# Organize data for plotting below     
   family.preds.limited = family.preds.limited[order(family.preds.limited$family_mu,decreasing = TRUE),]
   fam.post.mean = family.preds.limited 
     
################################
# plot family
    if(model.probs.uniform){filename="./plots/uniform_prior_family_accuracy.model.based2.tiff"}else{filename="./plots/informed_prior_family_accuracy.model.based2.tiff"}
     
     #windows(4.5,8)
     tiff(
       filename = filename,
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
         rev(fam.post.mean$family),
         at = 1:26,
         las = 1,
         side = 2,
         cex = 1.15,
         line = .5
       )
       bbplot::axis_text(
         rev(
           paste0("n = ", fam.post.mean$n_species)
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
     post.sp.mean <- post.sp.mean[post.sp.mean$family %in% fam.post.mean$family,]
     post.sp.mean$y <- as.numeric(
       factor(
         post.sp.mean$family,
         levels = rev(fam.post.mean$family)
       )
     )
     
     post.sp.mean$col <- my_cols[(post.sp.mean$y %%2)+1]
     set.seed(3)
     post.sp.mean$y <- jitter(post.sp.mean$y, 1.75)
     
     
     points(
       x = post.sp.mean$prob.mean,
       y = post.sp.mean$y,
       pch = 19,
       col =  scales::alpha(post.sp.mean$col, 0.5)
     )
     
     fam.post.mean$col <- rep(my_cols, 13)
     for(i in 1:26){
       j <- 27 - i
       my_lo <- fam.post.mean$family_mu[i] - fam.post.mean$family_sd[i]
       my_hi <- fam.post.mean$family_mu[i] + fam.post.mean$family_sd[i]
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
         col = fam.post.mean$col[i]
       )
       points(
         x = fam.post.mean$family_mu[i],
         y = j,
         pch = 19,
         col = "white",
         cex = 2.5
       )
       points(
         x = fam.post.mean$family_mu[i],
         y = j,
         pch = 19,
         col = fam.post.mean$col[i],
         cex = 1.5
       )
     }
 dev.off()
