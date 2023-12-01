################################################################
# Function for species diel niche modeling - analysis unit data
################################################################
# Install package (if needed)

# Install Kaleido
# install.packages('reticulate')
# reticulate::install_miniconda()
 reticulate::conda_install('r-reticulate', 'python-kaleido=')
 
 reticulate::conda_install(packages = "python-kaleido==0.1.0post1")
 
# reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
# reticulate::use_miniconda('r-reticulate')
################################################################


p <- plot_ly(x = 1:10)
tmp <- "test.png"
save_image(p, tmp)
file.show(tmp)

# Load libraries

library(Diel.Niche)
library(fs)
#library(plotly)
#library(reticulate)
#reticulate::use_miniconda('r-reticulate')
# Work around issue with sys not being imported 
# https://stackoverflow.com/questions/73604954/error-when-using-python-kaleido-from-r-to-convert-plotly-graph-to-static-image
#reticulate::py_run_string("import sys")


diel_hyp <- function(sci_name, hyp_set) {
  print(paste0("Processing for species: ", sci_name))
  index <- which(dat$scientificName==sci_name)
  sp.data <- dat[index,]
  # Number of analysis unit for this species  
  n_units <- nrow(sp.data)

  # This code will loop over the analysis unit data
  # and fit hypotheses and extract model support probabilities
  # and get posterior samples for the most supported hypothesis.
  
  # Extract twilight/day/night observations and force into matrix
  y <- t(matrix(c(sp.data$twilight,
               sp.data$day,
               sp.data$night),
             nrow=3,byrow = TRUE))
  # Create a matrix that records the model probabilities for 
  # the most supported hypothesis (from a given set)
  hyp.best <- as.data.frame(matrix(NA, nrow=nrow(y),ncol=2))
  colnames(hyp.best)=c("hyp","prob")
  
  # Create a list to save posterior samples of the most supported hypothesis
  # for each analysis unit and the convergence criteria (gelman-rubin diagnostic)
  hyp.model.fit <- vector("list", n_units)
  
  # Define which hypothesis set to use - Traditional or General
  # hyp.set <- hyp.sets("Traditional")
  # hyp.set <- hyp.sets("General")
  hyp.set <- hyp.sets(hyp_set)
  
  # Start for loop  
  for(i in 1:nrow(y)){
    # Fit hyps using diel.fit function
    # only a "small number of MCMC iterations are needed"
    out <- diel.fit(y[i,,drop=FALSE],hyp=hyp.set,prints=FALSE,
                 n.mcmc = 30000, burnin = 5000)
    
    # Get most supported hyp and prob
    index <- which.max(out$bf.table[,2])
    extract.hyp.and.prob <- c(rownames(out$bf.table)[index],
                           out$bf.table[index,2])
    # Print every iteration
    # print(c(i, extract.hyp.and.prob))
    # store object  
    hyp.best[i,] <- extract.hyp.and.prob
    
    # Use most supported model to get posterior samples using 3 chains
    out <- diel.fit(y[i,,drop=FALSE],hyp=hyp.best[i,1],
                 prints=FALSE, bf.fit=FALSE, post.fit=TRUE,
                 n.chains = 3,
                 n.mcmc = 50000, burnin = 10000)
    
    # Combine the three chains  
    post.samps <- do.call("rbind", out$post.samp[[1]])
    
    # Only random 10k samples
    post.samps <- post.samps[sample(1:nrow(post.samps), 10000),] 
    
    # Save the gelman-rubin diagnostics of convergence
    # and posterior samples
    hyp.model.fit[[i]][[1]] <- out$gelm.diag
    hyp.model.fit[[i]][[2]] <- post.samps
    # hyp.model.fit[[i]][[3]] <- "name" 
    # where name is the file name
    # this will help get it for each project 
    hyp.model.fit[[i]][[3]] <- sp.data$file_name[i]

    # Monitor loop/print every 10th iteration
    if(i%%10==0) cat("\n i =", i)
    
  } #End for loop
  
  dir_name <- path_join(c('.', 'Output', 'DielNiche', 'Species',
                          gsub(" ", "_", sci_name), hyp_set))
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
  }
  
  # Save the main objects
  # save(hyp.model.fit,file="../Output/DielNiche/Test/hyp.model.fit")
  save(hyp.model.fit,file=path_join(c(dir_name, "hyp.model.fit")))
  save(hyp.best,file=path_join(c(dir_name, "hyp.best")))  
  
  # Explore hyp support for the species
  table(hyp.best[,1])
  
  # Remove hyps that were not supported 
  # by more than >0.7
  table(hyp.best[which(hyp.best[,2]>0.7),1])
  
  # Plot posteriors on main hypotheses
  # need to extract and combine all posteriors, also subsample
  # to reduce plotting time    
  combine.samples <- lapply(hyp.model.fit,FUN=function(x){
    temp <- x[[2]]
    temp[sample(1:nrow(temp),1000),]
  })
  
  combine.samples <- do.call("rbind",combine.samples)
  
  # Plot posterior samples with full diel niche space for this
  # hypothesis set
  
  # plot_diel <- triplot(hyp=hyp.set,
  #                        posteriors=combine.samples,
  #                        more.points = TRUE)
  # diel_plot_file <- path_join(c(dir_name, "diel.png"))
  # Plotly causes an array buffer issue for some species
  # with a lot of points (eg. Canis latrans) 
  # Throw up a warning for the plot (so plot is empty/not saved and 
  # avoid terminating the run) but continue running for the remaining species
  # tryCatch(
  #   save_image(plot_diel, diel_plot_file),
  #   error=function(e) {print("Error saving plot")}
  # )
  # Plot posterior distributions for crepuscular/day/night
  fname <- path_join(c(dir_name, "hist.png"))
  png(fname)
  hist(hyp.model.fit[[1]][[2]][,1])
  hist(hyp.model.fit[[1]][[2]][,2],add=TRUE,col=2)
  hist(hyp.model.fit[[1]][[2]][,3],add=TRUE,col=3)
  dev.off()
}

# Load analysis unit data
#dat <- read.csv("../data/diel_data.txt")
#dat <- read.csv("../data/analysis_units/diel_data.csv")
dat <- read.csv("./data/analysis_units/diel_data.csv")

sci_names <- unique(dat$scientificName)

# Test output for fossa
# diel_hyp("Cryptoprocta ferox", "Traditional")
# diel_hyp("Cryptoprocta ferox", "General")
# diel_hyp("Cryptoprocta ferox", "Maximizing")

# Test output for Virginia opossum
# diel_hyp("Didelphis virginiana", "Traditional")
# diel_hyp("Didelphis virginiana", "General")

hyp_sets <- c("Traditional", "General", "Maximizing")

# Diel niche output for a subset of species
# Assign a and b values such that a>=1 and b<=total species
# for (i in sci_names[a:b]) {
#  for (hyp_set in hyp_sets) {
#      diel_hyp(i, hyp_set)
#  }
#}

# Diel niche output for all species
for (i in sci_names) {
 for (hyp_set in hyp_sets) {
     diel_hyp(i, hyp_set)
 }
}

###Issues####
# For sp 108: plotly error:
# Error in py_run_string_impl(code, local, convert) : 
# ValueError: Transform failed with error code 525: Array buffer allocation failed
#
# (1) Fixed with the trycatch statement
# This one fails with a plotly issue "Array buffer allocation failed"
# diel_hyp("Canis latrans", "Traditional")
# The diel.png file alone will not be generated - rather it'll be empty
# 
# (2) Alternatively use the range functionality for the remaining species:
# for (i in sci_names[109:length(sci_names)]) {
#   for (hyp_set in hyp_sets) {
#       diel_hyp(i, hyp_set)
#     }
#   }


