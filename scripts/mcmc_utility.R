
split_mcmc <- function(x){
  # get parameter names
  pars <- colnames(x)
  # unique parameters
  unq_pars <- unique(
    gsub(
      "\\[.*\\]",
      "",
      pars
    )
  )
  # make list object to store arrays in
  result_list <- vector(
    "list",
    length = length(unq_pars)
  )
  # name the list
  names(result_list) <- unq_pars
  # fill in the arrays
  for(i in 1:length(result_list)){
    # get just the parameters
    tmp <- pars[grep(
      paste0(
        "^",unq_pars[i], "\\["
      ),
      pars
    )]
    if(length(tmp) == 0){
      tmp <- pars[grep(
        paste0("^",unq_pars[i],"$"),
        pars
      )]
    }
    # and then the array dimensions
    arr_dim <- gsub(
      paste0(
        unq_pars[i],"\\[|\\]"
      ),
      "",
      tmp
    )
    arr_dim <- strsplit(
      arr_dim,
      ","
    )
    ndim <- length(arr_dim[[1]])
    npar <- length(arr_dim)
    # make a matrix
    arr_ind <- suppressWarnings(
      matrix(
        as.numeric(
          unlist(arr_dim)
        ),
        ncol = ndim,
        nrow = npar,
        byrow = TRUE
      )
    )
    if(nrow(arr_ind) == 1 & ncol(arr_ind) == 1){
      arr_ind[1,1] <- 1
    }
    # get max index for each
    max_ind <- apply(arr_ind, 2, max)
    # and then fill in the array
    result_list[[i]] <- array(
      x[,tmp],
      dim = c(nrow(x), max_ind)
    )
    
  }
  return(result_list)
}


qw <- function(x){
  quantile(x, probs = c(0.025,0.5,0.975))
}


gen_preds_unit_inxs <- function(
    mcmc, beta_mu_idx, beta_trait_idx, covar, inxs_covar
){
  noct_pred <- mcmc$noct_beta_mu[,beta_mu_idx] %*% rbind(1,covar)
  noct_pred_lo <- noct_pred +  
    mc$noct_trait_beta[,beta_trait_idx] %*% rbind(
      inxs_covar[1], inxs_covar[1] * covar)
  noct_pred_hi <- noct_pred +  
    mcmc$noct_trait_beta[,beta_trait_idx] %*% rbind(
      inxs_covar[2], inxs_covar[2] * covar)
  noct_pred_lo <- exp(noct_pred_lo)
  noct_pred_hi <- exp(noct_pred_hi)
  
  diur_pred <- mcmc$diur_beta_mu[,beta_mu_idx] %*% rbind(1,covar)
  diur_pred_lo <- diur_pred +  
    mc$diur_trait_beta[,beta_trait_idx] %*% rbind(
      inxs_covar[1], inxs_covar[1] * covar)
  diur_pred_hi <- diur_pred +  
    mcmc$diur_trait_beta[,beta_trait_idx] %*% rbind(
      inxs_covar[2], inxs_covar[2] * covar)
  diur_pred_lo <- exp(diur_pred_lo)
  diur_pred_hi <- exp(diur_pred_hi)
  
  denom_lo <- 1 + noct_pred_lo + diur_pred_lo
  denom_hi <- 1 + noct_pred_hi + diur_pred_hi
  
  # generate probability
  cath_lo <- 1 / denom_lo
  cath_hi <- 1 / denom_hi
  noct_lo <- noct_pred_lo / denom_lo
  noct_hi <- noct_pred_hi / denom_hi
  diur_lo <- diur_pred_lo / denom_lo
  diur_hi <- diur_pred_hi / denom_hi
  
  cath_lo <- apply(
    cath_lo,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  cath_hi <- apply(
    cath_hi,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  noct_lo <- apply(
    noct_lo,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  noct_hi <- apply(
    noct_hi,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  diur_lo <- apply(
    diur_lo,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  diur_hi <- apply(
    diur_hi,
    2,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  to_return <- list(
    cathemeral = list(
      lo = t(cath_lo),
      hi = t(cath_hi)
    ),
    nocturnal = list(
      lo = t(noct_lo),
      hi = t(noct_hi)
    ),
    diurnal = list(
      lo = t(diur_lo),
      hi = t(diur_hi)
    )
  )
  return(to_return)
}
