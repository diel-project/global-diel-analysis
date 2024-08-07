my_model <- nimble::nimbleCode(
  {
    for(i in 1:ndata){
      # cathemeral, reference
      mu[i,1] <- 1
      # diurnal
      mu[i,2] <- exp( 
        inprod(
          diur_beta_mu[1:ncov_unit],
          unit_dm[i, 1:ncov_unit]
        ) +
        inprod(
        diur_unit_beta[species_vec[i],1:ncov_unit],
        unit_dm[i, 1:ncov_unit]
        ) + inprod(
          diur_trait_beta[1:ncov_trait],
          trait_dm[i, 1:ncov_trait]
        )
        + diur_family_beta[family_vec[i]] +
          diur_project_beta[project_vec[i]]
      )
      # nocturnal
      mu[i,3] <- exp(
        inprod(
          noct_beta_mu[1:ncov_unit],
          unit_dm[i, 1:ncov_unit]
        ) +
        inprod(
          noct_unit_beta[species_vec[i],1:ncov_unit],
          unit_dm[i, 1:ncov_unit]
        ) + inprod(
          noct_trait_beta[1:ncov_trait],
          trait_dm[i, 1:ncov_trait]
        ) + noct_family_beta[family_vec[i]] +
            noct_project_beta[project_vec[i]]
      )
      mu_prob[i,1:3] <- mu[i,1:3] / sum(mu[i,1:3])
      y[i] ~ dcat(
        mu_prob[i,1:3]
      )
    }
    # phylo_diurnal[1:nspecies,1:nspecies] ~ dmnorm(
    #   zero_vec[1:nspecies],
    #   cov = Sigma[1:nspecies, 1:nspecies]
    # )
    # phylo_nocturnal[1:nspecies,1:nspecies] ~ dmnorm(
    #   zero_vec[1:nspecies],
    #   cov = Sigma[1:nspecies, 1:nspecies]
    # )
    # family-level random effects
    diur_family_sd ~ dgamma(1,1)
    noct_family_sd ~ dgamma(1,1)
    for(fam in 1:nfamily){
      diur_family_beta[fam] ~ dnorm(
        0,
        sd = diur_family_sd
      )
      noct_family_beta[fam] ~ dnorm(
        0,
        sd = noct_family_sd
      )
    }
    diur_project_sd ~ dgamma(1,1)
    noct_project_sd ~ dgamma(1,1)
    for(proj in 1:nproject){
      diur_project_beta[proj] ~ dnorm(
        0,
        sd = diur_project_sd
      )
      noct_project_beta[proj] ~ dnorm(
        0,
        sd = noct_project_sd
      )
    }
    for(k in 1:ncov_unit){
      diur_beta_mu[k] ~ dnorm(0, sd = 1)
      noct_beta_mu[k] ~ dnorm(0, sd = 1)
      diur_sd[k] ~ dgamma(1,1)
      noct_sd[k] ~ dgamma(1,1)
      for(sp in 1:nspecies){
        diur_unit_beta[sp, k] ~ dnorm(
          0,
          sd = diur_sd[k]
        )
        noct_unit_beta[sp, k] ~ dnorm(
          0,
          sd = noct_sd[k]
        )
      }
    }
    for(j in 1:ncov_trait){
      diur_trait_beta[j] ~ dnorm(0, sd = 1)
      noct_trait_beta[j] ~ dnorm(0, sd = 1)
    }
}
)



model_single_species <- nimble::nimbleCode(
  {
    for(i in 1:ndata){
      # cathemeral, reference
      mu[i,1] <- 1
      # diurnal
      mu[i,2] <- exp( 
        inprod(
          diur_unit_beta[1:ncov_unit],
          unit_dm[i, 1:ncov_unit]
        )
      )
      # nocturnal
      mu[i,3] <- exp(
        inprod(
          noct_unit_beta[1:ncov_unit],
          unit_dm[i, 1:ncov_unit]
        )
      )
      mu_prob[i,1:3] <- mu[i,1:3] / sum(mu[i,1:3])
      y[i] ~ dcat(
        mu_prob[i,1:3]
      )
    }
    for(k in 1:ncov_unit){
      diur_unit_beta[k] ~ dnorm(
        0,
        sd = 2
      )
      noct_unit_beta[k] ~ dnorm(
        0,
        sd = 2
      )
    }
  }
)

