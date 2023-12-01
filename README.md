# global-diel-analysis

### data (`./data/`)

*diel_data.csv*

*Traditional.species.hyps.with.analysis.units.csv*




### Scripts (`./scripts/`)

*diel_niche_all_species.R*: Fits the Traditional hypothesis set via the `Diel.Niche` package to all analysis units in the dataset.

**fit_hypothesis_multinomial_nimble.R:** Fits the Categorical regression to the most supported diel phenotype for each species and analysis unit that has substantial support (Pr(hypothesis) > 0.8). There were a few other species that were removed because they were either arboreal or too small.

**ghf.R:** Pulls the global human footprint data for each analysis unit and year. This script assumes you have already download the global human footprint data. As of 2023-12-1, this can be downloaded at https://wcshumanfootprint.org/data-access. 

**mcmc_utility.R:** Some helper functions to fit or parse the results of the `NIMBLE` model we fitted to the data.

**plot_constrained.R:** Plots out the diel-phenotype level and family-level results from the output of `Diel.Niche` (i.e., the constrained multinomial model).

**query_ranges.R:** Calculates the overall distributional extent for each species in our analysis. As of 2023-12-1, this can be downloaded at https://www.iucnredlist.org/resources/spatial-data-download. 

**summarise_multinomial.R:** Generate many of the supplemental figures and the like from the Categorical model we fitted via `NIMBLE`.






