# global-diel-analysis

### data (`./data/`)

**./analysis_units/diel_data.csv**: The summarised camera trap data to be used in the constrained multinomial model. The supporting information has all the information about how we created our analysis units, which are summarised to the level of a camera trapping project. Each project can have multiple analysis units if the project sampled for long enough (and species had sufficient records). This dataset has 14587 rows and 18 rows of data.


| Column           | Data type              | Explanation                                                                                                                                                                                                                                                                                                                |
|------------------|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `scientificName` | Categorical            | The binomial nomenclature for the species related to the given analysis unit                                                                                                                                                                                                                                               |
| `analysis_unit`  | Key                    | A unique identifier for each row. Used to join to other datasets as needed                                                                                                                                                                                                                                                 |
| `twilight`       | Integer                | The number of independent records taken during twilight                                                                                                                                                                                                                                                                    |
| `day`            | Integer                | The number of independent records taken during the day                                                                                                                                                                                                                                                                     |
| `night`          | Integer                | The number of independent records taken at night                                                                                                                                                                                                                                                                           |
| `trap_nights`    | Integer                | The number of trap nights (cameras by days of sampling)                                                                                                                                                                                                                                                                    |
| `nsite`          | Integer                | The number of camera trapping sites                                                                                                                                                                                                                                                                                        |
| `min_date`       | Date (m/d/yyyy)        | The start date for this analysis unit                                                                                                                                                                                                                                                                                      |
| `max_date`       | Date (m/d/yyyy)        | The end date for this analysis unit                                                                                                                                                                                                                                                                                        |
| `mean_lat`       | Spatial coordinate (y) | The average latitude (crs = 4326) of this project                                                                                                                                                                                                                                                                          |
| `mean_lon`       | Spatial coorindate (x) | The average longitude (crs = 4326) of this project                                                                                                                                                                                                                                                                         |
| `country`        | Categorical            | The name of the country where camera trapping occurred                                                                                                                                                                                                                                                                     |
| `phylum`         | Categorical            | The Phylum of the species related to the given analysis unit                                                                                                                                                                                                                                                               |
| `class`          | Categorical            | The Class of the species related to the given analysis unit                                                                                                                                                                                                                                                                |
| `order`          | Categorical            | The Order of the species related to the given analysis unit                                                                                                                                                                                                                                                                |
| `family`         | Categorical            | The Family of the species related to the given analysis unit                                                                                                                                                                                                                                                               |
| `file_name`      | Categorical            | The name of the file this analysis unit was generated from                                                                                                                                                                                                                                                                 |
| `unit_type`      | Categorical            | The type of analysis unit. Can be one of three values: `28day`, which are analysis units that are 28 days long. `56day`, which are analysis units that are 56 days long. `allday`, which use all the data from a camera trapping project. These different windows were used to increase the sample size for rarer species. |

*Traditional.species.hyps.with.analysis.units.csv*




### Scripts (`./scripts/`)

*diel_niche_all_species.R*: Fits the Traditional hypothesis set via the `Diel.Niche` package to all analysis units in the dataset.

**fit_hypothesis_multinomial_nimble.R:** Fits the Categorical regression to the most supported diel phenotype for each species and analysis unit that has substantial support (Pr(hypothesis) > 0.8). There were a few other species that were removed because they were either arboreal or too small.

**ghf.R:** Pulls the global human footprint data for each analysis unit and year. This script assumes you have already download the global human footprint data. As of 2023-12-1, this can be downloaded at https://wcshumanfootprint.org/data-access. 

**mcmc_utility.R:** Some helper functions to fit or parse the results of the `NIMBLE` model we fitted to the data.

**plot_constrained.R:** Plots out the diel-phenotype level and family-level results from the output of `Diel.Niche` (i.e., the constrained multinomial model).

**query_ranges.R:** Calculates the overall distributional extent for each species in our analysis. As of 2023-12-1, this can be downloaded at https://www.iucnredlist.org/resources/spatial-data-download. 

**summarise_multinomial.R:** Generate many of the supplemental figures and the like from the Categorical model we fitted via `NIMBLE`.






