# global-diel-analysis

### data (`./data`)

**./analysis_units/diel_data.csv**: The summarised camera trap data to be used in the constrained multinomial model. The supporting information has all the information about how we created our analysis units, which are summarised to the level of a camera trapping project. Each project can have multiple analysis units if the project sampled for long enough (and species had sufficient records). This dataset has 14587 rows and 18 columns of data.


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


**./analysis_units/global_human_footprint_2023_09_25.csv**: The global human footprint data quieried via `./scripts/ghf.R`. This dataset has 4 columns and 2580 rows. 



| Column           | Data type              | Explanation          |
|------------------|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`year` | Year |The year sampling occured for a given analysis unit |
|`ghf` | Proportion | The global human footprint summarised from the average location (the mean latitude and longitude) of a camera trap project. The average global human footprint was queried with a 20 km radius buffer around this point.  |
| `analysis_unit`  | Key     | A unique identifier for each row. Used to join to other datasets as needed    |          
| `file_name`      | Categorical     | The name of the file this analysis unit was generated from    |                                                                     

**./analysis_units/species_distrib_range.csv:** The area of the global extent of a species distribution, queried from IUCN range list data. This dataset has 2 columns and 448 rows.

| Column | Date type | Explanation |
|--------|-----------|-------------|  
|`scientificName`|  Categorical            | The binomial nomenclature for a species, used to join to other datasets that contain a species scientific name |
|`EOO` | Numeric | The area of a species global extent in square kilometers. EOO stands for extent of occurrence. |



**./diel.niche_results/informed/Generalspecies.hyp.ref.PRIOR90.csv:** The probability of the literature informed diel phenotype for species while using an informed prior of 0.9 for the reference. This is for the General hypothesis set, which further splits the Cathemeral hypothesis into bimodal hypothesis sets (e.g., Diurnal-Noctural for species that spend more time in two diel time periods). This dataset has 3 columns and 14587 rows. Analysis unit names are not included here, but this dataset is ordered the same as the other datasets of this size.

| Column | Date type | Explanation |
|--------|-----------|-------------|
| `species` | Categorical | The binomial nomenclature for a species, used to join to other datasets that contain a species scientific name. In other datasets this is `scientificName`.|
|`hyp.ref` | Categorical | The most supported diel phenotype based on the data. For the general hypothesis set this can include `Diurnal`, `Nocturnal`, `Crepuscular`, `Cathemeral`, `Diurnal-Crepuscular`, `Diurnal-Nocturnal`, and `Crepuscular-Nocturnal`.|
|`prob` | Probability | The probability of the reference phenotype.|

**./diel.niche_results/informed/Traditionalspecies.hyp.ref.PRIOR90.csv:** The probability of the literature informed diel phenotype for species while using an informed prior of 0.9 for the reference. This is for the Traditional hypothesis set. This dataset has 3 columns and 14587 rows. Analysis unit names are not included here, but this dataset is ordered the same as the other datasets of this size. 

| Column | Date type | Explanation |
|--------|-----------|-------------|
| `species` | Categorical | The binomial nomenclature for a species, used to join to other datasets that contain a species scientific name. In other datasets this is `scientificName`.|
|`hyp.ref` | Categorical | The most supported diel phenotype based on the data. For the traditional hypothesis set this can include `Diurnal`, `Nocturnal`, `Crepuscular`,and `Cathemeral`.|
|`prob` | Probability | The probability of the reference phenotype.|

**./diel.niche_results/uniform/Generalspecies.hyp.ref.csv:** The probability of the literature informed diel phenotype for species while using an uniform prior across all hypotheses (i.e., equal probability for each). This is for the General hypothesis set, which further splits the Cathemeral hypothesis into bimodal hypothesis sets (e.g., Diurnal-Noctural for species that spend more time in two diel time periods). This dataset has 3 columns and 14587 rows. Analysis unit names are not included here, but this dataset is ordered the same as the other datasets of this size.

| Column | Date type | Explanation |
|--------|-----------|-------------|
| `species` | Categorical | The binomial nomenclature for a species, used to join to other datasets that contain a species scientific name. In other datasets this is `scientificName`.|
|`hyp.ref` | Categorical | The most supported diel phenotype based on the data. For the general hypothesis set this can include `Diurnal`, `Nocturnal`, `Crepuscular`, `Cathemeral`, `Diurnal-Crepuscular`, `Diurnal-Nocturnal`, and `Crepuscular-Nocturnal`.|
|`prob` | Probability | The probability of the reference phenotype.|

**./diel.niche_results/uniform/Traditionalspecies.hyp.ref.csv:** The probability of the literature informed diel phenotype for species while using an uniform prior across all hypotheses (i.e., equal probability for each). This is for the Traditional hypothesis set. This dataset has 3 columns and 14587 rows. Analysis unit names are not included here, but this dataset is ordered the same as the other datasets of this size. 

| Column | Date type | Explanation |
|--------|-----------|-------------|
| `species` | Categorical | The binomial nomenclature for a species, used to join to other datasets that contain a species scientific name. In other datasets this is `scientificName`.|
|`hyp.ref` | Categorical | The most supported diel phenotype based on the data. For the traditional hypothesis set this can include `Diurnal`, `Nocturnal`, `Crepuscular`,and `Cathemeral`.|
|`prob` | Probability | The probability of the reference phenotype.|


**./diel.niche_results/total.sample.size.by.analysis.units.csv:** The number of independent records tied to a given analysis unit. This dataset has 2 columns and 14587 rows (there are actually three columns, one of them being row names numbered sequentially).As a reminder, sample size was calculated as the number of unique images at least 15 minutes apart from one another at unique camera trapping locations within a given analysis unit (which is a temporal subset of a projects data).


| Column | Date type | Explanation |
|--------|-----------|-------------|
| `analysis_unit`  | Key     | A unique identifier for each row. Used to join to other datasets as needed    |
| `ss`| Integer | Sample size |  


**./data/phylo/Phylacine_Trait_data_sp_corrected.csv:** This is the Phylacine trait database for the species used in our analysis. For information about what each column means, refer to their documentation. [Link to the paper here](https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.2443). The only difference is that we updated a few species names to line up with our own analysis.

**./hours_per_day_lat.csv:** The number of daylights hours per day given date and latitude. This dataset has 3 columns and 660665 rows.

| Column | Date type | Explanation |
|--------|-----------|-------------|
|`Julien Day` | Julien date | The ordinal day of the year, ranges from 1 to 365 |
|`latitude` | Spatial coordinate | Latitudinal band, rounded down to an integer|
|`hours.per.day`| Integer | The number of hours of daylight|

**./Traditional.species.hyps.with.analysis.units.csv:**  The probability of the most supported Traditional diel phenotype for each analysis unit. This dataset has 8 columns  and 14587 rows. 

| Column | Date type | Explanation |
|--------|-----------|-------------|
|`species`| Categorical|The binomial nomenclature for a species, used to join to other datasets that contain a species scientific name. In other datasets this is `scientificName`.|
|`family`| Categorical | The family of the associated species. |
|`hypothesis` | Categorical | The most supported diel hypothesis for that analysis unit. Can either be `Cathemeral`, `Crepuscular`, `Nocturnal`, or `Diurnal`. |
|`p_hypothesis` | Probability | The probability of the diel hypothesis. |
|`unit_type` | Categorial | The type of analysis unit. Can be one of three values: `28day`, which are analysis units that are 28 days long. `56day`, which are analysis units that are 56 days long. `allday`, which use all the data from a camera trapping project. These different windows were used to increase the sample size for rarer species. |
| `mass` | grams | The average weight of the species, in grams. Collected from the Phylacine trait database. |
|`foraging_strata` | Categorical | Where a species forages, collected from the Phylacine trait database. 1 = Terrestrial, 2 = Scansorial, 3 = Arboreal. |
|`analysis_unit` | Key     | A unique identifier for each row. Used to join to other datasets as needed  |  








### Scripts (`./scripts/`)

*diel_niche_all_species.R*: Fits the Traditional hypothesis set via the `Diel.Niche` package to all analysis units in the dataset.

**fit_hypothesis_multinomial_nimble.R:** Fits the Categorical regression to the most supported diel phenotype for each species and analysis unit that has substantial support (Pr(hypothesis) > 0.8). There were a few other species that were removed because they were either arboreal or too small.

**ghf.R:** Pulls the global human footprint data for each analysis unit and year. This script assumes you have already download the global human footprint data. As of 2023-12-1, this can be downloaded at https://wcshumanfootprint.org/data-access. 

**mcmc_utility.R:** Some helper functions to fit or parse the results of the `NIMBLE` model we fitted to the data.

**plot_constrained.R:** Plots out the diel-phenotype level and family-level results from the output of `Diel.Niche` (i.e., the constrained multinomial model).

**query_ranges.R:** Calculates the overall distributional extent for each species in our analysis. As of 2023-12-1, this can be downloaded at https://www.iucnredlist.org/resources/spatial-data-download. 

**summarise_multinomial.R:** Generate many of the supplemental figures and the like from the Categorical model we fitted via `NIMBLE`.






