# LC-SC Model: Replication Code

R code for Fosse and Winship, *Varieties of Cross-Cohort Differentiation: Generalizing the LC-SC Model for Cohort Analysis* (Sociological Science).

| Script | Example | Data file |
|---|---|---|
| `01_confidence_banks_gss.R` | Confidence in U.S. banks, GSS 1975-2022 | `gss2022.RData` |
| `02_mortality_france.R` | French male mortality, 1816-2020 | `mortality.RData` |
| `03_fertility_japan.R` | Japanese fertility, 1947-2020 | `fertility.RData` |

## Data

The data are not included. Each file holds a data frame `df` with these columns:

- `gss2022.RData`: `confinan`, `age`, `year`, from the [General Social Survey](https://gss.norc.org/)
- `mortality.RData`: `age`, `period`, `Male`, `MaleExposure`, from the [Human Mortality Database](https://www.mortality.org/) (France)
- `fertility.RData`: `age`, `period`, `Total`, `Exposure`, from the [Human Fertility Database](https://www.humanfertility.org/) (Japan)

Put the files in `Data/`.

## Running

```r
install.packages(c("mgcv", "plot3D", "RColorBrewer", "ggplot2", "ggridges", "dplyr", "tidyr"))
```

From this folder, run any script on its own:

```sh
Rscript --vanilla 03_fertility_japan.R
```

Each script writes the paper's figures to `Figures/` and its tables to `Output/`. Tested with R 4.5.0 and mgcv 1.9-3. A first run takes about 2 minutes for the banks example, about 15 for mortality and about 35 for fertility.
