# Multidimensional Skill-Adjusted Real Wage Cyclicality in the United States

This repository contains the R replication code for the paper's analysis of how changes in occupational skill composition affect measured U.S. real wage cyclicality.

The empirical analysis combines monthly Current Population Survey (CPS) Outgoing Rotation Group microdata with O*NET occupational skill measures, a SOC-to-IPUMS OCC2010 crosswalk, and macroeconomic series for unemployment, the output gap, recessions, and consumer prices. The main specification holds the joint distribution of eight occupational skill factors fixed through entropy calibration and compares observed aggregate wage changes with skill-standardized wage changes.

## Repository contents

- `skill_adjusted_wage_cyclicality.R` — complete low-memory analysis and replication script.
- `README.md` — data requirements, empirical design, and replication instructions.

Raw CPS microdata are not included in the repository. Users should obtain the required CPS data directly from IPUMS CPS and comply with the applicable data-use and citation requirements.

## Main empirical specification

The primary outcome is real usual weekly earnings. The analysis covers 1982–2025 and uses workers ages 16–64 in the CPS ORG sample.

The preferred adjustment is a continuous eight-factor occupational skill model derived from O*NET Importance measures. In each month, entropy calibration reweights the observed occupational distribution to match the baseline joint skill distribution. The target moments include factor means, factor squares, and pairwise cross-products.

The main decomposition is:

- **Aggregate wage change:** year-over-year percentage change in observed real weekly earnings.
- **Skill-standardized wage change:** year-over-year percentage change after holding the occupational skill distribution fixed.
- **Composition component:** the difference between the aggregate and skill-standardized wage changes.

Monthly cyclicality regressions relate each component to the 12-month change in the unemployment rate, include calendar-month fixed effects, and use Newey-West standard errors with a 12-month lag window.

The script also estimates quarterly specifications using four-quarter changes in unemployment and the output gap.

## Occupational skill measures

The preferred O*NET specification uses eight latent factors:

1. Physical and manual capability
2. Management and interpersonal coordination
3. Education and professional knowledge
4. Engineering and technical systems
5. Perceptual and cognitive vigilance
6. Quantitative and scientific capability
7. Information analysis and administrative execution
8. Transportation and public-safety operations

The main scalar benchmark uses `SCALAR_EQUAL_DOMAIN_Z`. The script also supports five-, ten-, and thirteen-factor robustness specifications when the corresponding candidate factor-score files are available.

## Required software

The script requires R and the following packages:

```r
install.packages(c(
  "data.table",
  "fixest",
  "ggplot2",
  "lubridate",
  "zoo",
  "scales"
))
```

The analysis was designed to run on a standard desktop computer without loading all CPS microdata into memory at once. It processes one CPS extract at a time, saves compact intermediate panels to disk, and then performs the decomposition and regression analysis using those panels.

## Required data

Set `DATA_DIR` in the script, or define the environment variable `SKILL_WAGE_DATA_DIR`, so that it points to a directory containing the files below.

### CPS microdata

The final analysis used nine CPS extracts spanning 1981–2025, named according to this pattern:

```text
ipums_cps_1981_1985.csv
ipums_cps_1986_1990.csv
ipums_cps_1991_1995.csv
ipums_cps_1996_2000.csv
ipums_cps_2001_2005.csv
ipums_cps_2006_2010.csv
ipums_cps_2011_2015.csv
ipums_cps_2016_2020.csv
ipums_cps_2021_2025.csv
```

The script discovers files matching `ipums_cps_YYYY_YYYY.csv`, so equivalent partitions can be used as long as they collectively cover the required years.

The CPS files must contain these core variables:

```text
YEAR MONTH EARNWT AGE SEX RACE HISPAN EDUC STATEFIP
OCC OCC2010 IND IND1990 EMPSTAT EARNWEEK2 HOURWAGE2 PAIDHOUR
```

The following variables are optional but are used when available for secondary outcomes and robustness checks:

```text
AHRSWORKT UHRSWORKORG WKSWORKORG UNION METFIPS METRO
CLASSWKR MISH ELIGORG
```

`EARNWT` is used as the CPS ORG analysis weight.

### CPS–O*NET crosswalk

```text
SOCXOCC2010.csv
```

The crosswalk must contain `SOC_CODE` and either `OCC2010` or `OCC`. If a mapping-weight column is available, the script recognizes the first available field among:

```text
CROSSWALK_WEIGHT
EMPLOYMENT_WEIGHT
CONVERSION_FACTOR
SHARE
WEIGHT
```

If no weight is supplied, SOC mappings are equally weighted within OCC2010.

### O*NET skill measures

The primary skill file is expected at:

```text
onet_skill_measure_results_v2/data/onet_skill_measures_by_soc.csv
```

It must contain:

```text
SOC_CODE
SCALAR_EQUAL_DOMAIN_Z
FACTOR_01 ... FACTOR_08
```

For the optional factor-count robustness checks, the script will also use these files if present:

```text
onet_skill_measure_results_v2/data/candidate_factor_scores_f5.csv
onet_skill_measure_results_v2/data/candidate_factor_scores_f10.csv
onet_skill_measure_results_v2/data/candidate_factor_scores_f13.csv
```

These skill measures are treated as fixed occupational attributes in the wage analysis. The construction and validation of the O*NET skill space are conceptually prior to this replication script.

### Macroeconomic data

The script expects four local files:

```text
Unemployment Rate.csv
Output Gap.csv
Recessions.csv
CPIAUCSL.csv
```

The readers accept conventional date fields and common alternative value-column names. The CPI series is used to express earnings in 2025 dollars.

## Running the analysis

The simplest approach is to place the R script in the data directory and run:

```r
setwd("/path/to/data")
source("skill_adjusted_wage_cyclicality.R", echo = FALSE)
```

Alternatively, keep the data outside the repository and define its location before sourcing the script:

```r
Sys.setenv(SKILL_WAGE_DATA_DIR = "/path/to/private/data")
source("skill_adjusted_wage_cyclicality.R", echo = FALSE)
```

An optional output directory can also be supplied:

```r
Sys.setenv(SKILL_WAGE_OUTPUT_DIR = "/path/to/results")
```

By default, results are written to:

```text
skill_wage_results/
```

### Cache behavior

`REBUILD_CACHE <- FALSE` is the recommended setting. On the first run, missing cache files are created automatically. On later runs, the compact cached panels are reused.

Set:

```r
REBUILD_CACHE <- TRUE
```

only when the underlying CPS extracts, crosswalk, CPI construction, or occupation-level skill inputs have changed.

## Main outputs

The output directory contains `data/`, `tables/`, `figures/`, `logs/`, `models/`, and `cache/` subdirectories.

Key files include:

```text
tables/main_cyclicality_coefficients.csv
tables/main_attribution_statistics.csv
tables/main_attribution_components_long.csv
tables/main_continuous_calibration_diagnostics.csv
tables/quarterly_cyclicality_coefficients.csv
tables/robustness_cyclicality_coefficients.csv
tables/main_block_bootstrap_summary.csv
tables/episode_specific_cyclicality_coefficients.csv
tables/recession_influence_windows.csv
tables/leave_one_recession_influence_window_out_coefficients.csv
tables/factor_employment_cyclicality.csv
tables/skill_crosswalk_coverage_by_year.csv
tables/weekly_earnings_upper_support_breaks.csv
```

Main figures are written as:

```text
figures/figure_1_weekly_earnings_decomposition.png
figures/figure_2_method_comparison.png
figures/figure_3_episode_specific_composition.png
```

A successful run ends by writing:

```text
logs/run_status.txt
```

with the value `SUCCESS`, and records the R environment in `logs/session_info.txt`.

## Robustness and secondary analyses

The script includes the following analyses used to evaluate the main result:

- real earnings per actual hour;
- reported hourly wages and usual-hour earnings;
- prime-age and full-time samples;
- common real upper-cap specifications;
- five-, ten-, and thirteen-factor skill spaces;
- alternative 1991–1995 calibration baseline;
- scalar skill groups using 5, 10, and 20 groups;
- conventional age × education × sex × race/ethnicity worker-composition cells;
- quarterly unemployment and output-gap specifications;
- occupational-factor employment responses;
- a 999-replication moving-block bootstrap using 24-month blocks.

The optional demographic × industry-cell robustness is disabled by default because it generates a large number of cells. Set `RUN_INDUSTRY_CELL_ROBUSTNESS <- TRUE` to run it.

## Historical-episode definitions

Because the dependent variables are 12-month changes, the paper distinguishes current-month dates from the dates entering the year-over-year comparison.

The reported episode samples are:

- **Clean pre-COVID:** current month before March 2020.
- **Pandemic/recovery transition:** March 2020 through December 2022.
- **Post-transition current date:** January 2023 onward; some 2023 observations still compare with 2022.
- **Strict clean post-transition:** both the current month and the 12-month reference month occur after December 2022; this begins in January 2024.
- **Strict non-transition sample:** clean pre-COVID observations combined with strict clean post-transition observations.

For recession-sensitivity tests, the omitted influence window begins 12 months before each NBER recession and ends 12 months after the recession. The COVID recession influence window is extended through December 2022 to remove the pandemic/reopening transition from that specification.

## Reference results from the final run

The final 1982–2025 run produced the following main weekly-earnings coefficients on the 12-month unemployment-rate change:

| Component | Coefficient | Standard error |
|---|---:|---:|
| Aggregate real weekly earnings | 0.623 | 0.168 |
| Skill-standardized weekly earnings | 0.398 | 0.104 |
| Occupational skill composition | 0.225 | 0.083 |

The continuous eight-factor composition component equals approximately **36.1%** of the full-sample aggregate coefficient. These values can be used as a replication check when the same source data and derived O*NET skill measures are used.

The moving-block bootstrap is intentionally more conservative than the Newey-West inference. In the final run, the bootstrap 95% interval for the composition coefficient included zero.

## Notes on interpretation

The skill-standardized series is a composition-adjusted aggregate wage index. It should not be interpreted as a causal estimate of the price of skill. The procedure holds the occupational skill distribution fixed while occupation-specific earnings and CPS employment weights vary over time.

Similarly, the conventional demographic-cell adjustment and the occupational-skill adjustment are alternative composition controls. Their estimated contributions should not be mechanically added because worker demographics, education, occupation, and skill composition overlap.

## Reproducibility and data sharing

The R code and nonrestricted derived inputs may be posted publicly. Raw IPUMS CPS extracts should not be uploaded to this repository. Researchers reproducing the analysis should obtain the CPS microdata directly from IPUMS CPS and cite IPUMS and O*NET according to their current documentation and data-use requirements.

For exact reproduction, retain the original crosswalk, the same O*NET skill-measure files, and the macroeconomic source files used in the paper.
