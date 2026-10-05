# CH4 Genetic Correlation and Methane-Efficiency Simulation Pipeline

This repository contains an R-based simulation and validation workflow developed to test the statistical procedures used in the CH4 manuscript before applying them to real genetic evaluation data.

The main script generates synthetic animal-level GEBV and reliability data, adjusts GEBV for genetic trends, estimates approximate genetic correlations (AGC), constructs methane-efficiency indices, and performs simulated daughter and LR validation. The generated files are test outputs and must not be interpreted as results from real animals.

## Main script

`ch4_AGC_ME_indices_simulation_test.R`

The script is organized into nine modules:

1. **Basic settings**
   - Defines the 37 comparison traits, PCH4, trait groups, output directory, random seed, and number of bootstrap replicates.

2. **Simulation of GEBV and GREL data**
   - Simulates 1,200 bulls and 5,000 cows by default.
   - Generates GEBV for PCH4 and 37 production, functional, conformation, and economic-index traits.
   - Adds trait-specific linear and quadratic birth-year trends.
   - Simulates trait reliabilities (`GREL`), PCH4 daughter counts, and PCH4-record indicators.

3. **Selection of highly reliable bulls**
   - Retains bulls with at least 30 PCH4 daughters.
   - Requires `GREL >= 0.50` for PCH4 and every comparison trait.
   - The default simulation intentionally produces 1,020 bulls that meet these criteria.

4. **Genetic-trend adjustment**
   - Calculates mean GEBV by birth year in the selected bulls.
   - Fits a weighted quadratic regression for each trait:

     `mean GEBV ~ birth year + birth year^2`

   - Subtracts the predicted birth-year trend from each animal's GEBV to obtain trend-adjusted GEBV (`*_GEBV_adj`).

5. **Approximate genetic correlations (AGC)**
   - Calculates correlations between trend-adjusted GEBV.
   - Applies a reliability-based correction factor.
   - Uses non-parametric bootstrap resampling to calculate the standard error and percentile 95% confidence interval.
   - Restricts AGC estimates to the valid interval from -1 to 1.

6. **Methane-efficiency indices**
   - Fits expected PCH4 using three alternative predictor sets:
     - **MEP:** milk yield (`MY`), fat yield (`FY`), and protein yield (`PY`).
     - **MEF:** functional index (`VEF`).
     - **MEG:** global index (`VEG`).
   - Defines methane efficiency as:

     `ME = expected trend-adjusted PCH4 - observed trend-adjusted PCH4`

   - Therefore, a larger ME value indicates lower PCH4 than expected and greater methane efficiency.
   - Standardizes MEP, MEF, and MEG to RMEP, RMEF, and RMEG using reference cows born in 2020 with a PCH4 record:

     `RME = 100 + 10 * (ME - reference mean) / reference SD`

   - The reference population consequently has a mean of 100 and a standard deviation of 10.

7. **AGC analyses involving RME indices**
   - Estimates AGC between each RME index and the 37 comparison traits.
   - Estimates AGC among RMEP, RMEF, and RMEG.
   - Reports bootstrap SE and 95% confidence intervals.

8. **Simulated daughter validation**
   - Generates sire-level mean daughter PCH4 values.
   - Tests the association between each RME index and daughter mean PCH4 using linear regression.
   - Summarizes daughter PCH4 across five RME classes: `<85`, `85-95`, `95-105`, `105-115`, and `>115`.
   - A negative regression slope is expected because higher RME denotes greater methane efficiency.

9. **Simulated LR validation**
   - Generates partial- and whole-data GEBV for 2,038 validation animals.
   - Calculates LR-style validation statistics, including accuracy, dispersion, bias, and ratio accuracy.

## Requirements

- R (version 4.1 or later recommended)
- R package: `data.table`

Install the required package in R if necessary:

```r
install.packages("data.table")
```

## How to run

Place the main R script in the working directory and run:

```bash
Rscript ch4_AGC_ME_indices_simulation_test.R
```

By default, the analysis uses 1,000 bootstrap replicates, as intended for manuscript-level analysis. Results are written to:

```text
CH4_simulation_test_results/
```

### Quick test

To reduce execution time while checking that the workflow runs correctly, set a smaller bootstrap count:

```bash
B_AGC=50 Rscript ch4_AGC_ME_indices_simulation_test.R
```

On Windows Command Prompt:

```bat
set B_AGC=50
Rscript ch4_AGC_ME_indices_simulation_test.R
```

The value of `B_AGC` must be a positive integer. Use `B_AGC=1000` for the final analysis.

## Output files

| Output file | Description |
|---|---|
| `sim_gebv_all_animals.csv` | Simulated animal information, GEBV, and GREL for PCH4 and 37 traits. |
| `selected_bulls.csv` | Bulls passing the daughter-count and reliability criteria. |
| `genetic_trend_by_birth_year.csv` | Observed and fitted trait-specific GEBV trends by birth year. |
| `adjusted_trend_check.csv` | Birth-year means of trend-adjusted GEBV for checking residual trends. |
| `AGC_PCH4_vs_37traits_adjusted_with_95CI.csv` | AGC of PCH4 with each of the 37 traits, including bootstrap SE and 95% CI. |
| `ME_expected_PCH4_regression_coefficients.csv` | Regression coefficients used to predict expected PCH4 for MEP, MEF, and MEG. |
| `ME_RME_standardization_check.csv` | Means and SD of ME and RME indices in the reference cows. |
| `all_animals_with_ME_RME_indices.csv` | All simulated animals with adjusted GEBV, expected PCH4, ME, and RME indices. |
| `selected_bulls_with_ME_RME_indices.csv` | Selected bulls with ME/RME indices and temporary simulated RME reliabilities. |
| `AGC_RME_indices_vs_37traits_with_95CI.csv` | AGC of RMEP, RMEF, and RMEG with the 37 traits. |
| `AGC_among_RME_indices_with_95CI.csv` | Pairwise AGC among RMEP, RMEF, and RMEG. |
| `sim_daughter_pch4_by_sire.csv` | Simulated sire-level RME indices and daughter mean PCH4. |
| `daughter_PCH4_regression_on_RME.csv` | Regression results for daughter PCH4 on each RME index. |
| `daughter_PCH4_by_RME_class_summary.csv` | Daughter PCH4 summary statistics by RME class. |
| `sim_lr_validation_animals.csv` | Simulated partial- and whole-data GEBV used for LR validation. |
| `LR_validation_statistics.csv` | LR validation accuracy, dispersion, bias, and ratio accuracy. |
| `test_summary.txt` | Compact summary of the main checks and results from the current run. |

## Important interpretation notes

- This code is intended to verify the analysis workflow and expected input/output structure. All supplied CSV results are based on simulated data.
- The simulated genetic relationships among traits were deliberately specified in the data-generation step; they are not empirical biological estimates.
- `RMEP_GREL`, `RMEF_GREL`, and `RMEG_GREL` are temporary approximations calculated from the component-trait reliabilities for simulation testing only. In the real analysis, replace them with official index reliabilities when available.
- The daughter PCH4 values are constructed to decrease as RMEG increases. The simulated daughter validation therefore checks whether the code recovers the intended direction; it is not independent biological validation.
- The LR module also uses simulated partial and whole evaluations. Real LR validation requires GEBV from actual reduced- and complete-data genetic evaluations, the appropriate additive genetic variance, and animal inbreeding coefficients.
- Trait abbreviations and economic-index definitions should be documented according to the official breeding-value system used in the manuscript.

## Adapting the workflow to real data

To apply the pipeline to real results:

1. Replace the output of `simulate_gebv_data()` with an animal-level input table containing animal ID, animal type, birth year, PCH4 daughter count, trait GEBV, and trait GREL.
2. Retain the naming convention `<TRAIT>_GEBV` and `<TRAIT>_GREL`, or modify the corresponding column definitions in the script.
3. Confirm that the selected-bull criteria match the Methods section.
4. Confirm whether a quadratic birth-year model is appropriate for each trait and inspect `adjusted_trend_check.csv` for residual trends.
5. Define the intended real reference population for RME standardization.
6. Replace approximate RME reliabilities, simulated daughter means, and simulated LR inputs with official data.
7. Use at least 1,000 bootstrap replicates for the final reported confidence intervals.

## Reproducibility

The script uses the fixed random seed `20260606`. Re-running the same version of the script with the same R environment and bootstrap count should reproduce the simulation results.

