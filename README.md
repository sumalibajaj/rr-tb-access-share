# Tuberculosis and rifampicin-resistant TB modelling: manuscript code

This repository contains R code for a tuberculosis (TB) transmission modelling analysis, with a focus on rifampicin-resistant TB (RR-TB), diagnosis, treatment outcomes, and counterfactual intervention scenarios. The code accompanies a research manuscript and is organised as a reproducible workflow using the [`targets`](https://docs.ropensci.org/targets/) package.

## Scope of this repository

**One illustrative country is included for each of four model structures**, rather than the full set of country-specific analyses. These examples demonstrate how the model structures are defined, calibrated, and used for counterfactual simulations. The examples should **not** be interpreted as the complete set of country runs or, by themselves, as reproducing all manuscript estimates.

| Model | Example country | Main structural feature |
| --- | --- | --- |
| M1 | India | Public, engaged-private and non-engaged-private treatment sectors |
| M2 | Philippines | Public and non-engaged-private treatment sectors |
| M3 | China | TB model without explicit private-sector or HIV stratification |
| M4 | South Africa | TB model stratified by HIV status and antiretroviral therapy (ART) |

The example countries are identified in `_targets.R`. The broader country list used for descriptive inputs and plots is separate from the four illustrative model runs.

## Repository files

| File | Purpose |
| --- | --- |
| `_targets.R` | Main pipeline: data preparation, model inputs, calibration, posterior summaries, counterfactual simulations, figures and tables |
| `data_exploration.R` | Extracts and processes country-level TB burden, notification and treatment-outcome inputs |
| `get_hiv_data.R` | Prepares HIV incidence, prevalence and ART coverage inputs |
| `get_addresses.R` | Constructs compartment indices and model state mappings |
| `allocate_parameters.R` | Maps sampled parameters to model rates and probabilities |
| `goveqs.R` | Differential equations for the four model structures and scale-up dynamics |
| `get_distribution_fns.R` | Constructs likelihood functions from uncertainty intervals |
| `get_objective.R` | Simulates model trajectories and calculates calibration objectives |
| `MCMC_adaptive.R` | Adaptive Metropolis-style MCMC calibration sampler |
| `get_intervention.R` | Historical counterfactual and prospective intervention simulations |
| `summarise_outputs.R` | Posterior, simulation and aggregated-output summaries |
| `misc_plotting.R` | Figures and diagnostic plots |
| `renv.lock` | Recorded R package environment |

## Requirements

- R (the supplied `renv.lock` records **R 4.5.1**).
- The [`renv`](https://rstudio.github.io/renv/) package for restoring the recorded dependency environment.
- The [`targets`](https://docs.ropensci.org/targets/) package for running the workflow.
- Input data files referenced by `_targets.R` and helper scripts, in their expected relative paths under `data/`.

The code uses packages including `dplyr`, `tidyr`, `purrr`, `ggplot2`, `deSolve`, `MASS`, `R.utils`, `future`, `lhs`, `readxl`, `cowplot`, and `patchwork`. Refer to `renv.lock` for the recorded package versions. Some packages are loaded directly in individual scripts rather than through `tar_option_set()`.

## Running the analysis

1. Clone or download this repository and set the working directory to its root.
2. Ensure the required input data are available in the `data/` directory, preserving the filenames and subdirectories expected by the scripts. The uploaded code references TB burden, RR-TB burden, notification and treatment-outcome files, HIV/ART input files, and population projections. These data files are **not included in the set of code files supplied for this README**.
3. In R, restore the package environment:

   ```r
   install.packages("renv") # if not already installed
   renv::restore()
   ```

4. Inspect `_targets.R` for data paths, iteration settings and the targets you want to run. Then run:

   ```r
   targets::tar_manifest()  # inspect the pipeline
   targets::tar_make()      # execute the pipeline
   ```

   For selected outputs, use, for example, to run 4 MCMC chains in parallel:

   ```r
   targets::tar_make_future(output_mcmc_in_1, workers = 4)
   targets::tar_read(posterior_summary_in_1)
   ```

### Calibration settings

The supplied `_targets.R` sets `iters_init_screen = 100`, `iters_mcmc_stage1 = 100`, `iters_mcmc_stage2 = 100`, and `intx_sams = 5`. **These are short example settings, not evidence of the iteration counts used for final manuscript results.** For substantive inference, configure appropriate chain lengths, posterior sampling and convergence diagnostics before running the analysis.

The calibration workflow uses initial parameter screening, adaptive MCMC, posterior summaries and diagnostic plots. The MCMC function accepts `hiv_model = TRUE` for the HIV-stratified model and `FALSE` for the other models, changing which TB mortality multiplier is checked against a model constraint.

## Counterfactual and intervention analysis

`get_intervention.R` evaluates historical trajectories and counterfactual scenarios, as well as forward projections. Scenarios in the code include a historical counterfactual without RR-TB recognition/services and prospective improvements in drug-susceptibility testing (DST) and second-line (SL) treatment outcomes. Outputs include RR-TB incidence and TB mortality trajectories, with posterior uncertainty summaries.

The pipeline also contains targets for aggregating example-country draws. **An aggregation of the four illustrative countries is not a global estimate.** Any global results reported in the manuscript require the full intended set of country-specific analyses and appropriate population-weighted aggregation.

## Outputs

Depending on which targets are run, the pipeline generates MCMC traces, posterior distributions, posterior predictive checks, intervention and counterfactual plots, and summary tables. Several figures are written to `outputs/`, and a summary table is written under `outputs/processed data/`. Create these directories if they do not already exist.

## Data availability and reproducibility

The scripts refer to external epidemiological and demographic data.

This repository provides **illustrative model implementations**, not end-to-end reproduction of all manuscript results.

## Citation

If you use this code, please cite the associated manuscript. It will be added here in due course.

