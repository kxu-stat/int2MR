---
Integrative Mendelian randomization for detecting exposure-by-group interactions using group-specific and combined summary statistics
---

This repository contains the int2MR R package, which implements an integrative Mendelian randomization (int2MR) method for detecting both the direct exposure–outcome effect within comparison and reference groups, as well as the exposure–group interaction effect.

# Installation
Before installing int2MR, ensure that you have the devtools package installed, along with rstan for Bayesian modeling. To install the development version of int2MR from GitHub, run the following commands in R:

```
# Load devtools package
library(devtools)

# Install the int2MR package from GitHub
install_github("kxu-stat/int2MR")

# Load the int2MR package
library(int2MR)
```

# Usage

Please refer to the ['int2MR' vignette](https://github.com/kxu-stat/int2MR/blob/main/vignettes/int2MR.pdf) for a tutorial to use the 'int2MR' package.

By default, `int2MR()` uses MCMC sampling through `rstan::sampling()`. For a purely optimization-based run, set `estimation = "optimizing"`:

```r
result_2sample_opt <- int2MR(
  data_list_2sample = example_2sample_data,
  model_type = "2sample",
  estimation = "optimizing",
  prior_inv_gamma_shape = 0.1,
  prior_inv_gamma_scale = 0.1
)

result_2sample_opt$result_2sample
```

The `estimation` argument controls how `est_beta` and `est_beta_int` are estimated: `estimation = "sampling"` uses MCMC posterior means, while `estimation = "optimizing"` uses the MAP estimate from `rstan::optimizing()`. In both modes, standard errors such as `se_beta` and `se_beta_int` are computed from the inverse negative Hessian.

The inverse-gamma hyperparameters control the prior scale of uncorrelated pleiotropic-effect variance components. As a rule of thumb, start with `prior_inv_gamma_shape = prior_inv_gamma_scale = 0.02` for weak or sparse uncorrelated pleiotropic effects, try `0.1` for moderate effects, and check stronger effects with a small sensitivity grid such as `c(0.02, 0.05, 0.1, 0.2, 0.5)`. Report whether `est_beta` and `est_beta_int` are stable across the grid rather than relying on a single setting.

# Other
The code for simulation and visualization is store in the folder [supp](https://github.com/kxu-stat/int2MR/blob/main/supp). The GWAS summary statistics from ROSMAP and primary data analysis results are availale at [Zenodo](https://doi.org/10.5281/zenodo.16341091).

# Reference
Ke Xu, Nathaniel Maydanchik, Bowei Kang, Jianhai Chen, Qixiang Chen, Gongyao Xu, Shinya Tasaki, David A. Bennett, Lin S. Chen. Integrative Mendelian Randomization for Detecting Exposure-by-group Interactions Using Group-Specific and Combined Summary Statistics. [doi.org/10.1371/journal.pgen.1011819](https://doi.org/10.1371/journal.pgen.1011819)

# Development
This package is maintained by [kxu6@nd.edu](kxu6@nd.edu).
