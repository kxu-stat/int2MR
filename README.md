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
install_github("Likeli-Ke/int2MR")

# Load the int2MR package
library(int2MR)
```

# Usage

Please refer to the ['int2MR' vignette](https://github.com/Likeli-Ke/int2MR/blob/main/vignettes/int2MR.Rmd) for a tutorial to use the 'int2MR' package.

# Reference
This package is maintained by [kxu6@nd.edu](kxu6@nd.edu).

