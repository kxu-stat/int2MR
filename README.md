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

# Other
The code for simulation and visualization is store in the folder [supp](https://github.com/kxu-stat/int2MR/blob/main/supp). The GWAS summary statistics from ROSMAP and primary data analysis results are availale at [Zenodo](https://doi.org/10.5281/zenodo.16341091).

# Reference
Ke Xu, Nathaniel Maydanchik, Bowei Kang, Jianhai Chen, Qixiang Chen, Gongyao Xu, Shinya Tasaki, David A. Bennett, Lin S. Chen. Integrative Mendelian Randomization for Detecting Exposure-by-group Interactions Using Group-Specific and Combined Summary Statistics. [doi.org/10.1101/2025.01.26.25321136](https://doi.org/10.1101/2025.01.26.25321136)

# Development
This package is maintained by [kxu6@nd.edu](kxu6@nd.edu).

