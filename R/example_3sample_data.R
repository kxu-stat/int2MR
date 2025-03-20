#' Three-Sample Example Data for the int2MR Package
#'
#' This dataset contains simulated summary statistics for a three-sample Mendelian Randomization analysis. It provides data from multiple GWAS sources for both the IV-to-exposure and IV-to-outcome associations, along with corresponding variances and comparison group proportions for interaction effect detection.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{p}{An integer representing the number of genetic instruments.}
#'   \item{hat_gamma}{A numeric vector of the estimated IV-to-exposure effects from the first IV-to-exposure GWAS summary data.}
#'   \item{hat_s_gamma_sq}{A numeric vector of variances (squared standard errors) for the IV-to-exposure effects from the second IV-to-exposure GWAS summary data.}
#'   \item{hat_s1_sq}{A numeric vector of variances (squared standard errors) for the IV-to-outcome effects from the first IV-to-outcome GWAS summary data.}
#'   \item{hat_s2_sq}{A numeric vector of variances for the IV-to-outcome effects from the second IV-to-outcome GWAS summary data.}
#'   \item{hat_s3_sq}{A numeric vector of variances for the IV-to-outcome effects from the third IV-to-outcome GWAS summary data.}
#'   \item{hat_Gamma1}{A numeric vector of the estimated IV-to-outcome effects from the first IV-to-outcome GWAS summary data.}
#'   \item{hat_Gamma2}{A numeric vector of the estimated IV-to-outcome effects from the second IV-to-outcome GWAS summary data.}
#'   \item{hat_Gamma3}{A numeric vector of the estimated IV-to-outcome effects from the third IV-to-outcome GWAS summary data.}
#'   \item{rho1}{A numeric value specifying the proportion of the comparison group for detecting the interaction effect from the first IV-to-outcome GWAS summary data.}
#'   \item{rho2}{A numeric value specifying the proportion of the comparison group for detecting the interaction effect from the second IV-to-outcome GWAS summary data.}
#'   \item{rho3}{A numeric value specifying the proportion of the comparison group for detecting the interaction effect from the third IV-to-outcome GWAS summary data.}
#' }
#' @source Simulated data for demonstration purposes.
"example_3sample_data"
