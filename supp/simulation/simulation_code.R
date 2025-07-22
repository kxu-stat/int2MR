# ------------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------------
library(lmerTest)          # Linear mixed models with Satterthwaite’s degrees of freedom
library(MendelianRandomization)  # Two-sample Mendelian randomization methods
library(rstan)             # Interface to Stan for Bayesian modeling
library(ggplot2)           # Grammar of graphics plotting
library(ggpubr)            # 'ggplot2' based publication ready plots
library(devtools)          # Tools to make developing R packages easier
library(rstanarm)          # Bayesian applied regression modeling via Stan
library(mvtnorm)           # Multivariate normal and t distributions
library(dplyr)             # Data manipulation
library(parallel)          # Parallel computing
library(Rcpp)              # Seamless R and C++ integration
library(RcppArmadillo)     # Rcpp integration with Armadillo C++ linear algebra library
library(pbapply)           # Progress bar for *apply functions
library(RcppDist)          # Rcpp integration for common distributions
library(mr.raps)           # Mendelian randomization via robust adjusted profile score

# Source custom data-generation function (joint model with UHP)
source("~/data_generation_func_joint_v1UHP.R")


# ------------------------------------------------------------------
# 2. Define first Stan model (3-sample joint analysis)
# ------------------------------------------------------------------
stan_model_code <- "
data {
  int<lower=1> p;                 // number of instruments
  vector[p] hat_s1_sq;            // sampling variances for sample 1
  vector[p] hat_s2_sq;            // sampling variances for sample 2
  vector[p] hat_s3_sq;            // sampling variances for sample 3
  vector[p] hat_gamma;            // estimated gamma (IV→exposure)
  vector[p] hat_Gamma1;           // estimated Gamma (IV→outcome) for sample 1
  vector[p] hat_Gamma2;           // ... sample 2
  vector[p] hat_Gamma3;           // ... sample 3
  vector[p] hat_s_gamma_sq;       // variance of hat_gamma
  real rho1; real rho2; real rho3; // proportions for reparameterization
}

parameters {
  real beta;                      // main causal effect
  real beta_int;                  // interaction effect
  vector[p] gamma;                // true gamma values
  vector[p] alpha1;               // pleiotropic effects (sample 1)
  vector[p] alpha2;               // ... sample 2
  vector[p] alpha3;               // ... sample 3
  real<lower=0> sigma_alpha1_sq;  // prior variances
  real<lower=0> sigma_alpha2_sq;
  real<lower=0> sigma_alpha3_sq;
  real<lower=0> sigma_gamma_sq;
}

transformed parameters {
  real<lower=0> sigma_alpha1 = sqrt(sigma_alpha1_sq);
  real<lower=0> sigma_alpha2 = sqrt(sigma_alpha2_sq);
  real<lower=0> sigma_alpha3 = sqrt(sigma_alpha3_sq);
  real<lower=0> sigma_gamma  = sqrt(sigma_gamma_sq);
}

model {
  // --- Priors on variance components ---
  sigma_gamma_sq  ~ inv_gamma(0.02, 0.02);
  sigma_alpha1_sq ~ inv_gamma(0.02, 0.02);
  sigma_alpha2_sq ~ inv_gamma(0.02, 0.02);
  sigma_alpha3_sq ~ inv_gamma(0.02, 0.02);

  // --- Priors on effects ---
  gamma  ~ normal(0, sigma_gamma);
  alpha1 ~ normal(0, sigma_alpha1);
  alpha2 ~ normal(0, sigma_alpha2);
  alpha3 ~ normal(0, sigma_alpha3);

  // --- Likelihoods ---
  hat_gamma  ~ normal(gamma, sqrt(hat_s_gamma_sq));
  hat_Gamma1 ~ normal((beta + rho1 * beta_int) .* gamma + alpha1, sqrt(hat_s1_sq));
  hat_Gamma2 ~ normal((beta + rho2 * beta_int) .* gamma + alpha2, sqrt(hat_s2_sq));
  hat_Gamma3 ~ normal((beta + rho3 * beta_int) .* gamma + alpha3, sqrt(hat_s3_sq));
}
"
stan_model <- stan_model(model_code = stan_model_code)


# ------------------------------------------------------------------
# 3. Define second Stan model (2-sample joint analysis)
# ------------------------------------------------------------------
stan_model_code2 <- "
data {
  int<lower=1> p;
  vector[p] hat_s1_sq;
  vector[p] hat_s2_sq;
  vector[p] hat_gamma;
  vector[p] hat_Gamma1;
  vector[p] hat_Gamma2;
  vector[p] hat_s_gamma_sq;
  real rho1; real rho2;
}

parameters {
  real beta_int;               // interaction effect
  real beta;                   // main effect
  vector[p] gamma;
  vector[p] alpha1;
  vector[p] alpha2;
  real<lower=0> sigma_alpha1_sq;
  real<lower=0> sigma_alpha2_sq;
  real<lower=0> sigma_gamma_sq;
}

transformed parameters {
  real<lower=0> sigma_alpha1 = sqrt(sigma_alpha1_sq);
  real<lower=0> sigma_alpha2 = sqrt(sigma_alpha2_sq);
  real<lower=0> sigma_gamma  = sqrt(sigma_gamma_sq);
}

model {
  sigma_gamma_sq  ~ inv_gamma(0.02, 0.02);
  sigma_alpha1_sq ~ inv_gamma(0.02, 0.02);
  sigma_alpha2_sq ~ inv_gamma(0.02, 0.02);

  gamma  ~ normal(0, sigma_gamma);
  alpha1 ~ normal(0, sigma_alpha1);
  alpha2 ~ normal(0, sigma_alpha2);

  hat_gamma  ~ normal(gamma, sqrt(hat_s_gamma_sq));
  hat_Gamma1 ~ normal((beta + rho1 * beta_int) .* gamma + alpha1, sqrt(hat_s1_sq));
  hat_Gamma2 ~ normal((beta + rho2 * beta_int) .* gamma + alpha2, sqrt(hat_s2_sq));
}
"
stan_model2 <- stan_model(model_code = stan_model_code2)


# ------------------------------------------------------------------
# 4. Data‐generation function (PART 1)
#
#    Simulates individual‐level and summary‐level GWAS
# ------------------------------------------------------------------
data_generation <- function(p, prob, 
                            true_beta, true_beta_SX, bs, 
                            n1, n2, n3, n4,
                            a_alpha, b_alpha3, b_alpha4,
                            beta_U3, beta_U4) 
{
  # Call the custom C++ function 'dgm_joint' to generate data
  GWAS_summary <- dgm_joint(
    m = p, n1x = n1, n2y = n2, n3y = n3, n4y = n4,
    a_f_1 = 0.1, b_f_1 = 0.3,
    a_alpha = a_alpha, b_alpha3 = b_alpha3, b_alpha4 = b_alpha4,
    beta = true_beta, beta_int = true_beta_SX, prob = prob,
    beta_U3 = beta_U3, beta_U4 = beta_U4
  )
  
  # Extract individual‐level outcomes, exposures, genotype indicators, etc.
  yall3 <- as.vector(GWAS_summary$X3y)
  yall4 <- as.vector(GWAS_summary$X4y)
  zall3 <- as.vector(GWAS_summary$Y3)
  zall4 <- as.vector(GWAS_summary$Y4)
  sex3  <- as.vector(GWAS_summary$I3y)
  sex4  <- as.vector(GWAS_summary$I4y)
  
  # Build data frames for individual‐level and summary‐level data
  df_individual <- data.frame(
    yall = c(yall3, yall4),
    zall = c(zall3, zall4),
    sex  = c(sex3, sex4)
  )
  df_summary <- data.frame(
    gammah = GWAS_summary$b_exp_1,
    se1    = GWAS_summary$se_exp_1,
    Gammah = GWAS_summary$b_out_2,
    se2    = GWAS_summary$se_out_2,
    Gammahs = GWAS_summary$b_out_3,
    se3     = GWAS_summary$se_out_3,
    Gammahs2 = GWAS_summary$b_out_4,
    se4      = GWAS_summary$se_out_4
  )
  df_geno <- t(cbind(GWAS_summary$g3y, GWAS_summary$g4y))
  
  # Return a list containing all relevant pieces
  return(list(
    df_individual = df_individual,
    df_summary    = df_summary,
    df_geno       = df_geno,
    varcomp_X1x   = GWAS_summary$varcomp_X1x,
    varcomp_U1x   = GWAS_summary$varcomp_U1x,
    varcomp_X2y   = GWAS_summary$varcomp_X2y,
    varcomp_U2y   = GWAS_summary$varcomp_U2y,
    varcomp_UHP2y = GWAS_summary$varcomp_UHP2y,
    varcomp_X3y   = GWAS_summary$varcomp_X3y,
    varcomp_U3y   = GWAS_summary$varcomp_U3y,
    varcomp_UHP3y = GWAS_summary$varcomp_UHP3y,
    varcomp_X4y   = GWAS_summary$varcomp_X4y,
    varcomp_U4y   = GWAS_summary$varcomp_U4y,
    varcomp_UHP4y = GWAS_summary$varcomp_UHP4y,
    varcomp_X2y_interact = GWAS_summary$varcomp_X2y_interact,
    gamma = GWAS_summary$gamma
  ))
}


# ------------------------------------------------------------------
# 5. Gibbs‐sampler power & type I error simulation (PART 2)
# ------------------------------------------------------------------
power_simulation <- function(n_sim = 20, n_iter = 1000, n_select = 100, 
                             prob1 = 0.5, prob2 = 0, prob3 = 1, 
                             true_beta = 0.1, true_beta_SX = -0.1, bs = 0,
                             n1 = 5000, n2 = 10000, n3 = 5000, n4 = 5000,
                             a_alpha, b_alpha3, b_alpha4,
                             beta_U3, beta_U4) 
{
  # Initialize result containers
  pval <- est_beta <- vector("list", length = 0)
  
  # Parallel loop over n_sim replications
  result <- pblapply(1:n_sim, function(sim_i) {
    # 5.1 Generate data
    data_mix     <- data_generation(p = n_select, prob = prob1,
                                    true_beta = true_beta, true_beta_SX = true_beta_SX,
                                    bs = bs, n1 = n1, n2 = n2, n3 = n3, n4 = n4,
                                    a_alpha = a_alpha, b_alpha3 = b_alpha3, b_alpha4 = b_alpha4,
                                    beta_U3 = beta_U3, beta_U4 = beta_U4)
    df_ind      <- data_mix$df_individual
    df_sum      <- data_mix$df_summary
    df_geno     <- data_mix$df_geno
    true_gamma  <- data_mix$gamma
    
    # 5.2 Prepare individual‐level data for 2SLS and interaction analyses
    colnames(df_geno) <- paste0("IV", seq_len(ncol(df_geno)))
    df_ind <- cbind(df_ind, df_geno)
    df_ind$inter <- df_ind$yall * df_ind$sex
    
    # — 2SLS estimations for main and interaction effects —
    #      first‐stage fits and predicted values
    iv_names <- paste(colnames(df_geno), collapse = " + ")
    first_stage_y    <- lm(as.formula(paste("yall ~", iv_names)), data = df_ind)
    first_stage_int  <- lm(as.formula(paste("inter ~", iv_names)), data = df_ind)
    df_ind$pred_yall <- predict(first_stage_y, newdata = df_ind)
    df_ind$pred_int  <- predict(first_stage_int, newdata = df_ind)
    
    #      second‐stage regressions
    tsls_fit <- summary(lm(zall ~ pred_yall + pred_int, data = df_ind))
    est_beta_tsls_inter     <- tsls_fit$coefficients["pred_yall", "Estimate"]
    est_beta_SX_tsls_inter  <- tsls_fit$coefficients["pred_int", "Estimate"]
    pval_beta_tsls_inter    <- tsls_fit$coefficients["pred_yall", "Pr(>|t|)"]
    pval_beta_SX_tsls_inter <- tsls_fit$coefficients["pred_int", "Pr(>|t|)"]
    
    # — Simple interaction analysis (no IV) —
    interact_fit <- summary(lm(zall ~ yall * sex, data = df_ind))
    est_beta_interact     <- interact_fit$coefficients["yall", "Estimate"]
    est_beta_SX_interact  <- interact_fit$coefficients["yall:sex", "Estimate"]
    pval_beta_interact    <- interact_fit$coefficients["yall", "Pr(>|t|)"]
    pval_beta_SX_interact <- interact_fit$coefficients["yall:sex", "Pr(>|t|)"]
    
    # 5.3 Prepare summary‐level inputs for various MR methods
    hat_gamma    <- df_sum$gammah
    Shat_gamma_sq<- df_sum$se1^2
    hat_G1       <- df_sum$Gammah; Shat_G1_sq <- df_sum$se2^2
    hat_G2       <- df_sum$Gammahs; Shat_G2_sq <- df_sum$se3^2
    hat_G3       <- df_sum$Gammahs2; Shat_G3_sq <- df_sum$se4^2
    p            <- length(hat_gamma)
    
    # 5.4 Run MR methods (IVW, Egger, Median, cML, RAPS)
    mr_in_single   <- mr_input(bx = hat_gamma, bxse = sqrt(Shat_gamma_sq),
                               by = hat_G2, byse = sqrt(Shat_G2_sq))
    mr_ivw_mix     <- mr_ivw(mr_input(bx = hat_gamma, bxse = sqrt(Shat_gamma_sq),
                                      by = hat_G3, byse = sqrt(Shat_G3_sq)), model = "fixed")
    est_beta_ivw   <- mr_ivw_mix$Estimate
    pval_beta_ivw  <- mr_ivw_mix$Pvalue
    egger_mix      <- mr_egger(mr_input_mix)
    est_beta_egger <- egger_mix$Estimate
    pval_beta_egger<- egger_mix$Pvalue.Est
    median_mix     <- mr_median(mr_input_mix)
    est_beta_median<- median_mix$Estimate
    pval_beta_median<- median_mix$Pvalue
    cml_mix        <- mr_cML(mr_input_mix, n = c(n3 + n4), DP = FALSE)
    est_beta_cml   <- cml_mix$Estimate
    pval_beta_cml  <- cml_mix$Pvalue
    raps_single    <- mr.raps(b_exp = hat_gamma, b_out = hat_G2,
                              se_exp = sqrt(Shat_gamma_sq), se_out = sqrt(Shat_G2_sq))
    est_beta_raps  <- raps_single$beta.hat
    pval_beta_raps <- raps_single$beta.p.value
    
    # 5.5 Fit Stan models via optimization (approximate posterior)
    data_list1 <- list(
      p = p,
      hat_s1_sq      = Shat_G1_sq,
      hat_s2_sq      = Shat_G2_sq,
      hat_s3_sq      = Shat_G3_sq,
      hat_gamma      = hat_gamma,
      hat_Gamma1     = hat_G1,
      hat_Gamma2     = hat_G2,
      hat_Gamma3     = hat_G3,
      hat_s_gamma_sq = Shat_gamma_sq,
      rho1 = prob1, rho2 = prob2, rho3 = prob3
    )
    opt1 <- optimizing(object = stan_model, data = data_list1, hessian = TRUE)
    fit_1 <- sampling(object = stan_model1, data = data_list_1, 
                      chains = 2, iter = 10000, warmup = 5000)
    samples_1 <- rstan::extract(fit_1)
    est1 <- c(mean(samples_1$beta), mean(samples_1$beta))
    names(est1) <- c("beta","beta_int")
    # Extract Hessian‐based covariance and compute standard errors
    cov1 <- MASS::ginv(-opt1$hessian)
    sd_beta1     <- sqrt(cov1["beta", "beta"])
    sd_beta_int1 <- sqrt(cov1["beta_int", "beta_int"])
    sd_total1    <- sqrt(cov1["beta", "beta"] + cov1["beta_int","beta_int"] + 
                           2 * cov1["beta","beta_int"])
    # est1 <- opt1$par[c("beta","beta_int")]
    
    # 5.6 Repeat for 2-sample Stan model
    data_list2 <- list(
      p = p,
      hat_s1_sq      = Shat_G3_sq,
      hat_s2_sq      = Shat_G2_sq,
      hat_gamma      = hat_gamma,
      hat_Gamma1     = hat_G3,
      hat_Gamma2     = hat_G2,
      hat_s_gamma_sq = Shat_gamma_sq,
      rho1 = prob3, rho2 = prob2
    )
    opt2 <- optimizing(object = stan_model2, data = data_list2, hessian = TRUE)
    fit_2 <- sampling(object = stan_model2, data = data_list_2, 
                      chains = 2, iter = 10000, warmup = 5000)
    est2 <- c(mean(samples_2$beta), mean(samples_2$beta))
    names(est2) <- c("beta","beta_int")
    est2 <- mean(samples_2$beta)
    cov2 <- MASS::ginv(-opt2$hessian)
    sd_beta2     <- sqrt(cov2["beta", "beta"])
    sd_beta_int2 <- sqrt(cov2["beta_int","beta_int"])
    sd_total2    <- sqrt(cov2["beta","beta"] + cov2["beta_int","beta_int"] +
                           2 * cov2["beta","beta_int"])
    # est2 <- opt2$par[c("beta","beta_int")]
    
    # 5.7 Compute p‐values for Stan‐based estimates
    pval_beta1     <- 2 * (1 - pnorm(abs(est1["beta"]    / sd_beta1)))
    pval_beta_int1 <- 2 * (1 - pnorm(abs(est1["beta_int"] / sd_beta_int1)))
    pval_total1    <- 2 * (1 - pnorm(abs(sum(est1) / sd_total1)))
    pval_beta2     <- 2 * (1 - pnorm(abs(est2["beta"]    / sd_beta2)))
    pval_beta_int2 <- 2 * (1 - pnorm(abs(est2["beta_int"] / sd_beta_int2)))
    pval_total2    <- 2 * (1 - pnorm(abs(sum(est2) / sd_total2)))
    
    # 5.8 Collect all estimates/p‐values into a data.frame for this replicate
    output <- data.frame(
      # true parameters
      true_beta       = true_beta,
      true_beta_SX    = true_beta_SX,
      # 3‐sample estimates
      est_beta        = est1["beta"],
      est_beta_SX     = est1["beta_int"],
      pval_beta       = pval_beta1,
      pval_beta_SX    = pval_beta_int1,
      pval_total      = pval_total1,
      # 2‐sample estimates
      est_beta_compare    = est2["beta"],
      est_beta_SX_compare = est2["beta_int"],
      pval_beta_compare   = pval_beta2,
      pval_beta_SX_compare= pval_beta_int2,
      pval_total_compare  = pval_total2,
      # 2SLS and interaction
      est_beta_tsls_inter     = est_beta_tsls_inter,
      est_beta_SX_tsls_inter  = est_beta_SX_tsls_inter,
      pval_beta_tsls_inter    = pval_beta_tsls_inter,
      pval_beta_SX_tsls_inter = pval_beta_SX_tsls_inter,
      est_beta_interact       = est_beta_interact,
      est_beta_SX_interact    = est_beta_SX_interact,
      pval_beta_interact      = pval_beta_interact,
      pval_beta_SX_interact   = pval_beta_SX_interact,
      # MR methods
      est_beta_ivw    = est_beta_ivw,    pval_beta_ivw    = pval_beta_ivw,
      est_beta_egger  = est_beta_egger,  pval_beta_egger  = pval_beta_egger,
      est_beta_median = est_beta_median, pval_beta_median = pval_beta_median,
      est_beta_cml    = est_beta_cml,    pval_beta_cml    = pval_beta_cml,
      est_beta_raps   = est_beta_raps,   pval_beta_raps   = pval_beta_raps
      # (可根据需要继续添加其他列)
    )
    return(output)
  }, cl = 8)  # Use 8 cores for parallelization
  
  # Combine all replicates into one data.frame
  result <- do.call(rbind, result)
  return(result)
}


# ------------------------------------------------------------------
# 6. Main simulation loops (PART 3)
#    Iterates over settings to compute type I error and power
# ------------------------------------------------------------------
# Define the grid of parameter settings
param_grid <- expand.grid(
  iter       = 2,                               # number of replicate iterations
  set_id     = 1:3,                             # index for different setting lists
  a_alpha    = c(1.5),                          # hyperparameter a_alpha
  n_select   = c(200),                          # number of selected IVs
  n3         = c(2000),                         # sample size of group-specific GWAS (female)
  n4         = c(1000),                         # sample size of group-specific GWAS (male)
  true_beta  = c(0.05),                         # main causal effect
  prob1      = c(0.5),                          # proportion of reference group
  n2         = c(10000, 20000, 50000),          # sample size of combined GWAS
  stringsAsFactors = FALSE
)

# Pre-allocate result lists
type1_results <- list()
power_results <- list()

# Iterate over each row of the grid
for (i in seq_len(nrow(param_grid))) {
  # Extract parameters for this run
  params <- param_grid[i, ]
  set_id <- params$set_id
  
  # Determine UHP and alpha settings based on set_id
  if (params$true_beta == 0) {
    setting_vals <- setting2[[set_id]]
  } else {
    setting_vals <- setting3[[set_id]]
    params$true_beta_SX <- 0
  }
  beta_U3  <- setting_vals[1]
  beta_U4  <- setting_vals[2]
  b_alpha3 <- setting_vals[3]
  b_alpha4 <- setting_vals[4]
  
  # Run power_simulation depending on true_beta
  if (params$true_beta == 0 && params$n2 == 20000) {
    result <- power_simulation(
      n_sim     = 1000,
      n_iter    = 20000,
      n_select  = params$n_select,
      prob1     = params$prob1,
      prob2     = 0,
      prob3     = 1,
      true_beta = 0,
      true_beta_SX = 0,
      bs        = 0,
      n1        = 20000,
      n2        = params$n2,
      n3        = params$n3,
      n4        = params$n4,
      a_alpha   = params$a_alpha,
      b_alpha3  = b_alpha3,
      b_alpha4  = b_alpha4,
      beta_U3   = beta_U3,
      beta_U4   = beta_U4
    )
    # Collect type I error results
    type1_results[[length(type1_results) + 1]] <- result
    
  } else if (params$true_beta > 0) {
    result <- power_simulation(
      n_sim     = 200,
      n_iter    = 20000,
      n_select  = params$n_select,
      prob1     = 0.5,
      prob2     = 0,
      prob3     = 1,
      true_beta = params$true_beta,
      true_beta_SX = 0,
      bs        = 0,
      n1        = 20000,
      n2        = params$n2,
      n3        = params$n3,
      n4        = params$n4,
      a_alpha   = params$a_alpha,
      b_alpha3  = b_alpha3,
      b_alpha4  = b_alpha4,
      beta_U3   = beta_U3,
      beta_U4   = beta_U4
    )
    # Collect power results
    power_results[[length(power_results) + 1]] <- result
  }
}

# Combine and save
result_type1error <- do.call(rbind, type1_results)
result_power      <- do.call(rbind, power_results)
save(
  result_power,
  result_type1error,
  file = "simulation.RData"
)

# end of script
