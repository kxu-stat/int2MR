#' Run One Round of int2MR Estimation with Sampling and Model Type Option
#'
#' This function takes input data lists for the three-sample and/or two-sample Stan models
#' and returns posterior summaries from running rstan's sampling function. The user can specify
#' which model(s) to run via the model_type parameter.
#'
#' @param data_list_3sample A named list containing data for the 3-sample Stan model.
#' @param data_list_2sample A named list containing data for the 2-sample Stan model.
#' @param model_type A character string specifying which model to run:
#'        "3sample" to run only the three-sample model,
#'        "2sample" to run only the two-sample model,
#'        or "both" (default) to run both.
#' @param prior_inv_gamma_shape Prior shape parameter for the inverse-gamma distributions (default 0.2).
#' @param prior_inv_gamma_scale Prior scale parameter for the inverse-gamma distributions (default 0.2).
#' @param chains Number of chains for MCMC sampling (default 4).
#' @param iter Total number of iterations per chain (default 2000).
#' @param warmup Number of warmup (burn-in) iterations (default 1000).
#' @param adapt_delta Target acceptance probability for the sampler (default 0.95).
#'
#' @return A list with elements depending on model_type. For each executed model, a data frame summarizing
#'         the posterior mean, standard deviation, and p-values for beta, beta_int, and the total effect is provided.
#' @export
int2MR <- function(data_list_3sample = NULL,
                              data_list_2sample = NULL,
                              model_type = c("3sample", "2sample"),
                              prior_inv_gamma_shape = 0.1,
                              prior_inv_gamma_scale = 0.1,
                              chains = 2, iter = 5000, warmup = 2500,
                              adapt_delta = 0.95) {
  # Match the model_type argument
  model_type <- match.arg(model_type)

  # Load required packages
  require(rstan)
  rstan_options(auto_write = TRUE)

  # Set common control parameters for sampling
  control_list <- list(adapt_delta = adapt_delta)

  results <- list()

  if (model_type %in% c("3sample")) {
    if (is.null(data_list_3sample))
      stop("data_list_3sample must be provided for the three-sample model.")

    # Compile the Stan model for the 3-sample setting.
    stan_model_code_3 <- sprintf("
data {
  int<lower=1> p;
  vector[p] hat_s1_sq;
  vector[p] hat_s2_sq;
  vector[p] hat_s3_sq;
  vector[p] hat_gamma;
  vector[p] hat_Gamma1;
  vector[p] hat_Gamma2;
  vector[p] hat_Gamma3;
  vector[p] hat_s_gamma_sq;
  real rho1;
  real rho2;
  real rho3;
}
parameters {
  real beta;
  real beta_int;
  vector[p] gamma;
  vector[p] alpha1;
  vector[p] alpha2;
  vector[p] alpha3;
  real<lower=0> sigma_alpha1_sq;
  real<lower=0> sigma_alpha2_sq;
  real<lower=0> sigma_alpha3_sq;
  real<lower=0> sigma_gamma_sq;
}
transformed parameters {
  real<lower=0> sigma_alpha1 = sqrt(sigma_alpha1_sq);
  real<lower=0> sigma_alpha2 = sqrt(sigma_alpha2_sq);
  real<lower=0> sigma_alpha3 = sqrt(sigma_alpha3_sq);
  real<lower=0> sigma_gamma = sqrt(sigma_gamma_sq);
}
model {
  // Priors for variance parameters
  sigma_gamma_sq ~ inv_gamma(%f, %f);
  sigma_alpha1_sq ~ inv_gamma(%f, %f);
  sigma_alpha2_sq ~ inv_gamma(%f, %f);
  sigma_alpha3_sq ~ inv_gamma(%f, %f);

  // Priors for gamma and alpha
  gamma ~ normal(0, sigma_gamma);
  alpha1 ~ normal(0, sigma_alpha1);
  alpha2 ~ normal(0, sigma_alpha2);
  alpha3 ~ normal(0, sigma_alpha3);

  // Observation models
  hat_gamma ~ normal(gamma, sqrt(hat_s_gamma_sq));
  hat_Gamma1 ~ normal((beta + rho1 * beta_int) .* gamma + alpha1, sqrt(hat_s1_sq));
  hat_Gamma2 ~ normal((beta + rho2 * beta_int) .* gamma + alpha2, sqrt(hat_s2_sq));
  hat_Gamma3 ~ normal((beta + rho3 * beta_int) .* gamma + alpha3, sqrt(hat_s3_sq));
}
", prior_inv_gamma_shape, prior_inv_gamma_scale,
prior_inv_gamma_shape, prior_inv_gamma_scale,
prior_inv_gamma_shape, prior_inv_gamma_scale,
prior_inv_gamma_shape, prior_inv_gamma_scale)

    stan_mod_3sample <- stan_model(model_code = stan_model_code_3, verbose = FALSE)

    # Run sampling for the 3-sample model.
    fit_3 <- sampling(object = stan_mod_3sample,
                      data = data_list_3sample,
                      chains = chains,
                      iter = iter,
                      warmup = warmup,
                      control = control_list,
                      refresh = 0)

    opt_3 <- optimizing(
      object = stan_mod_3sample,
      data = data_list_3sample,
      hessian = TRUE)

    hessian_reg_3 <- opt_3$hessian
    cov_3 <- MASS::ginv(-hessian_reg_3)
    rownames(cov_3) <- rownames(hessian_reg_3)
    colnames(cov_3) <- colnames(hessian_reg_3)

    # Extract posterior samples and compute summaries.
    samples_3 <- rstan::extract(fit_3)
    est_beta_3 <- mean(samples_3$beta)
    se_beta_3 <- sqrt(cov_3["beta", "beta"])
    est_beta_int_3 <- mean(samples_3$beta_int)
    se_beta_int_3 <-  sqrt(cov_3["beta_int", "beta_int"])
    total_3 <- mean(samples_3$beta+samples_3$beta_int)
    se_total_3 <- sqrt(cov_3["beta_int", "beta_int"]
                       + cov_3["beta", "beta"]
                       + 2 * cov_3["beta_int", "beta"])

    result_3sample <- data.frame(
      est_beta = est_beta_3,
      se_beta = se_beta_3,
      pval_beta = 2 * (1 - pnorm(abs(est_beta_3)/se_beta_3)),
      est_beta_int = est_beta_int_3,
      se_beta_int = se_beta_int_3,
      pval_beta_int = 2 * (1 - pnorm(abs(est_beta_int_3)/se_beta_int_3)),
      total_effect = total_3,
      pval_total = 2 * (1 - pnorm(abs(total_3)/se_total_3))
    )

    results$result_3sample <- result_3sample
  }

if (model_type %in% c("2sample")) {
  if (is.null(data_list_2sample))
    stop("data_list_2sample must be provided for the two-sample model.")

  # Compile the Stan model for the 2-sample setting.
  stan_model_code_2 <- sprintf("
data {
  int<lower=1> p;
  vector[p] hat_s1_sq;
  vector[p] hat_s2_sq;
  vector[p] hat_gamma;
  vector[p] hat_Gamma1;
  vector[p] hat_Gamma2;
  vector[p] hat_s_gamma_sq;
  real rho1;
  real rho2;
}
parameters {
  real beta_int;
  real beta;
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
  real<lower=0> sigma_gamma = sqrt(sigma_gamma_sq);
}
model {
  sigma_gamma_sq ~ inv_gamma(%f, %f);
  sigma_alpha1_sq ~ inv_gamma(%f, %f);
  sigma_alpha2_sq ~ inv_gamma(%f, %f);

  gamma ~ normal(0, sigma_gamma);
  alpha1 ~ normal(0, sigma_alpha1);
  alpha2 ~ normal(0, sigma_alpha2);

  hat_gamma ~ normal(gamma, sqrt(hat_s_gamma_sq));
  hat_Gamma1 ~ normal((beta + rho1 * beta_int) .* gamma + alpha1, sqrt(hat_s1_sq));
  hat_Gamma2 ~ normal((beta + rho2 * beta_int) .* gamma + alpha2, sqrt(hat_s2_sq));
}
", prior_inv_gamma_shape, prior_inv_gamma_scale,
prior_inv_gamma_shape, prior_inv_gamma_scale,
prior_inv_gamma_shape, prior_inv_gamma_scale)

  stan_mod_2sample <- stan_model(model_code = stan_model_code_2, verbose = FALSE)

  # Run sampling for the 2-sample model.
  fit_2 <- sampling(object = stan_mod_2sample,
                    data = data_list_2sample,
                    chains = chains,
                    iter = iter,
                    warmup = warmup,
                    control = control_list,
                    refresh = 0)

  opt_2 <- optimizing(
    object = stan_mod_2sample,
    data = data_list_2sample,
    hessian = TRUE)

  hessian_reg_2 <- opt_2$hessian
  cov_2 <- MASS::ginv(-hessian_reg_2)
  rownames(cov_2) <- rownames(hessian_reg_2)
  colnames(cov_2) <- colnames(hessian_reg_2)

  # Extract posterior samples and compute summaries.
  samples_2 <- rstan::extract(fit_2)
  est_beta_2 <- mean(samples_2$beta)
  se_beta_2 <- sqrt(cov_2["beta", "beta"])
  est_beta_int_2 <- mean(samples_2$beta_int)
  se_beta_int_2 <-  sqrt(cov_2["beta_int", "beta_int"])
  total_2 <- mean(samples_2$beta+samples_2$beta_int)
  se_total_2 <- sqrt(cov_2["beta_int", "beta_int"]
                     + cov_2["beta", "beta"]
                     + 2 * cov_2["beta_int", "beta"])

  result_2sample <- data.frame(
    est_beta = est_beta_2,
    se_beta = se_beta_2,
    pval_beta = 2 * (1 - pnorm(abs(est_beta_2)/se_beta_2)),
    est_beta_int = est_beta_int_2,
    se_beta_int = se_beta_int_2,
    pval_beta_int = 2 * (1 - pnorm(abs(est_beta_int_2)/se_beta_int_2)),
    total_effect = total_2,
    pval_total = 2 * (1 - pnorm(abs(total_2)/se_total_2))
  )

  results$result_2sample <- result_2sample
}

return(results)
}

#usethis::use_data(example_3sample_data, overwrite = TRUE)
#usethis::use_data(example_2sample_data, overwrite = TRUE)
