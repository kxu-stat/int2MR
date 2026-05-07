#' Run One Round of int2MR Estimation with Sampling or Optimization
#'
#' This function takes input data lists for the three-sample and/or two-sample Stan models
#' and returns effect estimates from either rstan's sampling function or rstan's
#' optimization function. Standard errors are computed from the inverse negative
#' Hessian in both modes. The user can specify which model(s) to run via the
#' model_type parameter.
#'
#' @param data_list_3sample A named list containing data for the 3-sample Stan model.
#' @param data_list_2sample A named list containing data for the 2-sample Stan model.
#' @param model_type A character string specifying which model to run:
#'        "3sample" to run only the three-sample model,
#'        "2sample" to run only the two-sample model,
#'        or "both" (default) to run both.
#' @param estimation A character string specifying the estimation method:
#'        "sampling" (default) to estimate effects using MCMC posterior means, or
#'        "optimizing" to use Stan's optimization-based MAP estimate with
#'        standard errors from the inverse negative Hessian in both modes.
#' @param prior_inv_gamma_shape Prior shape parameter for the inverse-gamma distributions (default 0.02).
#' @param prior_inv_gamma_scale Prior scale parameter for the inverse-gamma distributions (default 0.02).
#'        Smaller values put more prior mass near sparse/small uncorrelated pleiotropic-effect
#'        variances, while larger values allow more dispersion. As a rule of thumb, start with
#'        `prior_inv_gamma_shape = prior_inv_gamma_scale = 0.02` for weak or sparse uncorrelated
#'        pleiotropic effects, try `0.1` for moderate effects, and check stronger effects with a
#'        small sensitivity grid such as `c(0.02, 0.05, 0.1, 0.2, 0.5)`.
#' @param chains Number of chains for MCMC sampling (default 2).
#' @param iter Total number of iterations per chain (default 5000).
#' @param warmup Number of warmup (burn-in) iterations (default 2000).
#' @param adapt_delta Target acceptance probability for the sampler (default 0.95).
#'
#' @return A list with elements depending on model_type. For each executed model, a data frame summarizing
#'         the estimate, Hessian-based standard error, and p-values for beta, beta_int, and the total effect is provided.
#'
#' @note  
#' When you call `stan_model(model_code = ...)` with the same Stan code string **and** the same
#' `prior_inv_gamma_shape`/`prior_inv_gamma_scale`, rstan will look for a previously compiled
#' binary (in `~/.R/stan/`); if found, it reuses it, otherwise it recompiles.  
#' **To avoid slow repeated compilation** in large-scale runs, either  
#' 1. always pass identical hyperparameters so the cache hits; or  
#' 2. move the two `stan_model()` calls outside of this function (compile once),  
#'    then pass the precompiled model objects into `int2MR()`.  
#'
#' @export
int2MR <- function(data_list_3sample = NULL,
                              data_list_2sample = NULL,
                              model_type = c("both", "3sample", "2sample"),
                              estimation = c("sampling", "optimizing"),
                              prior_inv_gamma_shape = 0.02,
                              prior_inv_gamma_scale = 0.02,
                              chains = 2, iter = 5000, warmup = 2000,
                              adapt_delta = 0.95) {
  # Match the model_type argument
  model_type <- match.arg(model_type)
  estimation <- match.arg(estimation)

  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Package 'rstan' is required to run int2MR.")
  rstan::rstan_options(auto_write = TRUE)

  # Set common control parameters for sampling
  control_list <- list(adapt_delta = adapt_delta)

  results <- list()

  summarize_stan_estimate <- function(est_beta, se_beta, est_beta_int, se_beta_int,
                                      total_effect, se_total) {
    data.frame(
      est_beta = unname(est_beta),
      se_beta = unname(se_beta),
      pval_beta = unname(2 * (1 - stats::pnorm(abs(est_beta) / se_beta))),
      est_beta_int = unname(est_beta_int),
      se_beta_int = unname(se_beta_int),
      pval_beta_int = unname(2 * (1 - stats::pnorm(abs(est_beta_int) / se_beta_int))),
      total_effect = unname(total_effect),
      se_total = unname(se_total),
      pval_total = unname(2 * (1 - stats::pnorm(abs(total_effect) / se_total)))
    )
  }

  summarize_with_hessian <- function(est_beta, est_beta_int, cov_mat) {
    se_beta <- sqrt(cov_mat["beta", "beta"])
    se_beta_int <- sqrt(cov_mat["beta_int", "beta_int"])
    total_effect <- est_beta + est_beta_int
    se_total <- sqrt(cov_mat["beta", "beta"] +
                       cov_mat["beta_int", "beta_int"] +
                       2 * cov_mat["beta", "beta_int"])

    summarize_stan_estimate(
      est_beta = est_beta,
      se_beta = se_beta,
      est_beta_int = est_beta_int,
      se_beta_int = se_beta_int,
      total_effect = total_effect,
      se_total = se_total
    )
  }

  run_stan_estimation <- function(stan_model, data_list) {
    opt <- rstan::optimizing(
      object = stan_model,
      data = data_list,
      hessian = TRUE
    )

    hessian_reg <- opt$hessian
    cov_mat <- MASS::ginv(-hessian_reg)
    rownames(cov_mat) <- rownames(hessian_reg)
    colnames(cov_mat) <- colnames(hessian_reg)

    if (estimation == "optimizing") {
      return(
        summarize_with_hessian(
          est_beta = opt$par["beta"],
          est_beta_int = opt$par["beta_int"],
          cov_mat = cov_mat
        )
      )
    }

    fit <- rstan::sampling(
      object = stan_model,
      data = data_list,
      chains = chains,
      iter = iter,
      warmup = warmup,
      control = control_list,
      refresh = 0
    )
    samples <- rstan::extract(fit)
    summarize_with_hessian(
      est_beta = mean(samples$beta),
      est_beta_int = mean(samples$beta_int),
      cov_mat = cov_mat
    )
  }

  if (model_type %in% c("3sample", "both")) {
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

    stan_mod_3sample <- rstan::stan_model(model_code = stan_model_code_3, verbose = FALSE)

    result_3sample <- run_stan_estimation(stan_mod_3sample, data_list_3sample)

    results$result_3sample <- result_3sample
  }

if (model_type %in% c("2sample", "both")) {
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

  stan_mod_2sample <- rstan::stan_model(model_code = stan_model_code_2, verbose = FALSE)

  result_2sample <- run_stan_estimation(stan_mod_2sample, data_list_2sample)

  results$result_2sample <- result_2sample
}

return(results)
}

#usethis::use_data(example_3sample_data, overwrite = TRUE)
#usethis::use_data(example_2sample_data, overwrite = TRUE)
