#include <RcppArmadillo.h>
#include <random>
#include <math.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include <ctime>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// Performs a fast linear regression of y on X (intercept + single predictor)
void fastLm(const arma::vec & y, const arma::mat & X, double &coefb, double &stdb) {
  int n = X.n_rows;
  int k = X.n_cols;
  
  // Solve for coefficients
  arma::colvec coef = arma::solve(X, y);
  arma::colvec resid = y - X * coef;
  
  // Estimate variance of residuals
  double sig2 = as_scalar(trans(resid) * resid / (n - k));
  // Standard errors of coefficients
  arma::colvec stderrest = sqrt(sig2 * diagvec(inv(trans(X) * X)));
  
  // Return the slope coefficient and its standard error (ignore intercept)
  coefb = coef[1];
  stdb = stderrest[1];
#include <RcppArmadillo.h>
#include <random>
#include <math.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include <ctime>
#include <cmath>
  
  // [[Rcpp::depends(RcppArmadillo)]]
  using namespace Rcpp;
  using namespace arma;
  using namespace std;
  
  // Performs a fast linear regression of y on X (intercept + single predictor)
  void fastLm(const arma::vec & y, const arma::mat & X, double &coefb, double &stdb) {
    int n = X.n_rows;
    int k = X.n_cols;
    
    // Solve for coefficients
    arma::colvec coef = arma::solve(X, y);
    arma::colvec resid = y - X * coef;
    
    // Estimate variance of residuals
    double sig2 = as_scalar(trans(resid) * resid / (n - k));
    // Standard errors of coefficients
    arma::colvec stderrest = sqrt(sig2 * diagvec(inv(trans(X) * X)));
    
    // Return the slope coefficient and its standard error (ignore intercept)
    coefb = coef[1];
    stdb = stderrest[1];
  }
  
  // [[Rcpp::export]]
  // Computes slope coefficients and standard errors for each column of X against y
  List fastSigLm(const arma::vec & y, const arma::mat & X) {
    int p = X.n_cols;
    int n = X.n_rows;
    arma::mat xx = zeros(p, 2);
    double coefb = 0;
    double stdb = 0;
    vec coef = zeros(p, 1);
    vec std = zeros(p, 1);
    
    for (int j = 0; j < p; ++j) {
      // Build design matrix with intercept and j-th predictor
      xx = join_rows(ones(n, 1), X.col(j));
      fastLm(y, xx, coefb, stdb);
      coef[j] = coefb;
      std[j]  = stdb;
    }
    
    return List::create(
      Named("coef") = coef,
      Named("std")  = std
    );
  }
  
  // [[Rcpp::export]]
  // Generate a gamma vector of length n: values drawn uniformly from two intervals
  arma::vec generate_gamma(int n) {
    arma::vec gamma(n);
    for (int i = 0; i < n; ++i) {
      double u = R::runif(0, 1);
      if (u < 0.5) {
        // Sample from [0.1, 0.2]
        gamma[i] = R::runif(0.1, 0.2);
      } else {
        // Sample from [-0.2, -0.1]
        gamma[i] = R::runif(-0.2, -0.1);
      }
    }
    return gamma;
  }
  
  // Function to generate an m x n matrix of binomial(2, f_i) entries; f is a vector of probabilities
  arma::Mat<int> generate_binomial_matrix(int m, int n, const arma::vec & f) {
    arma::Mat<int> matrix(m, n);
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0; i < m; ++i) {
      std::binomial_distribution<> binom_dist(2, f[i]);
      for (int j = 0; j < n; ++j) {
        matrix(i, j) = binom_dist(gen);
      }
    }
    return matrix;
  }
  
  // [[Rcpp::export]]
  // Data-generating function: simulates individual-level genetic and covariate data across four studies
  List dgm_individual_cpp(
      int m, int n1x, int n2y, int n3y, int n4y,
      double a_f_1, double b_f_1,
      double a_alpha,
      double b_alpha3, double b_alpha4,
      double beta, double beta_int, double prob,
      double beta_U3 = 1, double beta_U4 = 1
  ) {
    // Generate gamma for m variants
    vec gamma = generate_gamma(m);
    
    // Generate allele frequency vectors f1..f4
    vec f1 = a_f_1 + (b_f_1 - a_f_1) * randu<vec>(m);
    vec f2 = a_f_1 + (b_f_1 - a_f_1) * randu<vec>(m);
    vec f3 = a_f_1 + (b_f_1 - a_f_1) * randu<vec>(m);
    vec f4 = a_f_1 + (b_f_1 - a_f_1) * randu<vec>(m);
    
    // Generate genotype matrices for each study
    Mat<int> g1x = generate_binomial_matrix(m, n1x, f1);
    Mat<int> g2y = generate_binomial_matrix(m, n2y, f2);
    Mat<int> g3y = generate_binomial_matrix(m, n3y, f3);
    Mat<int> g4y = generate_binomial_matrix(m, n4y, f4);
    
    // Genetic effect heterogeneity (alpha) vectors for studies 2-4
    vec alpha_3_2 = b_alpha3 * randn<vec>(m);
    vec alpha_4_2 = b_alpha4 * randn<vec>(m);
    vec alpha_3_3 = b_alpha3 * randn<vec>(m);
    vec alpha_4_3 = b_alpha4 * randn<vec>(m);
    vec alpha_3_4 = b_alpha3 * randn<vec>(m);
    vec alpha_4_4 = b_alpha4 * randn<vec>(m);
    
    // Study 1: generate X1
    vec U1x = randn<vec>(n1x);
    vec I1x = conv_to<vec>::from(randu<vec>(n1x) < prob);
    vec X1x = (gamma.t() * g1x).t()
      + beta_U3 * U1x
    + beta_U4 * U1x % I1x
    + a_alpha * randn<vec>(n1x);
    double varcomp_X1x = var((gamma.t() * g1x).t()) / var(X1x);
    
    // Study 2: generate X2 and Y2
    vec U2y = randn<vec>(n2y);
    vec I2y = conv_to<vec>::from(randu<vec>(n2y) < prob);
    vec X2y = (gamma.t() * g2y).t()
      + beta_U3 * U2y
    + beta_U4 * U2y % I2y
    + a_alpha * randn<vec>(n2y);
    vec Y2 = (alpha_3_2.t() * g2y).t()
      + ((alpha_4_2.t() * g2y).t()) % I2y
    + beta_U3 * U2y
    + beta_U4 * U2y % I2y
    + beta * X2y
    + beta_int * X2y % I2y
    + a_alpha * randn<vec>(n2y);
    double varcomp_X2y = var(beta * X2y) / var(Y2);
    double varcomp_X2y_interact = var(beta_int * X2y % I2y) / var(Y2);
    double varcomp_U2y = var(beta_U3 * U2y + beta_U4 * U2y % I2y) / var(Y2);
    double varcomp_UHP2y = var((alpha_3_2.t() * g2y).t() + ((alpha_4_2.t() * g2y).t()) % I2y) / var(Y2);
    
    // Study 3: generate X3 and Y3 (no interaction)
    vec U3y = randn<vec>(n3y);
    vec I3y = zeros<vec>(n3y);
    vec X3y = (gamma.t() * g3y).t()
      + beta_U3 * U3y
    + beta_U4 * U3y % I3y
    + a_alpha * randn<vec>(n3y);
    vec Y3 = (alpha_3_3.t() * g3y).t()
      + ((alpha_4_3.t() * g3y).t()) % I3y
    + beta_U3 * U3y
    + beta_U4 * U3y % I3y
    + beta * X3y
    + beta_int * X3y % I3y
    + a_alpha * randn<vec>(n3y);
    double varcomp_X3y = var(beta * X3y + beta_int * X3y % I3y) / var(Y3);
    double varcomp_U3y = var(beta_U3 * U3y + beta_U4 * U3y % I3y) / var(Y3);
    double varcomp_UHP3y = var((alpha_3_3.t() * g3y).t() + ((alpha_4_3.t() * g3y).t()) % I3y) / var(Y3);
    
    // Study 4: generate X4 and Y4 (full interaction)
    vec U4y = randn<vec>(n4y);
    vec I4y = ones<vec>(n4y);
    vec X4y = (gamma.t() * g4y).t()
      + beta_U3 * U4y
    + beta_U4 * U4y % I4y
    + a_alpha * randn<vec>(n4y);
    vec Y4 = (alpha_3_4.t() * g4y).t()
      + ((alpha_4_4.t() * g4y).t()) % I4y
    + beta_U4 * U4y % I4y
    + beta * X4y
    + beta_int * X4y % I4y
    + a_alpha * randn<vec>(n4y);
    double varcomp_X4y = var(beta * X4y) / var(Y4);
    double varcomp_U4y = var(beta_U3 * U4y + beta_U4 * U4y % I4y) / var(Y4);
    double varcomp_UHP4y = var((alpha_3_4.t() * g4y).t() + ((alpha_4_4.t() * g4y).t()) % I4y) / var(Y4);
    
    // Return all simulated variables and variance components
    List final_list;
    final_list["X1"] = X1x;
    final_list["X2y"] = X2y;
    final_list["Y2"] = Y2;
    final_list["X3y"] = X3y;
    final_list["Y3"] = Y3;
    final_list["X4y"] = X4y;
    final_list["Y4"] = Y4;
    final_list["g1x"] = g1x;
    final_list["g2y"] = g2y;
    final_list["g3y"] = g3y;
    final_list["g4y"] = g4y;
    final_list["I3y"] = I3y;
    final_list["I4y"] = I4y;
    
    final_list["varcomp_X1x"] = varcomp_X1x;
    final_list["varcomp_X2y"] = varcomp_X2y;
    final_list["varcomp_U2y"] = varcomp_U2y;
    final_list["varcomp_UHP2y"] = varcomp_UHP2y;
    final_list["varcomp_X3y"] = varcomp_X3y;
    final_list["varcomp_U3y"] = varcomp_U3y;
    final_list["varcomp_UHP3y"] = varcomp_UHP3y;
    final_list["varcomp_X4y"] = varcomp_X4y;
    final_list["varcomp_U4y"] = varcomp_U4y;
    final_list["varcomp_UHP4y"] = varcomp_UHP4y;
    final_list["varcomp_X2y_interact"] = varcomp_X2y_interact;
    final_list["gamma"] = gamma;
    
    return final_list;
  }
  
