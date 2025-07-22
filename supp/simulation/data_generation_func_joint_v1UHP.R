
# gibbs_v5/data_generation_func_single_cml_v2.R

# function to generate two-population GWAS summary statistics

# dgm = function (
    # K, sig_01, sig_02, rho_0, sig_x1, sig_x2, rho_1, 
# a_f_1, b_f_1, a_f_2, b_f_2,
# 
# mu_a1, mu_a2, mu_b1, mu_b2, sig_a1, sig_a2, sig_b1, sig_b2, 
# sig_t, rho_t, rho_e, q1, q2, sig_y1, sig_y2, rho_2
# ) {

# generate the same number of indep. IVs: m
# rho_gamma: Target (Spearman) correlation
sourceCpp("~/data_generation_func_joint_v1UHP.cpp")
# generate summary statistics

dgm_summary = function(GWAS_individual) {
  gammaSS1 = fastSigLm(GWAS_individual$X1, t(GWAS_individual$g1x))
  b_exp_1 = gammaSS1$coef; se_exp_1 = gammaSS1$std
  GammaSS2 = fastSigLm(GWAS_individual$Y2, t(GWAS_individual$g2y))
  b_out_2 = GammaSS2$coef; se_out_2 = GammaSS2$std  
  GammaSS3 = fastSigLm(GWAS_individual$Y3, t(GWAS_individual$g3y))
  b_out_3 = GammaSS3$coef; se_out_3 = GammaSS3$std  
  GammaSS4 = fastSigLm(GWAS_individual$Y4, t(GWAS_individual$g4y))
  b_out_4 = GammaSS4$coef; se_out_4 = GammaSS4$std  
  return(list(b_exp_1 = b_exp_1, se_exp_1 = se_exp_1, 
              b_out_2 = b_out_2, se_out_2 = se_out_2,
              b_out_3 = b_out_3, se_out_3 = se_out_3,
              b_out_4 = b_out_4, se_out_4 = se_out_4,
              g3y = GWAS_individual$g3y, X3y = GWAS_individual$X3y, Y3 = GWAS_individual$Y3, I3y = GWAS_individual$I3y, 
              g4y = GWAS_individual$g4y, X4y = GWAS_individual$X4y, Y4 = GWAS_individual$Y4, I4y = GWAS_individual$I4y,
              varcomp_X1x = GWAS_individual$varcomp_X1x, varcomp_U1x = GWAS_individual$varcomp_U1x, 
              varcomp_X2y = GWAS_individual$varcomp_X2y, 
              varcomp_U2y = GWAS_individual$varcomp_U2y, varcomp_UHP2y = GWAS_individual$varcomp_UHP2y, 
              varcomp_X3y = GWAS_individual$varcomp_X3y, 
              varcomp_U3y = GWAS_individual$varcomp_U3y, varcomp_UHP3y = GWAS_individual$varcomp_UHP3y, 
              varcomp_X4y = GWAS_individual$varcomp_X4y, 
              varcomp_U4y = GWAS_individual$varcomp_U4y, varcomp_UHP4y = GWAS_individual$varcomp_UHP4y,
              varcomp_X2y_interact = GWAS_individual$varcomp_X2y_interact,
              gamma = GWAS_individual$gamma))
}

# wrap up
dgm_joint = function(m, n1x, n2y, n3y, n4y,
                     a_f_1, b_f_1,
                     a_alpha, 
                     b_alpha3, b_alpha4,
                     beta, beta_int, prob,
                     beta_U3 = 1, beta_U4 = 1)
{
  GWAS_individual = dgm_individual_cpp(m, n1x, n2y, n3y, n4y,
                                   a_f_1, b_f_1,
                                   a_alpha, 
                                   b_alpha3, b_alpha4,
                                   beta, beta_int, prob,
                                   beta_U3, beta_U4)
  
  GWAS_summary = dgm_summary(GWAS_individual)
  return(GWAS_summary)
}

