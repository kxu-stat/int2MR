

devtools::document()   # To generate documentation and NAMESPACE
devtools::build()      # To build the package tarball
devtools::install()    # To install the package locally

library(int2MR)

result_3sample <- int2MR(data_list_3sample = example_3sample_data,
                 model_type = "3sample",
                 prior_inv_gamma_shape = 0.1,
                 prior_inv_gamma_scale = 0.1,
                 chains = 2, iter = 5000, warmup = 1000,
                 adapt_delta = 0.95)
result_3sample$result_3sample

result_2sample <- int2MR(data_list_2sample = example_2sample_data,
                 model_type = "2sample",
                 prior_inv_gamma_shape = 0.1,
                 prior_inv_gamma_scale = 0.1,
                 chains = 2, iter = 5000, warmup = 1000,
                 adapt_delta = 0.95)
result_2sample$result_2sample
