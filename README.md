
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nlss

<!-- badges: start -->
<!-- badges: end -->

Network Latent Source Separation (nlss), a blind source separation
algorithm designed for network data.

## Installation

Install the released version of nlss from Github with:

``` r
devtools::install_github("benwu233/nlss")
```

## Example

This is a basic example which shows you how to use nlss:

First we load the package and simulate networks with a nlss model with
three latent sources.

``` r
library(nlss)  
set.seed(612)

#simulate data with NLSS
sim = sim_NLSS2(n_node = 50, n = 50, alpha_0 = c(3,2,1.5), alpha_1 = 2)
```

Then, we solve the nlss model with the MCMC algorithm:

``` r
res = NLSS(data=sim$X, q=3, k0=0, total_iter = 5000, burn_in = 1000, thin = 10, show_step=1000, joint=TRUE, sprs = 1)
```

We summarize the results:

``` r
sum_res = NLSS_sum(res,th=0.95, k0 = 1, nstart = 1, nend = 400)
```
