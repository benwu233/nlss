# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Simu_DLSS <- function(S, A, K) {
    .Call('_nlss_Simu_DLSS', PACKAGE = 'nlss', S, A, K)
}

Simu_data_bi <- function(X, mu0, sd0, th) {
    .Call('_nlss_Simu_data_bi', PACKAGE = 'nlss', X, mu0, sd0, th)
}

DLSS_Deviance_c <- function(X, S, A, K) {
    .Call('_nlss_DLSS_Deviance_c', PACKAGE = 'nlss', X, S, A, K)
}

DLSS_logLik_noisefree <- function(X, S, A, K) {
    .Call('_nlss_DLSS_logLik_noisefree', PACKAGE = 'nlss', X, S, A, K)
}

NLSS_logLik_noise <- function(X, A, beta, group, K, G) {
    .Call('_nlss_NLSS_logLik_noise', PACKAGE = 'nlss', X, A, beta, group, K, G)
}

NLSS_logLik_noise_group <- function(X, A, beta, group, g0, K, G) {
    .Call('_nlss_NLSS_logLik_noise_group', PACKAGE = 'nlss', X, A, beta, group, g0, K, G)
}

DLSS_logLik_noise_1 <- function(X, S, A, beta, group, K, G) {
    .Call('_nlss_DLSS_logLik_noise_1', PACKAGE = 'nlss', X, S, A, beta, group, K, G)
}

DLSS_logLik_noise0_group <- function(X, S, A, group, g0, K, G) {
    .Call('_nlss_DLSS_logLik_noise0_group', PACKAGE = 'nlss', X, S, A, group, g0, K, G)
}

DLSS_logLik_noise0 <- function(X, S, A, beta, K) {
    .Call('_nlss_DLSS_logLik_noise0', PACKAGE = 'nlss', X, S, A, beta, K)
}

parallelDLSS_update_Y_n <- function(Y, X, A, S, seed, q, p, n, K) {
    invisible(.Call('_nlss_parallelDLSS_update_Y_n', PACKAGE = 'nlss', Y, X, A, S, seed, q, p, n, K))
}

parallelDLSS_update_Y_n_z <- function(Y, X, S, seed, q, p, n, K, alpha_1, alpha_0) {
    invisible(.Call('_nlss_parallelDLSS_update_Y_n_z', PACKAGE = 'nlss', Y, X, S, seed, q, p, n, K, alpha_1, alpha_0))
}

update_S_n <- function(S, X, A, beta, group, q, p, K, G, n) {
    invisible(.Call('_nlss_update_S_n', PACKAGE = 'nlss', S, X, A, beta, group, q, p, K, G, n))
}

update_S_n_z <- function(S, X, Y, beta, group, q, p, K, G, n) {
    invisible(.Call('_nlss_update_S_n_z', PACKAGE = 'nlss', S, X, Y, beta, group, q, p, K, G, n))
}

parallelDLSS_update_A_n <- function(A, X, Y, alpha_1, alpha_0, q, p, n, seed, seed2) {
    invisible(.Call('_nlss_parallelDLSS_update_A_n', PACKAGE = 'nlss', A, X, Y, alpha_1, alpha_0, q, p, n, seed, seed2))
}

update_beta_n <- function(beta, S, group, beta_0, q, p, K, G) {
    invisible(.Call('_nlss_update_beta_n', PACKAGE = 'nlss', beta, S, group, beta_0, q, p, K, G))
}

NLSS_gibbs_sampler_n <- function(X, A0, S0, Y0, beta0, group, gamma, alpha, kk, total_iter = 1000L, burn_in = 100L, thin = 10L, show_step = 50L) {
    .Call('_nlss_NLSS_gibbs_sampler_n', PACKAGE = 'nlss', X, A0, S0, Y0, beta0, group, gamma, alpha, kk, total_iter, burn_in, thin, show_step)
}

