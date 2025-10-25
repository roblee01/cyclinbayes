#' Bayesian Lingam Causal Discovery
#'
#' @description
#' A bayesian implementation of the lingam method
#'
#'
#' @param data_matrix numeric matrix (N x p) matrix where N is the sample size and p is the number of features
#' @param a_mu Hyperparameter for the mean of the normal prior used for each of the cluster means in the mixed normal distribution. Default value is 0.
#' @param b_mu Hyperparameter for the variance of the normal prior used for each of the cluster means in the mixed normal distribution. Default value is 2.
#' @param a_gamma Hyperparameter for the first shape of the beta prior for the probability there is a node from one feature to another.
#' @param b_gamma Hyperparameter for the second shape of the beta prior for the probability there is a node from one feature to another.
#' @param a_tao Hyperparameter for the shape of the inverse gamma prior for each of the cluster variance in the mixed normal distribution.
#' @param b_tao Hyperparameter for the scale of the inverse gamma prior for each of the cluster variance in the mixed normal distribution.
#' @param a_og_tao The shape parameter of the proposal when sampling tao with the adjacency matrix entries
#' @param b_og_tao The scale parameter of the proposal when sampling tao with the adjacency matrix entries
#' @param b_og_tao Hyperparameter for the scale of the inverse gamma prior for each of the cluster variance in the mixed normal distribution.
#' @param alpha Hyperparameter for the concentration for the dirichlet prior for the cluster probability of the corresponding feature.
#' @param M Number of possible clusters.
#' @param num_iter The number of iterations for the bayessclingam algorithm
#'
#' @return A list (from the C backend) typically containing adjacency and causal-effect matrices and, optionally, samples/diagnostics.
#'
#' @export

BayesSCLingam <- function(data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_og_tao, b_og_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter) {
  return(BayesSCLingam_cpp(data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_og_tao, b_og_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter))
  #.Call('_cyclinbayes_BayesSCLingam', PACKAGE = 'cyclinbayes', data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_og_tao, b_og_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter)
}



