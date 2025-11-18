#' Bayesian Cyclic Causal Discovery
#'
#' @description
#' BayesDCG fits a Bayesian linear causal model on a directed cyclic graph (DCG) with non-Gaussian errors, returning posterior samples of the adjacency and causal-effect matrices for systems with feedback.
#'
#' @param data_matrix Numeric matrix of dimension \eqn{N \times p}, where rows correspond to observations and columns correspond to variables (features) included in the causal graph.
#' @param a_mu Hyperparameter for the mean of the normal prior on each mixture component mean \eqn{\mu_k} in the error mixture model (location parameter). Default value is 0.
#' @param b_mu Hyperparameter for the variance of the normal prior on each mixture component mean \eqn{\mu_k} (controls how tightly the component means are shrunk toward \code{a_mu}).
#' Default value is 2.
#' @param a_gamma First shape parameter of the Beta prior on the edge-inclusion probability \eqn{\gamma} (probability that there is an edge from one node to another).
#' Default value is 2.
#' @param b_gamma Second shape parameter of the Beta prior on the edge-inclusion probability \eqn{\gamma}. Along with \code{a_gamma} this controls the expected sparsity of the graph.
#' Default value is 1.
#' @param a_tao Shape parameter of the inverse gamma prior on each mixture component variance in the error distribution (controls the prior tail heaviness for component variances).
#' Default value is 2.
#' @param b_tao Scale parameter of the inverse gamma prior on each mixture component variance in the error distribution (sets the typical size of the component variances).
#' Default value is 1.
#' @param alpha Concentration parameter of the Dirichlet prior on the mixture weights for the normal mixture error distribution (controls how evenly the mixture components are used).
#' Default value is 1
#' @param M Integer giving the maximum number of mixture components allowed in the normal mixture error model.
#' @param num_iter Integer giving the total number of MCMC iterations for the \code{BayesDCG} algorithm.
#' @return A list (from the C backend) typically containing adjacency and causal-effect matrices and samples/diagnostics.
#'
#' @export

BayesDCG <- function(data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter) {
  return(BCD_cpp(data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter))
  #.Call('_cyclinbayes_BCD', PACKAGE = 'cyclinbayes', data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_og_tao, b_og_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter)
}

