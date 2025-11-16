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
#'
#' @examples
#' # Run BayesSCLingam based on simulated example
#'
#' # Data Generation
#' set.seed(21)
#'
#'
#' N = 300 # Sample size of data
#' num_covariates = 10 # Number of features for test data
#' M = 2 # Number of finite clusters for mixed normal in likelihood
#' num_iter = 10000 # Number of iterations MCMC runs
#'
#' # function to generate DAG example
#' example_list = generates_examples_DAG(num_covariates, N, M, 0.9, 21)
#'
#' data_matrix = example_list$data_matrix # generated data
#' Adjacency_matrix_true = example_list$Adjacency_matrix_true # true graph structure
#'
#'
#' # Input parameters for the BayesSCLingam function
#'
#' params = list(
#' a_mu = 0,
#' b_mu = 2,
#' a_gamma = 0.5,
#' b_gamma = 0.5,
#' a_gamma_1 = 2,
#' b_gamma_1 = 1,
#' a_tao = 2,
#' b_tao = 1,
#' a_og_tao = 0.01,
#' b_og_tao = 0.01,
#' alpha = 1
#' )
#'
#'
#' results_lists = BayesSCLingam(
#' data_matrix,
#' params$a_mu,
#' params$b_mu,
#' params$a_gamma,
#' params$b_gamma,
#' params$a_tao,
#' params$b_tao,
#' params$a_og_tao,
#' params$b_og_tao,
#' params$a_gamma_1,
#' params$b_gamma_1,
#' params$alpha,
#' M,
#' num_iter
#' ) # Runs the Acyclic algorithm
#'
#'
#' Adjacency_matrix_means = results_lists$Adjacency_matrix_means
#' Adjacency_matrix_list = results_lists$Adjacency_matrix_list
#' Causal_effect_matrix_list = results_lists$Causal_effect_matrix_list
#' gamma_list = results_lists$gamma_list
#' gamma_1_list = results_lists$gamma_1_list
#' mu_matrix_list = results_lists$mu_matrix_list
#' tao_matrix_list = results_lists$tao_matrix_list
#' pi_matrix_list = results_lists$pi_matrix_list

BayesDAG <- function(data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_og_tao, b_og_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter) {
  return(BayesSCLingam_cpp(data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_og_tao, b_og_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter))
  #.Call('_cyclinbayes_BayesSCLingam', PACKAGE = 'cyclinbayes', data_matrix, a_mu, b_mu, a_gamma, b_gamma, a_tao, b_tao, a_og_tao, b_og_tao, a_gamma_1, b_gamma_1, alpha, M, num_iter)
}



