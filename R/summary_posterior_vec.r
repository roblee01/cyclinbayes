#' Generate summary based on posterior sample provided
#'
#' @description
#' Computes basic posterior summaries for a univariate parameter, including the posterior mean, a central equal-tailed credible interval, and a highest posterior density (HPD) interval.
#'
#' @param posterior_samples Posterior draws for a single parameter.
#' @param level Credible level for the intervals, given as a probability between 0 and 1. For example, \code{0.95} for a 95% interval (default), \code{0.90} for 90%.


#' @return A list with components:
#' \itemize{
#'   \item \code{posterior_mean:} posterior mean of \code{posterior_samples}.
#'   \item \code{credible_interval:} equal-tailed credible interval at the
#'         specified level (e.g., 2.5% and 97.5% for \code{level = 0.95}).
#'   \item \code{hpd_interval:} HPD interval at the same level, as returned by
#'         \code{HDInterval::hdi()}.
#' }
#'
#' @export
#' @examples
#' # Run BayesDAG simulated example
#'
#' N = 300
#' num_covariates = 10
#' M = 2
#' num_iter = 10000
#'
#'
#' example_list = generates_examples_DAG(num_covariates, N, M, 0.9, 21)
#'
#' data_matrix = example_list$data_matrix
#' Adjacency_matrix_true = example_list$Adjacency_matrix_true
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
#' results_lists = BayesDAG(
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
#' gamma_list = results_lists$gamma_list
#' gamma_posterior_summary = summary_posterior_vec(gamma_list,0.95)
#' gamma_posterior_summary$credible_interval
#' gamma_posterior_summary$hpd_interval[,1]



summary_posterior_vec = function(posterior_samples, level){


  credible_interval = quantile(posterior_samples,c((1-level)/2,1-(1-level)/2))
  posterior_mean = mean(posterior_samples)

  hpd_interval = HDInterval::hdi(posterior_samples, credMass = level)

  return(list(posterior_mean = posterior_mean, credible_interval = credible_interval, hpd_interval = hpd_interval))
}
