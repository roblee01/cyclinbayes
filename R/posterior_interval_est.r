#' Generate summary based on posterior sample provided
#'
#' @description
#' Compute posterior summaries from a matrix of MCMC samples. When \code{adjacency = TRUE}, the function treats each row as a flattened
#' adjacency matrix and returns edge wise posterior inclusion probabilities and posterior probabilities for unique graph structures. Otherwise,
#' the function treats each column as a parameter and generates a credible and hpd interval based on the level input.
#'
#' @param posterior_matrix numeric posterior sample matrix, where each row corresponds to the MCMC iteration
#' @param level Credible level for the intervals, given as a probability between 0 and 1. Only utilized for posterior analysis of non graph structure samples.
#' @return
#' If \code{adjacency = TRUE}, a list with components:
#' \itemize{
#'    \item \code{pip_matrix:} numeric matrix (P), with edge wise posterior inclusion probabilities, where entry P_{ij} is the proportion of sampled graphs contain an edge from i to j.
#'    \item \code{pip_graph_results:} list in which each element corresponds to a unique graph structure and the proportion of times that structure appears in the posterior graph samples.
#' }
#'
#'
#' If \code{adjacency = FALSE}, a list with components:
#' \itemize{
#'   \item \code{ci_matrix:} numeric matrix, where each row is the equal-tailed credible interval for each parameter provide in the posterior matrix at the
#'         specified level (e.g., 2.5\% and 97.5\% for \code{level = 0.95}).
#'   \item \code{hpd_matrix:} numeric matrix, where each row is the HPD interval for each parameter at the same level.
#' }
#' @export
#'
#' @examples
#'
#' # Run BayessDAG simulated example
#'
#' N = 300
#' num_covariates = 10
#' M = 2
#' num_iter = 5000
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
#'
#' # Extract posterior outputs
#' Adjacency_matrix_means = results_list$Adjacency_matrix_means
#' Adjacency_matrix_list = results_list$Adjacency_matrix_list
#' Causal_effect_matrix_list = results_list$Causal_effect_matrix_list
#' gamma_list = results_list$gamma_list
#' gamma_1_list = results_list$gamma_1_list
#' mu_matrix_list = results_list$mu_matrix_list
#' tao_matrix_list = results_list$tao_matrix_list
#' pi_matrix_list = results_list$pi_matrix_list
#'
#' # Compute HPD and CI summaries for causal effect matrix
#' Causal_effect_matrix_summary = posterior_interval_est(Causal_effect_matrix_list, level = 0.95)
#' hpd_matrix_cyclic = Causal_effect_matrix_summary$hpd_matrix
#' ci_matrix_cyclic = Causal_effect_matrix_summary$ci_matrix

posterior_interval_est = function(posterior_matrix, level, adjacency = FALSE){

  num_iter = nrow(posterior_matrix)

  if(adjacency){
    num_iter = nrow(posterior_matrix)
    num_features = sqrt(ncol(posterior_matrix))
    Adjacency_matrix_list_75 = posterior_matrix[(0.75*num_iter):num_iter,]
    pip_matrix = matrix(colMeans(Adjacency_matrix_list_75), num_features, num_features)


    denominator_val = nrow(Adjacency_matrix_list_75)
    probabilities_vec = rep(0,denominator_val)
    for(i in 1:denominator_val){
      probabilities_vec[i] = paste(Adjacency_matrix_list_75[i,],collapse='')
    }
    pip_graph_results = sort(table(probabilities_vec),decreasing=TRUE)/denominator_val

    return(list(pip_matrix = pip_matrix, pip_graph_results = pip_graph_results))
  }

  ci_matrix = matrix(0,ncol(posterior_matrix),3)

  if(ncol(posterior_matrix) == 1){
    posterior_75 = as.matrix(posterior_matrix[(0.75*num_iter):num_iter,])
  } else{
    posterior_75 = posterior_matrix[(0.75*num_iter):num_iter,]
  }

  hpd_matrix = HDInterval::hdi(posterior_75, credMass = level)

  for(i in 1:ncol(posterior_matrix)){
    ci_matrix[i,] =  quantile(posterior_75[,i],c((1-level)/2,0.5,1-(1-level)/2))
  }

  return(list(hpd_matrix = hpd_matrix, ci_matrix = ci_matrix))
}
