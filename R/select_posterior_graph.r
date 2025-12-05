#' Generate summary based on
#'
#' @description
#' Selects a representative graph structure from posterior samples by computing a **weighted medoid** under a chosen distance. Intended for lists of sampled adjacency matrices
#'
#' @param posterior_adjacency_matrices A matrix of dimension \code{num_iter × (p*p)} where each row is a flattened adjacency matrix sampled during the posterior (0/1 entries).
#' @param dist_type Character string specifying the distance metric.
#'   One of \code{"shd"}, \code{"sid"}, \code{"forb"}, or \code{"custom"}.
#' @param dist_fun Optional user-supplied distance function used when
#'   \code{dist_type = "custom"}. Must have signature \code{dist_fun(A, B)},
#'   where \code{A} and \code{B} are \eqn{p \times p} adjacency matrices and
#'   return a non negative scalar distance.
#' @param burn_in_frac Fraction of iterations to discard as burn in (default 0.75).
#' @return The best possible graph structure found through finding the smallest distance through finding the weighted medoid based on chosen distance.
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
#'
#' Adjacency_matrix_means = results_lists$Adjacency_matrix_means
#'
#' Adjacency_matrix = select_posterior_graph(Adjacency_matrix_list,dist_type = 'shd') # Best graph structure using shd
#' Adjacency_matrix
#'
#' Adjacency_matrix = select_posterior_graph(Adjacency_matrix_list,dist_type = 'sid') # Best graph structure using sid
#' Adjacency_matrix
#'
#' Adjacency_matrix = select_posterior_graph(Adjacency_matrix_list,dist_type = 'forb') # Best graph structure using forb
#' Adjacency_matrix



select_posterior_graph = function(Adjacency_matrix_list, dist_type = 'shd', dist_fun = NULL, burn_in_frac = 0.75){
  num_covariates = sqrt(ncol(Adjacency_matrix_list))
  num_iter = nrow(Adjacency_matrix_list)
  p = sqrt(ncol(Adjacency_matrix_list))

  keep_inds = seq(floor(burn_in_frac * num_iter), num_iter)
  A_keep = Adjacency_matrix_list[keep_inds, , drop = FALSE]

  structure_strings = apply(A_keep, 1, paste, collapse = ",")
  unique_strings = unique(structure_strings)

  tab = table(factor(structure_strings, levels = unique_strings))
  weights = as.numeric(tab)

  unique_graphs = lapply(unique_strings, function(s) {
    vec = as.numeric(strsplit(s, ",")[[1]])
    matrix(vec, p, p)
  })


  internal_dist = switch(
    dist_type,
    "shd" = function(A, B) {
      g1 <- as(A, "graphNEL")
      g2 <- as(B, "graphNEL")
      pcalg::shd(g1, g2)
    },
    "sid" = function(A, B) {
      # SID requires DAGs
      if (!gRbase::is.DAG(A) || !gRbase::is.DAG(B)) {
        stop("SID distance requires all graphs to be DAGs.")
      }
      SID::structIntervDist(A, B)$sid
    },
    "forb" = function(A, B) {
      base::norm(A - B, type = "F")
    },
    "custom" = {
      if (is.null(dist_fun)) {
        stop("dist_type = 'custom' requires a user-supplied dist_fun(A, B).")
      }
      dist_fun
    }
  )

  v = length(unique_graphs)

  total_distance = numeric(v)

  for (i in seq_len(v)) {
    d_sum = 0
    for (j in seq_len(v)) {
      d_ij = internal_dist(unique_graphs[[i]], unique_graphs[[j]])
      d_sum = d_sum + weights[j] * d_ij
    }
    total_distance[i] = d_sum
  }

  best_index = which.min(total_distance)
  best_adjacency_matrix= unique_graphs[[best_index]]


  return(best_adjacency_matrix)
}
