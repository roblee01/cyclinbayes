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
#' result_list = BayesDAG(
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
#' Adjacency_matrix_list = result_list$Adjacency_matrix_list
#'
#' Adjacency_matrix = point_est_graph(Adjacency_matrix_list,dist_type = 'shd') # Best graph structure using shd
#' Adjacency_matrix
#'
#' Adjacency_matrix = point_est_graph(Adjacency_matrix_list,dist_type = 'sid') # Best graph structure using sid
#' Adjacency_matrix
#'
#' custom_edge_mismatch = function(A, B) {
#' return(sum(abs(A - B)))
#' }
#'
#' Adjacency_matrix = point_est_graph(Adjacency_matrix_list, dist_type = 'custom', dist_fun = custom_edge_mismatch)
#' Adjacency_matrix



point_est_graph = function(Adjacency_matrix_list, dist_type = 'shd', dist_fun = NULL, burn_in_frac = 0.75){
  num_covariates = sqrt(ncol(Adjacency_matrix_list))
  num_iter = nrow(Adjacency_matrix_list)
  p = sqrt(ncol(Adjacency_matrix_list))


  start = floor(burn_in_frac * num_iter) + 1L
  keep_inds = start:num_iter
  A_keep = Adjacency_matrix_list[keep_inds, , drop = FALSE]





  structure_strings = apply(A_keep, 1, paste, collapse = ",")
  unique_strings = unique(structure_strings)

  tab = table(factor(structure_strings, levels = unique_strings))
  weights = as.numeric(tab)

  unique_graphs = lapply(unique_strings, function(s) {
    vec = as.numeric(strsplit(s, ",")[[1]])
    matrix(vec, p, p)
  })

  v = length(unique_graphs)

  D = matrix(0, nrow = v, ncol = v)

  if(dist_type == 'shd'){
    graph_list =  lapply(unique_graphs, function(A) as(A, "graphNEL"))

    for (i in seq_len(v - 1L)) {
      gi = graph_list[[i]]
      for (j in (i + 1L):v) {
        d_ij = pcalg::shd(gi, graph_list[[j]])
        D[i, j] = d_ij
        D[j, i] = d_ij
      }
    }
  } else if(dist_type == 'sid'){
    is_dag = vapply(unique_graphs, gRbase::is.DAG, logical(1L))
    if (!all(is_dag)) {
      stop("SID distance requires all graphs to be DAGs.")
    }

    for (i in seq_len(v - 1L)) {
      Ai <- unique_graphs[[i]]
      for (j in (i + 1L):v) {
        d_ij <- SID::structIntervDist(Ai, unique_graphs[[j]])$sid
        D[i, j] <- d_ij
        D[j, i] <- d_ij
      }
    }
  } else if (dist_type == "custom") {
    if (is.null(dist_fun)) {
      stop("dist_type = 'custom' requires a user-supplied dist_fun(A, B).")
    }

    for (i in seq_len(v - 1L)) {
      Ai <- unique_graphs[[i]]
      for (j in (i + 1L):v) {
        d_ij <- dist_fun(Ai, unique_graphs[[j]])
        D[i, j] <- d_ij
        D[j, i] <- d_ij
      }
    }
  } else {
    stop("Unknown dist_type: ", dist_type)
  }

  total_distance = as.vector(D %*% weights)
  best_index = which.min(total_distance)

  best_adjacency_matrix = unique_graphs[[best_index]]
  return(best_adjacency_matrix)


  return(best_adjacency_matrix)
}
