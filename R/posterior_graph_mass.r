#' Compute Posterior Mass of a Given Graph Structure
#'
#' @description
#' This function computes the posterior mass (relative frequency) of a specified graph structure among the collection of
#' sampled posterior adjacency matrices. The input graph is converted into an adjacency matrix vector, and function
#' counts how often all edges appear in sampled posterior graph.
#'
#'
#' @param graph_structure An \code{igraph} object representing the target graph.
#' @param posterior_graph_structures A matrix of dimension \code{num_iter × p^2}, where each row corresponds
#' to a flattened adjacency matrix from one MCMC iteration.
#'
#' @return A numeric value between 0 and 1 giving the posterior mass (frequency) with which the target
#' graph's edges appear in the posterior samples.
#'
#' @export

posterior_graph_mass = function(graph_structure, posterior_graph_structures){
  count = 0
  A = as.matrix(as_adjacency_matrix(graph_structure))
  graph_structure_vec = c(A)
  causal_effects_index = which(graph_structure_vec == 1)

  for(i in 1:nrow(posterior_graph_structures)){
    causal_effects_posterior_index = which(posterior_graph_structures[i,]==1)
    all(causal_effects_index %in% causal_effects_posterior_index)
    count = count + 1
  }
  return(count/nrow(posterior_graph_structures))
}

