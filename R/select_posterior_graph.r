#' Generate summary based on
#'
#' @description
#' Selects a representative graph structure from posterior samples by computing a **weighted medoid** under a chosen distance. Intended for lists of sampled adjacency matrices
#'
#' @param posterior_adjacency_matrices A matrix of dimension \code{num_iter × (p*p)} where each row is a flattened adjacency matrix sampled during the posterior (0/1 entries).
#' @param dist_type Type of distance to use when selecting the medoid (e.g., `"shd"`, `"sid"`,
#'   `"forb"`).

#' @return The best possible graph structure found through finding the smallest distance through finding the weighted medoid based on chosen distance.
#'
#' @export
#'
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



select_posterior_graph = function(Adjacency_matrix_list, dist_type){
  num_covariates = sqrt(ncol(Adjacency_matrix_list))
  num_iter = nrow(Adjacency_matrix_list)

  if(dist_type == 'shd'){

    Adjacency_matrix_75 = Adjacency_matrix_list[(0.75*num_iter):num_iter,]
    each_graph_structure = rep(0,nrow(Adjacency_matrix_75))
    for(i in 1:nrow(Adjacency_matrix_75)){
      each_graph_structure[i] = paste(Adjacency_matrix_75[i,], collapse=',')
    }

    unique_graph_structures = unique(each_graph_structure)


    tab = table(factor(each_graph_structure, levels = unique_graph_structures))

    possible_pairs = combn(1:length(unique_graph_structures),2)

    unique_graph_structures_mat = matrix(0,length(unique_graph_structures),ncol(Adjacency_matrix_75))


    for(i in 1:length(unique_graph_structures)){
      unique_graph_structures_mat[i,] = as.numeric(strsplit(unique_graph_structures[i][1],',')[[1]])
    }

    shd_vec = rep(0,ncol(possible_pairs))

    structure_matrix = matrix(0, num_covariates, num_covariates)

    for(i in 1:ncol(possible_pairs)){
      current_graph_pairs = t(possible_pairs)[i,]

      graph_1_vec = unique_graph_structures_mat[current_graph_pairs[1],]
      structure_matrix[] = graph_1_vec
      #graph_1 = matrix(graph_1_vec,num_covariates,num_covariates)
      graph_1 = as(structure_matrix, "graphNEL")

      graph_2_vec = unique_graph_structures_mat[current_graph_pairs[2],]
      structure_matrix[] = graph_2_vec
      #graph_2 = matrix(graph_2_vec,num_covariates,num_covariates)
      graph_2 = as(structure_matrix, "graphNEL")

      shd_vec[i] = pcalg::shd(graph_1,graph_2)
    }

    shd_per_unique = rep(0, length(unique_graph_structures))


    for(i1 in 1:length(unique_graph_structures)){

      sum_portion_1 = sum(shd_vec[which(t(possible_pairs)[,1]==i1)]*tab[t(possible_pairs)[which(t(possible_pairs)[,1]==i1),2]])

      sum_portion_2 = sum(shd_vec[which(t(possible_pairs)[,2]==i1)]*tab[t(possible_pairs)[which(t(possible_pairs)[,2]==i1),1]])

      shd_per_unique[i1] = sum_portion_1 + sum_portion_2
    }

    best_adjacency_mat = matrix(unique_graph_structures_mat[which.min(shd_per_unique),],num_covariates,num_covariates)
  } else if(dist_type == 'sid'){
    Adjacency_matrix_75 = Adjacency_matrix_list[(0.75*num_iter):num_iter,]
    each_graph_structure = rep(0,nrow(Adjacency_matrix_75))
    for(i in 1:nrow(Adjacency_matrix_75)){
      each_graph_structure[i] = paste(Adjacency_matrix_75[i,], collapse=',')
    }

    unique_graph_structures = unique(each_graph_structure)

    tab = table(factor(each_graph_structure, levels = unique_graph_structures))

    possible_pairs = combn(1:length(unique_graph_structures),2)

    unique_graph_structures_mat = matrix(0,length(unique_graph_structures),ncol(Adjacency_matrix_75))

    structure_matrix = matrix(0, num_covariates, num_covariates)

    for(i in 1:length(unique_graph_structures)){
      graph_structure = as.numeric(strsplit(unique_graph_structures[i][1],',')[[1]])
      structure_matrix[] = graph_structure
      if(gRbase::is.DAG(structure_matrix)){
        unique_graph_structures_mat[i,] = as.numeric(strsplit(unique_graph_structures[i][1],',')[[1]])
      } else{
        return('Graph structure needs to be DAG')
      }
    }

    sid_vec = rep(0,ncol(possible_pairs))

    for(i in 1:ncol(possible_pairs)){
      current_graph_pairs = t(possible_pairs)[i,]

      graph_1_vec = unique_graph_structures_mat[current_graph_pairs[1],]
      structure_matrix[] = graph_1_vec
      #graph_1 = matrix(graph_1_vec,num_covariates,num_covariates)
      graph_1 = structure_matrix

      graph_2_vec = unique_graph_structures_mat[current_graph_pairs[2],]
      structure_matrix[] = graph_2_vec
      #graph_2 = matrix(graph_2_vec,num_covariates,num_covariates)
      graph_2 = structure_matrix

      sid_result = SID::structIntervDist(graph_1,graph_2)

      sid_vec[i] = sid_result
    }

    sid_vec = rep(0,ncol(possible_pairs))

    for(i in 1:ncol(possible_pairs)){
      current_graph_pairs = t(possible_pairs)[i,]

      graph_1_vec = unique_graph_structures_mat[current_graph_pairs[1],]
      structure_matrix[] = graph_1_vec
      #graph_1 = matrix(graph_1_vec,num_covariates,num_covariates)
      graph_1 = structure_matrix

      graph_2_vec = unique_graph_structures_mat[current_graph_pairs[2],]
      structure_matrix[] = graph_2_vec
      #graph_2 = matrix(graph_2_vec,num_covariates,num_covariates)
      graph_2 = structure_matrix

      sid_result = SID::structIntervDist(graph_1,graph_2)

      sid_vec[i] = sid_result$sid
    }

    sid_per_unique = rep(0, length(unique_graph_structures))


    for(i1 in 1:length(unique_graph_structures)){

      sum_portion_1 = sum(sid_vec[which(t(possible_pairs)[,1]==i1)]*tab[t(possible_pairs)[which(t(possible_pairs)[,1]==i1),2]])

      sum_portion_2 = sum(sid_vec[which(t(possible_pairs)[,2]==i1)]*tab[t(possible_pairs)[which(t(possible_pairs)[,2]==i1),1]])

      sid_per_unique[i1] = sum_portion_1 + sum_portion_2
    }

    best_adjacency_mat = matrix(unique_graph_structures_mat[which.min(sid_per_unique),],num_covariates,num_covariates)
  } else if(dist_type == 'forb'){
    Adjacency_matrix_75 = Adjacency_matrix_list[(0.75*num_iter):num_iter,]
    each_graph_structure = rep(0,nrow(Adjacency_matrix_75))
    for(i in 1:nrow(Adjacency_matrix_75)){
      each_graph_structure[i] = paste(Adjacency_matrix_75[i,], collapse=',')
    }

    unique_graph_structures = unique(each_graph_structure)

    tab = table(factor(each_graph_structure, levels = unique_graph_structures))

    possible_pairs = combn(1:length(unique_graph_structures),2)

    unique_graph_structures_mat = matrix(0,length(unique_graph_structures),ncol(Adjacency_matrix_75))

    structure_matrix = matrix(0, num_covariates, num_covariates)

    for(i in 1:length(unique_graph_structures)){
      unique_graph_structures_mat[i,] = as.numeric(strsplit(unique_graph_structures[i][1],',')[[1]])
    }

    forb_vec = rep(0,ncol(possible_pairs))

    structure_matrix = matrix(0, num_covariates, num_covariates)

    for(i in 1:ncol(possible_pairs)){
      current_graph_pairs = t(possible_pairs)[i,]

      graph_1_vec = unique_graph_structures_mat[current_graph_pairs[1],]
      structure_matrix[] = graph_1_vec
      #graph_1 = matrix(graph_1_vec,num_covariates,num_covariates)
      graph_1 = structure_matrix

      graph_2_vec = unique_graph_structures_mat[current_graph_pairs[2],]
      structure_matrix[] = graph_2_vec
      #graph_2 = matrix(graph_2_vec,num_covariates,num_covariates)
      graph_2 = structure_matrix

      forb_vec[i] = norm(graph_1 - graph_2, type = "F")
    }

    forb_per_unique = rep(0, length(unique_graph_structures))


    for(i1 in 1:length(unique_graph_structures)){

      sum_portion_1 = sum(forb_vec[which(t(possible_pairs)[,1]==i1)]*tab[t(possible_pairs)[which(t(possible_pairs)[,1]==i1),2]])

      sum_portion_2 = sum(forb_vec[which(t(possible_pairs)[,2]==i1)]*tab[t(possible_pairs)[which(t(possible_pairs)[,2]==i1),1]])

      forb_per_unique[i1] = sum_portion_1 + sum_portion_2
    }

    best_adjacency_mat = matrix(unique_graph_structures_mat[which.min(forb_per_unique),],num_covariates,num_covariates)
  }
  return(best_adjacency_mat)
}
