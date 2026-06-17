#include <RcppArmadillo.h>
#include <queue>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// Distribution arma versions

arma::vec rbeta_cpp(int n, double alpha, double beta) {
  arma::vec g1 = arma::randg<arma::vec>(n, arma::distr_param(alpha, 1.0));
  arma::vec g2 = arma::randg<arma::vec>(n, arma::distr_param(beta, 1.0));
  arma::vec result = (g1 / (g1 + g2)).eval();
  return result;
}

arma::vec rdirichlet_cpp(const arma::vec& alpha) {
  int k = alpha.n_elem;
  arma::vec y(k);
  for (int i = 0; i < k; i++) {
    y(i) = randg(arma::distr_param(alpha(i), 1.0));
  }

  return y / arma::sum(y);
}

arma::vec rinvgamma_cpp(arma::uword n, double shape, double scale) {
  arma::vec gamma_samples = arma::randg<arma::vec>(n, arma::distr_param(shape, 1.0 / scale));
  return 1.0 / gamma_samples;
}


int sample_categorical_cpp(const arma::rowvec& probs) {
  arma::rowvec norm_probs = probs / arma::accu(probs);
  arma::rowvec cumprobs = arma::cumsum(norm_probs);
  double u = arma::randu();

  for (arma::uword i = 0; i < cumprobs.n_elem; ++i) {
    if (u < cumprobs[i]) {
      return i;
    }
  }
  return probs.n_elem-1;
}

double log_dmvn(const arma::rowvec& x, const arma::colvec& mu, const arma::mat& Sigma) {
  const double log2pi = std::log(2.0 * M_PI);
  arma::vec diff = arma::conv_to<arma::vec>::from(x.t()) - mu; // make column
  arma::mat L;
  bool status = arma::chol(L, Sigma, "lower"); // Sigma = L * L.t()
  if (!status) {
    return -arma::datum::inf; // singular; return very small log-prob
  }
  arma::vec sol = arma::solve(arma::trimatl(L), diff);          // L * y = diff
  arma::vec quad = arma::solve(arma::trimatu(L.t()), sol);      // L.t() * z = y -> z = Sigma^{-1} diff
  double maha = arma::dot(sol, sol);                           // diff' * inv(Sigma) * diff = ||sol||^2
  double logDet = 2.0 * arma::accu(arma::log(arma::diagvec(L))); // log |Sigma|
  int d = x.n_elem;
  return -0.5 * (d * log2pi + logDet + maha);
}


inline double fast_dnorm_log(double x, double mean, double sd) {
  const double log_sqrt_2pi = 0.9189385332046727;
  double z = (x - mean) / sd;
  return -log_sqrt_2pi - std::log(sd) - 0.5 * z * z;
}


inline arma::vec fast_dnorm_log_vec(const arma::vec& x,
                                    double mean,
                                    double sd) {
  const double log_sqrt_2pi = 0.9189385332046727;
  arma::vec z = (x - mean) / sd;              // elementwise
  return -log_sqrt_2pi
  - std::log(sd)
    - 0.5 * arma::square(z);            // elementwise square
}

// [[Rcpp::export]]
double logSumExp(const arma::rowvec& x) {
  double max_val = x.max();
  return max_val + log(sum(exp(x - max_val)));
}

// [[Rcpp::export]]
double log_dgamma(double x, double a, double b) {
  return a * std::log(b) - std::lgamma(a) - (a + 1.0) * std::log(x) - b / x;
}

// DAG Function

//
bool is_dag(const arma::mat& adj){
  arma::uword n = adj.n_rows;

  // Use an integer vector for in-degree counts:
  arma::Col<int> in_degree(n);
  for (arma::uword col = 0; col < n; ++col) {
    // sum of column = in-degree of node col
    in_degree(col) = static_cast<int>(arma::accu(adj.col(col)));
  }

  std::queue<arma::uword> q;
  for (arma::uword i = 0; i < n; ++i) {
    if (in_degree(i) == 0) {
      q.push(i);
    }
  }

  int count = 0;
  while (!q.empty()) {
    arma::uword u = q.front();
    q.pop();
    ++count;

    for (arma::uword v = 0; v < n; ++v) {
      if (adj(u, v) != 0) {
        in_degree(v)--;
        if (in_degree(v) == 0) {
          q.push(v);
        }
      }
    }
  }

  return count == static_cast<int>(n);
}

bool is_positive_definite(const arma::mat& A) {
  if (!arma::approx_equal(A, A.t(), "absdiff", 1e-10)) {
    return false;
  }
  arma::mat L;
  return arma::chol(L, A);
}

double log_integral_result_calculator_cpp(arma::mat Adjacency_matrix_enter, arma::mat Z_matrix, arma::mat data_matrix, double current_y, arma::vec tao_input, bool need_tao, int N, int M, double gamma_1, double b_mu, double a_tao, double b_tao, double which_iter, double num_iter){
  arma::mat Z_current = Z_matrix.rows(static_cast<arma::uword>(current_y-1)*static_cast<arma::uword>(N),(static_cast<arma::uword>(current_y-1)+1)*static_cast<arma::uword>(N)-1);
  //Rcpp::Rcout << "Z_current" << Z_current << std::endl;

  arma::vec Z_vec = arma::vectorise(Z_current.t());
  std::vector<int> current_assignments_vec;
  //Rcpp::IntegerVector tmp = Rcpp::wrap(Z_vec);
  //Rcpp::Rcout << "Z_vec" << Z_vec << std::endl;
  for(arma::uword i = 0; i < Z_vec.n_elem; i++){
    if(Z_vec[i] == 1){
      int index = (i+1) % M;

      if(index == 0){
        index = M;
      }
      current_assignments_vec.push_back(index);
    }
  }

  arma::uvec indices(static_cast<arma::uword>(current_assignments_vec.size()));
  for(arma::uword i = 0; i < current_assignments_vec.size(); i++){
    indices[i] = current_assignments_vec[i] - 1;
  }

  arma::vec current_tao_list = tao_input.elem(indices);


  arma::uvec selected_cols = arma::find(Adjacency_matrix_enter.row(static_cast<arma::uword>(current_y-1)) == 1);
  arma::mat used_ys = data_matrix.cols(selected_cols);

  double ncol_used_ys = used_ys.n_cols;


  arma::mat used_ys_plus = arma::join_rows(used_ys, Z_current);

  arma::mat weighted_input = used_ys_plus.each_col() / current_tao_list;
  arma::mat V_mat = trans(weighted_input) * used_ys_plus;

  arma::vec diag_correction_1(static_cast<arma::uword>(ncol_used_ys));
  diag_correction_1.fill(1/gamma_1);

  arma::vec diag_correction_2(M);
  diag_correction_2.fill(1/b_mu);

  arma::vec diag_correction = join_cols(diag_correction_1,diag_correction_2);

  V_mat.diag() += diag_correction;

  double tao_density_portion = 0;
  if(need_tao){
    for(arma::uword k = 0; k < tao_input.n_elem; k++){
      double x = tao_input(k);
      tao_density_portion += log_dgamma(x, a_tao, b_tao);
    }
  }
  double first_part = -(N + M + ncol_used_ys) * std::log(2 *arma::datum::pi) + arma::sum(arma::log(1/arma::sqrt(current_tao_list))) + arma::sum(arma::log(arma::sqrt(diag_correction)));
  double first_part_1 = (1 - std::exp(-20 * which_iter / num_iter)) * (N / 2.0 + M + selected_cols.n_elem) * std::log(2 * M_PI);
  //Rcpp::Rcout << "first_part" << first_part << std::endl;
  //Rcpp::Rcout << "first_part_1 " << first_part_1 << std::endl;

  double logdet;
  double sign;

  arma::log_det(logdet, sign, V_mat);
  double second_part = -0.5 * logdet - 0.5 * sum((1 / current_tao_list) % square(data_matrix.col(static_cast<arma::uword>(current_y-1))));
  //Rcpp::Rcout << "second part" << second_part << std::endl;

  arma::mat V_inv;

  //arma::mat V_inv = arma::inv_sympd(V_mat);

  if(is_positive_definite(V_mat)){
    V_inv = arma::inv_sympd(V_mat);
  } else{
    V_inv = arma::pinv(V_mat);
  }



  arma::vec weight_vec = data_matrix.col(current_y-1) % (1 / current_tao_list);
  arma::rowvec third_part_1 = trans(weight_vec) * used_ys_plus;

  double third_part = (1 + std::exp(-20 * which_iter / num_iter)) * as_scalar(third_part_1 * V_inv * third_part_1.t());
  //Rcpp::Rcout << "third part" << third_part << std::endl;
  if(need_tao){
    //Rcpp::Rcout << "tao density portion" << tao_density_portion << std::endl;
    return first_part + first_part_1 + second_part + 0.5 * third_part + tao_density_portion;
  } else{
    return first_part + first_part_1 + second_part + 0.5 * third_part;
  }
}

arma::vec Metropolis_hastings_portions_original_cpp(arma::mat data_matrix, arma::mat Adjacency_matrix_enter, arma::mat Causal_effect_matrix_enter, arma::mat Z_matrix_enter, arma::mat mu_mat, arma::mat tao_mat, double N, double M, double gamma_1, double gamma_result){
  arma::vec final_metropolis_vec = arma::zeros<arma::vec>(3);

  arma::mat I = arma::eye<arma::mat>(Causal_effect_matrix_enter.n_rows, Causal_effect_matrix_enter.n_cols);
  arma::mat IminusB = I - Causal_effect_matrix_enter;
  arma::mat IminusB_inv = arma::inv(IminusB);


  double num_features = Causal_effect_matrix_enter.n_cols;

  arma::vec result_vec = arma::zeros<arma::vec>(N);

  arma::vec mu_mat_vec = arma::vectorise(mu_mat.t());
  arma::vec tao_mat_vec = arma::vectorise(tao_mat.t());

  for(arma::uword i = 0; i < N; i++){
    arma::uvec idx = i + N * arma::regspace<arma::uvec>(0, num_features - 1);

    arma::mat Z_assignments_sample = Z_matrix_enter.rows(idx);

    arma::vec Z_vec = arma::vectorise(Z_assignments_sample.t());
    arma::uvec which1 = arma::find(Z_vec == 1);

    arma::colvec current_mu = IminusB_inv * (mu_mat_vec.elem(which1));
    arma::vec current_tao = tao_mat_vec.elem(which1);


    arma::mat temp = IminusB_inv;

    temp.each_row() %= current_tao.t();

    arma::mat variance_mat = temp * IminusB_inv.t();

    result_vec[i] = log_dmvn(data_matrix.row(i),current_mu,variance_mat);
  }

  double first_part = arma::accu(result_vec);
  final_metropolis_vec[0] = first_part;

  arma::vec Causal_effect_matrix_vec = arma::vectorise(Causal_effect_matrix_enter.t());
  arma::vec Causal_effect_matrix_no_zero_vec = Causal_effect_matrix_vec.elem(arma::find(Causal_effect_matrix_vec != 0));

  double second_part = arma::sum(fast_dnorm_log_vec(Causal_effect_matrix_no_zero_vec, 0, std::sqrt(gamma_1)));
  final_metropolis_vec[1] = second_part;

  arma::vec adj_all = arma::vectorise(Adjacency_matrix_enter);  // contains 0s and 1s

  // 2) count the 1's
  arma::uword count1 = arma::accu(adj_all);

  // 3) total entries (0's + 1's)
  arma::uword total = adj_all.n_elem;
  arma::uword count0 = total - count1;

  // 4) log‐pmf sum for Bernoulli(gamma_result):
  //    count1 * log(p) + count0 * log(1-p)
  double lp = std::log(gamma_result);
  double lq = std::log(1.0 - gamma_result);

  double third_part = count1 * lp + count0 * lq;
  final_metropolis_vec[2] = third_part;

  return final_metropolis_vec;
}

arma::vec Metropolis_hastings_portions_cpp(arma::mat data_matrix, arma::mat Adjacency_matrix_enter, arma::mat Causal_effect_matrix_enter, arma::mat Z_matrix_enter, arma::mat mu_mat, arma::mat tao_mat, double N, double M, double gamma_1, double gamma_result){
  //const arma::uword N = data_matrix.n_rows;
  const arma::uword p = Causal_effect_matrix_enter.n_cols;

  arma::vec out(3, arma::fill::zeros);

  arma::vec final_metropolis_vec = arma::zeros<arma::vec>(3);

  arma::mat IminusB = arma::eye<arma::mat>(p, p) - Causal_effect_matrix_enter;
  double log_det_IminusB;
  double sign;
  arma::log_det(log_det_IminusB, sign, IminusB);

  const double const_term = -0.5 * p * std::log(2.0 * M);

  arma::vec mu_vec  = arma::vectorise(mu_mat.t());   // length p * M
  arma::vec tao_vec = arma::vectorise(tao_mat.t());

  arma::vec result_vec(N);

  arma::uvec base_idx = N * arma::regspace<arma::uvec>(0, p - 1);

  for(arma::uword i = 0; i < N; i++){
    arma::uvec idx = i + base_idx;

    arma::mat Z_assignments_sample = Z_matrix_enter.rows(idx);
    arma::vec Z_vec = arma::vectorise(Z_assignments_sample.t());
    arma::uvec which1 = arma::find(Z_vec == 1);


    arma::vec mu_i  = mu_vec.elem(which1);   // length p
    arma::vec tao_i = tao_vec.elem(which1);  // length p

    arma::colvec y_i = data_matrix.row(i).t();
    arma::colvec eps_i = IminusB * y_i;

    arma::vec diff = eps_i - mu_i;
    double quad     = arma::sum( (diff % diff) / tao_i );
    double log_tau  = arma::sum( arma::log(tao_i) );

    result_vec[i] = const_term
    - 0.5 * (log_tau + quad)
      + log_det_IminusB;
  }

  double first_part = arma::accu(result_vec);
  out[0] = first_part;


  arma::uvec nz_idx = arma::find(Causal_effect_matrix_enter != 0.0);
  arma::vec  B_nz   = Causal_effect_matrix_enter.elem(nz_idx);
  out[1] = arma::sum(fast_dnorm_log_vec(B_nz, 0.0, std::sqrt(gamma_1)));

  arma::vec adj_all = arma::vectorise(Adjacency_matrix_enter);  // contains 0s and 1s
  arma::uword count1 = arma::accu(adj_all);
  arma::uword total = adj_all.n_elem;
  arma::uword count0 = total - count1;

  // 4) log‐pmf sum for Bernoulli(gamma_result):
  //    count1 * log(p) + count0 * log(1-p)
  double lp = std::log(gamma_result);
  double lq = std::log(1.0 - gamma_result);

  double third_part = count1 * lp + count0 * lq;
  out[2] = third_part;

  return out;
}

// [[Rcpp::export]]
List BayesSCLingam_cpp(arma::mat data_matrix, double a_mu, double b_mu,
                       double a_gamma, double b_gamma, double a_tao, double b_tao,
                       double a_og_tao, double b_og_tao, double a_gamma_1, double b_gamma_1,
                       double alpha, double M, double num_iter){

  const arma::uword p  = static_cast<arma::uword>(data_matrix.n_cols);
  const arma::uword uN = static_cast<arma::uword>(data_matrix.n_rows);
  const arma::uword uM = static_cast<arma::uword>(M);
  const double N       = data_matrix.n_rows;

  // --- Storage ---
  arma::vec log_likelihood_list(num_iter, arma::fill::zeros);
  arma::vec gamma_1_list(num_iter,        arma::fill::zeros);
  arma::vec gamma_list(num_iter,          arma::fill::zeros);
  arma::mat Adjacency_matrix_list(num_iter,     p * p,  arma::fill::zeros);
  arma::mat Causal_effect_matrix_list(num_iter, p * p,  arma::fill::zeros);
  arma::mat mu_matrix_list(num_iter,  uM * p, arma::fill::zeros);
  arma::mat tao_matrix_list(num_iter, uM * p, arma::fill::zeros);
  arma::mat pi_matrix_list(num_iter,  uM * p, arma::fill::zeros);

  // --- Initialization ---
  arma::mat Adjacency_matrix(p, p, arma::fill::zeros);
  arma::mat Causal_effect_matrix(p, p, arma::fill::zeros);

  double gamma_1      = rinvgamma_cpp(1, a_gamma, b_gamma)(0);
  double gamma_result = rbeta_cpp(1, a_gamma, b_gamma)(0);

  arma::vec rand_norm_vals_1 = arma::randn<arma::vec>(p * uM);
  arma::vec rand_norm_vals   = a_mu + b_mu * rand_norm_vals_1;
  arma::mat mu_mat = arma::reshape(rand_norm_vals, uM, p).t();

  arma::vec rand_inv_gamma_vals = rinvgamma_cpp(p * uM, a_tao, b_tao);
  arma::mat tao_mat = arma::reshape(rand_inv_gamma_vals, uM, p).t();

  arma::mat pi_mat(p, uM);
  arma::vec alpha_vec(uM, arma::fill::value(alpha));
  for(int i4 = 0; i4 < p; i4++)
    pi_mat.row(i4) = rdirichlet_cpp(alpha_vec).t();

  arma::mat Z_matrix = arma::zeros<arma::mat>(p * uN, uM);
  for(int i5 = 0; i5 < p; i5++){
    arma::rowvec probs = pi_mat.row(i5);
    for(int j5 = 0; j5 < uN; j5++){
      arma::uword ij = i5 * uN + j5;
      Z_matrix(ij, sample_categorical_cpp(probs)) = 1;
    }
  }

  arma::mat epsilon_mat = ((arma::eye(p, p) - Causal_effect_matrix)
                             * data_matrix.t()).t();

  arma::vec numerator_result(uM,      arma::fill::zeros);
  arma::vec numerator_portion_1(uN,   arma::fill::zeros);
  arma::vec denominator_portion_1(uN, arma::fill::zeros);
  arma::mat practice_prob_mat(uN, uM);
  arma::uvec all_indices = arma::regspace<arma::uvec>(0, p - 1);

  // --- Precompute exclude_j vectors once before the loop ---
  std::vector<arma::uvec> remaining_index_by_node(p);
  std::vector<arma::uvec> exclude_j_by_node(p);
  for(arma::uword k = 0; k < p; k++){
    remaining_index_by_node[k] = all_indices.elem(arma::find(all_indices != k));
    exclude_j_by_node[k]       = all_indices.elem(arma::find(all_indices != k));
  }

  // --- MH accept lambda ---
  auto mh_accept = [](double log_r) -> bool {
    if(log_r >= 0) return true;
    return sample_categorical_cpp(
      arma::rowvec{std::exp(log_r), 1.0 - std::exp(log_r)}) == 0;
  };

  // --- Main MCMC loop ---
  for(int i = 1; i <= num_iter; i++){

    // 1. gamma update
    double total_entries = std::pow(Adjacency_matrix.n_rows, 2);
    double a = a_gamma + arma::accu(Adjacency_matrix);
    double b = b_gamma + total_entries - arma::accu(Adjacency_matrix)
      - Adjacency_matrix.n_rows;
    gamma_result = rbeta_cpp(1, a, b)(0);
    gamma_list(i - 1) = gamma_result;

    // 2. Adjacency matrix + tao update
    for(int i1 = 0; i1 < p; i1++){
      const arma::uvec& remaining_index = remaining_index_by_node[i1];

      arma::vec og_tao = rinvgamma_cpp(uM, a_og_tao, b_og_tao);

      double log_numerator_portion_1 = log_integral_result_calculator_cpp(
        Adjacency_matrix, Z_matrix, data_matrix, i1+1, og_tao,
        true, N, M, gamma_1, b_mu, a_tao, b_tao, i, num_iter);
      double log_denominator_portion_1 = log_integral_result_calculator_cpp(
        Adjacency_matrix, Z_matrix, data_matrix, i1+1,
        tao_mat.row(i1).t(), true, N, M, gamma_1, b_mu, a_tao, b_tao,
        i, num_iter);

      if(mh_accept(log_numerator_portion_1 - log_denominator_portion_1))
        tao_mat.row(i1) = og_tao.t();

      for(int j1 = 0; j1 < remaining_index.n_elem; j1++){
        Adjacency_matrix(i1, remaining_index(j1)) = 1;
        if(!is_dag(Adjacency_matrix)){
          Adjacency_matrix(i1, remaining_index(j1)) = 0;
        } else {
          double log_numerator_portion_1 = log_integral_result_calculator_cpp(
            Adjacency_matrix, Z_matrix, data_matrix, i1+1,
            tao_mat.row(i1).t(), false, N, M, gamma_1, b_mu, a_tao,
            b_tao, i, num_iter);

          Adjacency_matrix(i1, remaining_index(j1)) = 0;
          double log_denominator_portion_1 = log_integral_result_calculator_cpp(
            Adjacency_matrix, Z_matrix, data_matrix, i1+1,
            tao_mat.row(i1).t(), false, N, M, gamma_1, b_mu, a_tao,
            b_tao, i, num_iter);

          double r_log_tmp = log_numerator_portion_1 - log_denominator_portion_1
          - (1 + std::exp(-35.0*i/num_iter)) * std::log(1 - gamma_result)
            + (1 + std::exp(-35.0*i/num_iter)) * std::log(gamma_result);

            if(mh_accept(r_log_tmp))
              Adjacency_matrix(i1, remaining_index(j1)) = 1;

            if(Adjacency_matrix(i1, remaining_index(j1)) == 1){
              Adjacency_matrix(i1, remaining_index(j1)) = 0;
              Adjacency_matrix(remaining_index(j1), i1) = 1;

              if(!is_dag(Adjacency_matrix)){
                Adjacency_matrix(i1, remaining_index(j1)) = 1;
                Adjacency_matrix(remaining_index(j1), i1) = 0;
              } else {
                double log_numerator_portion_1_1 =
                  log_integral_result_calculator_cpp(
                    Adjacency_matrix, Z_matrix, data_matrix, i1+1,
                    tao_mat.row(i1).t(), false, N, M, gamma_1, b_mu,
                    a_tao, b_tao, i, num_iter);
                double log_numerator_portion_1_2 =
                  log_integral_result_calculator_cpp(
                    Adjacency_matrix, Z_matrix, data_matrix,
                    remaining_index(j1)+1,
                    tao_mat.row(remaining_index(j1)).t(), false, N, M,
                    gamma_1, b_mu, a_tao, b_tao, i, num_iter);

                Adjacency_matrix(i1, remaining_index(j1)) = 1;
                Adjacency_matrix(remaining_index(j1), i1) = 0;

                double log_numerator_portion_2_1 = log_numerator_portion_1;
                double log_numerator_portion_2_2 =
                  log_integral_result_calculator_cpp(
                    Adjacency_matrix, Z_matrix, data_matrix,
                    remaining_index(j1)+1,
                    tao_mat.row(remaining_index(j1)).t(), false, N, M,
                    gamma_1, b_mu, a_tao, b_tao, i, num_iter);

                double r_log_tmp_1 =
                  (log_numerator_portion_1_1 + log_numerator_portion_1_2)
                  - (log_numerator_portion_2_1 + log_numerator_portion_2_2);

                if(mh_accept(r_log_tmp_1)){
                  Adjacency_matrix(i1, remaining_index(j1)) = 0;
                  Adjacency_matrix(remaining_index(j1), i1) = 1;
                }
              }
            }
        }
      }
    }
    Adjacency_matrix_list.row(i - 1) = arma::vectorise(Adjacency_matrix).t();

    // 3. mu update
    for(int i2 = 0; i2 < p; i2++){
      arma::mat first_part = Z_matrix.rows(
        i2 * uN, (i2 + 1) * uN - 1);
      first_part.each_row() /= tao_mat.row(i2);

      for(int j2 = 0; j2 < uM; j2++)
        numerator_result(j2) = (a_mu / b_mu)
        + arma::dot(first_part.col(j2), epsilon_mat.col(i2));

      arma::vec denominator_result = (1.0 / b_mu)
        + arma::sum(first_part, 0).t();
      arma::vec mean_vec = numerator_result / denominator_result;
      arma::vec sd_vec   = arma::sqrt(1.0 / denominator_result);
      arma::rowvec mu_row = mean_vec.t()
        + arma::randn<arma::rowvec>(uM) % sd_vec.t();
      arma::rowvec mu_row_ordered = arma::sort(mu_row, "ascend");
      mu_mat.row(i2) = mu_row_ordered;
    }
    mu_matrix_list.row(i - 1) = arma::vectorise(mu_mat).t();

    // 4. Causal effect update
    for(int i3 = 0; i3 < p; i3++){
      arma::mat category_mat = Z_matrix.rows(
        i3 * uN, (i3 + 1) * uN - 1);
      category_mat.each_row() /= tao_mat.row(i3);

      // Hoisted outside j3/z3 loops since it only depends on i3,
      // not j3 or z3 — Causal_effect_matrix(i3,*) doesn't change
      // within this i3 iteration until after all j3 are processed
      arma::rowvec mu_row_i3 = mu_mat.row(i3);

      for(int j3 = 0; j3 < p; j3++){
        if(Adjacency_matrix(i3, j3) == 0){
          Causal_effect_matrix(i3, j3) = 0;
        } else {
          const arma::uvec& exclude_j = exclude_j_by_node[j3];

          // Fix 1: hoisted out of z3 loop — these only depend on i3 and j3,
          // not z3, so they were being recomputed N times per edge for no reason
          arma::rowvec ce_row = arma::rowvec(Causal_effect_matrix.row(i3));
          arma::rowvec ce_portion =
            arma::conv_to<arma::rowvec>::from(ce_row.elem(exclude_j));

          for(int z3 = 0; z3 < uN; z3++){
            arma::rowvec Y = data_matrix.row(z3);

            arma::rowvec y_sub =
              arma::conv_to<arma::rowvec>::from(Y.elem(exclude_j));

            double dot_portion     = arma::dot(ce_portion, y_sub);
            arma::rowvec adjusted_mu = Y(i3) - (mu_row_i3 + dot_portion);
            arma::rowvec category_row = Y(j3) * category_mat.row(z3);

            numerator_portion_1(z3)   = arma::dot(category_row, adjusted_mu);
            denominator_portion_1(z3) = arma::accu(
              std::pow(Y(j3), 2) * category_mat.row(z3));
          }

          double numerator_final   = arma::sum(numerator_portion_1);
          double denominator_final = 1.0 / gamma_1
          + arma::sum(denominator_portion_1);
          double mean_result     = numerator_final / denominator_final;
          double variance_result = 1.0 / denominator_final;
          Causal_effect_matrix(i3, j3) = mean_result
          + std::sqrt(variance_result) * arma::randn();
        }
      }
    }
    Causal_effect_matrix_list.row(i - 1) =
      arma::vectorise(Causal_effect_matrix).t();

    // 5. Epsilon update
    epsilon_mat = ((arma::eye(p, p) - Causal_effect_matrix)
                     * data_matrix.t()).t();

    // 6. tao update
    for(int i4 = 0; i4 < p; i4++){
      arma::mat Z_portion = Z_matrix.rows(
        i4 * uN, (i4 + 1) * uN - 1);
      arma::rowvec a_tao_vec = a_tao + arma::sum(Z_portion, 0) / 2.0;
      for(int j4 = 0; j4 < uM; j4++){
        arma::vec eps_col      = epsilon_mat.col(i4);
        double mu_ij           = mu_mat(i4, j4);
        arma::vec squared_diff = arma::square(eps_col - mu_ij);
        arma::vec b_portion    = 0.5 * Z_portion.col(j4) % squared_diff;
        double b_val           = b_tao + arma::sum(b_portion);
        tao_mat(i4, j4)        = rinvgamma_cpp(1, a_tao_vec(j4), b_val)(0);
      }
    }
    tao_matrix_list.row(i - 1) = arma::vectorise(tao_mat).t();

    // 7. gamma_1 update
    double a_1 = a_gamma_1 + arma::accu(Adjacency_matrix) / 2.0;
    double b_1 = b_gamma_1 + arma::accu(
      Adjacency_matrix % (Causal_effect_matrix % Causal_effect_matrix))
      / 2.0;
    gamma_1 = rinvgamma_cpp(1, a_1, b_1)(0);
    gamma_1_list(i - 1) = gamma_1;

    // 8. Z matrix update
    for(int i5 = 0; i5 < p; i5++){
      for(int z5 = 0; z5 < uN; z5++){
        int iz = i5 * uN + z5;
        Z_matrix.row(iz).zeros();

        for(int j5 = 0; j5 < uM; j5++){
          double mu_1  = mu_mat(i5, j5);
          double tao_1 = tao_mat(i5, j5);
          double log_probability = fast_dnorm_log(
            epsilon_mat(z5, i5), mu_1, std::sqrt(tao_1));
          practice_prob_mat(z5, j5) = log_probability;
        }

        arma::rowvec log_probs = arma::log(pi_mat.row(i5))
          + practice_prob_mat.row(z5);
        arma::rowvec probs = arma::exp(log_probs - logSumExp(log_probs));

        double u       = arma::randu();
        double cum_prob = 0.0;
        int sampled    = uM - 1;
        for(int m = 0; m < uM; m++){
          cum_prob += probs[m];
          if(u < cum_prob){ sampled = m; break; }
        }
        Z_matrix(iz, sampled) = 1;
      }
    }

    // 9. pi update
    for(int i6 = 0; i6 < p; i6++){
      arma::mat Z_portion = Z_matrix.rows(
        i6 * uN, (i6 + 1) * uN - 1);
      arma::rowvec Z_colsum = arma::sum(Z_portion, 0);
      arma::vec alpha_vec_local = Z_colsum.t() + alpha;
      pi_mat.row(i6) = rdirichlet_cpp(alpha_vec_local).t();
    }
    pi_matrix_list.row(i - 1) = arma::vectorise(pi_mat).t();

    // 10. Log likelihood
    arma::vec log_parts = Metropolis_hastings_portions_cpp(
      data_matrix, Adjacency_matrix, Causal_effect_matrix,
      Z_matrix, mu_mat, tao_mat, N, M, gamma_1, gamma_result);
    log_likelihood_list(i - 1) = log_parts[0];
  }

  return List::create(
    Named("Adjacency_matrix_list")     = Adjacency_matrix_list,
    Named("Causal_effect_matrix_list") = Causal_effect_matrix_list,
    Named("gamma_list")                = gamma_list,
    Named("gamma_1_list")              = gamma_1_list,
    Named("mu_matrix_list")            = mu_matrix_list,
    Named("tao_matrix_list")           = tao_matrix_list,
    Named("pi_matrix_list")            = pi_matrix_list,
    Named("log_likelihood_list")       = log_likelihood_list
  );
}

// [[Rcpp::export]]
List BCD_cpp(arma::mat data_matrix, double a_mu, double b_mu,
             double a_gamma, double b_gamma, double a_tao, double b_tao,
             double a_gamma_1, double b_gamma_1, double alpha, double M, double num_iter){

  const arma::uword p  = static_cast<arma::uword>(data_matrix.n_cols);
  const arma::uword uN = static_cast<arma::uword>(data_matrix.n_rows);
  const arma::uword uM = static_cast<arma::uword>(M);
  const double N       = data_matrix.n_rows;

  // --- Storage ---
  arma::vec log_likelihood_list(num_iter, arma::fill::zeros);
  arma::vec gamma_1_list(num_iter,        arma::fill::zeros);
  arma::vec gamma_list(num_iter,          arma::fill::zeros);
  arma::mat Adjacency_matrix_list(num_iter,     p * p,  arma::fill::zeros);
  arma::mat Causal_effect_matrix_list(num_iter, p * p,  arma::fill::zeros);
  arma::mat mu_matrix_list(num_iter,  uM * p, arma::fill::zeros);
  arma::mat tao_matrix_list(num_iter, uM * p, arma::fill::zeros);
  arma::mat pi_matrix_list(num_iter,  uM * p, arma::fill::zeros);

  // --- Initialization ---
  arma::mat Adjacency_matrix(p, p, arma::fill::zeros);
  arma::mat Causal_effect_matrix(p, p, arma::fill::zeros);

  double gamma_1      = rinvgamma_cpp(1, a_gamma, b_gamma)(0);
  double gamma_result = rbeta_cpp(1, a_gamma, b_gamma)(0);

  arma::vec rand_norm_vals_1 = arma::randn<arma::vec>(p * uM);
  arma::vec rand_norm_vals   = a_mu + b_mu * rand_norm_vals_1;
  arma::mat mu_mat = arma::reshape(rand_norm_vals, uM, p).t();

  arma::vec rand_inv_gamma_vals = rinvgamma_cpp(p * uM, a_tao, b_tao);
  arma::mat tao_mat = arma::reshape(rand_inv_gamma_vals, uM, p).t();

  arma::mat pi_mat(p, uM);
  arma::vec alpha_vec(uM, arma::fill::value(alpha));
  for(int i4 = 0; i4 < p; i4++)
    pi_mat.row(i4) = rdirichlet_cpp(alpha_vec).t();

  arma::mat Z_matrix = arma::zeros<arma::mat>(p * uN, uM);
  for(int i5 = 0; i5 < p; i5++){
    arma::rowvec probs = pi_mat.row(i5);
    for(int j5 = 0; j5 < uN; j5++){
      arma::uword ij = i5 * uN + j5;
      Z_matrix(ij, sample_categorical_cpp(probs)) = 1;
    }
  }

  arma::mat epsilon_mat = ((arma::eye(p, p) - Causal_effect_matrix)
                             * data_matrix.t()).t();

  arma::vec numerator_result(uM,      arma::fill::zeros);
  arma::vec numerator_portion_1(uN,   arma::fill::zeros);
  arma::vec denominator_portion_1(uN, arma::fill::zeros);
  arma::mat practice_prob_mat(uN, uM);
  arma::uvec all_indices = arma::regspace<arma::uvec>(0, p - 1);

  // --- Precompute remaining_indices once before the loop ---
  std::vector<arma::uvec> remaining_indices(p);
  for(arma::uword k = 0; k < p; k++)
    remaining_indices[k] = all_indices.elem(arma::find(all_indices != k));

  // --- Stability check helper ---
  // Returns true if the matrix is stable (max eigenvalue modulus < 1)
  // Uses a cheap row-sum pre-filter before falling through to full eig_gen.
  // The Gershgorin circle theorem guarantees that if the maximum absolute
  // row sum is < 1 then all eigenvalues lie strictly inside the unit disk,
  // so eig_gen can be skipped entirely in that case. Only when the row sum
  // exceeds 1 do we need the full check.
  auto is_stable = [](const arma::mat& C) -> bool {
    // Gershgorin pre-filter: cheap O(p^2) check
    arma::vec row_sums = arma::sum(arma::abs(C), 1);
    if(row_sums.max() < 1.0) return true;
    // Fall through to full O(p^3) eigenvalue decomposition
    arma::cx_vec eigenvalues  = arma::eig_gen(C);
    arma::vec mod_eigenvalues = arma::abs(eigenvalues);
    return mod_eigenvalues.max() < 1.0;
  };

  // --- Initialize log_current_vec before loop ---
  arma::vec log_current_vec = Metropolis_hastings_portions_cpp(
    data_matrix, Adjacency_matrix, Causal_effect_matrix,
    Z_matrix, mu_mat, tao_mat, N, M, gamma_1, gamma_result);
  double log_current = arma::sum(log_current_vec);

  // --- Main MCMC loop ---
  for(int i = 1; i <= num_iter; i++){

    // 1. gamma update
    double total_entries = std::pow(Adjacency_matrix.n_rows, 2);
    double a = a_gamma + arma::accu(Adjacency_matrix);
    double b = b_gamma + total_entries - arma::accu(Adjacency_matrix)
      - Adjacency_matrix.n_rows;
    // NOTE: original shadows outer gamma_result here, preserved as-is
    gamma_result = rbeta_cpp(1, a, b)(0);
    gamma_list(i - 1) = gamma_result;

    // gamma_result just changed so one recompute is unavoidable
    log_current_vec = Metropolis_hastings_portions_cpp(
      data_matrix, Adjacency_matrix, Causal_effect_matrix,
      Z_matrix, mu_mat, tao_mat, N, M, gamma_1, gamma_result);
    log_current = arma::sum(log_current_vec);

    // Fix 1: record likelihood here reusing the unavoidable gamma recompute
    log_likelihood_list(i - 1) = log_current_vec(0);

    // 2. Structure + causal effect MH update (add/remove edge)
    double proposal_mean_add = 0;
    double burn_in           = 1000;
    double proposal_sd_add;

    if(i < burn_in){
      proposal_sd_add = 0.15;
    } else {
      proposal_sd_add = std::sqrt(gamma_1);
    }

    double log_q_fwd = 0;
    double log_q_rev = 0;

    for(arma::uword i6 = 0; i6 < p; i6++){
      // Fix 2: use precomputed remaining_indices
      arma::uvec entries_order = arma::shuffle(remaining_indices[i6]);

      for(arma::uword j6 = 0; j6 < entries_order.n_elem; j6++){
        log_q_fwd = 0;
        log_q_rev = 0;

        arma::mat Causal_prop    = Causal_effect_matrix;
        arma::mat Adjacency_prop = Adjacency_matrix;
        double current_column    = entries_order[j6];

        if(Adjacency_prop(i6, current_column) == 0){
          Adjacency_prop(i6, current_column) = 1;
          double proposal_coef = proposal_mean_add
          + proposal_sd_add * arma::randn();
          Causal_prop(i6, current_column) = proposal_coef;
          log_q_fwd = fast_dnorm_log(
            proposal_coef, proposal_mean_add, proposal_sd_add);
          log_q_rev = 0;
        } else {
          Adjacency_prop(i6, current_column) = 0;
          double b_old = Causal_effect_matrix(i6, current_column);
          Causal_prop(i6, current_column) = 0;
          log_q_fwd = 0;
          log_q_rev = fast_dnorm_log(
            b_old, proposal_mean_add, proposal_sd_add);
        }

        // Fix 4: use Gershgorin pre-filter before full eig_gen
        if(is_stable(Causal_prop)){
          arma::vec log_prop_vec = Metropolis_hastings_portions_cpp(
            data_matrix, Adjacency_prop, Causal_prop,
            Z_matrix, mu_mat, tao_mat, N, M, gamma_1, gamma_result);

          double log_proposal = arma::sum(log_prop_vec);
          double log_r        = log_proposal - log_current
          + log_q_rev - log_q_fwd;
          double prob         = std::exp(log_r);
          arma::rowvec bool_vec = {1, prob};

          double sampled_val = arma::randu();
          if(sampled_val < bool_vec.min()){
            Causal_effect_matrix = Causal_prop;
            Adjacency_matrix     = Adjacency_prop;
            log_current          = log_proposal;
            log_current_vec      = log_prop_vec;
          }
        }
      }
    }

    // 3. Causal effect weight MH update
    double burn_in_sd     = 15000;
    arma::vec v           = {static_cast<double>(i), burn_in_sd};
    double weight_prop_sd = 0.03 + 0.07 * (arma::min(v) / burn_in_sd);

    log_current = log_current_vec(0) + log_current_vec(1);

    for(arma::uword i6 = 0; i6 < p; i6++){
      // Fix 3: skip current_row copy
      arma::uvec present_rows = arma::find(Adjacency_matrix.row(i6) != 0);

      if(present_rows.n_elem != 0){
        for(arma::uword j6 = 0; j6 < present_rows.n_elem; j6++){
          arma::mat Causal_prop = Causal_effect_matrix;
          Causal_prop(i6, present_rows(j6)) =
            Causal_effect_matrix(i6, present_rows(j6))
            + weight_prop_sd * arma::randn();

          // Fix 4: use Gershgorin pre-filter before full eig_gen
          if(is_stable(Causal_prop)){
            arma::vec log_prop_vec = Metropolis_hastings_portions_cpp(
              data_matrix, Adjacency_matrix, Causal_prop,
              Z_matrix, mu_mat, tao_mat, N, M, gamma_1, gamma_result);

            double log_proposal = log_prop_vec(0) + log_prop_vec(1);
            double log_r        = log_proposal - log_current;
            double prob         = std::exp(log_r);
            arma::rowvec bool_vec = {1, prob};

            double sampled_val = arma::randu();
            if(sampled_val < bool_vec.min()){
              Causal_effect_matrix(i6, present_rows(j6)) =
                Causal_prop(i6, present_rows(j6));
              log_current_vec = log_prop_vec;
              log_current     = log_proposal;
            }
          }
        }
      }
    }

    Adjacency_matrix_list.row(i - 1) =
      arma::vectorise(Adjacency_matrix).t();
    Causal_effect_matrix_list.row(i - 1) =
      arma::vectorise(Causal_effect_matrix).t();

    // 4. Epsilon update
    epsilon_mat = ((arma::eye(p, p) - Causal_effect_matrix)
                     * data_matrix.t()).t();

    // 5. mu update
    for(int i2 = 0; i2 < p; i2++){
      arma::mat first_part = Z_matrix.rows(
        i2 * uN, (i2 + 1) * uN - 1);
      first_part.each_row() /= tao_mat.row(i2);

      for(int j2 = 0; j2 < uM; j2++)
        numerator_result(j2) = (a_mu / b_mu)
        + arma::dot(first_part.col(j2), epsilon_mat.col(i2));

      arma::vec denominator_result = (1.0 / b_mu)
        + arma::sum(first_part, 0).t();
      arma::vec mean_vec = numerator_result / denominator_result;
      arma::vec sd_vec   = arma::sqrt(1.0 / denominator_result);
      arma::rowvec mu_row = mean_vec.t()
        + arma::randn<arma::rowvec>(uM) % sd_vec.t();
      // NOTE: original does not sort mu_row in BCD, preserved as-is
      mu_mat.row(i2) = mu_row;
    }
    mu_matrix_list.row(i - 1) = arma::vectorise(mu_mat).t();

    // 6. tao update
    for(int i4 = 0; i4 < p; i4++){
      arma::mat Z_portion = Z_matrix.rows(
        i4 * uN, (i4 + 1) * uN - 1);
      arma::rowvec a_tao_vec = a_tao + arma::sum(Z_portion, 0) / 2.0;
      for(int j4 = 0; j4 < uM; j4++){
        arma::vec eps_col      = epsilon_mat.col(i4);
        double mu_ij           = mu_mat(i4, j4);
        arma::vec squared_diff = arma::square(eps_col - mu_ij);
        arma::vec b_portion    = 0.5 * Z_portion.col(j4) % squared_diff;
        double b_val           = b_tao + arma::sum(b_portion);
        tao_mat(i4, j4)        = rinvgamma_cpp(1, a_tao_vec(j4), b_val)(0);
      }
    }
    tao_matrix_list.row(i - 1) = arma::vectorise(tao_mat).t();

    // 7. gamma_1 update
    double a_1 = a_gamma_1 + arma::accu(Adjacency_matrix) / 2.0;
    double b_1 = b_gamma_1 + arma::accu(
      Adjacency_matrix % (Causal_effect_matrix % Causal_effect_matrix))
      / 2.0;
    gamma_1 = rinvgamma_cpp(1, a_1, b_1)(0);
    gamma_1_list(i - 1) = gamma_1;

    // 8. Z matrix update
    for(int i5 = 0; i5 < p; i5++){
      for(int z5 = 0; z5 < uN; z5++){
        int iz = i5 * uN + z5;
        Z_matrix.row(iz).zeros();

        for(int j5 = 0; j5 < uM; j5++){
          double mu_1  = mu_mat(i5, j5);
          double tao_1 = tao_mat(i5, j5);
          practice_prob_mat(z5, j5) = fast_dnorm_log(
            epsilon_mat(z5, i5), mu_1, std::sqrt(tao_1));
        }

        arma::rowvec log_probs = arma::log(pi_mat.row(i5))
          + practice_prob_mat.row(z5);
        arma::rowvec probs = arma::exp(log_probs - logSumExp(log_probs));

        double u        = arma::randu();
        double cum_prob = 0.0;
        int sampled     = uM - 1;
        for(int m = 0; m < uM; m++){
          cum_prob += probs[m];
          if(u < cum_prob){ sampled = m; break; }
        }
        Z_matrix(iz, sampled) = 1;
      }
    }

    // 9. pi update
    for(int i6 = 0; i6 < p; i6++){
      arma::mat Z_portion = Z_matrix.rows(
        i6 * uN, (i6 + 1) * uN - 1);
      arma::rowvec Z_colsum     = arma::sum(Z_portion, 0);
      arma::vec alpha_vec_local = Z_colsum.t() + alpha;
      pi_mat.row(i6) = rdirichlet_cpp(alpha_vec_local).t();
    }
    pi_matrix_list.row(i - 1) = arma::vectorise(pi_mat).t();

    // Step 10 removed — likelihood recorded at top of next iteration
  }

  // Record final iteration likelihood
  arma::vec final_log = Metropolis_hastings_portions_cpp(
    data_matrix, Adjacency_matrix, Causal_effect_matrix,
    Z_matrix, mu_mat, tao_mat, N, M, gamma_1, gamma_result);
  log_likelihood_list(num_iter - 1) = final_log(0);

  return List::create(
    Named("Adjacency_matrix_list")     = Adjacency_matrix_list,
    Named("Causal_effect_matrix_list") = Causal_effect_matrix_list,
    Named("gamma_list")                = gamma_list,
    Named("gamma_1_list")              = gamma_1_list,
    Named("mu_matrix_list")            = mu_matrix_list,
    Named("tao_matrix_list")           = tao_matrix_list,
    Named("pi_matrix_list")            = pi_matrix_list,
    Named("log_likelihood_list")       = log_likelihood_list
  );
}




