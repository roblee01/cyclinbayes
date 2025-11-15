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


summary_posterior_vec = function(posterior_samples, level){


  credible_interval = quantile(posterior_samples,c((1-level)/2,1-(1-level)/2))
  posterior_mean = mean(posterior_samples)

  hpd_interval = HDInterval::hdi(posterior_samples, credMass = level)

  return(list(posterior_mean = posterior_mean, credible_interval = credible_interval, hpd_interval = hpd_interval))
}
