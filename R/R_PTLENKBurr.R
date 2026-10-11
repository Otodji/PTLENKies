#' Random Generation for the PTLENKBurr Distribution
#'
#' Generates random samples from the PTLENKBurr distribution using inverse
#' transform sampling (quantile method).
#'
#' @param n_samples Number of random samples to generate (must be a positive integer).
#' @param m Parameter m of the distribution (m > 0).
#' @param n Parameter n of the distribution (n > 0).
#' @param p Parameter p of the distribution (p > 0).
#' @param lambda Parameter lambda of the distribution (lambda > 0).
#' @param beta Parameter beta of the distribution (beta > 0).
#' @param c Parameter c of the Burr XII baseline component (c > 0).
#' @param k Parameter k of the Burr XII baseline component (k > 0).
#'
#' @return A numeric vector of \code{n_samples} random values generated from
#'   the PTLENKBurr distribution.
#'
#' @details
#' The parameter for the sample size is named \code{n_samples} to avoid naming
#' conflicts with the distribution parameter \code{n}. Random variates are generated
#' by applying the quantile function \code{\link{Q_PTLENKBurr}} to uniform random
#' variables: \eqn{X = Q(U)} where \eqn{U \sim \text{Uniform}(0, 1)}.
#'
#' @importFrom stats runif
#' @export
#'
#' @examples
#' set.seed(123)
#' sim_data <- R_PTLENKBurr(
#'   n_samples = 100,
#'   m = 1, n = 1, p = 1,
#'   lambda = 0.5, beta = 1.5, c = 2, k = 1
#' )
#' head(sim_data)
R_PTLENKBurr <- function(n_samples, m, n, p, lambda, beta, c, k) {

  # 1. Validate sample size parameter
  if (!is.numeric(n_samples) || length(n_samples) != 1 || n_samples <= 0) {
    stop("'n_samples' must be a single strictly positive integer.")
  }
  n_samples <- as.integer(n_samples)

  # 2. Generate standard uniform random variates U ~ Uniform(0, 1)
  u <- runif(n_samples)

  # 3. Compute quantiles via inverse CDF transformation
  x_values <- Q_PTLENKBurr(u = u, m = m, n = n, p = p, lambda = lambda, beta = beta, c = c, k = k)

  return(x_values)
}
