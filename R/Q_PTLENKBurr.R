#' Quantile Function of the PTLENKBurr Distribution
#'
#' Calculates the quantile value (inverse CDF) of the PTLENKBurr distribution
#' for a given probability \code{u} or a vector of probabilities.
#'
#' @param u Numeric vector of probabilities (\eqn{0 \le u \le 1}).
#' @param m Parameter m of the distribution (m > 0).
#' @param n Parameter n of the distribution (n > 0).
#' @param p Parameter p of the distribution (p > 0).
#' @param lambda Parameter lambda of the distribution (lambda > 0).
#' @param beta Parameter beta of the distribution (beta > 0).
#' @param c Parameter c of the Burr XII baseline component (c > 0).
#' @param k Parameter k of the Burr XII baseline component (k > 0).
#' @param tol Tolerance level for numerical root finding via \code{\link[stats]{uniroot}}. Default is \code{1e-8}.
#'
#' @return A numeric vector of quantiles corresponding to the probabilities in \code{u}.
#'
#' @details
#' The quantile function \eqn{Q(u)} is obtained by numerically solving \eqn{F(x) - u = 0}
#' using \code{\link[stats]{uniroot}}. The probability argument is named \code{u} to avoid
#' name collision with the distribution parameter \code{p}.
#'
#' @importFrom stats uniroot
#' @export
#'
#' @examples
#' # Calculate quartiles (25%, 50%, 75%)
#' Q_PTLENKBurr(
#'   u = c(0.25, 0.5, 0.75),
#'   m = 1, n = 1, p = 1,
#'   lambda = 0.5, beta = 1.5, c = 2, k = 1
#' )
Q_PTLENKBurr <- function(u, m, n, p, lambda, beta, c, k, tol = 1e-8) {

  # Internal helper function to compute quantile for a single scalar probability
  get_single_quantile <- function(u_single) {
    # 1. Handle out-of-domain values and boundary conditions
    if (is.na(u_single)) return(NA)
    if (u_single < 0 || u_single > 1) return(NaN)
    if (u_single == 0) return(0)
    if (u_single == 1) return(Inf)

    # 2. Target objective function: F(x) - u = 0
    target_fun <- function(x) {
      C_PTLENKBurr(x = x, m = m, n = n, p = p, lambda = lambda, beta = beta, c = c, k = k) - u_single
    }

    # 3. Dynamic upper bound expansion for root finding
    upper_bound <- 10
    while (suppressWarnings(C_PTLENKBurr(upper_bound, m, n, p, lambda, beta, c, k)) < u_single && upper_bound < 1e7) {
      upper_bound <- upper_bound * 2
    }

    # 4. Numerical root-finding via uniroot
    res <- try(
      uniroot(target_fun, interval = c(1e-10, upper_bound), tol = tol)$root,
      silent = TRUE
    )

    if (inherits(res, "try-error")) {
      return(NA)
    }
    return(res)
  }

  # Vectorize evaluation across probability vector u with type safety
  quantiles <- vapply(u, get_single_quantile, FUN.VALUE = numeric(1))

  return(quantiles)
}
