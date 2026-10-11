#' Cumulative Distribution Function (CDF) of the PTLENKBurr Distribution
#'
#' This function calculates the Cumulative Distribution Function (CDF) of the PTLENKBurr distribution.
#'
#' @param x Value or vector of positive values up to which to calculate the CDF.
#' @param m Parameter m of the distribution (m > 0).
#' @param n Parameter n of the distribution (n > 0).
#' @param p Parameter p of the distribution (p > 0).
#' @param lambda Parameter lambda of the distribution (lambda > 0).
#' @param beta Parameter beta of the distribution (beta > 0).
#' @param c Parameter c of the Burr XII baseline component (c > 0).
#' @param k Parameter k of the Burr XII baseline component (k > 0).
#'
#' @return A numeric vector containing the values of the CDF evaluated at \code{x}.
#'
#' @details
#' It takes parameters x, m, n, p, lambda, beta, c, and k, and returns the CDF value at x based on these parameters.
#'
#' The Cumulative Distribution Function (CDF) of the PTLENKBurr distribution is defined as:
#' \deqn{
#'   F(x; m, n, p, \lambda, \beta, c, k) = \exp\left( m n \left(1 - \frac{1}{K(x)^p}\right) \right) \left[ 2 - \exp\left( n \left(1 - \frac{1}{K(x)^p}\right) \right) \right]^m
#' }
#' where
#' \deqn{
#'   K(x) = 1 - \exp\left( -\lambda \left( \frac{1 - (1+x^c)^{-k}}{(1+x^c)^{-k}} \right)^\beta \right)
#' }
#' and \eqn{m, n, p, \lambda, \beta, c, k > 0}, \eqn{x > 0}.
#'
#' @export
#'
#' @examples
#' C_PTLENKBurr(
#'   x = c(0.5, 1.0, 1.5),
#'   m = 1, n = 1, p = 1,
#'   lambda = 0.5, beta = 1.5, c = 2, k = 1
#' )
C_PTLENKBurr <- function(x, m, n, p, lambda, beta, c, k) {
  K <- 1 - exp(-lambda * ((1 - (1 + x^c)^(-k)) / (1 + x^c)^(-k))^beta)
  cdf <- exp(m * n * (1 - 1 / (K^p))) * (2 - exp(n * (1 - 1 / (K^p))))^m
  return(cdf)
}
