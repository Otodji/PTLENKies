#' Probability Density Function (PDF) of the PTLENKBurr Distribution
#'
#' This function calculates the Probability Density Function (PDF) of the PTLENKBurr distribution.
#'
#' \deqn{
#'   f(x; m, n, p, \lambda, \beta, c, k) = 2 m n p \lambda \beta c k x^{c-1} (1+x^c)^{-(k+1)} \frac{\left(1 - (1+x^c)^{-k}\right)^{\beta-1}}{\left((1+x^c)^{-k}\right)^{\beta+1}} \frac{1 - K(x)}{K(x)^{p+1}} \exp\left( m n \left(1 - \frac{1}{K(x)^p}\right) \right) \left( 1 - \exp\left( n \left(1 - \frac{1}{K(x)^p}\right) \right) \right) \left( 2 - \exp\left( n \left(1 - \frac{1}{K(x)^p}\right) \right) \right)^{m-1}
#' }
#' where
#' \deqn{
#'   K(x) = 1 - \exp\left( -\lambda \left( \frac{1 - (1+x^c)^{-k}}{(1+x^c)^{-k}} \right)^\beta \right)
#' }
#' and \eqn{m, n, p, \lambda, \beta, c, k > 0}, \eqn{x > 0}.
#'
#' @param x Value or vector of positive values to evaluate the PDF at.
#' @param m Parameter m of the distribution (m > 0).
#' @param n Parameter n of the distribution (n > 0).
#' @param p Parameter p of the distribution (p > 0).
#' @param lambda Parameter lambda of the distribution (lambda > 0).
#' @param beta Parameter beta of the distribution (beta > 0).
#' @param c Parameter c of the Burr XII baseline component (c > 0).
#' @param k Parameter k of the Burr XII baseline component (k > 0).
#'
#' @return A numeric vector containing the values of the PDF evaluated at \code{x}.
#' @export
#'
#' @examples
#' P_PTLENKBurr(
#'   x = c(0.5, 1.0, 1.5),
#'   m = 1, n = 1, p = 1,
#'   lambda = 0.5, beta = 1.5, c = 2, k = 1
#' )
P_PTLENKBurr <- function(x, m, n, p, lambda, beta, c, k) {
  K <- 1 - exp(-lambda * ((1 - (1 + x^c)^(-k)) / (1 + x^c)^(-k))^beta)
  term_bur1 <- 2 * m * n * p * lambda * beta * c * k * x^(c - 1) * (1 + x^c)^(-k - 1)
  term_bur2 <- (1 - (1 + x^c)^(-k))^(beta - 1) / ((1 + x^c)^(-k))^(beta + 1)
  term_tl1 <- (1 - K) / K^(p + 1) * exp(m * n * (1 - 1 / (K^p)))
  term_tl2 <- (1 - exp(n * (1 - 1 / (K^p)))) * (2 - exp(n * (1 - 1 / (K^p))))^(m - 1)

  result <- term_bur1 * term_bur2 * term_tl1 * term_tl2
  return(result)
}
