#' Maximum Likelihood Estimation (MLE) of the PTLENKBurr Distribution
#'
#' Estimates the parameters of the Power Transformed Lindley Extended Modified
#' Kies Burr XII (PTLENKBurr) distribution using maximum likelihood estimation (MLE)
#' via constrained optimization with \code{\link[stats]{nlminb}}.
#'
#' @param data Numeric vector of positive observations.
#' @param start_params Optional numeric vector of length 7 containing initial parameter values for
#'   \code{c(m, n, p, lambda, beta, c, k)}. Default is \code{c(1, 1, 1, 1, 1, 1, 1)}.
#' @param lower Lower bounds for parameters. Default is \code{1e-5} for all 7 parameters.
#'
#' @return A named numeric vector of estimated parameters: \code{m}, \code{n}, \code{p},
#'   \code{lambda}, \code{beta}, \code{c}, and \code{k}.
#'
#' @details
#' Optimization is performed using \code{\link[stats]{nlminb}}. Infeasible iterations
#' returning non-positive densities, \code{NA}, or \code{NaN} values are assigned a high penalty.
#'
#' @importFrom stats nlminb
#' @export
#'
#' @examples
#' data_sample <- c(0.5, 1.2, 0.8, 2.1, 0.3, 1.5, 0.9)
#' E_PTLENKBurr(data_sample)
E_PTLENKBurr <- function(data, start_params = NULL, lower = 1e-5) {

  # 1. Validate input data
  if (!is.numeric(data) || length(data) == 0) {
    stop("'data' must be a non-empty numeric vector.")
  }
  data <- data[!is.na(data)]
  if (any(data <= 0)) {
    stop("All elements of 'data' must be strictly positive (> 0).")
  }

  # 2. Handle initial parameter values
  if (is.null(start_params)) {
    start_params <- c(m = 1, n = 1, p = 1, lambda = 1, beta = 1, c = 1, k = 1)
  } else {
    if (length(start_params) != 7 || any(start_params <= 0)) {
      stop("'start_params' must be a numeric vector of length 7 with strictly positive values.")
    }
  }

  # 3. Define negative log-likelihood objective function
  neg_log_likelihood <- function(parameters) {
    m      <- parameters[1]
    n      <- parameters[2]
    p      <- parameters[3]
    lambda <- parameters[4]
    beta   <- parameters[5]
    c      <- parameters[6]
    k      <- parameters[7]

    # Strict parameter domain constraint check
    if (any(parameters <= 0)) {
      return(Inf)
    }

    # Evaluate density for each data point
    pdf_vals <- P_PTLENKBurr(x = data, m = m, n = n, p = p, lambda = lambda, beta = beta, c = c, k = k)

    # Numerical safety guard: return penalty if density returns NA, NaN, 0, or negative values
    if (any(is.na(pdf_vals)) || any(is.nan(pdf_vals)) || any(pdf_vals <= 0)) {
      return(1e10) # High penalty to guide optimization algorithm
    }

    log_lik <- sum(log(pdf_vals))
    return(-log_lik)
  }

  # 4. Define parameter lower bounds
  lower_bounds <- if (length(lower) == 1) rep(lower, 7) else lower

  # 5. Execute constrained optimization
  result <- suppressWarnings(
    nlminb(
      start = start_params,
      objective = neg_log_likelihood,
      lower = lower_bounds,
      upper = rep(Inf, 7),
      control = list(trace = 0, iter.max = 500, eval.max = 1000)
    )
  )

  # 6. Extract and label estimated parameters
  estimated_params <- result$par
  names(estimated_params) <- c("m", "n", "p", "lambda", "beta", "c", "k")

  return(estimated_params)
}
