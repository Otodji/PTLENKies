#' Histogram and Fitted Density Function (PDF) of the PTLENKBurr Distribution
#'
#' Generates a histogram of the provided dataset and overlays the fitted Probability
#' Density Function (PDF) of the PTLENKBurr distribution with parameters estimated
#' via Maximum Likelihood Estimation using \code{\link{E_PTLENKBurr}}.
#'
#' @param data Numeric vector of positive data values.
#' @param breaks Number of breaks or break specification for the histogram. Default is \code{15}.
#' @param col_hist Fill color for the histogram bars. Default is \code{"lightblue"}.
#' @param col_line Color for the fitted PDF curve. Default is \code{"firebrick3"}.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{hist}}.
#'
#' @return An invisible named numeric vector containing the 7 estimated parameters
#'   \code{c(m, n, p, lambda, beta, c, k)}.
#'
#' @details
#' The function uses \code{\link{E_PTLENKBurr}} to estimate distribution parameters
#' and evaluates \code{\link{P_PTLENKBurr}} over a sequence from near zero to 1.2 times
#' the maximum data value.
#'
#' @importFrom graphics grid hist legend lines
#' @export
#'
#' @examples
#' sample_data <- c(0.5, 1.2, 0.8, 2.1, 0.3, 1.5, 0.9)
#' Sim_PTLENKBurr(sample_data)
Sim_PTLENKBurr <- function(data, breaks = 15, col_hist = "lightblue", col_line = "firebrick3", ...) {

  # 1. Validate input data
  if (!is.numeric(data) || length(data) == 0) {
    stop("'data' must be a non-empty numeric vector.")
  }
  data <- data[!is.na(data)]
  if (any(data <= 0)) {
    stop("All elements of 'data' must be strictly positive (> 0).")
  }

  # 2. Estimate 7 parameters via E_PTLENKBurr
  estimated_params <- E_PTLENKBurr(data)

  # 3. Grid of x values for density curve
  x_seq <- seq(1e-5, max(data) * 1.2, length.out = 1000)

  # 4. Compute fitted probability density
  density_values <- P_PTLENKBurr(
    x = x_seq,
    m = estimated_params["m"],
    n = estimated_params["n"],
    p = estimated_params["p"],
    lambda = estimated_params["lambda"],
    beta = estimated_params["beta"],
    c = estimated_params["c"],
    k = estimated_params["k"]
  )
  density_values[is.nan(density_values) | is.na(density_values)] <- 0

  # 5. Adjust Y-axis scale to avoid clipping peak density
  h <- hist(data, plot = FALSE, breaks = breaks)
  max_y <- max(max(h$density), max(density_values)) * 1.15

  # 6. Plot empirical histogram
  hist(
    data,
    freq = FALSE,
    breaks = breaks,
    col = col_hist,
    border = "white",
    ylim = c(0, max_y),
    main = "Histogram with Fitted PTLENKBurr Density",
    xlab = "Values",
    ylab = "Density",
    ...
  )
  grid(col = "gray90", lty = "dotted")

  # 7. Overlay fitted theoretical density curve
  lines(x_seq, density_values, col = col_line, lwd = 2.5)

  # 8. Add legend
  legend(
    "topright",
    legend = c("Empirical Data", "Fitted PTLENKBurr PDF"),
    fill = c(col_hist, NA),
    border = c("white", NA),
    col = c(NA, col_line),
    lty = c(NA, 1),
    lwd = c(NA, 2.5),
    bg = "white",
    box.col = "gray80"
  )

  # 9. Return estimated parameters invisibly
  invisible(estimated_params)
}
