#' Plot the Cumulative Distribution Function (CDF) of the PTLENKBurr Distribution
#'
#' Generates an enhanced plot of the Cumulative Distribution Function (CDF) of the
#' PTLENKBurr distribution over [0, +Inf), automatically calculating and
#' displaying the numerical median.
#'
#' @param m Parameter m of the distribution (m > 0).
#' @param n Parameter n of the distribution (n > 0).
#' @param p Parameter p of the distribution (p > 0).
#' @param lambda Parameter lambda of the distribution (lambda > 0).
#' @param beta Parameter beta of the distribution (beta > 0).
#' @param c Parameter c of the Burr XII baseline component (c > 0).
#' @param k Parameter k of the Burr XII baseline component (k > 0).
#' @param max_x Optional maximum x-axis value. If NULL (default), it is automatically computed using the CDF.
#' @param show_median Logical. If TRUE (default), calculates and plots the median.
#' @param col_line Color of the CDF curve. Default is \code{"deepskyblue4"}.
#' @param col_median Color of the median lines and point. Default is \code{"firebrick3"}.
#' @param legend_pos Position of the legend. Default is \code{"right"}.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#'
#' @return A plot of the CDF and an invisible \code{list} containing the evaluation grid
#'   data frame and the numerical value of the median.
#'
#' @details
#' The median is numerically calculated by finding the root of \eqn{F(x) - 0.5 = 0}
#' using \code{\link[stats]{uniroot}}.
#'
#' @importFrom graphics grid legend plot points segments
#' @importFrom stats uniroot
#' @export
#'
#' @examples
#' Plot_CPTLENKBurr(
#'   m = 1, n = 1, p = 1, lambda = 0.5,
#'   beta = 1.5, c = 2, k = 1
#' )
Plot_CPTLENKBurr <- function(m, n, p, lambda, beta, c, k,
                             max_x = NULL,
                             show_median = TRUE,
                             col_line = "deepskyblue4",
                             col_median = "firebrick3",
                             legend_pos = "right", ...) {

  # 1. Automatic upper bound estimation via CDF (up to F(x) = 0.999)
  if (is.null(max_x)) {
    max_x <- 5
    while (suppressWarnings(C_PTLENKBurr(max_x, m, n, p, lambda, beta, c, k)) < 0.999 && max_x < 1e5) {
      max_x <- max_x * 2
    }
  }

  # 2. Grid evaluation on [0, max_x]
  n_points <- 1000
  x_seq <- seq(0, max_x, length.out = n_points)

  # 3. CDF computation and handling of boundary/undefined values
  cdf_values <- C_PTLENKBurr(x = x_seq, m = m, n = n, p = p, lambda = lambda, beta = beta, c = c, k = k)
  cdf_values[is.nan(cdf_values) | is.na(cdf_values)] <- 0

  # 4. Numerical estimation of Median (F(x) = 0.5)
  median_val <- NA
  if (show_median) {
    target_fun <- function(x) C_PTLENKBurr(x, m, n, p, lambda, beta, c, k) - 0.5
    try({
      res <- uniroot(target_fun, interval = c(1e-8, max_x), extendInt = "upX")
      median_val <- res$root
    }, silent = TRUE)
  }

  # 5. Graphic configuration
  plot_args <- list(
    x = x_seq,
    y = cdf_values,
    type = "l",
    col = col_line,
    lwd = 2.5,
    xlab = "x",
    ylab = "F(x)",
    main = "CDF of the PTLENKBurr Distribution",
    xlim = c(0, max_x),
    ylim = c(0, 1)
  )

  # User arguments override
  user_args <- list(...)
  plot_args[names(user_args)] <- user_args

  # Render primary curve and grid
  do.call(plot, plot_args)
  grid(col = "gray90", lty = "dotted")

  # 6. Plot median projections and legend
  if (show_median && !is.na(median_val)) {
    # Orthogonal projections to axes
    segments(x0 = median_val, y0 = 0, x1 = median_val, y1 = 0.5, col = col_median, lty = 2, lwd = 1.8)
    segments(x0 = 0, y0 = 0.5, x1 = median_val, y1 = 0.5, col = col_median, lty = 2, lwd = 1.8)

    # Point marker at median (x_0.5, 0.5)
    points(median_val, 0.5, col = col_median, pch = 19, cex = 1.3)

    # Legend rendering
    legend(legend_pos,
           inset = c(0.02, 0.05),
           legend = c("CDF F(x)", sprintf("Median = %.4f", median_val)),
           col = c(plot_args$col, col_median),
           lty = c(1, 2),
           lwd = c(plot_args$lwd, 1.8),
           pch = c(NA, 19),
           bg = "white", box.col = "gray80")
  }

  # 7. Return evaluation data frame and median value invisibly
  invisible(list(
    data = data.frame(x = x_seq, cdf = cdf_values),
    median = median_val
  ))
}
