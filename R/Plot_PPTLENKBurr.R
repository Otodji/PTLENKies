#' Plot the Probability Density Function (PDF) of the PTLENKBurr Distribution
#'
#' Generates an enhanced plot of the Probability Density Function (PDF) of the
#' PTLENKBurr distribution over [0, +Inf), automatically calculating and
#' displaying the numerical mode and median.
#'
#' @param m Parameter m of the distribution (m > 0).
#' @param n Parameter n of the distribution (n > 0).
#' @param p Parameter p of the distribution (p > 0).
#' @param lambda Parameter lambda of the distribution (lambda > 0).
#' @param beta Parameter beta of the distribution (beta > 0).
#' @param c Parameter c of the Burr XII baseline component (c > 0).
#' @param k Parameter k of the Burr XII baseline component (k > 0).
#' @param max_x Optional maximum x-axis value. If NULL (default), it is automatically computed using the CDF.
#' @param show_mode Logical. If TRUE (default), calculates and plots the mode.
#' @param show_median Logical. If TRUE (default), calculates and plots the median.
#' @param col_line Color of the PDF curve. Default is \code{"deepskyblue4"}.
#' @param col_mode Color of the mode marker. Default is \code{"darkorange2"}.
#' @param col_median Color of the median marker. Default is \code{"firebrick3"}.
#' @param legend_pos Position of the legend. Default is \code{"topright"}.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#'
#' @return A plot of the PDF and an invisible \code{list} containing the evaluation grid
#'   data frame, calculated mode, and median.
#'
#' @details
#' The upper bound \code{max_x} is determined automatically by finding where \eqn{F(x) \ge 0.999}.
#' The mode is calculated numerically by maximizing the PDF via \code{\link[stats]{optimize}},
#' and the median is found via \code{\link[stats]{uniroot}} on the CDF (\code{\link{C_PTLENKBurr}}).
#'
#' @importFrom graphics grid legend plot points segments
#' @importFrom stats optimize uniroot
#' @export
#'
#' @examples
#' Plot_PPTLENKBurr(
#'   m = 1, n = 1, p = 1, lambda = 0.5,
#'   beta = 1.5, c = 2, k = 1
#' )
Plot_PPTLENKBurr <- function(m, n, p, lambda, beta, c, k,
                             max_x = NULL,
                             show_mode = TRUE,
                             show_median = TRUE,
                             col_line = "deepskyblue4",
                             col_mode = "darkorange2",
                             col_median = "firebrick3",
                             legend_pos = "topright", ...) {

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

  # 3. PDF computation and handling of boundary/undefined values
  pdf_values <- P_PTLENKBurr(x = x_seq, m = m, n = n, p = p, lambda = lambda, beta = beta, c = c, k = k)
  pdf_values[is.nan(pdf_values) | is.na(pdf_values) | is.infinite(pdf_values)] <- 0

  # 4. Numerical estimation of Mode (maximizing f(x))
  mode_x <- NA
  mode_y <- NA
  if (show_mode) {
    pdf_fun <- function(x) P_PTLENKBurr(x, m, n, p, lambda, beta, c, k)
    try({
      opt <- optimize(pdf_fun, interval = c(1e-8, max_x), maximum = TRUE)
      mode_x <- opt$maximum
      mode_y <- opt$objective
    }, silent = TRUE)
  }

  # 5. Numerical estimation of Median (F(x) = 0.5)
  median_x <- NA
  median_y <- NA
  if (show_median) {
    target_fun <- function(x) C_PTLENKBurr(x, m, n, p, lambda, beta, c, k) - 0.5
    try({
      res <- uniroot(target_fun, interval = c(1e-8, max_x), extendInt = "upX")
      median_x <- res$root
      median_y <- P_PTLENKBurr(median_x, m, n, p, lambda, beta, c, k)
    }, silent = TRUE)
  }

  # 6. Y-axis scaling for legend headroom
  max_y_plot <- max(pdf_values, na.rm = TRUE) * 1.25

  plot_args <- list(
    x = x_seq,
    y = pdf_values,
    type = "l",
    col = col_line,
    lwd = 2.5,
    xlab = "x",
    ylab = "f(x)",
    main = "PDF of the PTLENKBurr Distribution",
    xlim = c(0, max_x),
    ylim = c(0, max_y_plot)
  )

  # User arguments override
  user_args <- list(...)
  plot_args[names(user_args)] <- user_args

  # Render primary curve and grid
  do.call(plot, plot_args)
  grid(col = "gray90", lty = "dotted")

  # 7. Plot Mode and Median markers
  legend_labels <- c("PDF f(x)")
  legend_cols <- c(plot_args$col)
  legend_ltys <- c(1)
  legend_pchs <- c(NA)

  # Mode marker
  if (show_mode && !is.na(mode_x)) {
    segments(x0 = mode_x, y0 = 0, x1 = mode_x, y1 = mode_y, col = col_mode, lty = 2, lwd = 1.8)
    points(mode_x, mode_y, col = col_mode, pch = 19, cex = 1.3)

    legend_labels <- c(legend_labels, sprintf("Mode = %.4f", mode_x))
    legend_cols <- c(legend_cols, col_mode)
    legend_ltys <- c(legend_ltys, 2)
    legend_pchs <- c(legend_pchs, 19)
  }

  # Median marker
  if (show_median && !is.na(median_x)) {
    segments(x0 = median_x, y0 = 0, x1 = median_x, y1 = median_y, col = col_median, lty = 3, lwd = 1.8)
    points(median_x, median_y, col = col_median, pch = 17, cex = 1.3)

    legend_labels <- c(legend_labels, sprintf("Median = %.4f", median_x))
    legend_cols <- c(legend_cols, col_median)
    legend_ltys <- c(legend_ltys, 3)
    legend_pchs <- c(legend_pchs, 17)
  }

  # Legend rendering
  legend(legend_pos,
         inset = c(0.02, 0.03),
         legend = legend_labels,
         col = legend_cols,
         lty = legend_ltys,
         lwd = c(plot_args$lwd, rep(1.8, length(legend_labels) - 1)),
         pch = legend_pchs,
         bg = "white", box.col = "gray80")

  # 8. Return evaluation data frame and point estimates invisibly
  invisible(list(
    data = data.frame(x = x_seq, pdf = pdf_values),
    mode = mode_x,
    median = median_x
  ))
}
