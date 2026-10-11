#' Dataset: GlassFibers
#'
#' This dataset contains measurements corresponding to the breaking strengths of 1.5 cm glass fibers.
#' The data are presented as a numeric vector.
#'
#' @name GlassFibers
#' @docType data
#' @format A numeric vector of breaking strengths.
#' @description The dataset contains observations collected at the UK National Physical Laboratory
#' regarding the breaking strengths of 1.5 cm glass fibers, previously reported in Atchade (2023).
#'
#' @references
#' Atchade, N. (2023). Statistical modeling of lifetime data and probability distributions.
#'
#' @usage data(GlassFibers)
#'
#' @keywords datasets
GlassFibers <- c(
  3.200, 3.400, 3.500, 3.500, 3.600, 3.600, 3.900, 4.150, 4.300, 4.400, 4.400, 4.400, 4.540, 4.600, 4.700, 4.700, 4.730, 4.750, 4.790, 4.800, 4.800, 4.820,
  4.900, 4.950, 5.100, 5.100, 5.150, 5.200, 5.200, 5.300, 5.400, 5.400, 5.400, 5.400, 5.400, 5.600, 5.600, 5.700, 5.700, 5.700, 5.800,
  5.800, 5.800, 5.800, 5.900, 5.900, 5.900, 5.900, 5.900, 5.900, 5.900, 5.900, 5.900, 5.900, 6.000, 6.000, 6.000, 6.000, 6.000, 6.000,
  6.000, 6.000, 6.200, 6.300, 6.400, 6.500, 6.500, 6.500, 6.600, 6.700, 6.800, 6.800, 6.900, 7.000, 7.000, 7.000, 7.100, 7.300, 7.400,
  7.500, 7.900, 8.000, 8.200, 8.500, 8.600
)

# Enregistrement du jeu de données dans le dossier data/ du package
# usethis::use_data(GlassFibers, overwrite = TRUE)
# ou alternativement :
# save(GlassFibers, file = "data/GlassFibers.rda")
