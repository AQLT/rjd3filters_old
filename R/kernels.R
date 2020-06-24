#' Get the coefficients of a kernel
#'
#' Function to get the coefficient associated to a kernel. Those
#' coefficients are then used to compute the different filters.
#'
#' @inheritParams localpolynomials
#'
#' @return \code{tskernel} object (see \link[stats]{kernel}).
#' @export
#'
#' @examples
#' get_kernel("Henderson", horizon = 3)
get_kernel <- function(kernel = c("Henderson","Uniform", "Triangular",
                                  "Epanechnikov","Parabolic","Biweight", "Triweight","Tricube",
                                  "Trapezoidal", "Gaussian"),
                       horizon){
  kernel = match.arg(kernel)
  if(kernel == "Parabolic")
    kernel = "Epanechnikov"
  h <- horizon
  if(kernel == "Gaussian")
    h <- 4 * horizon
  jkernel <- .jcall("jdplus/data/analysis/DiscreteKernel",
                    "Ljava/util/function/IntToDoubleFunction;",
                    tolower(kernel), as.integer(h))
  coef = sapply(as.integer(seq.int(from = -horizon, to = horizon, by = 1)),
                jkernel$applyAsDouble)
  m = horizon
  result <- list(coef = coef, m = m)
  attr(result, "name") <- kernel
  attr(result, "class") <- "tskernel"
  result
}
