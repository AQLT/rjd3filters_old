#' @import rJava
NULL

#' Apply Henderson Linear Filter
#'
#' @param x input time-series.
#' @param length length of the Henderson filter.
#' @param musgrave boolean indicating if musgrave asymmetric filters should be used.
#' @param ic ic ratio.
#'
#' @return the target signal.
#' @examples
#' x <- retailsa$AllOtherGenMerchandiseStores
#' trend <- henderson(x, length = 13)
#' plot(x)
#' lines(trend, col = "red")
#' @importFrom stats is.ts na.omit start
#' @export
henderson<-function(x, length, musgrave=TRUE, ic=4.5){
  result <- .jcall("jdplus/experimentalsa/base/r/X11Decomposition", "[D", "henderson",
                   as.numeric(x), as.integer(length), musgrave, ic)

  if(is.ts(x))
    result <- ts(result,start = start(x), frequency = frequency(x))
  result
}

#' Apply Local Polynomials Filters
#'
#' @inheritParams henderson
#' @param horizon horizon (bandwidth) of the symmetric filter.
#' @param degree degree of polynomial.
#' @param kernel kernel uses.
#' @param endpoints methode for endpoints.
#' @param tweight timeliness weight.
#' @param passband passband threshold.
#'
#' @return the target signal
#' @examples
#' x <- retailsa$AllOtherGenMerchandiseStores
#' trend <- localpolynomials(x, horizon = 6)
#' plot(x)
#' lines(trend, col = "red")
#' @references Proietti, Tommaso and Alessandra Luati (2008). “Real time estimation in local polynomial regression, with application to trend-cycle analysis”.
#' @seealso \code{\link{lp_filter}}
#' @export
localpolynomials<-function(x,
                           horizon = 6,
                           degree = 3,
                           kernel = c("Henderson", "Uniform", "Biweight", "Trapezoidal", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"),
                           endpoints = c("LC", "QL", "CQ", "CC", "DAF"),
                           ic = 4.5,
                           tweight = 0, passband = pi/12){
  if(2*horizon < degree)
    stop("You need more observation (2 * horizon + 1) than variables (degree + 1) to estimate the filter.")

  d<-2/(sqrt(pi)*ic)
  kernel=match.arg(kernel)
  endpoints=match.arg(endpoints)
  result <- .jcall("jdplus/experimentalsa/base/r/LocalPolynomialFilters", "[D", "filter",
                   as.numeric(x), as.integer(horizon), as.integer(degree), kernel, endpoints, d,
                   tweight, passband)
  if(is.ts(x))
    result <- ts(result,start = start(x), frequency = frequency(x))
  result
}

#' Local Polynomials Filters
#'
#' @inheritParams localpolynomials
#' @details
#' * "LC": Linear-Constant filter
#' * "QL": Quadratic-Linear filter
#' * "CQ": Cubic-Quadratic filter
#' * "CC": Constant-Constant filter
#' * "DAF": Direct Asymmetric filter
#' * "CN": Cut and Normalized Filter
#'
#' @return \code{list} with coefficients, gain and phase values
#'
#' @examples
#' henderson_f <- lp_filter(horizon = 6, kernel = "Henderson")
#' plot_coef(henderson_f)
#' @references Proietti, Tommaso and Alessandra Luati (2008). “Real time estimation in local polynomial regression, with application to trend-cycle analysis”.
#' @seealso \code{\link{lp_filter}}
#' @return An object of class \code{"rkhs_filter"}, which is a list of 4 elements:\itemize{
#' \item{\code{"internal"}}{Java object used for internal functions}
#' \item{\code{"filters.coef"}}{The coefficients of the selected filter}
#' \item{\code{"filters.gain"}}{The gain function between 0 and pi (601 observations)}
#' \item{\code{"filters.phase"}}{The phase function between 0 and pi (601 observations)}
#' }
#' @export
lp_filter <- function(horizon = 6, degree = 3,
                      kernel = c("Henderson", "Uniform", "Biweight", "Trapezoidal", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"),
                      endpoints = c("LC", "QL", "CQ", "CC", "DAF", "CN"),
                      ic = 4.5,
                      tweight = 0, passband = pi/12){
  if(2*horizon < degree)
    stop("You need more observation (2 * horizon + 1) than variables (degree + 1) to estimate the filter.")
  d<-2/(sqrt(pi)*ic)
  kernel=match.arg(kernel)
  endpoints=match.arg(endpoints)
  jprops<-.jcall("jdplus/experimentalsa/base/r/LocalPolynomialFilters",
                 "Ljdplus/experimentalsa/base/r/FiltersToolkit$FiniteFilters;",
                 "filterProperties", as.integer(horizon),
                 as.integer(degree), kernel, endpoints, d,
                 tweight, passband)

  return(.jd2r_finitefilters(jprops, first_to_last = FALSE))
}
coefficients_names <- function(lb, ub){
  x <- sprintf("t%+i", seq(lb,ub))
  x <- sub("+0", "", x, fixed = TRUE)
  x
}


