#' Direct Filter Approach
#'
#' @inheritParams localpolynomials
#' @inheritParams mse
#' @param targetfilter the weights of the symmetric target filters.
#' @param accuracy.weight,smoothness.weight,timeliness.weight the weight used for the
#' optimisation. The weight associated to the residual is derived so that the sum of
#' the four weights are equal to 1.
#' @export
dfa_filter <- function(horizon = 6, degree = 0,
                       density = c("uniform", "rw"),
                       targetfilter = lp_filter(horizon = horizon)$filters.coef[,horizon+1],
                       passband = 2*pi/12,
                       accuracy.weight = 1/3,
                       smoothness.weight = 1/3,
                       timeliness.weight = 1/3){
  density = match.arg(density)
  if(length(targetfilter) != 2*horizon + 1)
    stop("The symmetric targetfilter must be of length 2*horizon+1")
  dfa_filter = J("demetra/saexperimental/r/DFAFilters")$filterProperties(
    rev(targetfilter[seq_len(horizon+1)]),
    as.integer(horizon), as.integer(degree), density=="rw",
    passband,
    accuracy.weight, smoothness.weight, timeliness.weight
  )
  return(structure(FiniteFilters2R(dfa_filter, horizon, TRUE),
                   class=c("dfa_filter", "FiniteFilters")))
}
#

