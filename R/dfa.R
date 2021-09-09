#' @export
dfa_filter <- function(horizon = 6, degree = 0,
                       density = c("uniform", "rw"),
                       targetfilter = lp_filter(horizon = horizon)$filters.coef[,horizon+1],
                       passband = 2*pi/12,
                       accuracy.weight = 1/3,
                       smoothness.weight = 1/3,
                       timeliness.weight = 1/3){
  density = match.arg(density)
  rkhs_filter = J("demetra/saexperimental/r/DFAFilters")$filterProperties(
    targetfilter[1:7],
    as.integer(horizon), as.integer(degree), density=="rw",
    passband,
    accuracy.weight, smoothness.weight, timeliness.weight
  )
  return(structure(FiniteFilters2R(rkhs_filter, horizon, TRUE),
                   class="dfa_filter"))
}
#

