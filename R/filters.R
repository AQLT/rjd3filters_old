#' @include ts.R procresults.R
#' @import rJava
NULL

#' Apply Henderson linear filter
#'
#' @param y input time-series
#' @param length length of the Henderson filter
#' @param musgrave boolean indicating if musgrave asymmetric filters should be used
#' @param ic ic ratio
#'
#' @return the target signal
#' @export
#'
#' @examples
henderson<-function(y, length, musgrave=TRUE, ic=4.5){
  result <- .jcall("demetra/saexperimental/r/X11Decomposition", "[D", "henderson",
                   as.numeric(y), as.integer(length), musgrave, ic)

  if(is.ts(y))
    result <- ts(result,start = start(y), frequency = frequency(y))
  result
}

#' Apply local polynomials filters
#'
#' @inheritParams henderson
#' @param horizon horizon of the filter
#' @param degree degree of polynomial
#' @param kernel kernel uses
#' @param endpoints methode for endpoints
#'
#' @return the target signal
#' @export
#'
#' @examples
localpolynomials<-function(y,
                           horizon,
                           degree = 3,
                           kernel = c("Henderson", "Uniform", "Biweight", "Trapezoidal", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"),
                           endpoints = c("LC", "QL", "CQ", "CC", "DAF"),
                           ic = 4.5){
  if(2*horizon < degree)
    stop("You need more observation (2 * horizon + 1) than variables (degree + 1) to estimate the filter.")

  d<-2/(sqrt(pi)*ic)
  kernel=match.arg(kernel)
  endpoints=match.arg(endpoints)
  result <- .jcall("demetra/saexperimental/r/LocalPolynomialFilters", "[D", "filter",
                   as.numeric(y), as.integer(horizon), as.integer(degree), kernel, endpoints, d)
  if(is.ts(y))
    result <- ts(result,start = start(y), frequency = frequency(y))
  result
}

#' Apply asymmetric polynomials filters
#'
#' @inheritParams localpolynomials
#' @param q number of available data for the assymetric filter.
#'  By default (\code{q = 0}) a real-time filter is used.
#'
#' @return the target signal
#' @export
#'
#' @examples
asymmetric_lp<-function(y,
                        horizon,
                        degree = 3,
                        kernel = c("Henderson", "Uniform", "Biweight", "Trapezoidal", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"),
                        endpoints = c("LC", "QL", "CQ", "CC", "DAF"),
                        ic = 4.5,
                        q = 0){
  first_date <- time(y)[1] + (horizon*2+q)/frequency(y)
  last_date <- time(y)[length(y)]-horizon/frequency(y)
  available_span <- time(ts(0, start =first_date,
                            end = last_date,
                            frequency = frequency(y)))
  ts(sapply(as.numeric(available_span), function(date_fin){
    res <- localpolynomials(window(y, end = date_fin),
                            horizon = horizon,
                            degree = degree,
                            kernel = kernel,
                            endpoints = endpoints,
                            ic = ic)
    res[length(res) - q]
  }),start = available_span[1],
  frequency = frequency(y))
}


#' Get properties of local polynomials filters
#'
#' @inheritParams localpolynomials
#' @details
#' * "LC": Linear-Constant filter $y_t$
#' * "QL": Quadratic-Linear filter
#' * "CQ": Cubic-Quadratic filter
#' * "CC": Constant-Constant filter
#' * "DAF": Direct Asymmetric filter
#' * "CN": Cut and Normalized Filter
#'
#' @return \code{list} with coefficients, gain and phase values
#' @export
#'
#' @examples
#' filterproperties(horizon = 3)
filterproperties <- function(horizon, degree = 3,
                           kernel = c("Henderson", "Uniform", "Biweight", "Trapezoidal", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"),
                           endpoints = c("LC", "QL", "CQ", "CC", "DAF", "CN"),
                           ic = 4.5){
  if(2*horizon < degree)
    stop("You need more observation (2 * horizon + 1) than variables (degree + 1) to estimate the filter.")
  d<-2/(sqrt(pi)*ic)
  kernel=match.arg(kernel)
  endpoints=match.arg(endpoints)
  jprops<-.jcall("demetra/saexperimental/r/LocalPolynomialFilters",
                 "Ldemetra/saexperimental/r/FiltersToolkit$FiniteFilters;",
                 "filterProperties", as.integer(horizon),
                 as.integer(degree), kernel, endpoints, d)
  sw<-proc_data(jprops, "sweights")
  swg<-proc_data(jprops, "sgain")
  aw<-sapply(0:(horizon-1), function(h){return(proc_data(jprops, paste0("aweights(", h,')')))})
  awg<-sapply(0:(horizon-1), function(h){return(proc_data(jprops, paste0("again(", h,')')))})
  awp<-sapply(0:(horizon-1), function(h){return(proc_data(jprops, paste0("aphase(", h,')')))})

  svariancereduction <- proc_data(jprops, "svariancereduction")
  avariancereduction <- sapply(0:(horizon-1), function(h){proc_data(jprops, paste0("avariancereduction(", h,')'))})
  abias0 <- sapply(0:(horizon-1), function(h){proc_data(jprops, paste0("abias0(", h,')'))})
  abias1 <- sapply(0:(horizon-1), function(h){proc_data(jprops, paste0("abias1(", h,')'))})
  abias2 <- sapply(0:(horizon-1), function(h){proc_data(jprops, paste0("abias2(", h,')'))})


  coefs = c(aw,list(sw))
  nbpoints = horizon*2+1
  coefs = sapply(coefs, function(x){
    c(x,rep(0,nbpoints-length(x)))
  })

  gain = cbind(awg,swg)
  phase = cbind(awp, 0)

  bias <- rbind(abias0, abias1, abias2)
  bias <- cbind(bias,
                c(sum(sw),
                  sum(seq(-horizon,horizon) * sw),
                  sum(seq(-horizon,horizon)^2 * sw)))
  variancereduction <- c(avariancereduction, svariancereduction)
  diagnostics <- rbind(variancereduction, bias)
  rownames(diagnostics) <- c("Variance reduction",
                             "Constant bias",
                             "Linear bias",
                             "Quadratic bias")

  filternames <- sprintf("q=%i", 0:(horizon))
  rownames(coefs) <- sprintf("t%+i", seq(-horizon,horizon))
  rownames(coefs) <- sub("+0", "", rownames(coefs), fixed = TRUE)
  colnames(gain) <- colnames(coefs) <-
    colnames(diagnostics) <- colnames(phase) <-
    filternames

  return(structure(list(
    internal = jprops,
    filters.coef = coefs,
    filters.gain = gain,
    filters.phase= phase,
    filters.diagnostics = diagnostics
  ),
  class="JD.Filters"))
}


#' Estimation of a filter using the Fidelity-Smoothness-Timeliness criterions
#'
#' @param lags Lags of the filter (should be positive)
#' @param leads Leads of the filter (should be positive or 0)
#' @param pdegree Local polynomials preservation: max degree
#' @param passband Passband of the filter (only used if tweight is )
#' @param smoothness.weight Weight for the smoothness criterion (in [0, 1])
#' @param smoothness.degree Degree of the smoothness criterion (3 for Henderson)
#' @param timeliness.weight Weight for the Timeliness criterion (in [0, 1[). sweight+tweight should be in [0,1]
#' @param timeliness.passband Passband for the timeliness criterion (in radians). The phase effect is computed in [0, passband]
#' @param timeliness.antiphase See Guggemos
#'
#' @return Contains the selected filter, its gain and phase, and the values of the 3 criterions
#' @export
#'
#' @examples
fstfilter<-function(lags, leads, pdegree=2, smoothness.weight=1, smoothness.degree=3, timeliness.weight=0, timeliness.passband=pi/6, timeliness.antiphase=T){
  jobj<-.jcall("demetra/saexperimental/r/FiltersToolkit", "Ldemetra/saexperimental/r/FiltersToolkit$FSTResult;",
               "fstfilter", as.integer(lags), as.integer(leads), as.integer(pdegree), smoothness.weight, as.integer(smoothness.degree),
              timeliness.weight, timeliness.passband, as.logical(timeliness.antiphase))
  return(fstresult(jobj))

}

fstresult<-function(jobj){
  filter<-.jcall(jobj,"[D", "weights")
  gain<-.jcall(jobj, "[D", "getGain")
  phase<-.jcall(jobj, "[D", "getPhase")
  criterions<-.jcall(jobj, "[D", "getCriterions")

  return(list(
    filter=filter,
    gain=gain,
    phase=phase,
    criterions=criterions
  ))
}

#' FST criterions
#'
#' @param weights Weights of the filter (from lower bound to upper bound)
#' @param lb Lower bound (usually negative)
#' @param passband Passband threshold for timeliness criterion
#'
#' @return
#' @export
#'
#' @examples
fst<-function(weights, lb, passband=pi/6){
  jobj<-.jcall("demetra/saexperimental/r/FiltersToolkit", "Ldemetra/saexperimental/r/FiltersToolkit$FSTResult;", "fst",
               weights, as.integer(lb), passband)
               return(fstresult(jobj))

}

#' Accuracy/smoothness/timeliness criterions through spectral decomposition
#' See Wildi/McElroy
#'
#' @param sweights Weights of the symmetric filter (from 0 to n only).
#' @param aweights Weights of the asymmetric filter (from -n to m)
#' @param density hypothesis on the spectral density
#' @param passband passband threshold
#'
#' @return The criterions
#' @export
#'
#' @examples
mse<-function(sweights, aweights, density=c("rw", "uniform"), passband=pi/6 ){
  spectral=match.arg(spectral)
  rslt<-.jcall("demetra/saexperimental/r/FiltersToolkit", "[D", "mseDecomposition",
                 sweights, aweights, spectral, passband)
  return (list(accuracy=rslt[1], smoothness=rslt[2], timeliness=rslt[3], residual=rslt[4]))
}

