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
  return (.jcall("demetra/saexperimental/r/X11Decomposition", "[D", "henderson", as.numeric(y), as.integer(length), musgrave, ic))
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
  d<-2/(sqrt(pi)*ic)
  kernel=match.arg(kernel)
  endpoints=match.arg(endpoints)
  return (.jcall("demetra/saexperimental/r/LocalPolynomialFilters", "[D", "filter",
                 as.numeric(y), as.integer(horizon), as.integer(degree), kernel, endpoints, d))
}


#' Get properties of local polynomials filters
#'
#' @inheritParams localpolynomials
#'
#' @return \code{list} with coefficients, gain and phase values
#' @export
#'
#' @examples
#' filterproperties(horizon = 3)
filterproperties <- function(horizon, degree = 3,
                           kernel = c("Henderson", "Uniform", "Biweight", "Trapezoidal", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"),
                           endpoints = c("LC", "QL", "CQ", "CC", "DAF"),
                           ic = 4.5){
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
  coefs = c(aw,list(sw))
  nbpoints = horizon*2+1
  coefs = sapply(coefs, function(x){
    c(x,rep(0,nbpoints-length(x)))
  })
  gain = cbind(awg,swg)
  filternames <- sprintf("q=%i", 0:(horizon))
  rownames(coefs) <- sprintf("t%+i", seq(-horizon,horizon))
  rownames(coefs) <- sub("+0", "", rownames(coefs), fixed = TRUE)
  colnames(gain) <- colnames(coefs) <- filternames
  colnames(awp) <- filternames[-length(filternames)]

  return(structure(list(
    internal=jprops,
    filters.coef = coefs,
    filters.gain = gain,
    asymmetricfilter.phase=awp
  ),
  class="JD.Filters"))
}
# filterproperties<-function(horizon, degree=3,
#                            kernel=c("Henderson", "Uniform", "Biweight", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"),
#                            endpoints=c("LC", "QL", "CQ", "CC", "DAF"),
#                            ic=4.5){
#   d<-2/(sqrt(pi)*ic)
#   kernel=match.arg(kernel)
#   endpoints=match.arg(endpoints)
#   jprops<-.jcall("demetra/saexperimental/r/LocalPolynomialFilters",
#                  "Ldemetra/saexperimental/r/FiltersToolkit$FiniteFilters;",
#                  "filterProperties", as.integer(horizon),
#                  as.integer(degree), kernel, endpoints, d)
#   sw<-proc_data(jprops, "sweights")
#   swg<-proc_data(jprops, "sgain")
#   asymetric.names <- sprintf("q=%i", 0:(horizon-1))
#   aw<-sapply(0:(horizon-1), function(h){return(proc_data(jprops, paste0("aweights(", h,')')))})
#   awg<-sapply(0:(horizon-1), function(h){return(proc_data(jprops, paste0("again(", h,')')))})
#   awp<-sapply(0:(horizon-1), function(h){return(proc_data(jprops, paste0("aphase(", h,')')))})
#   names(aw) <- colnames(awg) <- colnames(awp) <- asymetric.names
#   return(structure(list(
#     internal=jprops,
#     symmetricFilter.coef=sw,
#     symmetricfilter.gain=swg,
#     asymmetricFilters.coef=aw,
#     asymmetricfilter.gain=awg,
#     asymmetricfilter.phase=awp
#   ),
#   class="JD.Filters"))
# }
