#' @include ts.R procresults.R
#' @import rJava
NULL

#' Title
#'
#' @param y
#' @param length
#' @param musgrave
#' @param ic
#'
#' @return
#' @export
#'
#' @examples
henderson<-function(y, length, musgrave=TRUE, ic=4.5){
  return (.jcall("demetra/saexperimental/r/X11Decomposition", "[D", "henderson", as.numeric(y), as.integer(length), musgrave, ic))
}

#' Title
#'
#' @param y
#' @param horizon
#' @param degree
#' @param kernel
#' @param endpoints
#' @param ic
#'
#' @return
#' @export
#'
#' @examples
localpolynomials<-function(y, horizon, degree=3, kernel=c("Henderson", "Uniform", "Biweight", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"), endpoints=c("LC", "QL", "CQ", "CC", "DAF"), ic=4.5){
  d<-2/(sqrt(pi)*ic)
  kernel=match.arg(kernel)
  endpoints=match.arg(endpoints)
  return (.jcall("demetra/saexperimental/r/LocalPolynomialFilters", "[D", "filter", as.numeric(y), as.integer(horizon), as.integer(degree), kernel, endpoints, d))
}


#' Title
#'
#' @param horizon
#' @param degree
#' @param kernel
#' @param endpoints
#' @param ic
#'
#' @return
#' @export
#'
#' @examples
filterproperties<-function(horizon, degree=3, kernel=c("Henderson", "Uniform", "Biweight", "Triweight", "Tricube", "Gaussian", "Triangular", "Parabolic"), endpoints=c("LC", "QL", "CQ", "CC", "DAF"), ic=4.5){
  d<-2/(sqrt(pi)*ic)
  kernel=match.arg(kernel)
  endpoints=match.arg(endpoints)
  jprops<-.jcall("demetra/saexperimental/r/LocalPolynomialFilters", "Ldemetra/saexperimental/r/FiltersToolkit$FiniteFilters;", "filterProperties", as.integer(horizon), as.integer(degree), kernel, endpoints, d)
  sw<-proc_data(jprops, "sweights")
  swg<-proc_data(jprops, "sgain")
  aw<-sapply(1:horizon, function(h){return(proc_data(jprops, paste0("aweights(", horizon-h,')')))})
  awg<-sapply(1:horizon, function(h){return(proc_data(jprops, paste0("again(", horizon-h,')')))})
  awp<-sapply(1:horizon, function(h){return(proc_data(jprops, paste0("aphase(", horizon-h,')')))})
  return(structure(list(
    internal=jprops,
    symmetricFilter=sw,
    symmetricfilter.gain=swg,
    asymmetricFilters=aw,
    asymmetricfilter.gain=awg,
    asymmetricfilter.phase=awp
    ),
    class="JD.Filters"))
}
