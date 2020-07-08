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


