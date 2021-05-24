#' Reproducing Kernel Hilbert Space (RKHS) filters
#'
#' Estimation of a filter using Reproducing Kernel Hilbert Space (RKHS)
#' @inheritParams localpolynomials
#' @inheritParams mse
#' @param asymmetricCriterion the criteria used to compute the optimal bandwidth. If \code{"Undefined"}, \eqn{m+1} is used.
#' @references Dagum, Estela Bee and Silvia Bianconcini (2008). “The Henderson Smoother in Reproducing Kernel Hilbert Space”. In: Journal of Business & Economic Statistics 26, pp. 536–545. URL: \url{https://ideas.repec.org/a/bes/jnlbes/v26y2008p536-545.html}.
#' @examples
#' rkhs <- rkhs_filter(horizon = 6, asymmetricCriterion = "Timeliness")
#' @export
rkhs_filter <- function(horizon = 6, degree = 2,
                        kernel = c("BiWeight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "TriWeight", "TriCube"),
                        # optimalbw = FALSE,
                        asymmetricCriterion = c("Undefined", "FrequencyResponse", "Accuracy", "Smoothness", "Timeliness"),
                        density = c("uniform", "rw"),
                        passband = 2*pi/12){
  tspec = .jnew("jdplus/rkhs/RKHSFilterSpec")
  kernel = match.arg(kernel)
  asymmetricCriterion = match.arg(asymmetricCriterion)
  spectral = match.arg(density)
  tspec$setFilterLength(as.integer(horizon))
  tspec$setPolynomialDegree(as.integer(degree))
  tspec$setKernel(.jcall("jdplus/filters/KernelOption",
                         "Ljdplus/filters/KernelOption;",
                         "valueOf", kernel))
  # tspec$setOptimalBandWidth(optimalbw); # Ne semble rien changer
  tspec$setAsymmetricBandWith(.jcall("jdplus/filters/AsymmetricCriterion", "Ljdplus/filters/AsymmetricCriterion;",
                                     "valueOf",asymmetricCriterion));
  if(spectral=="rw"){
    tspec$setDensity(.jcall("jdplus/filters/SpectralDensity", "Ljdplus/filters/SpectralDensity;",
                            "valueOf","RandomWalk"));
  } else{
    tspec$setDensity(.jcall("jdplus/filters/SpectralDensity", "Ljdplus/filters/SpectralDensity;",
                            "valueOf","WhiteNoise"))
  }
  rkhs_filter <- .jcall("jdplus/rkhs/RKHSFilterFactory","Ljdplus/filters/ISymmetricFiltering;",
                        "of",tspec)
  sfilter = rkhs_filter$symmetricFilter()
  afilter = rkhs_filter$endPointsFilters()

  builder <- .jcall("demetra/saexperimental/r/FiltersToolkit$FiniteFilters",
                    "Ldemetra/saexperimental/r/FiltersToolkit$FiniteFilters$FiniteFiltersBuilder;",
                    "builder")
  builder$filter(sfilter)
  builder$afilters(afilter)
  jprops <- builder$build()

  sw<-proc_data(jprops, "sweights")
  swg<-proc_data(jprops, "sgain")
  aw<-sapply((horizon-1):0, function(h){return(proc_data(jprops, paste0("aweights(", h,')')))})
  awg<-sapply((horizon-1):0, function(h){return(proc_data(jprops, paste0("again(", h,')')))})
  awp<-sapply((horizon-1):0, function(h){return(proc_data(jprops, paste0("aphase(", h,')')))})

  coefs = c(aw,list(sw))
  nbpoints = horizon*2+1
  coefs = sapply(coefs, function(x){
    c(x,rep(0,nbpoints-length(x)))
  })

  gain = cbind(awg,swg)
  phase = cbind(awp, 0)

  filternames <- sprintf("q=%i", 0:(horizon))
  rownames(coefs) <- coefficients_names(-horizon, horizon)
  colnames(gain) <- colnames(coefs) <- colnames(phase) <-
    filternames

  return(structure(list(
    internal = jprops,
    filters.coef = coefs,
    filters.gain = gain,
    filters.phase= phase
  ),
  class="rkhs_filter"))
}
