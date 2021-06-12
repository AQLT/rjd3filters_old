#' Reproducing Kernel Hilbert Space (RKHS) filters
#'
#' Estimation of a filter using Reproducing Kernel Hilbert Space (RKHS)
#' @inheritParams localpolynomials
#' @inheritParams mse
#' @param asymmetricCriterion the criteria used to compute the optimal bandwidth. If \code{"Undefined"}, \eqn{m+1} is used.
#' @param optimalbw boolean indicating if the bandwith should be choosen by optimisation (with if \code{asymmetricCriterion} parameter).
#' If \code{optimalbw = FALSE} then the bandwith specified in \code{bw} will be  used.
#' @param bw the bandwidth to use if \code{optimalbw = FALSE}.
#' @references Dagum, Estela Bee and Silvia Bianconcini (2008). “The Henderson Smoother in Reproducing Kernel Hilbert Space”. In: Journal of Business & Economic Statistics 26, pp. 536–545. URL: \url{https://ideas.repec.org/a/bes/jnlbes/v26y2008p536-545.html}.
#' @examples
#' rkhs <- rkhs_filter(horizon = 6, asymmetricCriterion = "Timeliness")
#' plot_coef(rkhs)
#' @return An object of class \code{"rkhs_filter"}, which is a list of 4 elements:\itemize{
#' \item{\code{"internal"}}{Java object used for internal functions}
#' \item{\code{"filters.coef"}}{The coefficients of the selected filter}
#' \item{\code{"filters.gain"}}{The gain function between 0 and pi (601 observations)}
#' \item{\code{"filters.phase"}}{The phase function between 0 and pi (601 observations)}
#' }
#' @export
rkhs_filter <- function(horizon = 6, degree = 2,
                        kernel = c("BiWeight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "TriWeight", "TriCube"),
                        asymmetricCriterion = c("Timeliness", "FrequencyResponse", "Accuracy", "Smoothness", "Undefined"),
                        density = c("uniform", "rw"),
                        passband = 2*pi/12,
                        optimalbw = TRUE,
                        bw = horizon + 1){
  if(optimalbw){
    jprops <- rkhs_filter_optimalbw(horizon = horizon,
                                    degree = degree,
                                    kernel = kernel,
                                    asymmetricCriterion = asymmetricCriterion,
                                    density = density,
                                    passband = passband)
  }else{
    jprops <- rkhs_filter_bw(horizon = horizon,
                             degree = degree,
                             kernel = kernel,
                             bandWidth = bw)
  }

  return(structure(FiniteFilters2R(jprops, horizon),
                   class="rkhs_filter"))
}

rkhs_filter_bw <- function(horizon = 6,
                           degree = 2,
                           kernel = c("BiWeight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "TriWeight", "TriCube"),
                           bandWidth = horizon + 1){
  kernel = match.arg(kernel)
  jkernel = kernelToDoubleUnaryOperator(kernel, degree, horizon)
  afilter <- lapply(0:(horizon-1),function(i){
    .jcall("jdplus/rkhs/CutAndNormalizeFilters",
           "Ljdplus/math/linearfilters/FiniteFilter;",
           "of",
           jkernel,
           bandWidth, as.integer(horizon),
           as.integer(i))
  })
  afilter <- .jarray(afilter,
                     "jdplus/math/linearfilters/FiniteFilter")
  sfilter <- .jcall("jdplus/rkhs/KernelsUtility",
                    "Ljdplus/math/linearfilters/SymmetricFilter;",
                    "symmetricFilter",
                    jkernel,
                    horizon + 1,
                    as.integer(horizon))
  builder <- .jcall("demetra/saexperimental/r/FiltersToolkit$FiniteFilters",
                    "Ldemetra/saexperimental/r/FiltersToolkit$FiniteFilters$FiniteFiltersBuilder;",
                    "builder")
  builder$filter(sfilter)
  builder$afilters(afilter)
  jprops <- builder$build()
  jprops
}

rkhs_filter_optimalbw <- function(horizon = 6, degree = 2,
                                  kernel = c("BiWeight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "TriWeight", "TriCube"),
                                  asymmetricCriterion = c("Undefined", "FrequencyResponse", "Accuracy", "Smoothness", "Timeliness"),
                                  density = c("uniform", "rw"),
                                  passband = 2*pi/12){

  kernel = match.arg(kernel)
  asymmetricCriterion = match.arg(asymmetricCriterion)
  spectral = match.arg(density)

  KernelOption = J("jdplus.filters.KernelOption")
  AsymmetricCriterion= J("jdplus.filters.AsymmetricCriterion")
  SpectralDensity = J("jdplus/filters/SpectralDensity")

  tspec = new(J("jdplus.rkhs.RKHSFilterSpec"))

  tspec$setFilterLength(as.integer(horizon))
  tspec$setPolynomialDegree(as.integer(degree))
  tspec$setKernel(KernelOption$valueOf(kernel))

  # tspec$setOptimalBandWidth(optimalbw); # Ne semble rien changer
  tspec$setAsymmetricBandWith(AsymmetricCriterion$valueOf(asymmetricCriterion));
  tspec$setDensity(SpectralDensity$valueOf(switch(spectral,
                                                  rw="RandomWalk",
                                                  "WhiteNoise")))
  rkhs_filter <- J("jdplus.rkhs.RKHSFilterFactory")$of(tspec)

  sfilter = rkhs_filter$symmetricFilter()
  afilter = rkhs_filter$endPointsFilters()
  # reverse the order of the asymmetric filters
  afilter <- lapply(rev(seq_along(afilter)),
                    function(i){
                      afilter[[i]]
                    }
  )
  afilter <- .jarray(afilter,
                     "jdplus/math/linearfilters/FiniteFilter")
  builder <- .jcall("demetra/saexperimental/r/FiltersToolkit$FiniteFilters",
                    "Ldemetra/saexperimental/r/FiltersToolkit$FiniteFilters$FiniteFiltersBuilder;",
                    "builder")
  builder$filter(sfilter)
  builder$afilters(afilter)
  jprops <- builder$build()

  return(jprops)
}

kernelToDoubleUnaryOperator <- function(kernel, degree, horizon){
  jkernel <- J("jdplus/stats/Kernels")
  jHighOrderKernels <- J("jdplus/rkhs/HighOrderKernels")
  res = switch(kernel,
               BiWeight=jHighOrderKernels$kernel(jkernel$BIWEIGHT,
                                                 as.integer(degree)),
               TriWeight = jHighOrderKernels$kernel(jkernel$TRIWEIGHT,
                                                    as.integer(degree)),
               Uniform = jHighOrderKernels$kernel(jkernel$UNIFORM,
                                                  as.integer(degree)),
               Triangular = jHighOrderKernels$kernel(jkernel$TRIANGULAR,
                                                     as.integer(degree)),
               Epanechnikov = jHighOrderKernels$kernel(jkernel$EPANECHNIKOV,
                                                       as.integer(degree)),
               Henderson = jHighOrderKernels$kernel(jkernel$henderson(as.integer(horizon)),
                                                    as.integer(degree))
  )
  .jcast(res,"java.util.function.DoubleUnaryOperator")
}
