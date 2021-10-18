# x11_lp <- function(y, period = frequency(y),
#                    decomposition = c("Additive", "Multiplicative", "LogAdditive", "PseudoAdditive"),
#                    horizon=6, degree=2,
#                    kernel=c("Henderson", "BiWeight", "TriWeight", "TriCube", "Uniform", "Triangular", "Epanechnikov", "Trapezoidal"),
#                    endpoints=c("LC", "QL", "CQ", "CC", "DAF", "CN"),
#                    ic = 4.5,
#                    tweight = 0, passband = pi/12,
#                    seas.s0 = c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"),
#                    seas.s1 = c("S3X5", "S3X3", "S3X1", "S3X9", "S3X15"),
#                    extreme.lsig=1.5, extreme.usig=2.5){
#   iperiod <- as.integer(period)
#   if(abs(period - iperiod) < 1e-9){
#     P = iperiod
#     P = new( J("java.lang.Integer"), as.character(P))
#   }else{
#     P = period
#     P = new( J("java.lang.Double"), as.character(P))
#   }
#   P = .jcast(P,"java.lang.Number")
#   aparams<-2/(sqrt(pi)*ic)
#   kernel=match.arg(kernel)
#   endpoints=match.arg(endpoints)
#   seas0=match.arg(seas.s0)
#   seas1=match.arg(seas.s1)
#   endpoints=match.arg(endpoints)
#   decomposition = match.arg(decomposition)
#   asym <- switch(endpoints,
#                  DAF="Direct",
#                  CN="CutAndNormalize",
#                  "MMSRE")
#
#   KernelOption = J("jdplus.filters.KernelOption")
#   AsymmetricFiltersOption = J("jdplus.math.linearfilters.AsymmetricFilters$Option")
#
#   tspec = new(J("jdplus.filters.LocalPolynomialFilterSpec"))
#   tspec$setFilterLength(as.integer(horizon))
#   tspec$setPolynomialDegree(as.integer(degree))
#   tspec$setKernel(KernelOption$valueOf(kernel))
#   tspec$setAsymmetricFilters(AsymmetricFiltersOption$valueOf(asym));
#   if(asym == "MMSRE"){
#     adegree = switch (endpoints,
#       LC = 0L,
#       QL = 1L,
#       CQ = 2L
#     )
#     tspec$setAsymmetricPolynomialDegree(adegree);
#   }
#
#   tspec$setLinearModelCoefficients(.jarray(aparams));
#   tspec$setTimelinessWeight(tweight);
#   tspec$setPassBand(passband);
#
#   X11SeasonalFiltersFactory = J("jdplus.x11plus.X11SeasonalFiltersFactory")
#   SeasonalFilterOption = J("jdplus.x11plus.SeasonalFilterOption")
#
#   DecompositionMode = J("demetra.sa.DecompositionMode")
#   context = J("jdplus.x11plus.X11Context")$
#     builder()$
#     mode(DecompositionMode$valueOf(decomposition))$
#     period(P)$
#     trendFiltering(J("jdplus.filters.LocalPolynomialFilterFactory")$of(tspec))$
#     initialSeasonalFiltering(X11SeasonalFiltersFactory$filter(P, SeasonalFilterOption$valueOf(seas0)))$
#     finalSeasonalFiltering(X11SeasonalFiltersFactory$filter(P, SeasonalFilterOption$valueOf(seas1)))$
#     lowerSigma(extreme.lsig)$
#   upperSigma(extreme.usig)$
#   build()
#
#   mul = decomposition =="Multiplicative"
#   x11_res(y, context, mul)
# }
# x11_fst <- function(y, period = frequency(y),
#                    decomposition = c("Additive", "Multiplicative", "LogAdditive", "PseudoAdditive"),
#                    horizon=6, degree=2,
#                    smoothness.weight = 1, smoothness.degree = 3, timeliness.weight = 0,
#                    timeliness.passband = pi/6, timeliness.antiphase = TRUE,
#                    seas.s0 = c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"),
#                    seas.s1 = c("S3X5", "S3X3", "S3X1", "S3X9", "S3X15"),
#                    extreme.lsig=1.5, extreme.usig=2.5){
#   iperiod <- as.integer(period)
#   if(abs(period - iperiod) < 1e-9){
#     P = iperiod
#     P = new( J("java.lang.Integer"), as.character(P))
#   }else{
#     P = period
#     P = new( J("java.lang.Double"), as.character(P))
#   }
#   P = .jcast(P,"java.lang.Number")
#
#   seas0=match.arg(seas.s0)
#   seas1=match.arg(seas.s1)
#   decomposition = match.arg(decomposition)
#
#   tspec = new(J("jdplus.filters.FSTFilterSpec"))
#   tspec$setLags(as.integer(horizon))
#   tspec$setLeads(as.integer(horizon))
#   tspec$setPolynomialPreservationDegree(as.integer(degree))
#   tspec$setSmoothnessDegree(as.integer(smoothness.degree))
#   tspec$setAntiphase(timeliness.antiphase)
#   tspec$setSmoothnessWeight(smoothness.weight)
#   tspec$setTimelinessWeight(timeliness.weight)
#   tspec$setW0(0)
#   tspec$setW1(timeliness.passband)
#
#   X11SeasonalFiltersFactory = J("jdplus.x11plus.X11SeasonalFiltersFactory")
#   SeasonalFilterOption = J("jdplus.x11plus.SeasonalFilterOption")
#
#   DecompositionMode = J("demetra.sa.DecompositionMode")
#   context = J("jdplus.x11plus.X11Context")$
#     builder()$
#     mode(DecompositionMode$valueOf(decomposition))$
#     period(P)$
#     trendFiltering(J("jdplus.filters.FSTFilterFactory")$of(tspec))$
#     initialSeasonalFiltering(X11SeasonalFiltersFactory$filter(P, SeasonalFilterOption$valueOf(seas0)))$
#     finalSeasonalFiltering(X11SeasonalFiltersFactory$filter(P, SeasonalFilterOption$valueOf(seas1)))$
#     lowerSigma(extreme.lsig)$
#     upperSigma(extreme.usig)$
#     build()
#
#   mul = decomposition =="Multiplicative"
#   x11_res(y, context, mul)
# }
# x11_rkhs <- function(y, period = frequency(y),
#                      decomposition = c("Additive", "Multiplicative", "LogAdditive", "PseudoAdditive"),
#                      horizon = 6,
#                      kernel = c("BiWeight", "Henderson", "Epanechnikov", "Triangular", "Uniform", "TriWeight", "TriCube"),
#                      asymmetricCriterion = c("Timeliness", "FrequencyResponse", "Accuracy", "Smoothness", "Undefined"),
#                      density = c("uniform", "rw"),
#                      degree = 3,
#                      passband = 2*pi/12,
#                      seas.s0 = c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"),
#                      seas.s1 = c("S3X5", "S3X3", "S3X1", "S3X9", "S3X15"),
#                      extreme.lsig=1.5, extreme.usig=2.5){
#   iperiod <- as.integer(period)
#   if(abs(period - iperiod) < 1e-9){
#     P = iperiod
#     P = new( J("java.lang.Integer"), as.character(P))
#   }else{
#     P = period
#     P = new( J("java.lang.Double"), as.character(P))
#   }
#   P = .jcast(P,"java.lang.Number")
#   kernel=match.arg(kernel)
#   asymmetricCriterion=match.arg(asymmetricCriterion)
#   seas0=match.arg(seas.s0)
#   seas1=match.arg(seas.s1)
#   decomposition = match.arg(decomposition)
#   spectral = match.arg(density)
#
#   KernelOption = J("jdplus.filters.KernelOption")
#   AsymmetricCriterion= J("jdplus.filters.AsymmetricCriterion")
#   SpectralDensity = J("jdplus/filters/SpectralDensity")
#
#   tspec = new(J("jdplus.rkhs.RKHSFilterSpec"))
#
#   tspec$setFilterLength(as.integer(horizon))
#   tspec$setPolynomialDegree(as.integer(degree))
#   tspec$setKernel(KernelOption$valueOf(kernel))
#
#   # tspec$setOptimalBandWidth(optimalbw); # Ne semble rien changer
#   tspec$setAsymmetricBandWith(AsymmetricCriterion$valueOf(asymmetricCriterion));
#   tspec$setDensity(SpectralDensity$valueOf(switch(spectral,
#                                                   rw="RandomWalk",
#                                                   "WhiteNoise")))
#
#   X11SeasonalFiltersFactory = J("jdplus.x11plus.X11SeasonalFiltersFactory")
#   SeasonalFilterOption = J("jdplus.x11plus.SeasonalFilterOption")
#
#   DecompositionMode = J("demetra.sa.DecompositionMode")
#   context = J("jdplus.x11plus.X11Context")$
#     builder()$
#     mode(DecompositionMode$valueOf(decomposition))$
#     period(P)$
#     trendFiltering(J("jdplus.rkhs.RKHSFilterFactory")$of(tspec))$
#     initialSeasonalFiltering(X11SeasonalFiltersFactory$filter(P, SeasonalFilterOption$valueOf(seas0)))$
#     finalSeasonalFiltering(X11SeasonalFiltersFactory$filter(P, SeasonalFilterOption$valueOf(seas1)))$
#     lowerSigma(extreme.lsig)$
#     upperSigma(extreme.usig)$
#     build()
#
#   mul = decomposition =="Multiplicative"
#   x11_res(y, context, mul)
# }
# x11_res <- function(y, context, mul){
#   X11Kernel = new(J("jdplus.x11plus.X11Kernel"))
#
#   data = J("demetra.data.DoubleSeq")$of(as.numeric(y))
#   X11Kernel$process(data, context)
#
#   jrslt = J("demetra.saexperimental.r.X11Decomposition$Results")$
#     builder()$
#     y(data)$
#     kernel(X11Kernel)$
#     multiplicative(mul)$
#     build()
#
#   res <- ts.union(y,
#                   cbind(proc_vector(jrslt, "d11"),
#                         proc_vector(jrslt, "d12"),
#                         proc_vector(jrslt, "d10"),
#                         proc_vector(jrslt, "d13")))
#   colnames(res) <- c("y", "sa", "t", "s", "i")
#   res
# }
#' X-11 Decomposition With Custom Trend Filters
#'
#' Perform the X-11 decomposition using custom trend filter
#' @param y input time-series.
#' @param period period.
#' @param trend.coefs coefficients of the filters used for the trend-cycle extraction from
#' the real-time asymmetric filter to the symmetric filter. Can be a, object of class \code{"list"},
#' \code{"matrix"}, \code{"lp_filter"} or \code{"rkhs_filter"}.
#' @param mul boolean indicating if the decomposition mode is multiplicative.
#' @param seas.s0,seas.s1 seasonal filters.
#' @param extreme.lsig,extreme.usig boundaries used for outlier correction in irregular.
#' @examples
#' y <- retailsa$AllOtherGenMerchandiseStores
#' decomposition_lp <- x11(y, trend.coefs = lp_filter())
#' decomposition_rkhs <- x11(y, trend.coefs = rkhs_filter())
#' plot(y)
#' lines(decomposition_lp$decomposition[,"t"], col = "red")
#' lines(decomposition_rkhs$decomposition[,"t"], col = "green")
#' @importFrom stats ts.union
#' @export
x11 <- function(y, period = frequency(y),
                trend.coefs,  mul=TRUE,
                seas.s0=c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"),
                seas.s1=c("S3X5", "S3X3", "S3X1", "S3X9", "S3X15"),
                extreme.lsig=1.5, extreme.usig=2.5){
  seas0=match.arg(seas.s0)
  seas1=match.arg(seas.s1)

  if(inherits(trend.coefs, c("lp_filter","rkhs_filter"))){
    trend.coefs <- trend.coefs$filters.coef
  }

  if(is.matrix(trend.coefs)){
    trend.coefs <- lapply(1:ncol(trend.coefs),
                          function(i) trend.coefs[,i])
  }
  trend.coefs <- lapply(trend.coefs, removeTrailingZeroOrNA)
  sym_filter <- trend.coefs[[length(trend.coefs)]]
  asy_filter <- trend.coefs[-length(trend.coefs)]
  leftTrendFilter <- lapply(rev(asy_filter), rev)
  if(length(sym_filter) != 2*length(leftTrendFilter)+1){
    stop(sprintf("The symmetric filter is of length %i but only %i asymmetric filters provided",
                 length(sym_filter),
                 2*length(leftTrendFilter)+1))
  }
  mt <- J("jdplus.math.matrices.Matrix")
  ltrendf = mt$make(as.integer(length(sym_filter)-1), # nrows
                    as.integer(length(leftTrendFilter))#ncols
  )
  # for(i in seq_along(asy_filter)){
  #   ltrendf$column(i-1L)$drop(i-1L,0L)$
  #     copyFrom(asy_filter[[i]],0L)
  # }
  for(i in seq_along(leftTrendFilter)){
    ltrendf$column(i-1L)$drop(0L, i-1L)$
      copyFrom(leftTrendFilter[[i]],0L)
  }
  ctrendf = J("demetra.data.DoubleSeq")$of(sym_filter)
  x11decomp = J("demetra/saexperimental/r/X11Decomposition")
  jrslt = x11decomp$trendX11(as.numeric(y), period, mul,
                             ctrendf, ltrendf,
                             seas0, seas1, extreme.lsig, extreme.usig)
  decomposition <- ts.union(y,
                            cbind(proc_vector(jrslt, "d11"),
                                  proc_vector(jrslt, "d12"),
                                  proc_vector(jrslt, "d10"),
                                  proc_vector(jrslt, "d13")))
  colnames(decomposition) <- c("y", "sa", "t", "s", "i")
  decomposition
  return(structure(list(
    java=jrslt,
    decomposition=decomposition),
    class="X11_RSLT"))
}

