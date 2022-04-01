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
  leftTrendFilter <- do.call(cbind, lapply(leftTrendFilter, function(x){
    c(x, rep(0, length(sym_filter) - 1 - length(x)))
  }))

  if(length(sym_filter) != 2*ncol(leftTrendFilter)+1){
    stop(sprintf("The symmetric filter is of length %i but only %i asymmetric filters provided",
                 length(sym_filter),
                 2*ncol(leftTrendFilter)+1))
  }
  # mt <- J("jdplus.math.matrices.FastMatrix")
  # ltrendf = mt$make(as.integer(length(sym_filter)-1), # nrows
  #                   as.integer(length(leftTrendFilter))#ncols
  # )
  # for(i in seq_along(leftTrendFilter)){
  #   ltrendf$column(i-1L)$drop(0L, i-1L)$
  #     copyFrom(leftTrendFilter[[i]],0L)
  # }
  ltrendf = matrix_r2jd(leftTrendFilter)

  ctrendf = J("demetra.data.DoubleSeq")$of(sym_filter)
  x11decomp = J("demetra/saexperimental/r/X11Decomposition")
  jrslt = x11decomp$trendX11(as.numeric(y), period, mul,
                             ctrendf, ltrendf,
                             seas0, seas1, extreme.lsig, extreme.usig)
  decomposition <- cbind(y,
                         proc_vector(jrslt, "d11"),
                         proc_vector(jrslt, "d12"),
                         proc_vector(jrslt, "d10"),
                         proc_vector(jrslt, "d13"))
  colnames(decomposition) <- c("y", "sa", "t", "s", "i")
  return(structure(list(
    java=jrslt,
    decomposition=decomposition),
    class="X11_RSLT"))
}

