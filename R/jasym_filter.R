#' Linear Filtering on a Time Series
#'
#' Applies linear filtering to a univariate time series or to each series separately of a multivariate time series using either a moving average (symmetric or asymmetric) or a combination of
#' symmetric moving average at the center and asymmetric moving averages at the bounds.
#'
#' @param y a univariate or multivariate time series.
#' @param coefs an object of class \code{"moving_average"} or the coefficients
#' of a moving average. In that case the argument \code{lags} must be provided.
#' @param coefs a \code{matrix} or a \code{list} that contains all the coefficients of the asymmetric and symmetric filters.
#'  (from the symmetric filter to the shortest). See details.
#'
#' @param remove_missing if `TRUE` (default) leading and trailing NA are removed before filtering.
#'
#' @details
#'
#'
#' The functions \code{jfilter} extends \code{\link[stats]{filter}} allowing to apply every kind of moving averages
#' (symmetric and asymmetric filters) or to apply aset multiple moving averages
#' to deal with the boundaries.
#'
#' Let \eqn{x_t} be the input time series to filter.
#'
#' * If `coef` is an object [moving_average()], of length \eqn{q}, the result \eqn{y} is equal at time \eqn{t} to:
#' \deqn{y[t] = x[t-lags] * coef[1] + x[t-lags+1] * coef[1] + ... + x[t-lags+q] * coef[q]}.
#' It extends the function \code{\link[stats]{filter}} that would add \code{NA} at the end of the time series.
#'
#' * If `coef` is a `matrix`, `list` or [finite_filters()] object,  at the center,
#' the symmetric moving average is used (first column/element of \code{coefs}).
#' At the boundaries, the last moving average of \code{coefs} is used to compute the filtered
#' time series \eqn{y[n]} (no future point known), the second to last to compute the filtered
#' time series \eqn{y[n-1]} (one future point known)...
#'
#' @examples
#' y <- retailsa$DrinkingPlaces
#'
#' lags <- 6
#' leads <- 2
#' fst_coef <- fst_filter(lags = lags, leads = leads, smoothness.weight = 0.3, timeliness.weight = 0.3)
#' lpp_coef <- lp_filter(horizon = lags, kernel = "Henderson", endpoints = "LC")
#'
#' fst_ma <- jfilter(y, fst_coef)
#' lpp_ma <- jfilter(y, lpp_coef[,"q=2"])
#'
#' plot(ts.union(y, fst_ma, lpp_ma), plot.type = "single", col = c("black","red","blue"))
#'
#' trend <- jfilter(y, lpp_coef)
#' # This is equivalent to:
#' trend <- localpolynomials(y, horizon = 6)
#' @export
jfilter <- function(y, coefs, remove_missing = TRUE){
  UseMethod("jfilter", y)
}
#' @export
jfilter.default <- function(y, coefs, remove_missing = TRUE){
  if (is.moving_average(coefs)) {
    filter_ma(y, coefs)
  } else {
    ff_ma(y, coefs = coefs, remove_missing = remove_missing)
  }
}
#' @export
jfilter.matrix <- function(y, coefs, remove_missing = TRUE){
  result <- y
  for (i in seq_len(ncol(y))){
    result[, i] <- jfilter(y[,i], coefs = coefs, remove_missing = remove_missing)
  }
  result
}

filter_ma <- function(y, coefs){
  # if (!is.moving_average(coefs)) {
  #   coefs <- moving_average(coefs, -abs(lags))
  # }
  lb = lower_bound(coefs)
  ub = upper_bound(coefs)

  if (length(y) <= length(coefs))
    return(y * NA)

  DataBlock = J("jdplus.data.DataBlock")
  jy = DataBlock$of(as.numeric(y))
  out = DataBlock$of(as.numeric(rep(NA, length(y) - length(coefs)+1)))
  .ma2jd(coefs)$apply(jy,
                     out)
  result = out$toArray()
  result <- c(rep(NA, abs(min(lb, 0))),
              result,
              rep(NA, abs(max(ub, 0))))

  if (is.ts(y))
    result <- ts(result,start = start(y), frequency = frequency(y))
  result
}

ff_ma <- function(y, coefs, remove_missing = TRUE) {
  if (!inherits(coefs, "finite_filters")) {
    coefs <- finite_filters(coefs)
  }
  jffilters <- .finite_filters2jd(coefs)

  if (remove_missing) {
    data_clean = remove_bound_NA(y)
    y2 = data_clean$data
  } else {
    y2 = y
  }

  jy <- .jcall("demetra/data/DoubleSeq",
               "Ldemetra/data/DoubleSeq;", "of", as.numeric(y2))

  result <- .jcall("jdplus/math/linearfilters/FilterUtility",
                   "Ldemetra/data/DoubleSeq;", "filter",
                   jy,
                   jffilters$jsymf,
                   jffilters$jlasym,
                   jffilters$jrasym
  )

  result <- .jcall(result, "[D", "toArray")

  if (remove_missing){
    result = c(rep(NA, data_clean$leading), result,
               rep(NA, data_clean$trailing))
  }
  if(is.ts(y))
    result <- ts(result,start = start(y), frequency = frequency(y))
  result
}

.finite_filters2jd <- function(ff) {
  jsymf <- .ma2jd(ff@sfilter)
  rfilters <- ff@rfilters
  lfilters <- ff@lfilters
  if (length(rfilters) < upper_bound(ff@sfilter)) {
    # last points as NA
    rfilters <- c(rfilters,
                  rep(list(moving_average(NA, lags = 0)), upper_bound(ff@sfilter) - length(rfilters)))
  }
  if (length(lfilters) < -lower_bound(ff@sfilter)) {
    # first points as NA
    lfilters <- c(rep(list(moving_average(NA, lags = 0)), -lower_bound(ff@sfilter) - length(lfilters)),
                  lfilters)
  }
  jrasym <- lapply(rfilters, .ma2jd)
  jlasym <- rev(lapply(lfilters, .ma2jd))

  jsymf <- .jcast(jsymf,
                  "jdplus/math/linearfilters/IFiniteFilter")
  if (length(jrasym) == 0) {
    jrasym <- .jnull("[Ljdplus/math/linearfilters/IFiniteFilter;")
  } else {
    jrasym <- .jarray(jrasym,
                      "jdplus/math/linearfilters/IFiniteFilter")
  }
  if (length(jlasym) == 0) {
    jlasym <- .jnull("[Ljdplus/math/linearfilters/IFiniteFilter;")
  } else {
    jlasym <- .jarray(jlasym,
                      "jdplus/math/linearfilters/IFiniteFilter")
  }
  list(jsymf = jsymf,
       jrasym = jrasym,
       jlasym = jlasym)
}
