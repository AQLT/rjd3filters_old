#' Retrieve implicit forecasts corresponding to the asymmetric filters
#'
#' Function to retrieve the implicit forecasts corresponding to the asymmetric filters
#'
#' @details Let \eqn{h} be the bandwidth of the symmetric filter,
#' \eqn{v_{-h}, \ldots, v_h} the coefficients of the symmetric filter and
#' \eqn{w_{-h}^q, \ldots, w_h^q} the coefficients of the asymmetric filter used to estimate
#' the trend when \eqn{q} future values are known (with the convention \eqn{w_{q+1}^q=\ldots=w_h^q=0}).
#' Let denote \eqn{y_{-h},\ldots, y_0} the las \eqn{h} available values of the input times series.
#' The implicit forecasts, \eqn{y_{1}*,\ldots, y_h*} solve:
#' \deqn{
#' \forall q, \quad \sum_{i=-h}^0 v_iy_i + \sum_{i=1}^h v_iy_i*
#' =\sum_{i=-h}^0 w_i^qy_i + \sum_{i=1}^h w_i^qy_i*
#' }
#' which is equivalent to
#' \deqn{
#' \forall q, sum_{i=1}^h (v_i- w_i^q) y_i*
#' =\sum_{i=-h}^0 (w_i^q-v_i)y_i.
#' }
#' Note that this is solved numerically: the solution isn't exact.
#' @inheritParams jfilter
#' @examples
#' y <- retailsa$AllOtherGenMerchandiseStores
#' ql <- lpp_properties(horizon = 6, kernel = "Henderson", endpoints = "QL")$filters.coef
#' lc <- lpp_properties(horizon = 6, kernel = "Henderson", endpoints = "LC")$filters.coef
#' f_ql <- implicit_forecast(y, ql)
#' f_lc <- implicit_forecast(y, lc)
#'
#' plot(window(y, start = 2007),
#'      xlim = c(2007,2012))
#' lines(ts(c(tail(y,1), f_ql), frequency = frequency(y), start = end(y)),
#'       col = "red", lty = 2)
#' lines(ts(c(tail(y,1), f_lc), frequency = frequency(y), start = end(y)),
#'       col = "blue", lty = 2)
#' @importFrom stats time
#' @export
implicit_forecast <- function(y, coefs){
  UseMethod("implicit_forecast", y)
}
#' @export
implicit_forecast.default <- function(y, coefs){
  if(is.matrix(coefs)){
    coefs <- lapply(1:ncol(coefs), function(i) coefs[,i])
  }
  lags <- length(coefs)-1
  scoef <- coefs[[lags+1]]

  jsymf <- .jcall("jdplus/math/linearfilters/FiniteFilter",
                  "Ljdplus/math/linearfilters/FiniteFilter;",
                  "of", scoef, as.integer(-lags))
  jsymf <- .jcast(jsymf, "jdplus/math/linearfilters/IFiniteFilter")
  jasym <- lapply(lags:1,
                  function(i){
                    .jcall("jdplus/math/linearfilters/FiniteFilter",
                           "Ljdplus/math/linearfilters/FiniteFilter;",
                           "of", coefs[[i]],as.integer(-lags))
                  }
  )
  jasym <- .jarray(jasym,
                   "jdplus/math/linearfilters/IFiniteFilter")
  jy <- .jcall("demetra/data/DoubleSeq",
               "Ldemetra/data/DoubleSeq;", "of",as.numeric(tail(y,lags+1)))
  prev <- .jcall("jdplus/math/linearfilters/AsymmetricFilters",
         "[D","implicitForecasts",
         jsymf, jasym, jy)
  ts(prev,
     frequency = frequency(y),
     start = tail(time(y),1)+1/frequency(y))
}
#' @export
implicit_forecast.matrix <- function(y, coefs){
  result <- do.call(cbind, lapply(seq_len(ncol(y)), function (i) implicit_forecast(y[,i], coefs = coefs)))
  colnames(result) <- colnames(y)
  result
}
