#' Linear Asymmetric Filtering on a Time Series
#'
#' Applies linear filtering to a univariate time series or to each series separately of a multivariate time series using an asymmetric moving average
#'
#' @param y a univariate or multivariate time series.
#' @param coefs the coefficients of the asymmetric moving average.
#'
#' @details
#' For a time series \eqn{x} the asymmetric moving average is applied.
#' If \eqn{coef} is of length \eqn{q}, the result \eqn{y} is equal at time \eqn{t} to:
#' \deqn{y[t] = x[t-lags] * coef[1] + x[t-lags+1] * coef[1] + ... + x[t-lags+q] * coef[q]}
#'
#' @examples
#' lags <- 6
#' leads <- 2
#' y <- retailsa$DrinkingPlaces
#' fst_coef <- fstfilter(lags = lags, leads = leads, smoothness.weight = 0.3, timeliness.weight = 0.3)
#' fst_coef <- fst_coef$filter
#' lpp_coef <- lpp_properties(horizon = lags, kernel = "Henderson", endpoints = "LC")
#' lpp_coef <- lpp_coef$filters.coef[,"q=2"]
#'
#' fst_ma <- jasym_filter(y, fst_coef, lags)
#' lpp_ma <- jasym_filter(y, lpp_coef, lags)
#'
#' plot(ts.union(y, fst_ma, lpp_ma), plot.type = "single", col = c("black","red","blue"))
#' @inheritParams fstfilter
#' @export
jasym_filter <- function(y, coefs, lags){
  UseMethod("jasym_filter", y)
}
#' @export
jasym_filter.default <- function(y, coefs, lags){
  coefs <- removeTrailingZeroOrNA(coefs)
  leads <- length(coefs)-(lags+1)
  jsymf <- .jcall("jdplus/math/linearfilters/FiniteFilter",
                  "Ljdplus/math/linearfilters/FiniteFilter;",
                  "of", coefs,as.integer(-lags))
  jsymf <- .jcast(jsymf, "jdplus/math/linearfilters/IFiniteFilter")
  jsymf_none <- .jcall("jdplus/math/linearfilters/FiniteFilter",
                       "Ljdplus/math/linearfilters/FiniteFilter;",
                       "of", rep(0, lags+1),as.integer(-lags))
  jsymf_none <- .jcast(jsymf_none, "jdplus/math/linearfilters/IFiniteFilter")

  jy <- .jcall("demetra/data/DoubleSeq",
               "Ldemetra/data/DoubleSeq;", "of",as.numeric(y))

  jasym <- lapply(0:(lags-1),
                  function(q){
                    if(q<leads){
                      jsymf_none#Those values will me removed
                    }else{
                      jsymf
                    }
                  }
  )
  jasym <- .jarray(jasym,
                   "jdplus/math/linearfilters/IFiniteFilter")

  result <- .jcall("jdplus/math/linearfilters/FilterUtility",
                   "Ldemetra/data/DoubleSeq;", "filter",
                   jy,
                   jsymf,
                   .jnull("[Ljdplus/math/linearfilters/IFiniteFilter;"),
                   jasym
  )

  result <- result$toArray()
  if(leads>0){
    result[(length(result)+1-leads):length(result)] <- NA
  }
  if(is.ts(y))
    result <- ts(result,start = start(y), frequency = frequency(y))
  result
}
#' @export
jasym_filter.matrix <- function(y, coefs, lags){
  result <- y
  for (i in seq_len(ncol(y))){
    result[, i] <- jasym_filter(y[,i], coefs = coefs, lags = lags)
  }
  result
}
