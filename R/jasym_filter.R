#' Linear Asymmetric Filtering on a Time Series
#'
#' Applies linear filtering to a univariate time series or to each series separately of a multivariate time series using an asymmetric moving average
#'
#' @param y a univariate or multivariate time series.
#' @param coefs an object of class \code{"moving_average"} or the coefficients
#' of a moving average. In that case the argument \code{lags} must be provided.
#' @inheritParams fst_filter
#'
#' @details
#' For a time series \eqn{x} the asymmetric moving average is applied.
#' If \eqn{coef} is of length \eqn{q}, the result \eqn{y} is equal at time \eqn{t} to:
#' \deqn{y[t] = x[t-lags] * coef[1] + x[t-lags+1] * coef[1] + ... + x[t-lags+q] * coef[q]}
#' It extends the function \code{\link[stats]{filter}} that would add \code{NA} at the end of the time series.
#'
#' @examples
#' lags <- 6
#' leads <- 2
#' y <- retailsa$DrinkingPlaces
#' fst_coef <- fst_filter(lags = lags, leads = leads, smoothness.weight = 0.3, timeliness.weight = 0.3)
#' fst_coef <- fst_coef$filters.coef
#' lpp_coef <- lp_filter(horizon = lags, kernel = "Henderson", endpoints = "LC")
#' lpp_coef <- lpp_coef$filters.coef[,"q=2"]
#'
#' fst_ma <- jasym_filter(y, fst_coef, lags)
#' lpp_ma <- jasym_filter(y, lpp_coef, lags)
#'
#' plot(ts.union(y, fst_ma, lpp_ma), plot.type = "single", col = c("black","red","blue"))
#' @seealso \code{\link{jfilter}}
#' @export
jasym_filter <- function(y, coefs, lags){
  UseMethod("jasym_filter", y)
}
#' @export
jasym_filter.default <- function(y, coefs, lags){
  if(!is.moving_average(coefs)){
    coefs <- moving_average(coefs, -abs(lags))
  }
  lb = get_lower_bound(coefs)
  ub = get_upper_bound(coefs)

  DataBlock = J("jdplus.data.DataBlock")
  jy = DataBlock$of(as.numeric(y))
  out = DataBlock$of(as.numeric(rep(NA, jy$length() - length(coefs)+1)))
  ma2jd(coefs)$apply(jy,
                     out)
  result = out$toArray()
  result <- c(rep(NA, abs(min(lb, 0))),
              result,
              rep(NA, abs(max(ub, 0))))
  if(is.ts(y))
    result <- ts(result,start = start(y), frequency = frequency(y))
  result
}
# jasym_filter.default <- function(y, coefs, lags){
#   coefs <- removeTrailingZeroOrNA(coefs)
#   leads <- length(coefs)-(lags+1)
#   jsymf <- .jcall("jdplus/math/linearfilters/FiniteFilter",
#                   "Ljdplus/math/linearfilters/FiniteFilter;",
#                   "of", coefs,as.integer(-lags))
#   jsymf <- .jcast(jsymf, "jdplus/math/linearfilters/IFiniteFilter")
#   jsymf_none <- .jcall("jdplus/math/linearfilters/FiniteFilter",
#                        "Ljdplus/math/linearfilters/FiniteFilter;",
#                        "of", rep(0, lags+1),as.integer(-lags))
#   jsymf_none <- .jcast(jsymf_none, "jdplus/math/linearfilters/IFiniteFilter")
#
#   jy <- .jcall("demetra/data/DoubleSeq",
#                "Ldemetra/data/DoubleSeq;", "of",as.numeric(y))
#
#   jasym <- lapply(0:(lags-1),
#                   function(q){
#                     if(q<leads){
#                       jsymf_none#Those values will be removed
#                     }else{
#                       jsymf
#                     }
#                   }
#   )
#   jasym <- .jarray(jasym,
#                    "jdplus/math/linearfilters/IFiniteFilter")
#
#   result <- .jcall("jdplus/math/linearfilters/FilterUtility",
#                    "Ldemetra/data/DoubleSeq;", "filter",
#                    jy,
#                    jsymf,
#                    .jnull("[Ljdplus/math/linearfilters/IFiniteFilter;"),
#                    jasym
#   )
#
#   result <- result$toArray()
#   if(leads>0){
#     result[(length(result)+1-leads):length(result)] <- NA
#   }
#   if(is.ts(y))
#     result <- ts(result,start = start(y), frequency = frequency(y))
#   result
# }
#' @export
jasym_filter.matrix <- function(y, coefs, lags){
  result <- y
  for (i in seq_len(ncol(y))){
    result[, i] <- jasym_filter(y[,i], coefs = coefs, lags = lags)
  }
  result
}


#' Linear Filtering on a Time Series
#'
#' Applies linear filtering to a univariate time series or to each series separately of a multivariate time series using a
#' symmetric moving average at the center and an asymmetric moving average at the bounds.
#'
#' @param coefs a \code{matrix} or a \code{list} that contains all the coefficients of the asymmetric and symmetric filters.
#'  (from the shortest to the symmetric filter). See details.
#'
#' @param remove_missing if `TRUE` (default) leading and trailing NA are removed before filtering.
#'
#' @details
#' The functions \code{jfilter} extends \code{\link[stats]{filter}} allowing to set multiple moving averages
#' to deal with the boundaries.
#'
#' At the center, the symmetric moving average is used (last column/element of \code{coefs}).
#' At the boundaries, the first moving average of \code{coefs} is used for computing the filtered
#' time series \eqn{y[n]} (no future point known), the second for computes the filtered
#' time series \eqn{y[n-1]} (one future point known)...
#'
#' @examples
#' y <- retailsa$AllOtherGenMerchandiseStores
#'
#' lpp_ma <- lp_filter(horizon = 6)$filters.coef
#' trend <- jfilter(y, lpp_ma)
#' # This is equivalent to:
#' trend <- localpolynomials(y, horizon = 6)
#' @inheritParams jasym_filter
#' @seealso \code{\link{jasym_filter}}
#' @export
jfilter <- function(y, coefs, remove_missing = TRUE){
  UseMethod("jfilter", y)
}
#' @export
jfilter.default <- function(y, coefs, remove_missing = TRUE){
  if(is.matrix(coefs)){
    coefs <- lapply(1:ncol(coefs), function(i) coefs[,i])
  }
  lags <- length(coefs)-1
  scoef <- coefs[[lags+1]]

  jsymf <- .jcast(ma2jd(moving_average(scoef, -lags)), "jdplus/math/linearfilters/IFiniteFilter")


  jrasym <- lapply(lags:1,
                   function(i){
                     ma2jd(moving_average(removeTrailingZeroOrNA(coefs[[i]]), -lags))
                   }
  )
  jrasym <- .jarray(jrasym,
                    "jdplus/math/linearfilters/IFiniteFilter")
  jlasym <- lapply(lags:1,
                   function(i){
                     x = removeTrailingZeroOrNA(coefs[[i]])
                     p = get_upper_bound(moving_average(x, -lags))
                     ma2jd(moving_average(rev(x), -p))
                   }
  )
  jlasym <- .jarray(jlasym,
                    "jdplus/math/linearfilters/IFiniteFilter")
  if (remove_missing){
    data_clean = remove_bound_NA(y)
    y2 = data_clean$data
  } else {
    y2 = y
  }


  jy <- .jcall("demetra/data/DoubleSeq",
               "Ldemetra/data/DoubleSeq;", "of",as.numeric(y2))

  result <- .jcall("jdplus/math/linearfilters/FilterUtility",
                   "Ldemetra/data/DoubleSeq;", "filter",
                   jy,
                   jsymf,
                   jlasym,
                   jrasym
  )
  result <- result$toArray()

  if (remove_missing){
    result = c(rep(NA, data_clean$leading), result,
               rep(NA, data_clean$trailing))
  }
  if(is.ts(y))
    result <- ts(result,start = start(y), frequency = frequency(y))
  result
}
#' @export
jfilter.matrix <- function(y, coefs){
  result <- y
  for (i in seq_len(ncol(y))){
    result[, i] <- jfilter(y[,i], coefs = coefs)
  }
  result
}

