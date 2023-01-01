#'@importFrom methods is
setClass("moving_average",
         slots = c(coefficients = "numeric",
                   lower_bound = "numeric",
                   upper_bound = "numeric")
)
#' Manipulation of moving averages
#'
#' @param x vector of coefficients
#' @param lags integer indicating the number of lags of the moving average.
#' @param trailing_zero,leading_zero boolean indicating wheter to remove leading/trailing zero and NA.
#' @param s seasonal period for the \code{to_seasonal()} function.
#' @param ...,na.rm see \link[base]{mean}
#' @param i,value indices specifying elements to extract or replace and the new value
#' @param object,e1,e2 moving averages
#'
#' @examples
#' y <- retailsa$AllOtherGenMerchandiseStores
#' e1 <- moving_average(rep(1,12), lags = -6)
#' e1 <- e1/sum(e1)
#' e2 <- moving_average(rep(1/12, 12), lags = -5)
#' M2X12 <- (e1 + e2)/2
#' coef(M2X12)
#' M3 <- moving_average(rep(1/3, 3), lags = -1)
#' M3X3 <- M3 * M3
#' # M3X3 moving average applied to each month
#' M3X3
#' M3X3_seasonal <- to_seasonal(M3X3, 12)
#' # M3X3_seasonal moving average applied to the global series
#' M3X3_seasonal
#'
#' def.par <- par(no.readonly = TRUE)
#' par(mai = c(0.5, 0.8, 0.3, 0))
#' layout(matrix(c(1,2), nrow = 1))
#' plot_gain(M3X3, main = "M3X3 applied to each month")
#' plot_gain(M3X3_seasonal, main = "M3X3 applied to the global series")
#' par(def.par)
#'
#' # To apply the moving average
#' t <- y * M2X12
#' # Or use the jasym_filter() function:
#' t <- jasym_filter(y, M2X12)
#' si <- y - t
#' s <- si * M3X3_seasonal
#' # or equivalently:
#' s_mm <- M3X3_seasonal * (1 - M2X12)
#' s <- y * s_mm
#' plot(s)
#' @export
moving_average <- function(x, lags = -length(x), trailing_zero = FALSE, leading_zero = FALSE){
  if (inherits(x, "moving_average"))
    return (x)
  x <- as.numeric(x)
  if (trailing_zero)
    x <- rm_trailing_zero_or_na(x)
  if (leading_zero){
    new_x <- rm_leading_zero_or_na(x)
    lags <- lags - (length(new_x) - length(x))
    x <- new_x
  }
  upper_bound = lags + length(x) -1
  # remove 1 if it is >= 0 (central term)
  # upper_bound = upper_bound - (upper_bound >= 0)

  names(x) <- coefficients_names(lags,
                                 upper_bound)
  res <- new("moving_average",
             coefficients = x, lower_bound = lags,
             upper_bound = upper_bound)
  res
}
jd2ma <- function(jobj, trailing_zero = FALSE){
  x <- .jcall(jobj, "[D", "weightsToArray")
  lags <- .jcall(jobj, "I", "getLowerBound")
  moving_average(x, lags, trailing_zero = trailing_zero)
}
ma2jd <- function(x){
  lags <- lower_bound(x)
  coefs = as.numeric(coef(x))
  if (length(x) == 1){
    coefs <- .jarray(coefs)
  }
  .jcall("jdplus/math/linearfilters/FiniteFilter",
         "Ljdplus/math/linearfilters/FiniteFilter;",
         "of", coefs,
         as.integer(lags))
}
#' @rdname moving_average
#' @export
is.moving_average <- function(x){
  is(x, "moving_average")
}
#' @importFrom stats coef
#' @export
coef.moving_average <- function(object, ...){
  coefs = object@coefficients
  return(coefs)
}
#' @rdname moving_average
#' @export
is_symmetric <- function(x){
  # .jcall(ma2jd(x), "Z", "isSymmetric")
  (upper_bound(x) == (-lower_bound(x))) &&
    isTRUE(all.equal(coef(x), rev(coef(x)), check.attributes = FALSE))
}
#' @rdname moving_average
#' @export
upper_bound <- function(x){
  x@upper_bound
}
#' @rdname moving_average
#' @export
lower_bound <- function(x){
  x@lower_bound
}
#' @rdname moving_average
#' @export
mirror <- function(x){
  jd2ma(.jcall(ma2jd(x), "Ljdplus/math/linearfilters/FiniteFilter;", "mirror"))
}
#' @method rev moving_average
#' @rdname moving_average
#' @export
rev.moving_average <- function(x){
  mirror(x)
}
#' @rdname moving_average
#' @export
length.moving_average <- function(x){
  length(coef(x))
}
#' @rdname moving_average
#' @export
to_seasonal <- function(x, s){
  UseMethod("to_seasonal", x)
}
#' @export
to_seasonal.default <- function(x, s){
  lb <- lower_bound(x)
  up <- upper_bound(x)
  coefs <- coef(x)
  new_coefs <- c(unlist(lapply(coefs[-length(x)],
                               function(x){
                                 c(x, rep(0, s - 1))
                               })),
                 coefs[length(x)])
  moving_average(new_coefs, lb * s)
}

#' @rdname moving_average
#' @export
sum.moving_average <- function(..., na.rm = FALSE){
  sum(
    unlist(lapply(list(...),
                  function(x) sum(coef(x),na.rm = na.rm)
    )
    )
  )
}
#' @rdname moving_average
#' @export
setMethod("[",
          signature(x = "moving_average",
                    i = "numeric"),
          function(x, i) {
            coefs <- coef(x)
            indices <- seq_along(coefs)[i]
            coefs[-indices] <- 0
            if (all(coefs == 0))
              return(moving_average(0, lags = 0))

            moving_average(coefs, lags = lower_bound(x),
                           leading_zero = TRUE, trailing_zero = TRUE)
          })
#' @rdname moving_average
#' @export
setMethod("[",
          signature(x = "moving_average",
                    i = "logical"),
          function(x, i) {
            coefs <- coef(x)
            indices <- seq_along(coefs)[i]
            coefs[!indices] <- 0
            moving_average(coefs, lags = lower_bound(x),
                           leading_zero = TRUE, trailing_zero = TRUE)
          })
#' @rdname moving_average
#' @export
setReplaceMethod("[",
                 signature(x = "moving_average",
                           i = "ANY",
                           j = "missing",
                           value = "numeric"),
                 function(x, i, value) {
                   x@coefficients[i] <- value
                   x
                 })
#' @rdname moving_average
#' @export
cbind.moving_average <- function(...){
  all_mm <- list(...)
  new_lb = min(sapply(all_mm, lower_bound))
  new_ub = max(sapply(all_mm, upper_bound))
  nb_uterms = max(sapply(all_mm, function(x) lower_bound(x) + length(x)))
  new_mm <- lapply(all_mm, function(x){
    c(rep(0, abs(new_lb - lower_bound(x))),
      coef(x),
      rep(0, abs(nb_uterms - (lower_bound(x) + length(x))))
    )
  })
  new_mm <- do.call(cbind, new_mm)
  rownames(new_mm) <- coefficients_names(new_lb, new_ub)
  new_mm
}
#' @rdname moving_average
#' @export
setMethod("+",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.math.linearfilters.FiniteFilter")
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/math/linearfilters/FiniteFilter;",
                           "add",
                           .jcast(ma2jd(e1), "jdplus/math/linearfilters/IFiniteFilter"),
                           .jcast(ma2jd(e2), "jdplus/math/linearfilters/IFiniteFilter"))

            jd2ma(jobj)
          })
#' @rdname moving_average
#' @export
setMethod("+",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 + moving_average(e2,0)
          })
#' @rdname moving_average
#' @export
setMethod("+",
          signature(e1 = "numeric",
                    e2 = "moving_average"),
          function(e1, e2) {
            e2 + e1
          })
#' @rdname moving_average
#' @export
setMethod("+", signature(e1 = "moving_average", e2 = "missing"), function(e1,e2) e1)
#' @rdname moving_average
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "missing"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.math.linearfilters.FiniteFilter")
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/math/linearfilters/FiniteFilter;",
                           "negate",
                           .jcast(ma2jd(e1), "jdplus/math/linearfilters/IFiniteFilter"))
            jd2ma(jobj)
          })
#' @rdname moving_average
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.math.linearfilters.FiniteFilter")
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/math/linearfilters/FiniteFilter;",
                           "subtract",
                           .jcast(ma2jd(e1), "jdplus/math/linearfilters/IFiniteFilter"),
                           .jcast(ma2jd(e2), "jdplus/math/linearfilters/IFiniteFilter"))
            jd2ma(jobj)
          })
#' @rdname moving_average
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 + (- e2)
          })
#' @rdname moving_average
#' @export
setMethod("-",
          signature(e1 = "numeric",
                    e2 = "moving_average"),
          function(e1, e2) {
            e1 + (- e2)
          })
#' @rdname moving_average
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.math.linearfilters.FiniteFilter")
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/math/linearfilters/FiniteFilter;",
                           "multiply",
                           .jcast(ma2jd(e1), "jdplus/math/linearfilters/IFiniteFilter"),
                           .jcast(ma2jd(e2), "jdplus/math/linearfilters/IFiniteFilter"))
            jd2ma(jobj)
          })
#' @rdname moving_average
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            if (length(e2) == 1) {
              e1 * moving_average(e2,0)
            } else {
              jasym_filter(e2, e1)
            }
          })

#' @rdname moving_average
#' @export
setMethod("*",
          signature(e1 = "numeric",
                    e2 = "moving_average"),
          function(e1, e2) {
            if (length(e1) == 1) {
              moving_average(e1,0) * e2
            } else {
              jasym_filter(e1, e2)
            }
          })
#' @rdname moving_average
#' @export
setMethod("*",
          signature(e2 = "moving_average"),
          function(e1, e2) {
            jasym_filter(e1,e2)
          })
#' @rdname moving_average
#' @export
setMethod("*",
          signature(e1 = "moving_average"),
          function(e1, e2) {
            jasym_filter(e2, e1)
          })
#' @rdname moving_average
#' @export
setMethod("/",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 * moving_average(1/e2,0)
          })
#' @rdname moving_average
#' @export
setMethod("^",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            if (e2 == 0) {
              moving_average(1, 0)
            } else {
              Reduce(`*`, rep(list(e1), e2))
            }
          })
#' @rdname plot_filters
#' @export
plot_coef.moving_average <- function(x, nxlab = 7, add = FALSE,
                                     zeroAsNa = FALSE, ...){
  n <- max(abs(c(upper_bound(x), lower_bound(x))))
  x_plot <- vector(mode = "double", length = 2*n+1)
  names(x_plot) <- coefficients_names(-n, n)
  coefs <- coef(x)
  x_plot[names(coefs)] <- coefs
  if(zeroAsNa)
    x_plot <- trailingZeroAsNa(x_plot)
  matplot(seq(-n, n, by = 1), x_plot,
          xaxt = "n", xlab = "", type = "o", pch = 20,
          ylab = "coefficient", add = add, ...)
  if(!add)
    axis(1, at=seq(-n, n, by = 1), labels = names(x_plot))
}
#' @rdname plot_filters
#' @export
plot_gain.moving_average<- function(x, nxlab = 7, add = FALSE,
                                    xlim = c(0, pi), ...){
  g = get_properties_function(x, "Symmetric Gain")
  plot(g, type = "l",
       xaxt = "n", xlab = "",
       ylab = "gain", add = add, xlim = xlim, ...)
  if(!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}
#' @rdname plot_filters
#' @export
plot_phase.moving_average<- function(x, nxlab = 7, add = FALSE,
                                     xlim = c(0, pi), normalized = FALSE, ...){
  p = get_properties_function(x, "Symmetric Phase")

  # if(normalized){
  #   p <- Vectorize(function(x){
  #     if(x == 0){
  #       p(0)
  #     }else{
  #       p(x)/x
  #     }
  #   })
  # }
  plot(p, type = "l",
       xaxt = "n", xlab = "",
       ylab = "phase", add = add, xlim = xlim, ...)
  if(!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}
#' @export
get_properties_function.moving_average <- function(x,
                                                   component = c("Symmetric Gain",
                                                                 "Symmetric Phase",
                                                                 "Symmetric transfer",
                                                                 "Asymmetric Gain",
                                                                 "Asymmetric Phase",
                                                                 "Asymmetric transfer"), ...){
  x = ma2jd(x)
  component = match.arg(component)
  switch(component,
         "Symmetric Gain" = {
           get_gain_function(x)
         },
         "Asymmetric Gain" = {
           get_gain_function(x)
         },
         "Symmetric Phase" = {
           get_phase_function(x)
         },
         "Asymmetric Phase" = {
           get_phase_function(x)
         },
         "Symmetric transfer" = {
           get_frequencyResponse_function(x)
         },
         "Asymmetric transfer" = {
           get_frequencyResponse_function(x)
         })
}
#'@export
simple_ma <- function(order, lags = - trunc((order-1)/2)) {
  moving_average(rep(1, order), lags = lags) / order
}
#'@export
as.list.moving_average <- function(x, ...) {
  lapply(seq_along(x), function(i) x[i])
}
