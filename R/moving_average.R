# library(rJava)
# library(rjdfilters)
setClass("moving_average",
         representation = representation(internal = "jobjRef" )
)
#' Manipulation of moving averages
#'
#' @param x vector of coefficients
#' @param lags integer indicating the number of lags of the moving average.
#' @param s seasonal period for the \code{to_seasonal()} function.
#' @param ...,na.rm see \link[base]{mean}
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
moving_average <- function(x, lags = -length(x)){
  x <- removeTrailingZeroOrNA(as.numeric(x))
  if (length(x) == 1){
    x <- .jarray(x)
  }
  jobj <-  .jcall("jdplus/math/linearfilters/FiniteFilter",
                  "Ljdplus/math/linearfilters/FiniteFilter;",
                  "of", x,
                  as.integer(lags))
  new(Class = "moving_average", internal = jobj)
}
#' @rdname moving_average
#' @export
is.moving_average <- function(x){
  inherits(x, "moving_average")
}
#' @importFrom stats coef
#' @export
coef.moving_average <- function(object, ...){
  coefs <- .jcall(object@internal, "[D", "weightsToArray")
  names(coefs) <- coefficients_names(get_lower_bound(object),
                                  get_upper_bound(object))
  coefs
}
#' @rdname moving_average
#' @export
is_symmetric <- function(x){
  .jcall(x@internal, "Z", "isSymmetric")
}
#' @rdname moving_average
#' @export
get_upper_bound <- function(x){
  .jcall(x@internal, "I", "getUpperBound")
}
#' @rdname moving_average
#' @export
get_lower_bound <- function(x){
  .jcall(x@internal, "I", "getLowerBound")
}
#' @rdname moving_average
#' @export
mirror <- function(x){
  .jcall(x@internal, "Ljdplus/math/linearfilters/FiniteFilter;", "mirror")
}
#' @rdname moving_average
#' @export
length.moving_average <- function(x){
  .jcall(x@internal, "I", "length")
}
#' @rdname moving_average
#' @export
to_seasonal <- function(x, s){
  lb <- get_lower_bound(x)
  up <- get_upper_bound(x)
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
setMethod(f = "show",
          signature = "moving_average",
          definition = function(object){
            print(.jcall(object@internal, "S", "toString"))
            invisible(object)
          })
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
setMethod("+",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.math.linearfilters.FiniteFilter")
            # jobj <- finiteFilter$add(e1@internal, e2@internal)
            jobj <- .jcall(finiteFilter,
                   "Ljdplus/math/linearfilters/FiniteFilter;",
                   "add",
                   .jcast(e1@internal, "jdplus/math/linearfilters/IFiniteFilter"),
                   .jcast(e2@internal, "jdplus/math/linearfilters/IFiniteFilter"))

            new(Class = "moving_average", internal = jobj)
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
            moving_average(e1,0) + e2
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
            # jobj <- finiteFilter$negate(e1@internal)
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/math/linearfilters/FiniteFilter;",
                           "negate",
                           .jcast(e1@internal, "jdplus/math/linearfilters/IFiniteFilter"))
            new(Class = "moving_average", internal = jobj)
          })
#' @rdname moving_average
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.math.linearfilters.FiniteFilter")
            # jobj <- finiteFilter$subtract(e1@internal, e2@internal)
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/math/linearfilters/FiniteFilter;",
                           "subtract",
                           .jcast(e1@internal, "jdplus/math/linearfilters/IFiniteFilter"),
                           .jcast(e2@internal, "jdplus/math/linearfilters/IFiniteFilter"))
            new(Class = "moving_average", internal = jobj)
          })
#' @rdname moving_average
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 - moving_average(e2,0)
          })
#' @rdname moving_average
#' @export
setMethod("-",
          signature(e1 = "numeric",
                    e2 = "moving_average"),
          function(e1, e2) {
            moving_average(e1,0) - e2
          })
#' @rdname moving_average
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "moving_average"),
          function(e1, e2) {
            finiteFilter <- J("jdplus.math.linearfilters.FiniteFilter")
            # jobj <- finiteFilter$multiply(e1@internal, e2@internal)
            jobj <- .jcall(finiteFilter,
                           "Ljdplus/math/linearfilters/FiniteFilter;",
                           "multiply",
                           .jcast(e1@internal, "jdplus/math/linearfilters/IFiniteFilter"),
                           .jcast(e2@internal, "jdplus/math/linearfilters/IFiniteFilter"))
            new(Class = "moving_average", internal = jobj)
          })
#' @rdname moving_average
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 * moving_average(e2,0)
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
            jasym_filter(e1,e2)
          })
#' @rdname moving_average
#' @export
setMethod("*",
          signature(e1 = "numeric",
                    e2 = "moving_average"),
          function(e1, e2) {
            moving_average(e1,0) * e2
          })
#' @rdname moving_average
#' @export
setMethod("/",
          signature(e1 = "moving_average",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 * moving_average(1/e2,0)
          })
#' @rdname plot_filters
#' @export
plot_coef.moving_average <- function(x, nxlab = 7, add = FALSE,
                                 zeroAsNa = FALSE, ...){
  n <- max(abs(c(get_upper_bound(x), get_lower_bound(x))))
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
get_properties_function.moving_average <-
  function(x,
           component = c("Symmetric Gain",
                         "Symmetric Phase",
                         "Symmetric transfer",
                         "Asymmetric Gain",
                         "Asymmetric Phase",
                         "Asymmetric transfer"), ...){
    component = match.arg(component)
  switch(component,
         "Symmetric Gain" = {
           get_gain_function(x@internal)
         },
         "Asymmetric Gain" = {
           get_gain_function(x@internal)
         },
         "Symmetric Phase" = {
           get_phase_function(x@internal)
         },
         "Asymmetric Phase" = {
           get_phase_function(x@internal)
         },
         "Symmetric transfer" = {
           get_frequencyResponse_function(x@internal)
         },
         "Asymmetric transfer" = {
           get_frequencyResponse_function(x@internal)
         })
  }

# sum(e1)
# sum(e2)
#
# y <- retailsa$DrinkingPlaces
# library(rJava)
# e1 <- moving_average(rep(1,12), lags =-5)
# e1 <- e1/sum(e1)
# length(e1)
# e1
# e2 <- moving_average(rep(1,12), lags =-6)/12
# x <- (e1+e2)/12
# x
# e1 <- moving_average(1:5, lags =-3)
# x = e1
#
# jy <- J("demetra/data/DoubleSeq")$of(as.numeric(y))
# jres <- J("demetra/data/DoubleSeq$Mutable")$of(as.numeric(y))
# tmp$of
# tmp$of()
# x@internal$apply(jy,jres)
#
#
# lb = get_lower_bound(x)
# ub = get_upper_bound(x)
# DataBlock = J("jdplus.data.DataBlock")
# out = DataBlock$of(as.numeric(rep(NA, jy$length() - length(x)+1)))
# out$length()
# x@internal$apply(DataBlock$of(as.numeric(y)),
#                  out)
# out$length()
# result = out$toArray()
# result <- c(rep(NA, abs(min(lb, 0))),
#             result,
#          rep(NA, abs(max(ub, 0))))
# if(is.ts(y))
#   result <- ts(result,start = start(y), frequency = frequency(y))
#
# ts.union(jasym_filter(y, coef(x), abs(get_lower_bound(x))),
#          result)
# all.equal(result, jasym_filter(y, coef(x), abs(get_lower_bound(x))))

# e1+e2
# e1-e2
# e1*e2
# e2/2
