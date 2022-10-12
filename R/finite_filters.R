library(rjdfilters)
setClass("finite_filters",
         slots = c(sfilter = "moving_average",
                   lfilters = "list",
                   rfilters = "list")
)
finite_filters <- function(sfilter,
                           rfilters = NULL,
                           lfilters = NULL){
  if (is.null(lfilters) & !is.null(rfilters)) {
    lfilters = rev(lapply(rfilters, rev.moving_average))
  } else if (!is.null(lfilters) & is.null(rfilters)) {
    rfilters = rev(lapply(lfilters, rev.moving_average))
  } else if (is.null(lfilters) & is.null(rfilters)) {
    return(sfilter)
  }
  res <- new("finite_filters",
             sfilter = sfilter, lfilters = lfilters,
             rfilters = rfilters)
  res
}
setMethod(f = "show",
          signature = "finite_filters",
          definition = function(object){
            x = do.call(cbind, c(list(object@sfilter), object@rfilters))
            # if (isTRUE(all.equal(object@rfilters, object@lfilters))) {
            #   x = do.call(cbind, c(list(object@sfilter), object@rfilters))
            # } else {
            #   x = do.call(cbind, c(object@lfilters, list(object@sfilter), object@rfilters))
            # }
            print(x)
          })

#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "moving_average"),
          function(e1, e2) {
            e1@sfilter = e1@sfilter * e2
            e1@lfilters = lapply(e1@lfilters, `*`, e2)
            e1@rfilters = lapply(e1@rfilters, `*`, e2)
            e1
          })
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 * moving_average(e2,0)
          })
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "ANY",
                    e2 = "finite_filters"),
          function(e1, e2) {
            jfilter(e1,
                    c(object@rfilters, list(object@sfilter)))
          })

#' @rdname finite_filters
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
#' @rdname finite_filters
#' @export
setMethod("+",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 + moving_average(e2,0)
          })
#' @rdname finite_filters
#' @export
setMethod("+",
          signature(e1 = "numeric",
                    e2 = "finite_filters"),
          function(e1, e2) {
            e2 + e1
          })
#' @rdname finite_filters
#' @export
setMethod("+", signature(e1 = "finite_filters", e2 = "missing"), function(e1,e2) e1)
#' @rdname finite_filters
#' @export
setMethod("-",
          signature(e1 = "finite_filters",
                    e2 = "missing"),
          function(e1, e2) {
            e1@sfilter = - e1@sfilter
            e1@lfilters = lapply(e1@lfilters, `-`)
            e1@rfilters = lapply(e1@rfilters, `-`)
            e1
          })
#' @rdname finite_filters
#' @export
setMethod("-",
          signature(e1 = "finite_filters",
                    e2 = "moving_average"),
          function(e1, e2) {
            e1@sfilter = e1@sfilter - e2
            e1@lfilters = lapply(e1@lfilters, `-`, e2)
            e1@rfilters = lapply(e1@rfilters, `-`, e2)
            e1
          })
#' @rdname finite_filters
#' @export
setMethod("-",
          signature(e1 = "moving_average",
                    e2 = "finite_filters"),
          function(e1, e2) {
            e1 + (- e2)
          })
#' @rdname finite_filters
#' @export
setMethod("-",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            e1 - moving_average(e2,0)
          })
#' @rdname finite_filters
#' @export
setMethod("-",
          signature(e1 = "numeric",
                    e2 = "finite_filters"),
          function(e1, e2) {
            moving_average(e1,0) - e2
          })

#' @rdname finite_filters
#' @export
setMethod("/",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            e1@sfilter = e1@sfilter / e2
            e1@lfilters = lapply(e1@lfilters, `/`, e2)
            e1@rfilters = lapply(e1@rfilters, `/`, e2)
            e1
          })

#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "finite_filters"),
          function(e1, e2) {
            nb_ma_e1 <- length(e1@lfilters)
            nb_ma_e2 <- length(e2@lfilters)
            nb_ma <- max(nb_ma_e1, nb_ma_e2)

            new_e1 <- c(e1@lfilters,
                       rep(list(e1@sfilter), 2 * (nb_ma - nb_ma_e1) + 1),
                       e1@rfilters)
            new_e2 <- c(e2@lfilters,
                        rep(list(e2@sfilter), 2 * (nb_ma - nb_ma_e2)+ 1),
                        e2@rfilters)

            if (nb_ma_e1 < nb_ma_e2) {
              matrix_e2 <- t(do.call(cbind, lapply(1:length(new_e2), function(i){
                new_e2[[i]] * moving_average(1, lags = -nb_ma+(i-1))
              })))

              matrix_e1 <- t(do.call(cbind, lapply(1:length(new_e1), function(i){
                new_e1[[i]] * moving_average(1, lags = -abs(nb_ma_e1 - nb_ma_e2)*2 +i)
              })))
            } else if (nb_ma_e1 > nb_ma_e2) {
              matrix_e1 <- t(do.call(cbind, lapply(1:length(new_e1), function(i){
                new_e1[[i]] * moving_average(1, lags = -nb_ma+(i-1))
              })))

              matrix_e2 <- t(do.call(cbind, lapply(1:length(new_e2), function(i){
                new_e2[[i]] * moving_average(1, lags = -abs(nb_ma_e1 - nb_ma_e2)*2 +i-1)
              })))
            } else {
              matrix_e1 <- t(do.call(cbind, lapply(1:length(new_e1), function(i){
                new_e1[[i]] * moving_average(1, lags = -nb_ma+(i-1))
              })))
              matrix_e2 <- t(do.call(cbind, lapply(1:length(new_e2), function(i){
                new_e2[[i]] * moving_average(1, lags = -nb_ma+(i-1))
              })))
            }

            new_mat <- matrix_e1 %*% matrix_e2
            lags = -(ncol(new_mat)-1)/2
            sym = moving_average(new_mat[(nrow(new_mat)+1)/2,],
                                 lags = lags)
            rfilters = new_mat[-(1:((nrow(new_mat)+1)/2)),]
            rfilters = lapply(1:nrow(rfilters),function(i){
              moving_average(rfilters[i,],
                             lags = lags - i +1)
            })
            # rfilters = rev(rfilters)

            lfilters = new_mat[(1:((nrow(new_mat)-1)/2)),]
            lfilters = lapply(1:nrow(lfilters),function(i){
              moving_average(lfilters[i,],
                             lags = lags + i -1)
            })
            # lfilters = rev(lfilters)
            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)
          })
all_mm <- apply(lp_filter()$filters.coef,2,
                moving_average, lags = -6)

ff = finite_filters(all_mm[["q=6"]], rfilters = rev(all_mm[-7]))


ff*2/10
all_mm <- apply(lp_filter(horizon = 2)$filters.coef,2,
                moving_average, lags = -2)

# ff = finite_filters(all_mm[["q=2"]], rfilters = all_mm[-3])
ff2 = finite_filters(all_mm[["q=2"]], rfilters = all_mm[-3])
ff = finite_filters(moving_average(c(1/3,1/3,1/3), -1),
                    rfilters = list(moving_average(c(1/2,1/2), -1)))

e1 = finite_filters(all_mm[["q=2"]], rfilters = rev(all_mm[-3]))
e2 = finite_filters(moving_average(c(1/3,1/3,1/3), -1),
                    rfilters = list(moving_average(c(1/2,1/2), -1)))
ff2 * ff
class(ff2)
class(ff)
tmp = print(ff)
t(tmp)
object = ff
e1 = ff2
e2 = ff
rev()
