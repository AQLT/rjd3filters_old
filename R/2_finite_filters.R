setClass("finite_filters",
         slots = c(sfilter = "moving_average",
                   lfilters = "list",
                   rfilters = "list")
)
#' Manipulating Finite Filters
#'
#' @param sfilter the symmetric filter or a `"FiniteFilter` object.
#' @param rfilters the right filters (used on the last points).
#' @param lfilters the left filters (used on the first points).
#' @param first_to_last boolean indicating if the first element of `rfilters` is the
#' first asymmetric filter (when only one observation is missing) or the last one (real-time estimates).
#'
#' @examples
#' ff_lp = finite_filters(lp_filter())
#' ff_simple_ma = finite_filters(moving_average(c(1, 1, 1), lags = -1)/3,
#'                rfilters = list(moving_average(c(1, 1), lags = -1)/2))
#' ff_lp
#' ff_simple_ma
#' ff_lp * ff_simple_ma
#'
#' @export
finite_filters <- function(sfilter,
                           rfilters = NULL,
                           lfilters = NULL,
                           first_to_last = FALSE){
  UseMethod("finite_filters", sfilter)
}
#' @export
finite_filters.moving_average <- function(sfilter,
                                          rfilters = NULL,
                                          lfilters = NULL,
                                          first_to_last = FALSE){
  if (is.null(lfilters) & !is.null(rfilters)) {
    if (first_to_last) {
      rfilters = rev(rfilters)
    }
    lfilters = rev(lapply(rfilters, rev.moving_average))
  } else if (!is.null(lfilters) & is.null(rfilters)) {
    if (!first_to_last) {
      lfilters = rev(lfilters)
    }
    rfilters = rev(lapply(lfilters, rev.moving_average))
  } else if (is.null(lfilters) & is.null(rfilters)) {
    rfilters = lfilters = list()

  }
  res <- new("finite_filters",
             sfilter = sfilter, lfilters = lfilters,
             rfilters = rfilters)
  res
}
#' @export
finite_filters.FiniteFilters <- function(sfilter,
                                         rfilters = NULL,
                                         lfilters = NULL,
                                         first_to_last = FALSE){
  all_f = sfilter$filters.coef
  lags = -(nrow(all_f) - 1)/2
  all_f = lapply(seq_len(ncol(all_f)), function(i) {
    moving_average(all_f[,i], lags, trailing_zero = TRUE)
  })
  sfilter = all_f[[1]]
  rfilters = all_f[-1]
  res <- new("finite_filters",
             sfilter = sfilter,
             rfilters = rfilters)
  res
}
#' @export
finite_filters.list <- function(sfilter,
                                rfilters = NULL,
                                lfilters = NULL,
                                first_to_last = FALSE){
  lags <- length(sfilter)-1

  all_f <- lapply(sfilter,
                   function(x){
                     moving_average(rm_trailing_zero_or_na(x), -lags)
                   }
  )
  if (first_to_last)
    all_f <- rev(all_f)
  sfilter <- all_f[[1]]
  rfilters <- all_f[-1]
  finite_filters(sfilter = sfilter, rfilters = rfilters)
}
#' @export
finite_filters.matrix <- function(sfilter,
                                  rfilters = NULL,
                                  lfilters = NULL,
                                  first_to_last = FALSE){
  if (first_to_last) {
    index <- ncol(sfilter):1
  } else {
    index <- 1:ncol(sfilter)
  }
  coefs <- lapply(index, function(i) sfilter[,i])
  finite_filters(coefs, first_to_last = first_to_last)
}
.jd2finitefilters <- function(jf){
  jsfilter <- .jcall(jf, "Ljdplus/math/linearfilters/SymmetricFilter;", "symmetricFilter")
  jlfilter <- .jcall(jf, "[Ljdplus/math/linearfilters/IFiniteFilter;", "leftEndPointsFilters")
  jrfilter <- .jcall(jf, "[Ljdplus/math/linearfilters/IFiniteFilter;", "rightEndPointsFilters")
  finite_filters(.jd2ma(jsfilter),
                 rfilter = lapply(jrfilter, .jd2ma),
                 lfilter = rev(lapply(jlfilter, .jd2ma)))
}
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "moving_average"),
          function(e1, e2) {
            sym <- e1@sfilter * e2
            rfilters <- lapply(e1@rfilters, `*`, e2)
            lfilters <- lapply(e1@lfilters, `*`, e2)

            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)
          })
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "finite_filters"),
          function(e1, e2) {

            new_ub <- e1@upper_bound + e2@sfilter@upper_bound
            new_lb_sym <- e1@lower_bound + e2@sfilter@lower_bound
            new_lb <- e1@lower_bound + -length(e2@rfilters)

            new_e2 <- c(e2@lfilters,
                        rep(list(e2@sfilter), length(e1)),
                        e2@rfilters)
            new_e1 <- rep(list(e1), length(new_e2))
            new_e2 <- lapply(1:length(new_e2), function(i){
              new_e2[[i]] * moving_average(1,
                                           lags = (new_lb + (i - 1) + e1@lower_bound * (new_lb == new_lb_sym)))
            })
            new_e1 <- lapply(1:length(new_e1), function(i){
              new_e1[[i]] * moving_average(1, lags = (new_lb + (i - 1)))
            })
            all_f <- t(do.call(cbind,c(new_e1, new_e2)))
            mat_e1 <- all_f[seq_along(new_e1),]
            mat_e2 <- all_f[-seq_along(new_e1),]
            new_mat = (mat_e1[, seq_along(new_e2)] %*% mat_e2)[seq_len(1 + length(e2@lfilters) + length(e2@rfilters)),]

            max_lags <- min(sapply(new_e1, lower_bound), sapply(new_e2, lower_bound))

            # i_to_remove = seq_len(-(max_lags - new_lb))
            sym_mat <- new_mat[(nrow(new_mat)+1)/2, ]
            sym <- moving_average(sym_mat,
                                  lags = new_lb_sym, leading_zero = FALSE, trailing_zero = TRUE)
            rfilters <- new_mat[-(1:((nrow(new_mat)+1)/2)),, drop = FALSE]
            rfilters <- lapply(1:nrow(rfilters),function(i){
              moving_average(rfilters[i,-seq_len(i)],
                             lags = new_lb_sym, leading_zero = FALSE, trailing_zero = TRUE)
            })

            lfilters <- new_mat[(1:((nrow(new_mat)-1)/2)), , drop = FALSE]
            lfilters <- lapply(1:nrow(lfilters),function(i){
              moving_average(lfilters[i,],
                             lags = new_lb_sym + (nrow(lfilters) - i) + 1, leading_zero = FALSE, trailing_zero = TRUE)#why -1 ?
            })
            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)

          })
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            if (length(e2) == 1) {
              e1 * moving_average(e2,0)
            } else {
              jfilter(e2, e1)
            }
          })
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "ANY",
                    e2 = "finite_filters"),
          function(e1, e2) {
            jfilter(e1, e2)
          })
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "ANY"),
          function(e1, e2) {
            e2 * e1
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
setMethod("+",
          signature(e1 = "finite_filters",
                    e2 = "moving_average"),
          function(e1, e2) {
            e1@sfilter = e1@sfilter + e2
            e1@lfilters = lapply(e1@lfilters, `+`, e2)
            e1@rfilters = lapply(e1@rfilters, `+`, e2)
            e1
          })
#' @rdname finite_filters
#' @export
setMethod("+",
          signature(e1 = "moving_average",
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
#' @method as.matrix finite_filters
#' @export
as.matrix.finite_filters <- function(x, rfilters = TRUE, lfilters = FALSE, ...) {
  rfilters_s <- lfilters_s <-
    index_r <- index_l <- NULL
  index_sym = length(x@rfilters)
  if (lfilters) {
    lfilters_s <- x@lfilters
    index_l = seq(0, -(length(x@lfilters) - 1))
  }
  if (rfilters) {
    rfilters_s <- x@rfilters
    index_r = seq(length(x@rfilters) - 1, 0)
  }
  mat = do.call(cbind, c(lfilters_s, list(x@sfilter), rfilters_s))
  colnames(mat) <- sprintf("q=%i", c(index_l, index_sym, index_r))
  mat
}
#' @export
setMethod("^",
          signature(e1 = "finite_filters",
                    e2 = "numeric"),
          function(e1, e2) {
            Reduce(`*`, rep(list(e1), e2))
          })
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "finite_filters"),
          function(e1, e2) {
            new_ub <- length(e1@rfilters) + length(e2@rfilters)
            new_lb <- length(e1@lfilters) + length(e2@lfilters)
            new_lb_sym <- e1@sfilter@lower_bound + e2@sfilter@lower_bound

            new_e1 <- c(e1@lfilters,
                        rep(list(e1@sfilter), 1 + (new_ub - length(e1@rfilters)) + (new_lb - length(e1@lfilters))),
                        e1@rfilters)
            new_e2 <- c(e2@lfilters,
                        rep(list(e2@sfilter), 1 + (new_ub - length(e2@rfilters)) + (new_lb - length(e2@lfilters))),
                        e2@rfilters)

            new_e1 <- lapply(1:length(new_e1), function(i){
              new_e1[[i]] * moving_average(1, lags = (-new_lb + (i - 1)))
            })
            new_e2 <- lapply(1:length(new_e2), function(i){
              new_e2[[i]] * moving_average(1, lags = (-new_lb  + (i - 1)))
            })
            all_f <- t(do.call(cbind,c(new_e1, new_e2)))

            mat_e1 <- all_f[seq_along(new_e1),]
            mat_e2 <- all_f[-seq_along(new_e1),]
            new_mat <- (mat_e1[, seq_along(new_e2)] %*% mat_e2)

            max_lags <- min(sapply(new_e1, lower_bound), sapply(new_e2, lower_bound))

            sym_mat <- new_mat[(nrow(new_mat)+1)/2,]
            sym <- moving_average(sym_mat,
                                  lags = new_lb_sym, leading_zero = TRUE, trailing_zero = TRUE)
            rfilters <- new_mat[-(1:((nrow(new_mat)+1)/2)), , drop = FALSE]
            rfilters <- lapply(1:nrow(rfilters),function(i){
              moving_average(rfilters[i,],
                             lags = new_lb_sym - i, leading_zero = TRUE, trailing_zero = TRUE)
            })
            #
            lfilters <- new_mat[(1:((nrow(new_mat)-1)/2)), , drop = FALSE]
            lfilters <- lapply(1:nrow(lfilters),function(i){
              moving_average(lfilters[i,],
                             lags = new_lb_sym + (nrow(lfilters) - i) + 1, leading_zero = TRUE, trailing_zero = TRUE)#why -1 ?
            })
            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)
          })
#' @export
setMethod("+",
          signature(e1 = "finite_filters",
                    e2 = "finite_filters"),
          function(e1, e2) {

            sfilter <- e1@sfilter + e2@sfilter

            n_rfilter <- upper_bound(e1@sfilter) + upper_bound(e2@sfilter)
            n_lfilter <- lower_bound(e1@sfilter) + lower_bound(e2@sfilter)
            n_rfilter <- max(n_rfilter, 0)
            n_lfilter <- abs(min(n_lfilter, 0))

            e1_lfilters <- c(e1@lfilters,
                             rep(list(e1@sfilter),
                                 max(-(lower_bound(e2@sfilter) - lower_bound(e1@sfilter)), 0))
            )
            e2_lfilters <- c(e2@lfilters,
                             rep(list(e2@sfilter),
                                 max(-(lower_bound(e1@sfilter) - lower_bound(e2@sfilter)), 0))
            )
            e1_rfilters <- c(rep(list(e1@sfilter),
                                 max(upper_bound(e2@sfilter) - upper_bound(e1@sfilter), 0)),
                             e1@rfilters
            )
            e2_rfilters <- c(rep(list(e2@sfilter),
                                 max(upper_bound(e1@sfilter) - upper_bound(e2@sfilter), 0)),
                             e2@rfilters
            )
            e1_lfilters_f <- c(e1_lfilters, rep(list(0),
                                                max(length(e2_lfilters) - length(e1_lfilters), 0)))
            e2_lfilters_f <- c(e2_lfilters, rep(list(0),
                                                max(length(e1_lfilters) - length(e2_lfilters), 0)))
            e1_rfilters_f <- c(e1_rfilters, rep(list(0),
                                                max(length(e2_rfilters) - length(e1_rfilters), 0)))
            e2_rfilters_f <- c(e2_rfilters, rep(list(0),
                                                max(length(e1_rfilters) - length(e2_rfilters), 0)))

            lfilters <- mapply(`+`, e1_lfilters_f, e2_lfilters_f)
            rfilters <- mapply(`+`, e1_rfilters_f, e2_rfilters_f)
            finite_filters(sfilter = sfilter, rfilters = rfilters, lfilters = lfilters)
          })
#' @export
setMethod("-",
          signature(e1 = "finite_filters",
                    e2 = "finite_filters"),
          function(e1, e2) {
            e1 + (-e2)
          })
#' @export
to_seasonal.finite_filters <- function(x, s){
  x@sfilter = to_seasonal(x@sfilter, s)
  x@rfilters = unlist(lapply(x@rfilters, function(x){
    new_mm = to_seasonal(x, s)
    rep(list(new_mm), s)
  }))
  x@lfilters = unlist(lapply(x@lfilters, function(x){
    new_mm = to_seasonal(x, s)
    rep(list(new_mm), s)
  }))
  x
}
#'@export
imput_last_obs <- function(x, n, nperiod = 1) {
  nrfilters <- length(x@rfilters)
  nlfilters <- length(x@lfilters)
  if (missing(n))
    n <- max(nrfilters, nlfilters)
  n_r <- min(upper_bound(x@sfilter) - nrfilters, n)
  n_l <- min(abs(lower_bound(x@sfilter)) - nlfilters, n)
  new_rfilters <- c(x@rfilters, vector("list", n_r))
  new_lfilters <- c(vector("list", n_l), x@lfilters)
  for (i in seq_len(n_r)) {
    new_rfilters[[nrfilters + i]] <-
      new_rfilters[[nrfilters + i - nperiod]] *
      moving_average(1, lags = -nperiod)
  }
  for (i in rev(seq_len(n_l))) {
    new_lfilters[[i]] <-
      new_lfilters[[i + nperiod]] *
      moving_average(1, lags = nperiod)
  }
  finite_filters(x@sfilter, rfilters = new_rfilters, lfilters = new_lfilters)
}
