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
    rfilters = rep(list(sfilter), abs(get_upper_bound(sfilter)))
    # rfilters = lapply(seq_along(rfilters), function(i) {
    #   rfilters[[i]] * moving_average(1, lags = -i)
    # })
    lfilters = rep(list(sfilter), abs(get_lower_bound(sfilter)))
    # lfilters = lapply(seq_along(lfilters), function(i) {
    #   lfilters[[i]] * moving_average(1, lags = i)
    # })

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
  all_f = lapply(rev(seq_len(ncol(all_f))), function(i) {
    moving_average(all_f[,i], lags)
  })
  sfilter = all_f[[1]]
  rfilters = all_f[-1]
  lfilters = rev(lapply(rfilters, rev.moving_average))
  res <- new("finite_filters",
             sfilter = sfilter, lfilters = lfilters,
             rfilters = rfilters)
  res
}
jd2finitefilters <- function(jf){
  jsfilter <- .jcall(jf, "Ljdplus/math/linearfilters/SymmetricFilter;", "symmetricFilter")
  jlfilter <- .jcall(jf, "[Ljdplus/math/linearfilters/IFiniteFilter;", "leftEndPointsFilters")
  jrfilter <- .jcall(jf, "[Ljdplus/math/linearfilters/IFiniteFilter;", "rightEndPointsFilters")
  finite_filters(jd2ma(jsfilter),
                 rfilter = lapply(jrfilter, jd2ma),
                 lfilter = rev(lapply(jlfilter, jd2ma)))
}
#' @rdname finite_filters
#' @export
# e1 = finite_filters(moving_average(c(1/3,1/3, 1/3), -1), rfilters = list(moving_average(c(1/2, 1/2), -1)))
# e2 = moving_average(c(1/3,1/3, 1/3), -1)
# e1 = finite_filters(moving_average(rep(1,5), -2)/5,
#                     rfilters = list(moving_average(rep(1,4), -2)/4,
#                                     moving_average(rep(1,3), -2)/3))
setMethod("*",
          signature(e1 = "finite_filters",
                    e2 = "moving_average"),
          function(e1, e2) {
            sym <- e1@sfilter * e2
            rfilters <- lapply(e1@rfilters, `*`, e2)
            # rfilters <- c(rfilters,
            #               lapply(rev(seq_len(get_upper_bound(sym) - length(rfilters))), function(i){
            #                 moving_average(NA, lags = i - 1)
            #               })
            # )
            lfilters <- lapply(e1@lfilters, `*`, e2)
            # lfilters <- c(lapply(rev(seq_len(abs(length(lfilters) + get_lower_bound(sym)))),
            #                      function(i){
            #                 moving_average(NA, lags = -(i - 1))
            #               }),
            #               lfilters
            #               )

            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)
          })
#' @rdname finite_filters
#' @export
setMethod("*",
          signature(e1 = "moving_average",
                    e2 = "finite_filters"),
          function(e1, e2) {
            sym <- e1 * e2@sfilter
            rfilters <- lapply(e2@rfilters, `*`, e1)
            lfilters <- lapply(e2@lfilters, `*`, e1)

            finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)
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
                    c(rev(e2@rfilters), list(e2@sfilter)))
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
            # nb_ma_e1 <- length(e1@lfilters)
            # nb_ma_e2 <- length(e2@lfilters)
            # # nb_ma <- max(nb_ma_e1, nb_ma_e2)
            # # nb_ma <- nb_ma_e1 + nb_ma_e2
            # nb_ma <- e1@sfilter@upper_bound + e2@sfilter@upper_bound

            # new_e1 <- c(e1@lfilters,
            #             rep(list(e1@sfilter), 2 * (nb_ma - nb_ma_e1) + 1),
            #             e1@rfilters)
            # new_e2 <- c(e2@lfilters,
            #             rep(list(e2@sfilter), 2 * (nb_ma - nb_ma_e2)+ 1),
            #             e2@rfilters)
            #
            # matrix_e1 <- t(do.call(cbind, lapply(1:length(new_e1), function(i){
            #   new_e1[[i]] * moving_average(1, lags = -nb_ma+(i-1))
            # })))
            # matrix_e2 <- t(do.call(cbind, lapply(1:length(new_e2), function(i){
            #   new_e2[[i]] * moving_average(1, lags = -nb_ma+(i-1))
            # })))
            # # as.matrix(matrix_e2) |> round(3) |> View()
            # new_mat <- matrix_e1 %*% matrix_e2
            # lags = - nb_ma
            # sym_mat = new_mat[(nrow(new_mat)+1)/2,]
            # sym = moving_average(sym_mat,
            #                      lags = lags)
            # rfilters = new_mat[-(1:((nrow(new_mat)+1)/2)), , drop = FALSE]
            # rfilters = lapply(1:nrow(rfilters),function(i){
            #   moving_average(rfilters[i,-(1:i)],
            #                  lags = lags)
            # })
            #
            # lfilters = new_mat[(1:((nrow(new_mat)-1)/2)), , drop = FALSE]
            # lfilters = lapply(1:nrow(lfilters),function(i){
            #   moving_average(lfilters[i,],
            #                  lags = -lags - i - 1)#why -1 ?
            # })
            #
            # finite_filters(sfilter = sym, rfilters = rfilters, lfilters = lfilters)

            sfilter <- e1@sfilter * e2@sfilter

            n_rfilter <- get_upper_bound(e1@sfilter) + get_upper_bound(e2@sfilter)
            n_lfilter <- get_lower_bound(e1@sfilter) + get_lower_bound(e2@sfilter)
            n_rfilter <- max(n_rfilter, 0)
            n_lfilter <- abs(min(n_lfilter, 0))

            e1_lfilters <- c(e1@lfilters,
                             rep(list(e1@sfilter),
                                 max(n_lfilter + get_lower_bound(e1@sfilter), 0))
            )
            e2_lfilters <- c(e2@lfilters,
                             rep(list(e2@sfilter),
                                 max(n_lfilter + get_lower_bound(e2@sfilter), 0))
            )

            e1_rfilters <- c(rep(list(e1@sfilter),
                                 max(n_rfilter - get_upper_bound(e1@rfilters), 0)),
                             e1@rfilters
            )
            e2_rfilters <- c(rep(list(e2@sfilter),
                                 max(n_rfilter - get_upper_bound(e2@rfilters), 0)),
                             e2@rfilters
            )

            lfilters <- mapply(`*`, e1_lfilters, e2_lfilters)
            rfilters <- mapply(`*`, e1_rfilters, e2_rfilters)
            finite_filters(sfilter = sfilter, rfilters = rfilters, lfilters = lfilters)
          })
#' @export
setMethod("+",
          signature(e1 = "finite_filters",
                    e2 = "finite_filters"),
          function(e1, e2) {

            sfilter <- e1@sfilter + e2@sfilter

            n_rfilter <- get_upper_bound(e1@sfilter) + get_upper_bound(e2@sfilter)
            n_lfilter <- get_lower_bound(e1@sfilter) + get_lower_bound(e2@sfilter)
            n_rfilter <- max(n_rfilter, 0)
            n_lfilter <- abs(min(n_lfilter, 0))

            e1_lfilters <- c(e1@lfilters,
                             rep(list(e1@sfilter),
                                 max(-(get_lower_bound(e2@sfilter) - get_lower_bound(e1@sfilter)), 0))
            )
            e2_lfilters <- c(e2@lfilters,
                             rep(list(e2@sfilter),
                                 max(-(get_lower_bound(e1@sfilter) - get_lower_bound(e2@sfilter)), 0))
            )
            e1_rfilters <- c(rep(list(e1@sfilter),
                                 max(get_upper_bound(e2@sfilter) - get_upper_bound(e1@sfilter), 0)),
                             e1@rfilters
            )
            e2_rfilters <- c(rep(list(e2@sfilter),
                                 max(get_upper_bound(e1@sfilter) - get_upper_bound(e2@sfilter), 0)),
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
    # lapply(s:1, function(i){
    #   new_mm * moving_average(1, lags = (i-1))
    # })
  }))
  x@lfilters = unlist(lapply(x@lfilters, function(x){
    new_mm = to_seasonal(x, s)
    # lapply(1:s, function(i){
    #   new_mm * moving_average(1, lags = -(i-1))
    # })
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
  n_r <- min(get_upper_bound(x@sfilter) - nrfilters, n)
  n_l <- min(abs(get_lower_bound(x@sfilter)) - nlfilters, n)
  new_rfilters <- c(x@rfilters, vector("list", n_r))
  new_lfilters <- c(vector("list", n_l), x@lfilters)
  for (i in seq_len(n_r)) {
    new_rfilters[[nrfilters + i]] <-
      new_rfilters[[nrfilters + i - nperiod]] *
      moving_average(1, lags = -nperiod)
  }
  for (i in rev(seq_len(n_l))) {
    new_lfilters[[i]] <-
      new_lfilters[[n_l - i + nperiod]] *
      moving_average(1, lags = nperiod)
  }
  finite_filters(x@sfilter, rfilters = new_rfilters, lfilters = new_lfilters)
}
