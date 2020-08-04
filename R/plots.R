#' Plots
#'
#' @param x coefficients, gain or phase
#' @param q q
#' @param nxlab number of xlab
#' @param ... other arguments to \code{matplot}
#' @param add boolean indicating if the new plot is added to the previous one.
#' @param xlim vector containing x limits.
#' @param legend boolean indicating if the legend is printed.
#' @param legend.pos position of the legend.
#' @param normalized boolean indicatif if the phase function is normalized by the frequency.
#'
#' @examples
#' filter <- filterproperties(6, endpoints = "DAF", kernel = "Henderson")
#' plot_coef(filter, q = c(0,3), legend = TRUE)
#' plot_gain(filter, q = c(0,3), legend = TRUE)
#' plot_phase(filter, q = c(0,3), legend = TRUE)
#' @name plot_filters
#' @rdname plot_filters
#' @importFrom MASS fractions
#' @export
plot_coef <- function(x, nxlab = 7, add = FALSE,
                      zeroAsNa = FALSE, ...){
  UseMethod("plot_coef", x)
}
#' @export
plot_coef.lp_filter<- function(x, q = 0, nxlab = 7, add = FALSE, legend = FALSE,
                                 legend.pos = "topright",
                               zeroAsNa = FALSE, ...){
  x = x$filters.coef
  if(zeroAsNa)
    x  <- apply(x,2, trailingZeroAsNa)
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% colnames(x)]
  horizon <- (nrow(x)-1)/2
  if (length(col_to_plot) == 0) {
    if(!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "coefficient", xlim=c(-horizon, horizon), ylim=c(0, 1),
           ...)
      axis(1, at=seq(-horizon, horizon, by = 1), labels = rownames(x))
    }
    return(invisible(0))
  }
  matplot(seq(-horizon, horizon, by = 1),x[,col_to_plot],
          xaxt = "n", xlab = "", type = "o", pch = 20,
          ylab = "coefficient", add = add, ...)
  if(legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if(!add)
    axis(1, at=seq(-horizon, horizon, by = 1), labels = rownames(x))
}
#' @export
plot_coef.fst_filter <- function(x, nxlab = 7, add = FALSE,
                                 zeroAsNa = FALSE, ...){
  n <- max(abs(c(x$internal$ub(), x$internal$lb())))
  x_plot <- vector(mode = "double", length = 2*n+1)
  names(x_plot) <- coefficients_names(-n, n)
  x_plot[names(x$filter)] <- x$filter
  if(zeroAsNa)
    x_plot <- trailingZeroAsNa(x_plot)
  matplot(seq(-n, n, by = 1), x_plot,
          xaxt = "n", xlab = "", type = "o", pch = 20,
          ylab = "coefficient", add = add, ...)
  if(!add)
    axis(1, at=seq(-n, n, by = 1), labels = names(x_plot))
}
#' @name plot_filters
#' @rdname plot_filters
#' @export
plot_gain <- function(x, nxlab = 7, add = FALSE,
                      xlim = c(0, pi), ...){
  UseMethod("plot_gain", x)
}
#' @export
plot_gain.lp_filter<- function(x, q = 0, nxlab = 7, add = FALSE, legend = FALSE,
                      legend.pos = "topright", xlim = c(0, pi), ...){
  x = x$filters.gain
  x_values = seq(0, pi, length.out = nrow(x))
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% colnames(x)]
  if (length(col_to_plot) == 0) {
    if(!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "gain", xlim=xlim, ylim=c(0, 1),
           ...)
      x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
      axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
    }
    return(invisible(0))
  }
  matplot(x_values,x[, col_to_plot], type = "l",
          xaxt = "n", xlab = "",
          ylab = "gain", add = add, xlim = xlim, ...)

  if(legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if(!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}
#' @export
plot_gain.fst_filter<- function(x, nxlab = 7, add = FALSE,
                                xlim = c(0, pi), ...){
  x = x$gain
  x_values = seq(0, pi, length.out = length(x))

  matplot(x_values, x, type = "l",
          xaxt = "n", xlab = "",
          ylab = "gain", add = add, xlim = xlim, ...)
  if(!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}
#' @name plot_filters
#' @rdname plot_filters
#' @export
plot_phase <- function(x, nxlab = 7, add = FALSE,
                       xlim = c(0, pi), normalized = TRUE, ...){
  UseMethod("plot_phase", x)
}
#' @export
plot_phase.lp_filter<- function(x, q = 0, nxlab = 7, add = FALSE, legend = FALSE,
                       legend.pos = "topright", xlim = c(0, pi), normalized = TRUE, ...){
  x = x$filters.phase
  x_values = seq(0, pi, length.out = nrow(x))
  if(normalized){
    x[-1,] <- x[-1,] / x_values[-1]
  }
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% colnames(x)]
  if (length(col_to_plot) == 0) {
    if(!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "phase", xlim = xlim, ylim=c(0, 1),
           ...)
      x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
      axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
    }
    return(invisible(0))
  }
  matplot(x_values,x[, col_to_plot], type = "l",
          xaxt = "n", xlab = "",
          ylab = "phase", add = add, xlim = xlim, ...)
  if(legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if(!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }

}
#' @export
plot_phase.fst_filter<- function(x, nxlab = 7, add = FALSE,
                                 xlim = c(0, pi), normalized = TRUE, ...){
  x = x$phase
  x_values = seq(0, pi, length.out = length(x))
  if(normalized){
    x[-1] <- x[-1] / x_values[-1]
  }
  matplot(x_values, x, type = "l",
          xaxt = "n", xlab = "",
          ylab = "phase", add = add, xlim = xlim, ...)
  if(!add){
    x_lab_at <- seq(xlim[1]/pi, xlim[2]/pi, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}
xlabel <- function(x, symbol = "pi"){
  fracs <- strsplit(attr(MASS::fractions(x), "fracs"), "/")  # convert to fractions
  labels <- sapply(fracs, function(i)
    if (length(i) > 1) { paste(i[1], "*", symbol, "/", i[2]) }
    else { paste(i, "*", symbol) })
  labels <- sub("0 * pi", "0", labels, fixed = TRUE)
  labels <- sub("1 * pi", " pi", labels, fixed = TRUE)
  parse(text = labels)
}

trailingZeroAsNa <- function(x){
  if(x[length(x)]==0)
    x [seq(from = tail(which(!sapply(x, function(y) isTRUE(all.equal(y,0)))),1)+1,
           to = length(x),
           by = 1)] <- NA
  x
}
