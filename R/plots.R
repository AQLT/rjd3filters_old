#' Plots
#'
#' @param x coefficients, gain or phase
#' @param q q
#' @param nxlab number of xlab
#' @param ... other arguments to \code{matplot}
#' @param add boolean indicating if the new plot is added to the previous one
#' @param legend boolean indicating if the legend is printed
#' @param legend.pos position of the legend
#'
#'
#' @examples
#' filter <- filterproperties(3, endpoints = "DAF", kernel = "Henderson")
#' plot_coef(filter, q = c(0,3), legend = TRUE)
#' plot_gain(filter, q = c(0,3), legend = TRUE)
#' plot_phase(filter, q = c(0,3), legend = TRUE)
#' @name plot_filters
#' @rdname plot_filters
#' @importFrom MASS fractions
#' @export
plot_coef <- function(x, q = 0, add = FALSE, legend = FALSE,
                      legend.pos = "topright", ...){
  UseMethod("plot_coef", x)
}
#' @export
plot_coef.JD.Filters <- function(x, q = 0, nxlab = 7, add = FALSE, legend = FALSE,
                                 legend.pos = "topright", ...){
  x = x$filters.coef
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% colnames(x)]
  if (length(col_to_plot) == 0) {
    if(!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "coefficient", xlim=c(1, nrow(x)), ylim=c(0, 1),
           ...)
      axis(1, at=seq_len(nrow(x)), labels = rownames(x))
    }
    return(invisible(0))
  }
  matplot(x[, col_to_plot], type = "l",
          xaxt = "n", xlab = "",
          ylab = "coefficient", add = add, ...)
  if(legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if(!add)
    axis(1, at=seq_len(nrow(x)), labels = rownames(x))
}
#' @name plot_filters
#' @rdname plot_filters
#' @export
plot_gain <- function(x, q = 0, nxlab = 7, add = FALSE, legend = FALSE,
                      legend.pos = "topright", ...){
  UseMethod("plot_gain", x)
}
#' @export
plot_gain.JD.Filters <- function(x, q = 0, nxlab = 7, add = FALSE, legend = FALSE,
                      legend.pos = "topright", ...){
  x = x$filters.gain
  x_values = seq(0, pi, length.out = nrow(x))
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% colnames(x)]
  if (length(col_to_plot) == 0) {
    if(!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "gain", xlim=c(0, 2*pi), ylim=c(0, 1),
           ...)
      x_lab_at <- seq(0, 1, length.out = nxlab)
      axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
    }
    return(invisible(0))
  }
  matplot(x_values,x[, col_to_plot], type = "l",
          xaxt = "n", xlab = "",
          ylab = "gain", add = add)

  if(legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if(!add){
    x_lab_at <- seq(0, 1, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }
}
#' @name plot_filters
#' @rdname plot_filters
#' @export
plot_phase <- function(x, q = 0, nxlab = 7, add = FALSE, legend = FALSE,
                       legend.pos = "topright", ...){
  UseMethod("plot_phase", x)
}
#' @export
plot_phase.JD.Filters <- function(x, q = 0, nxlab = 7, add = FALSE, legend = FALSE,
                       legend.pos = "topright", ...){
  x = x$asymmetricfilter.phase
  x_values = seq(0, pi, length.out = nrow(x))
  col_to_plot <- sprintf("q=%i",q)
  col_to_plot <- col_to_plot[col_to_plot %in% colnames(x)]
  if (length(col_to_plot) == 0) {
    if(!add){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "phase", xlim=c(0, 2*pi), ylim=c(0, 1),
           ...)
      x_lab_at <- seq(0, 1, length.out = nxlab)
      axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
    }
    return(invisible(0))
  }
  matplot(x_values,x[, col_to_plot], type = "l",
          xaxt = "n", xlab = "",
          ylab = "phase", add = add)
  if(legend)
    legend(legend.pos,col_to_plot,
           col = seq_along(col_to_plot), lty=seq_along(col_to_plot), lwd=2)
  if(!add){
    x_lab_at <- seq(0, 1, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }

}

xlabel <- function(x, symbol = "pi"){
  fracs <- strsplit(attr(MASS::fractions(x), "fracs"), "/")  # convert to fractions
  labels <- sapply(fracs, function(i)
    if (length(i) > 1) { paste(i[1], "*", symbol, "/", i[2]) }
    else { paste(i, "*", symbol) })
  labels[1] <- "0"
  labels <- sub("1 * pi", " pi", labels, fixed = TRUE)
  parse(text = labels)
}
