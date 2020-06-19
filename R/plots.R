#' Plots
#'
#' @param x coefficients, gain or phase
#' @param q q
#' @param nxlab number of xlab
#' @param ... other arguments to \code{plot}
#' @param add add to previous plot
#'
#' @return
#'
#' @examples
#' @name plot_filters
#' @rdname plot_filters
#' @importFrom MASS fractions
#' @export
plot_coef <- function(x, q = 0, add = FALSE,...){
  x_values = seq(0, nrow(x)) - (nrow(x) - 1)/2
  y = x[, sprintf("q=%i",q)]
  if(!add){
    plot(y, type = "l",
         xaxt = "n", xlab = "",
         ylab = "coefficient", ...)
    axis(1, at=seq_len(nrow(x)), labels = rownames(x))
  }else{
    lines(y, ...)
  }
}
#' @name plot_filters
#' @rdname plot_filters
#' @export
plot_gain <- function(x, q = 0, nxlab = 6, add = FALSE,...){
  x_values = seq(0, 2 * pi, length.out = nrow(x))
  y = x[,sprintf("q=%i", q)]
  if(!add){
    plot(x_values, y, type = "l",
         xaxt = "n", xlab = "",
         ylab = "gain", ...)
    x_lab_at <- seq(0, 2, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }else{
    lines(x_values, y, ...)
  }
}
#' @name plot_filters
#' @rdname plot_filters
#' @export
plot_phase <- function(x, q = 0, nxlab = 6, add = FALSE,...){
  x_values = seq(0, 2 * pi, length.out = nrow(x))

  if(q < ncol(x))
    y = x[,sprintf("q=%i",q)]

  if(!add){
    if(q >= ncol(x)){
      plot(1, type="n",xaxt = "n", xlab = "",
           ylab = "phase", xlim=c(0, 2*pi), ylim=c(0, 1),
           ...)
    }else{
      plot(x_values, y, type = "l",
           xaxt = "n", xlab = "",
           ylab = "phase", ...)
    }

    x_lab_at <- seq(0, 2, length.out = nxlab)
    axis(1, at = x_lab_at * pi, labels = xlabel(x_lab_at))
  }else{
    lines(x_values, y, ...)
  }
}

xlabel <- function(x, symbol = "pi"){
  fracs <- strsplit(attr(MASS::fractions(x), "fracs"), "/")  # convert to fractions
  labels <- sapply(fracs, function(i)
    if (length(i) > 1) { paste(i[1], "*", symbol, "/", i[2]) }
    else { paste(i, "*", symbol) })
  labels[1] <- "0"
  parse(text = labels)
}
