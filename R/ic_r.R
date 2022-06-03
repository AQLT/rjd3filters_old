#' Compute IC-Ratio
#'
#' @param x input time series.
#' @param sc trend-cycle component.
#' @param boolean indication if the decomposition is multiplicative or additive.
#'
#' @examples
#' x = retailsa$AllOtherGenMerchandiseStores
#' sc <- henderson(x, length = 13, musgrave = FALSE)
#' ic_ratio(x, sc)
#'
#' @export
ic_ratio <- function(x, sc, mul = FALSE){
  remove_na <- is.na(x) | is.na(sc)
  x = as.numeric(x)[!remove_na]
  sc = as.numeric(sc)[!remove_na]
  result <- .jcall("demetra/saexperimental/r/X11Decomposition",
                   "D", "icratio",
                   x, sc, mul)
  result
  result
}
