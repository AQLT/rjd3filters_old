#' @export
print.fst_filter <- function(x, ...){
  print("Coefficients of FST filter:\n")
  print(x$filters.coef, ...)
}
#' @export
print.lp_filter <- function(x, ...){
  print("Coefficients of Local Polynomial filter:\n")
  print(x$filters.coef, ...)
}
#' @export
print.rkhs_filter <- function(x, ...){
  cat("Coefficients of the RKHS filter:\n")
  print(x$filters.coef, ...)
}
