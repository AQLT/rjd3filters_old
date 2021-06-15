#' @export
print.fst_filter <- function(x, ...){
  cat("$filters.coef\n")
  print(x$filters.coef, ...)
}
#' @export
print.lp_filter <- function(x, ...){
  cat("$filters.coef\n")
  print(x$filters.coef, ...)
}
#' @export
print.rkhs_filter <- function(x, ...){
  cat("$filters.coef\n")
  print(x$filters.coef, ...)
}
