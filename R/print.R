#' @export
print.fst_filter <- function(x, ...){
  cat("$filters.coef\n")
  print(x$filters.coef, ...)
}
#' @export
print.FiniteFilters <- function(x, ...){
  cat("$filters.coef\n")
  print(x$filters.coef, ...)
}
