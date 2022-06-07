#' Get Macurves Filters
#'
#' @param seas_filter the filter to extract
#'
#' @export
macurves <- function(seas_filter = c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15")){
  seas_filter = match.arg(toupper(seas_filter),
                          c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"))
  seas_opt <- .jcall("jdplus/x11plus/SeasonalFilterOption",
                     "Ljdplus/x11plus/SeasonalFilterOption;",
                     "valueOf",
                     seas_filter)
  P = .jcast(new( J("java.lang.Double"), as.character(12)),
             "java.lang.Number")
  seasFilter <- .jcall("jdplus/x11plus/X11SeasonalFiltersFactory",
                       "Ljdplus/filters/ISymmetricFiltering;",
                       "filter",
                       P, seas_opt)
  symFilter <- .jcall(seasFilter, "Ljdplus/math/linearfilters/SymmetricFilter;", "symmetricFilter")
  asymFilter <- .jcall(seasFilter, "[Ljdplus/math/linearfilters/IFiniteFilter;", "endPointsFilters")

  coefs = c(rev(asymFilter),list(symFilter))
  nbpoints = symFilter$length()
  horizon = (nbpoints - 1)/2
  coefs = sapply(coefs, function(x){
    x = x$weightsToArray()
    c(x,rep(0,nbpoints-length(x)))
  })
  rownames(coefs) <- coefficients_names(-horizon, horizon)
  colnames(coefs) <- sprintf("q=%i", 0:(horizon))

  coefs
}
