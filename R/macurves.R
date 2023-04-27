#' Get Macurves Filters
#'
#' @param seas_filter the filter to extract.
#' @param period period of the filter.
#' @examples
#' macurves("S3X3")
#' @export
macurves <- function(seas_filter = c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"), period = 12){
  seas_filter <- match.arg(toupper(seas_filter),
                          c("S3X3", "S3X1", "S3X5", "S3X9", "S3X15"))
  seas_opt <- .jcall("jdplus/experimentalsa/base/core/x11plus/SeasonalFilterOption",
                     "Ljdplus/experimentalsa/base/core/x11plus/SeasonalFilterOption;",
                     "valueOf",
                     seas_filter)
  P <- .jcast(new( J("java.lang.Double"), as.character(period)),
             "java.lang.Number")
  seasFilter <- .jcall("jdplus/experimentalsa/base/core/x11plus/X11SeasonalFiltersFactory",
                       "Ljdplus/experimentalsa/base/core/filters/ISymmetricFiltering;",
                       "filter",
                       P, seas_opt)
  .jd2r_finitefilters(seasFilter)
}
