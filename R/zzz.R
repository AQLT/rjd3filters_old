#' @import rJava
#' @importFrom graphics axis lines plot matplot
#' @importFrom stats frequency ts
#' @importFrom rjd3toolkit proc_data proc_dictionary proc_vector
#'@importFrom rjd3toolkit matrix_jd2r matrix_r2jd ts_jd2r ts_r2jd tsdomain_r2jd
NULL

.onLoad <- function(libname, pkgname) {
  # For debugts_ging: to see if Jars are effectively loaded
  # options(java.parameters = "-verbose:class")

  # TODO : devtools will look only in RJDemetra3\java for JAR files so copied them there too
  result <- .jpackage(pkgname, lib.loc=libname)
  if (!result) stop("Loading java packages failed")

  # what's your java  version?  Need > 1.5.0.
  jversion <- .jcall('java.lang.System','S','getProperty','java.version')
  if (jversion < "1.8.0") {
    stop(paste("Your java version is ", jversion,
               ".  Need 1.8.0 or higher.", sep=""))
  }

}

.onAttach <- function( libname , pkgname ){
  jversion <- .jcall('java.lang.System','S','getProperty','java.version')
  packageStartupMessage("Java requirements fullfilled, found version ",jversion)
}


#' Retail sales data
#'
#' A dataset containing monthly retailed sales
#'
#' @docType data
#' @format A \code{list} of \code{ts} objects from january 1992 to december 2010.
"retailsa"
