#' Get properties of local polynomials filters
#'
#' @param x a \code{"JD.Filters"} object.
#' @param component the component to extract.
#' @param ... unused other arguments.
#'
#' @examples
#' filter <- filterproperties(3, kernel = "Henderson")
#' sgain <- get_properties_function(filter, "Symmetric Gain")
#' plot(sgain, xlim= c(0, pi/12))
#' @export
get_properties_function <- function(x,
                                    component = c("Symmetric Gain",
                                                  "Symmetric Phase",
                                                  "Symmetric transfer",
                                                  "Asymmetric Gain",
                                                  "Asymmetric Phase",
                                                  "Asymmetric transfer"),
                                    ...){
  UseMethod("get_properties_function", x)
}
#' @export
get_properties_function.JD.Filters <- function(x,
                                               component = c("Symmetric Gain",
                                                             "Symmetric Phase",
                                                             "Symmetric transfer",
                                                             "Asymmetric Gain",
                                                             "Asymmetric Phase",
                                                             "Asymmetric transfer"), ...){
  component = match.arg(component)
  switch(component,
         "Symmetric Gain" = {
           get_gain_function(x$internal$getFilter())
         },
         "Asymmetric Gain" = {
           afunction = .jevalArray(x$internal$getAfilters())
           afunction = sapply(afunction, get_gain_function)
           names(afunction) <- sprintf("q=%s", seq(0, length(afunction)-1))
           afunction
         },
         "Symmetric Phase" = {
           get_phase_function(x$internal$getFilter())
         },
         "Asymmetric Phase" = {
           afunction = .jevalArray(x$internal$getAfilters())
           afunction = sapply(afunction, get_phase_function)
           names(afunction) <- sprintf("q=%s", seq(0, length(afunction)-1))
           afunction
         },
         "Symmetric transfer" = {
           get_frequencyResponse_function(x$internal$getFilter())
         },
         "Asymmetric transfer" = {
           afunction = .jevalArray(x$internal$getAfilters())
           afunction = sapply(afunction, get_frequencyResponse_function)
           names(afunction) <- sprintf("q=%s", seq(0, length(afunction)-1))
           afunction
         })
}
get_gain_function <- function(x){
  jgain <- x$gainFunction()$applyAsDouble
  Vectorize(function(x){
    jgain(x)
  })
}
get_phase_function <- function(x){
  jphase <- x$phaseFunction()$applyAsDouble
  Vectorize(function(x){
    jphase(x)
  })
}
get_frequencyResponse_function <- function(x){
  jfrf <- x$frequencyResponseFunction()$apply
  Vectorize(function(x){
    res <- jfrf(x)
    res$getRe()+res$getIm()*i
  })
}
