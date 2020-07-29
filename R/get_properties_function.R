#' Get properties of local polynomials filters
#'
#' @param x a \code{"lp_filter"} object.
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
get_properties_function.lp_filter <- function(x,
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
#' @export
get_properties_function.fst_filter <- function(x,
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
           get_gain_function(x$internal$getFilter())
         },
         "Symmetric Phase" = {
           get_phase_function(x$internal$getFilter())
         },
         "Asymmetric Phase" = {
           get_phase_function(x$internal$getFilter())
         },
         "Symmetric transfer" = {
           get_frequencyResponse_function(x$internal$getFilter())
         },
         "Asymmetric transfer" = {
           get_frequencyResponse_function(x$internal$getFilter())
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
    complex(real = res$getRe(), imaginary = res$getIm())
  })
}

#' @export
diagnostics_matrix <- function(x, lb, passband = pi/6){
  if (lb >=0)
    lb <- -lb

  results <- c(sum(x)-1, sum(x*seq(lb,length(x) + lb-1, by = 1)),
               sum(x*seq(lb,length(x) + lb-1, by = 1)^2),
               fst(x, lb, passband = passband)$criterions)
  names(results) <- c("b_c", "b_l", "b_q",
                      "F_g", "S_g", "T_g")
  results
}
