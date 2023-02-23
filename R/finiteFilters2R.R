FiniteFilters2R <- function(jprops, horizon, rev = FALSE){
  sw<-.proc_data(jprops, "sweights")
  swg<-.proc_data(jprops, "sgain")

  if (rev) {
    asym_h <- 0:(horizon-1)
  } else {
    asym_h <- (horizon-1):0
  }

  aw<-sapply(asym_h, function(h){return(.proc_data(jprops, paste0("aweights(", h,')')))})
  awg<-sapply(asym_h, function(h){return(.proc_data(jprops, paste0("again(", h,')')))})
  awp<-sapply(asym_h, function(h){return(.proc_data(jprops, paste0("aphase(", h,')')))})

  coefs = c(list(sw), aw)
  nbpoints = horizon*2+1
  coefs = sapply(coefs, function(x){
    c(x,rep(0,nbpoints-length(x)))
  })

  gain = cbind(swg, awg)
  phase = cbind(0,awp)

  filternames <- sprintf("q=%i", horizon:0)
  rownames(coefs) <- coefficients_names(-horizon, horizon)
  colnames(gain) <- colnames(coefs) <- colnames(phase) <-
    filternames
  return(list(
    internal = jprops,
    filters.coef = coefs,
    filters.gain = gain,
    filters.phase= phase
  ))
}
