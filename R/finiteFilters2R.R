FiniteFilters2R <- function(jprops, horizon, rev = FALSE){
  sw<-proc_data(jprops, "sweights")
  swg<-proc_data(jprops, "sgain")
  if(rev){
    asym_h <- (horizon-1):0
  }else{
    asym_h <- 0:(horizon-1)
  }

  aw<-sapply(asym_h, function(h){return(proc_data(jprops, paste0("aweights(", h,')')))})
  awg<-sapply(asym_h, function(h){return(proc_data(jprops, paste0("again(", h,')')))})
  awp<-sapply(asym_h, function(h){return(proc_data(jprops, paste0("aphase(", h,')')))})

  coefs = c(aw,list(sw))
  nbpoints = horizon*2+1
  coefs = sapply(coefs, function(x){
    c(x,rep(0,nbpoints-length(x)))
  })

  gain = cbind(awg,swg)
  phase = cbind(awp, 0)

  filternames <- sprintf("q=%i", 0:(horizon))
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
