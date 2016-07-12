##
## Code for R wrappers that call C functions to simulate from the multiplicate binomial distribution.
##  - includes code for conditioning  values to be within a certain range
##
##
##
##


## evaluates the density of a multiplicative binomial
dmultbinom <- function(x=NULL, size, p, psi, log=FALSE){
  
  out <- .C("dMultbinomC", size_R=as.integer(size),
       p_R=as.double(p), 
       psi_R=as.double(psi),  
       logden=as.double(rep(0, size+1)), 
       den=as.double(rep(0,size+1)))
  
  if(!is.null(x)) {
    if(log) return(out$logden[x+1])
    else return(out$den[x+1])
  }
  else{
    if(log) return(out$logden)
    else return(out$den)      
  }
}
## checked - seems fine


rmultbinom <- function(n=1, size, prob, psi){

  out <- .C("rMultbinomC", 
              n=as.integer(n), 
              size_R=as.integer(size),
              p_R=as.double(prob), 
              psi_R_=as.double(psi),  
              rndm=as.integer(rep(0,n)),
              seed=as.integer(round(10^8*runif(1))),
              den=as.double(rep(0,size+1))
    )
  return(out$rndm)
  
}

