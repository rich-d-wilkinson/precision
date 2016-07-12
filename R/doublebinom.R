##
## Code for R wrappers that call C functions to simulate from the double binomial distribution.
##
##
##
##


## evaluates the density of a multiplicative binomial
ddoublebinom <- function(x=NULL, size, p, psi, log=FALSE){
  
  out <- .C("dDoublebinomC", size_R=as.integer(size),
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


rdoublebinom <- function(n=1, size, prob, psi){
  
  out <- .C("rDoublebinomC", 
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