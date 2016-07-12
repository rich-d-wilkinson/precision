##
## Code for R wrappers that call C functions to simulate from distributions conditioned to be within a certain range
##  - for separate rmultbinom functions see the multbinom package, although it shouldn't be needed I think.
##
## Main function is rNM_given_nm - simulates proposals of N and M from the required distributions subject to the 
##  constraints that N>=n, M>=m, M-m<=N-n
##  - Distribtions it simulates from are N ~ Po(lambda) and then either i) M~ Bin(N, p) or M ~ multiplicateBinom(N, p, psi)
##
## C is slower than R for pure binomial random number generation, but C is much quicker for the greater than simulation

##

rNM_given_nm <- function(n, lambda, n.data, p, psi=NULL, m, family="binomial"){
  if (family=="binomial"){
    stopifnot(is.null(psi))
    if (length(n.data)==length(m) & length(m)==1){
     out <- .C("rNM_given_nm_binomial",
       n = as.integer(n), 
       lambda = as.double(lambda),  
       n_data = as.integer(n.data), 
       p = as.double(p), 
       m = as.integer(m), 
       N = as.integer(rep(0, n)), 
       M = as.integer(rep(0, n)),     
       seed = as.integer(round(10^8)*runif(1))
      )
      return(list(N=out$N, M=out$M))
    }
    else{
      stopifnot(length(n.data)==length(m))
      tmp <- cbind(n.data, m)
      out<-apply(tmp, 1, function(x){
            tmp2 <- .C("rNM_given_nm_binomial",
                  n = as.integer(n), 
                  lambda = as.double(lambda),  
                  n_data = as.integer(x[1]), 
                  p = as.double(p), 
                  m = as.integer(x[2]), 
                  N = as.integer(rep(0, n)), 
                  M = as.integer(rep(0, n)),     
                  seed = as.integer(round(10^8)*runif(1))
        )
        return(c(tmp2$N, tmp2$M))
      })
       N <- out[ 1:n, ]
       M <- out[ (n+1):(2*n), ]
       out2 <- list(N = N, M =M)
      return(out2)
    }
  }
  
  ## multiplicate binomial case
  if (family=="multbinom"){
    stopifnot(!is.null(psi))
    if (length(n.data)==length(m)& length(m)==1){
      out <- .C("rNM_given_nm_multbinom",
                n = as.integer(n), 
                lambda = as.double(lambda),  
                n_data = as.integer(n.data), 
                p = as.double(p), 
                psi = as.double(psi),
                m = as.integer(m), 
                N = as.integer(rep(0, n)), 
                M = as.integer(rep(0, n)),     
                seed = as.integer(round(10^8)*runif(1))
      )
      return(list(N=out$N, M=out$M))
    }
    else{
      stopifnot(length(n.data)==length(m))
      tmp <- cbind(n.data, m)
      out <- apply(tmp, 1, function(x){
        tmp2 <- .C("rNM_given_nm_multbinom",
                   n = as.integer(n), 
                   lambda = as.double(lambda),  
                   n_data = as.integer(x[1]), 
                   p = as.double(p), 
                   psi = as.double(psi),
                   m = as.integer(x[2]), 
                   N = as.integer(rep(0, n)), 
                   M = as.integer(rep(0, n)),     
                   seed = as.integer(round(10^8)*runif(1))
        )
        return(c(tmp2$N, tmp2$M))
      })
      N <- out[ 1:n, ]
      M <- out[ (n+1):(2*n), ]
      out2 <- list(N = N, M =M)
      return(out2)
    }
  }
  ## double binomial case
  if (family=="doublebinom"){
    stopifnot(!is.null(psi))
    if (length(n.data)==length(m)& length(m)==1){
      out <- .C("rNM_given_nm_doublebinom",
                n = as.integer(n), 
                lambda = as.double(lambda),  
                n_data = as.integer(n.data), 
                p = as.double(p), 
                psi = as.double(psi),
                m = as.integer(m), 
                N = as.integer(rep(0, n)), 
                M = as.integer(rep(0, n)),     
                seed = as.integer(round(10^8)*runif(1))
      )
      return(list(N=out$N, M=out$M))
    }
    else{
      stopifnot(length(n.data)==length(m))
      tmp <- cbind(n.data, m)
      out <- apply(tmp, 1, function(x){
        tmp2 <- .C("rNM_given_nm_doublebinom",
                   n = as.integer(n), 
                   lambda = as.double(lambda),  
                   n_data = as.integer(x[1]), 
                   p = as.double(p), 
                   psi = as.double(psi),
                   m = as.integer(x[2]), 
                   N = as.integer(rep(0, n)), 
                   M = as.integer(rep(0, n)),     
                   seed = as.integer(round(10^8)*runif(1))
        )
        return(c(tmp2$N, tmp2$M))
      })
      N <- out[ 1:n, ]
      M <- out[ (n+1):(2*n), ]
      out2 <- list(N = N, M =M)
      return(out2)
    }
  }
  
}
  
  
  



## checked seems to work fine - NOT SURE THIS IS NEEDED ANY MORE - DELETE IF NOT
rpoisb <- function(n, lambda, lb, up=NULL){
  
  ## ndata can be a vector
  stopifnot(!is.null(lb))
  if(length(lb)==1){
    out <- .C("rpois_lb",
              n = as.integer(n),
              lambda = as.double(lambda),
              lb = as.integer(lb),
              rndm = as.integer(rep(0,n)),
              seed = as.integer(round(10^8)*runif(1))
              )
    return(out$rndm)
    }
  else{
   out<- sapply(lb, function(x){tmp<- .C("rpois_lb",
                  n = as.integer(n),
                  lambda = as.double(lambda),
                  lb = as.integer(x),
                  rndm = as.integer(rep(0,n)),
                  seed = as.integer(round(10^8)*runif(1))
              )
              return(tmp$rndm)
    })
  return(out)
  }
}



##########################################################
###
###
### Bounded binomial random variables
### - NOT SURE THIS IS NEEDED ANY MORE - DELETE IF NOT
###
##########################################################

           
rbinomb <-function(n, size, p, lb, ub=NULL){
  
  stopifnot(!is.null(lb))
  stopifnot(length(lb)==length(size))
  ## no upper bound
  if(is.null(ub)){
    if(length(lb)==1){
      out <- .C("rbinom_lb",
                n = as.integer(n),
                size = as.integer(size),
                p = as.double(p),
                lb = as.integer(lb),
                rndm = as.integer(rep(0,n)),
                seed = as.integer(round(10^8)*runif(1))
                )
      return(out$rndm)
    }
    else{
      tmp1<-cbind(size, lb)
      out<- apply(tmp1, 1, function(x){tmp<- .C("rbinom_lb",
                                            n = as.integer(n),
                                            size = as.integer(x[1]),
                                            p = as.double(p),
                                            lb = as.integer(x[2]),
                                            rndm = as.integer(rep(0,n)),
                                            seed = as.integer(round(10^8)*runif(1))
                                            )
                                   return(tmp$rndm)
            })
      return(out)
    }
  }
  else{   ### lb and ub
    stopifnot(length(lb)==length(ub), length(lb)==length(size))
    if(length(lb)==1){
      out <- .C("rbinom_lb_ub",
                n = as.integer(n),
                size = as.integer(size),
                p = as.double(p),
                lb = as.integer(lb),
                ub = as.integer(ub),
                rndm = as.integer(rep(0,n)),
                seed = as.integer(round(10^8)*runif(1))
      )
      return(out$rndm)      
    }
    else{
      tmp1 <-cbind(size,lb, ub)
      out<- apply(tmp1, 1, function(x){tmp<- .C("rbinom_lb_ub",
                                            n = as.integer(n),
                                            size = as.integer(x[1]),
                                            p = as.double(p),
                                            lb = as.integer(x[2]),
                                            ub = as.integer(x[3]),
                                            rndm = as.integer(rep(0,n)),
                                            seed = as.integer(round(10^8)*runif(1))
                                            )
                                   return(tmp$rndm)
      })
      return(out)
    }
  }
}
      
      
      
      
    
  
  