###
### dm_given_nNM.R
###  - commands for the density p(m | n,N,M)
###  - and p(m,n | N, M)
###


dm_given_nNM <-function(n,m, N,M, log=FALSE){
  k=M - m  ## number of males who died
  r= N- n  ## number of deaths
  
  if((k<0 | r<0) | k>r) { print("Warning - we shouldn't be here - dm_given_nNM") ## we should never enter this loop
    if(log) return(-Inf)  ## this has prob zero
    else return(0) 
  }               
  out<-choose(M, k)*choose(N-M, r-k)/choose(N, r)  
#  print(paste("dm_given_nNM",n,m,N,M,log(out)))
  if(log) return(log(out))
  else return(out)
}


log_dmn_given_NM <- function(n,m,N,M, mort){
  tmp1 <- dm_given_nNM(n,m, N,M, log=TRUE)
  tmp2 <- dbinom(x=N-n, size=N, prob=mort, log=TRUE)   ## note N-n died out of N
  return(tmp1+tmp2)
}

log_dmn_given_NM_All <-  function(n.vec,m.vec,N.vec,M.vec, mort){
  tmp <- apply(cbind(n.vec, m.vec, N.vec, M.vec), 1, function(x) log_dmn_given_NM(n=x[1], m=x[2], N=x[3], M=x[4], mort=mort))
  return(sum(tmp))
}
