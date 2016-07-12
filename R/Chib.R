#################################################################################
### Part of wasps package
###
### Chib.R
###  - additional commands needed to implement Chib's method for calculating the model evidence
###
###
#################################################################################


## TEST ALL THESE FUNCTIONS

## function to calculate log(mean(exp(x_i))) when exp(x_i) very small
MeanLog <- function(log.vals, trim=0){
  
  min.val <- min(log.vals)
  if(is.infinite(min.val)){ return(log(mean(exp(log.vals))))}
  shifted <- log.vals - min.val
  logmean <- log(mean(exp(shifted), trim=trim))+min.val
  if(!is.infinite(logmean))  return(logmean)
  else return(log(mean(exp(log.vals))))
}

#vals <- -100:-90
# log(mean(exp(vals)))
# MeanLog(vals)

############################
### logprior() estimate the log of the prior at theta*
############################






## checked
LogPrior <- function(theta.star, hyper, family){
  
  if(family=="binomial") stopifnot(length(theta.star)==3)
  tmp<-c()
  tmp[1] <- dgamma(theta.star["lambda"], shape = hyper$alpha.lambda,   rate =  hyper$beta.lambda, log=TRUE)
  tmp[2] <- dbeta(x=theta.star["mort"], shape1 = hyper$a.m,   shape2 = hyper$b.m, log=TRUE)
  tmp[3] <- dbeta(x=theta.star["p"], shape1=hyper$a.p, shape2=hyper$b.p, log=TRUE)
  
  if(family=="multbinom" | family=="doublebinom") {
    stopifnot(length(theta.star)==4)
    tmp[4] <- dnorm(x=theta.star["psi"], mean=hyper$mu.psi, sd=hyper$sd.psi, log=TRUE)
  }
 return(sum(tmp))  
  
}





############################
### EstimateLogLikelihood() estimate the log of the likelihood at theta*
############################



## Given a value of N,M,n,m and theta, it returns log[ pi(nm|NM,theta)pi(N|theta)pi(M|N,theta) ]
LogLikeSummand <- function(N,M, n, m, theta.star, family){
  if(N < n) return(-Inf)
  if(M < m) return(-Inf)
  if(N - n < M - m) return(-Inf)
  
  out <- log_dmn_given_NM(n=n, 
                          m = m, 
                          N =N, 
                          M =M, 
                          mort=theta.star["mort"])
  out2 <- dpois(x=N, lambda=theta.star["lambda"], log=TRUE)
  if(family=="binomial")  out3 <- dbinom(x=M, size = N, prob= theta.star["p"], log=TRUE )
  if(family=="multbinom") out3 <- dmultbinom(x=M, size=N, p=theta.star["p"], psi=theta.star["psi"], log=TRUE)
  if(family=="doublebinom") out3 <- ddoublebinom(x=M, size=N, p=theta.star["p"], psi=theta.star["psi"], log=TRUE)
  return(out+out2+out3)
}






## Given a value of N,M,n,m and theta, it returns log E(pi(nm|NM,theta)) wrt prior pi(N|theta)pi(M|N,theta)
LogSumOneExp <- function(n, m, theta.star, family=family, max.n){
  
  
  val <-1
  total <- 0
  N=n
  while(val >10^-10| N<=max.n*5){
    #print(paste("*****N=",N))
    for(M in m:(N-n+m)){
      val <- exp(LogLikeSummand(N,M, n, m, theta.star, family=family))
      total <- total + val
     # print(paste('N=', N, 'M=', M, 'val=', val))
      # print(paste("M=",M))
      #  print(total)
    }
    N <- N+1
  }
  return(log(total))
}



## it gives log pi(n,m | theta) = log E(pi(n,m|NM,theta)) for the data vectors n,m - this is the product of the individual terms
LogLike <- function(data, theta.star, family="binomial"){
  max.n <- max(data[,"n"])
  out <- apply(data, 1, function(x) LogSumOneExp(n = x[1], m = x[2], theta.star, family=family, max.n =max.n))
  return(sum(out))
}
















####################################################################################
###
### EstimateLogPosterior() estimate the log of the posterior at theta*
###
### - use the MCMC output from the full chain and use the sum pi(theta*|NM^i, n,m)
###
####################################################################################






log_dtheta_NMnm <- function(theta.star, NM, data,  hyper, n.sum, family){
  ## N, M are vectors of length exp.size
  
  if(family=="binomial"){
    stopifnot(dim(data)==dim(NM))
    tmp<-c()
    sumNM <- colSums(NM)
    tmp["p"] <- dbeta(x=theta.star["p"], shape1 = sumNM["M"] + hyper$a.p,   shape2 = sumNM["N"] - sumNM["M"] + hyper$b.p, log=TRUE)   	
    tmp["mort"] <- dbeta(x=theta.star["mort"], shape1 = sumNM["N"] - n.sum + hyper$a.m,   shape2 = n.sum + hyper$b.m,  log=TRUE) 		
    tmp["lambda"] <- dgamma(x=theta.star["lambda"], shape = hyper$alpha.lambda+ sumNM["N"],   rate =  hyper$beta.lambda + dim(NM)[1], log=TRUE) 		
  }
  return(sum(tmp))    
}

### send mcmc.out and mcmc.out$fixed.p.psi ??
## split as two function EstimateLogPosterioBinom and MultBinom?

EstimateLogPosteriorBinomial <- function(data, theta.star, mcmc.out, family){
  
  n.sum = sum(data[,"n"])
 
  if(family=="binomial"){

    stopifnot(length(theta.star)==3)
    tmp <-cbind(mcmc.out$N.chain, mcmc.out$M.chain)
    out<-apply(tmp, 1, function(x){
                              NM<-matrix(x, ncol=2, byrow=FALSE)
                              colnames(NM)<-c("N", "M")
                              return(log_dtheta_NMnm(theta.star=theta.star, NM=NM, data=data,  hyper=hyper, n.sum=n.sum, family=family))
                              })  
  }

  return(list(estimate=MeanLog(out), logvals=out)) ## because the numbers are all small
}



## pi(\theta2* | n,m, theta1*)
## uses the mcmc fun with fixed p* and psi*
EstimateLogPosteriorMultbinomTerm1 <- function(data, theta.star, mcmc.out, family){
  
  n.sum = sum(data[,"n"])
  
  if(family=="multbinom" | family=="doublebinom"){
    stopifnot(length(theta.star)==4, var(mcmc.out$chain[,"p"])==0, var(mcmc.out$chain[,"psi"])==0)## use the fixed p and psi MCMC results
    
    ## pi(theta_2*|n,m, \theta_1*)
    out<-apply(mcmc.out$N.chain, 1, function(x){
      return(log_dtheta1_NMnmtheta2(theta.star=theta.star, N=x, data=data,  hyper=hyper, n.sum=n.sum, family=family))
    })
    
    ## Now do  pi(theta_1*|n,m)
    
  }

  return(list(estimate=MeanLog(out), logvals=out)) ## because the numbers are all small
}








log_dtheta1_NMnmtheta2 <- function(theta.star, N, data,  hyper, n.sum, family){
  ## N, M are vectors of length exp.size
  
  if(family=="multbinom" | family=="doublebinom"){
    stopifnot(dim(data)[1]==length(N))
    tmp<-c()
    sumN <- sum(N)
    tmp["mort"] <- dbeta(x=theta.star["mort"], shape1 = sumN - n.sum + hyper$a.m,   shape2 = n.sum + hyper$b.m,  log=TRUE) 		
    tmp["lambda"] <- dgamma(x=theta.star["lambda"], shape = hyper$alpha.lambda+ sumN,   rate =  hyper$beta.lambda + dim(data)[1], log=TRUE) 		
  }
  return(sum(tmp))    
}



LogQDensity <- function(theta.cur, theta.prop, step.size){
  stopifnot(step.size!=c(0,0))
  tmp1 <- dnorm(x=theta.prop["psi"], mean=theta.cur["psi"], sd=step.size["psi"], log=TRUE)
  tmp2 <- dnorm(x=Logit(theta.cur["p"]), mean=Logit(theta.prop["p"]), sd=step.size["p.logit"], log=TRUE)
  return(tmp1 + tmp2)  
}


#### check AcceptPPsi
LogPMHSubKernel <- function(theta.cur, theta.prop, NM, hyper, step.size, family){
  tmp1 <- AcceptPPsi(theta.cur=theta.cur, theta.prop=theta.prop, NM=NM, hyper=hyper, family=family)$log.prob
  tmp2 <- LogQDensity(theta.cur=theta.cur, theta.prop = theta.prop, step.size = step.size)
 # print(c(exp(tmp1), exp(tmp2)))
  return(tmp1+tmp2)
}


### mcmc.out should be from the full chain
ChibLogNumerator <- function(mcmc.out, theta.star, hyper, data, step.size, family){
  
  stopifnot(!is.null(mcmc.out$N.chain), !is.null(mcmc.out$M.chain))
  stopifnot(var(mcmc.out$chain[,"p"])>0, step.size[1]!=0, step.size[2]!=0)  ## uses the full MCMC chain
  chain <- mcmc.out$chain
  chainN <- mcmc.out$N.chain 
  chainM <- mcmc.out$M.chain
  
  nparam = dim(chain)[2]
  dimN <- dim(chainN)[2]
  
  tmp <- cbind( chainN, chainM,chain)
  
  out <- apply(tmp, 1, function(x){
                            param <- x[(2*dimN+1):(2*dimN+nparam)]
                            NM <- cbind( x[ 1:dimN ], x[ (dimN+1):(2*dimN)] )
                            colnames(NM) <- c("N", "M")                            
                            return( LogPMHSubKernel(theta.cur=param, theta.prop=theta.star, NM=NM, hyper=hyper, step.size=step.size,family=family))  
  })
  ## use log(mean(exp)) rather than LogMean here as the range of values is too large
  return(list(estimate=log(mean(exp(out))), logvals=out)) ## because the numbers are all small  
}



## mcmc.out should be from fixed p.psi
ChibLogDenom <- function(mcmc.out, theta.star, hyper, data, step.size, family){
  ## uses the  MCMC with p and psi fixed
  stopifnot(var(mcmc.out$chain[,"p"])==0, var(mcmc.out$chain[,"psi"])==0, step.size[1]!=0, step.size[2]!=0)  
  
  chain <- mcmc.out$chain
  chainN <- mcmc.out$N.chain 
  chainM <- mcmc.out$M.chain
  
  chain[,"p"] <- InvLogit(rnorm(n=dim(chain)[1], mean= Logit(theta.star["p"]), sd=step.size["p.logit"]))  ## simulate values of p and psi
  chain[,"psi"] <- rnorm(n=dim(chain)[1], mean= theta.star["psi"], sd=step.size["psi"])
  
  nparam = dim(chain)[2]
  dimN <- dim(chainN)[2]
  
  tmp <- cbind(chainN, chainM, chain)
  out <- apply(tmp, 1, function(x){
                             param <- x[(2*dimN+1):(2*dimN+nparam)]
                             NM <- cbind( x[ 1:dimN ], x[ (dimN+1):(2*dimN)] )
                             colnames(NM) <- c("N", "M")                            
                             return(AcceptPPsi(theta.cur=theta.star, theta.prop=param, NM=NM, hyper=hyper, family=family)$log.prob)
                           })
 ## use log(mean(exp)) rather than LogMean here as the range of values is too large
  return(list(estimate=log(mean(exp(out))), logvals=out)) ## because the numbers are all small  

}


### mcmc.out should be from the full chain
ChibLogNumeratorNoMort <- function(mcmc.out, theta.star, hyper, data, step.size, family){
  
  stopifnot(var(mcmc.out$chain[,"p"])>0, step.size[1]!=0, step.size[2]!=0)  ## uses the full MCMC chain
  chain <- mcmc.out$chain
  
  nparam = dim(chain)[2]
  
  NM <- data
  colnames(NM) <- c("N", "M")                            
  
  out <- apply(chain, 1, function(x){
    return( LogPMHSubKernel(theta.cur=x, theta.prop=theta.star, NM=NM, hyper=hyper, step.size=step.size, family=family))  
  })
  ## use log(mean(exp)) rather than LogMean here as the range of values is too large
  return(list(estimate=log(mean(exp(out))), logvals=out)) ## because the numbers are all small  
}


## mcmc.out should be from fixed p.psi
ChibLogDenomNoMort <- function( theta.star, hyper, data, step.size, n.mc, family){
  ## uses the  MCMC with p and psi fixed
  stopifnot(step.size[1]!=0, step.size[2]!=0)  
  
  chain <- matrix(nrow=n.mc, ncol=2)
  colnames(chain)<-c("p", "psi")
  chain[,"p"] <- InvLogit(rnorm(n=n.mc, mean= Logit(theta.star["p"]), sd=step.size["p.logit"]))  ## simulate values of p and psi
  chain[,"psi"] <- rnorm(n=n.mc, mean= theta.star["psi"], sd=step.size["psi"])
  
  
  out <- apply(chain, 1, function(x){
    param <- x
    NM <- data
    colnames(NM) <- c("N", "M")                            
    return(AcceptPPsi(theta.cur=theta.star, theta.prop=param, NM=NM, hyper=hyper, family=family)$log.prob)
  })
  ## use log(mean(exp)) rather than LogMean here as the range of values is too large
  return(list(estimate=log(mean(exp(out))), logvals=out)) ## because the numbers are all small  
  
}



####################################################################################
###
### Bootstrap.Var() estimate the variance of the estimators using bootstrap replicates
###
####################################################################################

Bootstrap.Var <- function(vals, B=1000){
  boot.rep <- c()
  for(i in 1:B){
    smp <- sample(vals, size=length(vals), replace = TRUE)
    boot.rep[i]<- MeanLog(smp)
    if(is.infinite(boot.rep[i])) boot.rep[i] <- log(mean(exp(smp)))
  }
  return(list(var=var(boot.rep), vals=boot.rep))
}




####################################################################################
###
### CalculateEvidence() - wraps all the other functions into a single easily callable function.
###
####################################################################################




CalculateEvidence <- function(mcmc.out, data,  hyper, family="binomial", sd=FALSE, nbatch=NULL, step.size=NULL, thin=FALSE, thinby=NULL, burnin=NULL){
  
  stopifnot(family == "binomial" | family == "multbinom" | family == "doublebinom")
  stopifnot(!is.null(mcmc.out$N.chain), !is.null(mcmc.out$N.chain)) 
  if(family!="binomial"){
    stopifnot(!is.null(nbatch), !is.null(step.size))
    if(thin) stopifnot(!is.null(thinby), !is.null(burnin))
  }
  ## nbatch only required for double and mult
  ## thin, thinby, burnin
  
  details <- list()
  theta.star <- colMeans(mcmc.out$chain)
  details$theta.star = theta.star
  log.prior <- LogPrior(theta.star = theta.star, hyper = hyper, family = family )
  
  if(family!="binomial"){
    
    log.post.term2.num.out <- ChibLogNumerator(mcmc.out = mcmc.out, theta.star = theta.star, hyper = hyper, data = data, step.size = step.size, family=family)
    log.post.term2.num <- log.post.term2.num.out$estimate
    if(sd) log.post.term2.num.sd <- Bootstrap.Var(log.post.term2.num.out$logvals, B=1000)$var
    
    details$log.post.term2.num = log.post.term2.num
    
    print("Running MCMC chain - may take some time...")
    step.size.fixed<-c(0,0)   ## 0.3 and 0.2 aree 
    names(step.size.fixed) <- c("p.logit", "psi")
    mcmc.fixed <- MCMCWithinGibbs( theta0 = theta.star,  data=data, hyper=hyper, nbatch=nbatch,  family=family, step.size=step.size.fixed, keepNM=TRUE)
    
    if(thin){
      mcmc.fixed.t  <- ThinChain(mcmc.fixed, thinby=thinby, burnin=burnin)
    }
    if(!thin){
      mcmc.fixed.t <- mcmc.fixed
    }
    
    log.post.term1.out <- EstimateLogPosteriorMultbinomTerm1(data=data, theta.star=theta.star, mcmc.out=mcmc.fixed.t, family=family)
    log.post.term1 <- log.post.term1.out$estimate
    details$log.post.term1 <- log.post.term1
     
    if(sd) log.post.term1.sd <- Bootstrap.Var(log.post.term1.out$logvals, B=1000)$var
    
    log.post.term2.denom.out <- ChibLogDenom(mcmc.out=mcmc.fixed.t, theta.star=theta.star, hyper=hyper, data=data, step.size=step.size, family=family)
    log.post.term2.denom <- log.post.term2.denom.out$estimate  
    if(sd) log.post.term2.denom.sd <- Bootstrap.Var(log.post.term2.denom.out$logvals, B=1000)$var
    
    details$log.post.term2.denom <- log.post.term2.denom
    
    log.post.term2 <- log.post.term2.num - log.post.term2.denom
    log.post <- log.post.term1 + log.post.term2
    
    details$log.post.term2 = log.post.term2
    details$log.post <- log.post  
    
    if(sd)  log.evidence.sd <-  log.post.term1.sd + log.post.term2.num.sd + log.post.term2.denom.sd
  }
  
  if(family=="binomial"){
    log.post.out <- EstimateLogPosteriorBinomial(data, theta.star, mcmc.out, family=family)
    log.post <- log.post.out$estimate
    if(sd) log.evidence.sd <- Bootstrap.Var(log.post.out$logvals, B=1000)$var
    details$log.post = log.post
  }
  
  log.like <- LogLike(data, theta.star, family=family)
  log.evidence <- log.like + log.prior - log.post
  
  details$log.like = log.like
  details$log.prior = log.prior
  details$log.evidence = log.evidence
  
  
  if(sd) return(list(log.evidence=log.evidence, sd=log.post.sd))
  else return(list(log.evidence=log.evidence, details=details))
}



CalcBF <-function(log.evidence){
  # a vector of  length 3 containing the log evidences.
  log.BF <- vector(length=choose(length(log.evidence), 2))
  diffs <- outer(log.evidence, log.evidence, FUN="-")
  BF.mat <- exp(diffs)
  BF <- c("mb"= BF.mat[2,1], "db"= BF.mat[3,1], "dm"= BF.mat[3,2])
  probH0 <-c("binomial"=1/(1+BF["mb"]+BF["db"]), "multbinom"=BF["mb"]/(1+BF["mb"]+BF["db"]), "doublebinom"=BF["db"]/(1+BF["mb"]+BF["db"])) 
  names(probH0) = c("binomial", "multbinom", "doublebinom")
  return(list(BF=BF, probH0=probH0))
}





