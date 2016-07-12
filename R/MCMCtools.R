###################################################################################
###
###
### Functions needed for the MCMC. Contains
###   - ProposeNM
###   - AcceptNM
###   - UpdateNM
###
###################################################################################












#########################################################################
### 
### ProposeNM()
###
### Given the data, it simulates a matrix of matching dimension with proposed values of NM according to the specified 
### probability model. It simulates proposals conditional on N>=n, M>=m and M-m<=N-n, so all proposals should be feasible.
###
###
#########################################################################

ProposeNM <- function(theta, data, family){
  ## we simulate from the prior pi(M|N, p)pi(N|lambda) conditioned on M>=m, N>=n, M<=N-n+m
  # add an option for the family, and make one function
  if(family=="binomial")  out <- rNM_given_nm(n=1, lambda=theta["lambda"], n.data=data[,"n"], p=theta["p"], m=data[,"m"], family=family)
  if(family=="multbinom" | family=="doublebinom") out <- rNM_given_nm(n=1, lambda=theta["lambda"], n.data=data[,"n"], p=theta["p"], psi=theta["psi"], m=data[,"m"], family=family)
  N <- out$N
  M <- out$M
  
  return(cbind("N"=N,"M"=M))
}


Logit <-function(x){
  log(x/{1-x})
}


InvLogit <-function(x){
  y<- exp(x)
  y/{1+y}
}



ProposePPsi <- function(theta, data, step.size){
  p.logit <- Logit(theta["p"])
  p.logit.new <- p.logit + rnorm(n=1, mean=0, sd=step.size["p.logit"])
  psi.new <- theta["psi"] + rnorm(n=1, mean=0, sd=step.size["psi"])
  theta.new <-theta
  theta.new["p"] <- InvLogit(p.logit.new)
  theta.new["psi"] <- psi.new	
  return(theta.new)	
}







#########################################################################
### 
### AcceptNM()
###
### Given three vectors of matching dimension, it gives calculates the probability of accepting each NM proposal, 
###  and does the updating. It returns a list with the updated NM matrix, the acceptance probabilities, and whether the move was accepted or not
###
###  WARNING - It currently has the acceptance probability for identical proposals artificially set to 0, when it should be 1
###          - this is so that the acceptance rate isn't inflated and so that difficult data points can be spotted
###          - it makes no difference to the MCMC algorithm, but it will affect Chib's method where we need to calculate 
###            the acceptance rate. It will need to be changed then I think (add a T/F option for Chib or just remove)
###
#########################################################################

AcceptNM <-function(theta, data, NM.cur, NM.prop){
  ## is the same for all versions of the model (at least with binomial mortality)
  ## move into generic position
  
  # we don't need checks that NM are feasible as the proposal should only provide feasible values
  
  stopifnot(dim(data) == dim(NM.cur), dim(NM.cur) == dim(NM.prop))
  tmp <- cbind(data, NM.cur, NM.prop)
  colnames(tmp) <- c("n", "m", "N.cur", "M.cur", "N.prop", "M.prop")
  
  acc.prob <- apply(tmp, 1, function(x){  
    if(x["N.prop"]==x["N.cur"] & x["M.prop"]==x["M.cur"]) return(0)   ## i.e. don't count identical moves as acceptances
    prob.m <- dm_given_nNM(n=x["n"], m=x["m"], N=x["N.prop"], M=x["M.prop"]) / dm_given_nNM(n=x["n"], m=x["m"], N=x["N.cur"], M=x["M.cur"]) 
    prob.n <- dbinom(x=x["n"], size=x["N.prop"], prob= 1-theta["mort"]) / dbinom(x=x["n"], size=x["N.cur"], prob= 1-theta["mort"])
    return (prob.m * prob.n)
  })
  
  U <- runif(dim(data)[1])
  acc <- (U < acc.prob)
  update <- which(acc)
  
  NM.cur[update, ] <- NM.prop[update, ]   ## just update successful proposals
  
  return(list(NM=NM.cur, prob=acc.prob, acc=acc))
}

#################################################################################################
###
###
### CreateMultbinomLookup
### - given a value of p and psi, and a max.size, it returns a matrix with dimension max.size * max.size+1
### - row n gives the complete density corresponding to size n for x=0, 1, ... n
###
#################################################################################################



CreateMultbinomLookup <-function(max.size, p, psi, family){
  
  # Each row contains all the possible log density values 
  out <- matrix(nrow=max.size, ncol = max.size+1)
  colnames(out)<-paste("x",0:max.size, sep="")
  rownames(out) <- paste("size",1:max.size, sep="")
  if(family =="multbinom"){
    for(i in 1:max.size){
      out[i,1:(i+1)] <- dmultbinom(size=i, p=p, psi=psi, log=TRUE)
    }}
  if(family =="doublebinom"){
    for(i in 1:max.size){
      out[i,1:(i+1)] <- ddoublebinom(size=i, p=p, psi=psi, log=TRUE)
    }}
  return(out)
}




AcceptPPsi <-function(theta.cur, theta.prop, NM, hyper, family){
  # work with the log
  ## note - this does not depend upon the data
  max.size <-max(NM[,"N"])
  logden.prop <- CreateMultbinomLookup(max.size=max.size, theta.prop["p"], theta.prop["psi"], family=family)
  
  logden.cur <- CreateMultbinomLookup(max.size=max.size, theta.cur["p"], theta.cur["psi"], family=family)
  
  log.prior.ratio.p <- dbeta(theta.prop["p"], shape1=hyper$a.p, shape2=hyper$b.p, log=TRUE) - dbeta(theta.cur["p"], shape1=hyper$a.p, shape2=hyper$b.p, log=TRUE)
  
  log.prior.ratio.psi <- dnorm(theta.prop["psi"], mean=hyper$mu.psi, sd=hyper$sd.psi, log=TRUE) - dnorm(theta.cur["psi"], mean=hyper$mu.psi, sd=hyper$sd.psi, log=TRUE)
  
  
  ## use the look up tables to assign the correct log likelihood
  loglik.num <- sum(apply(NM,1, function(x) logden.prop[x[1], x[2]+1]))    ## +1 needed in the column because of R's indexing from 1
  loglik.denom <- sum(apply(NM,1, function(x) logden.cur[x[1], x[2]+1])) 
  
  log.acc.prob <- log.prior.ratio.p + log.prior.ratio.psi + loglik.num - loglik.denom
  
  log.U <- log(runif(1))
  acc <- as.numeric(log.U < log.acc.prob)
  if(acc==1) theta=theta.prop
  else theta=theta.cur
  return(list(theta=theta, log.prob=min(0,log.acc.prob),acc=acc))
}


#########################################################################
### 
### UpdateNM()
###
### Given two matrices of matching dimension containing the data and current NM value, the function updates NM
### according to the specified probability model and returns a list containing the new NM value, the acceptance rate, and a list
### of which values were accepted and which weren't.
###
###
### Other updates are Gibbs updates for the parameters
###
###
#########################################################################

UpdateNM <-function(theta, data, NM.cur,   family)
{
  NM.prop <- ProposeNM(theta, data, family=family)
  acc.out <- AcceptNM(theta, data, NM.cur, NM.prop)
  
  NM <- acc.out$NM
  return(list(NM=NM, accept=sum(acc.out$acc)/dim(data)[1], acc.vec=acc.out$acc))
}


Update.p.binomial<-function(current, hyper){
  rbeta(1, shape1 = sum(current$NM[,"M"]) + hyper$a.p,   shape2 = sum(current$NM[,"N"]-current$NM[,"M"]) + hyper$b.p) 		
}

Update.mort<-function(current, hyper, n.sum){
  rbeta(1, shape1 = sum(current$NM[,"N"]) - n.sum + hyper$a.m,   shape2 = n.sum + hyper$b.m) 		
}

Update.lambda<-function(current, hyper, exp.size){
  rgamma(1, shape = hyper$alpha.lambda+ sum(current$NM[,"N"]),   rate =  hyper$beta.lambda+exp.size) 		
}



UpdatePPsi <-function(current, data,  hyper, step.size=step.size, family){
 
  theta.prop <- ProposePPsi(theta=current$theta, data=data, step.size=step.size)
  acc.out <- AcceptPPsi(theta.cur=current$theta, theta.prop=theta.prop, NM=current$NM, hyper=hyper,  family)
  
  return(acc.out$theta)
}














UpdateTheta <- function(current, data, hyper, n.sum, exp.size, step.size=NULL, family){
  
  if(family=="binomial"){
 
    current$theta["p"] <- Update.p.binomial(current, hyper)
    current$theta["mort"] <- Update.mort(current, hyper, n.sum)
    current$theta["lambda"] <- Update.lambda(current, hyper, exp.size)    
    }
  
  if(family=="multbinom" | family=="doublebinom"){
    current$theta["mort"] <- Update.mort(current, hyper, n.sum)
    current$theta["lambda"] <- Update.lambda(current, hyper, exp.size)
    if(step.size[1]>0) current$theta <- UpdatePPsi(current, data,  hyper, step.size=step.size, family)
  }
  
  return(current$theta)
}


#########################################################################
### 
### MCMCWithinGibbs()
###
### - start should be a list, with an element theta with parameter names, 
### and a matrix NM with dimension matching data with labelled columns.
### - The function does MH within Gibbs, with Gibbs steps where possible for the parameters
### - It returns three matrices, with the theta chains, the N chain, and the M chain, plus some information on the acceptance rate
### - To run the chain with p and psi fixed (needed for the Chib method), set step.size=c(0,0)
### MCMCWrapper should be used to call MCMCWithinGibbs
###
#########################################################################




### start$theta must be named correctly for the function to work
MCMCWithinGibbs<-function( theta0,  data, hyper, nbatch, family, step.size=NULL, keepNM=FALSE){
  
  stopifnot(family=="binomial" | family=="multbinom" | family=="doublebinom", dim(data)[2]==2)
  if(family=="binomial") stopifnot(length(theta0)==3, is.null(step.size))
  if(family=="multbinom" | family =="doublebinom") stopifnot(length(theta0)==4, !is.null(step.size))
  
  exp.size <- nrow(data)
  n.sum <- sum(data[,1])
  
  NM <- round(data/(1-theta0["mort"]))   ## choose a starting point for the latent variables
  colnames(NM) <-c("N", "M")
  
  param.chain <- matrix(nrow=nbatch, ncol=length(theta0))
  param.chain[1,] <- theta0
  colnames(param.chain) <- names(theta0)
  
  current <-list(theta=theta0, NM=NM)   ## current is a list with a theta value and a NM matrix
  accep<-c()
  acc.vec <- rep(0, exp.size)
  
    if(keepNM){
    N.chain <- matrix(nrow=nbatch, ncol=dim(data)[1])
    M.chain <- matrix(nrow=nbatch, ncol=dim(data)[1])
    N.chain[1,] <- NM[,"N"]
    M.chain[1,] <- NM[,"M"]
  }
  
  for(i in 1:(nbatch-1)){
    param.chain[i+1,] <- UpdateTheta(current=current, data=data, hyper=hyper, n.sum=n.sum, exp.size=exp.size, step.size=step.size, family=family)
    current$theta=param.chain[i+1,]
    NM.new <- UpdateNM(theta=current$theta, data=data, NM.cur=current$NM,   family=family) ## use the new theta.
    current$NM <- NM.new$NM
    
#    print(i)
    if(keepNM){
      N.chain[i+1,] <- current$NM[,"N"]
      M.chain[i+1,] <- current$NM[,"M"]
    }
    
    accep[i+1] <- NM.new$accept
    acc.vec <- acc.vec + NM.new$acc.vec
    if(floor(i/1000)==i/1000) print(paste("i=",i))
  }	
  if(keepNM)  { return(list(chain=param.chain, accep=accep, acc.vec=acc.vec/nbatch, N.chain=N.chain, M.chain=M.chain))	}
  else{	return(list(chain=param.chain, accep=accep, acc.vec=acc.vec/nbatch))	}
}



### start$theta must be named correctly for the function to work
MCMCWithinGibbsNoMort<-function( start,  data, hyper, nbatch, family, step.size=NULL){
  
  stopifnot(family=="binomial" | family=="multbinom" | family=="doublebinom", dim(data)[2]==2)
  if(family=="binomial") stopifnot(length(start$theta)==1, is.null(step.size))
  if(family=="multbinom" | family =="doublebinom") stopifnot(length(start$theta)==2, !is.null(step.size))
  
  if(family=="binomial") print("WARNING: we can do everything exactly for the binomial model so we don't need to run MCMC")
  stopifnot(start$NM==data)
            
  exp.size <- nrow(data)
  n.sum <- sum(data[,1])
  
  param.chain <- matrix(nrow=nbatch, ncol=length(start$theta))
  param.chain[1,] <- start$theta
  colnames(param.chain) <- names(start$theta)
  
  current <- start   ## current is a list with a theta value and a NM matrix
  accep<-c()
  acc.vec <- rep(0, exp.size)
  
  
  for(i in 1:(nbatch-1)){
  
    if(family=="binomial") param.chain[i+1,] <- Update.p.binomial(current=current, hyper=hyper)
    if(family=="multbinom" | family=="doublebinom") param.chain[i+1,] <- UpdatePPsi(current=current, data=data, hyper=hyper, step.size=step.size, family=family)
    current$theta=param.chain[i+1,]
    
    if(floor(i/1000)==i/1000) print(paste("i=",i))
  }  
  return(list(chain=param.chain, accep=accep))
}





ThinChain<-function(mcmc.out, thinby, burnin){

    thin.smp <- seq(from=burnin, to=length(mcmc.out$chain[,1]), by=thinby)
    mcmc.out$chain <- mcmc.out$chain[thin.smp,]
    mcmc.out$N.chain <-  mcmc.out$N.chain[thin.smp,]
    mcmc.out$M.chain <- mcmc.out$M.chain[thin.smp,]

    
    return(mcmc.out)
}


