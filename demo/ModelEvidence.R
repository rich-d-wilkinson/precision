#########################################################################
###
### Demonstration of how to calculate the model evidence for secondary data
###
###
########################################################################

library(precision)

#Load some data
data(GlegneriSecondary)
meelis.out <- Meelis.test(GlegneriSecondary, TwoSided=TRUE)
james.out = James.test(GlegneriSecondary, TwoSided=TRUE)

###########################################################
#### Define hyperparameters/priors
####  - the four built in secondary datasets have priors for d and lambda built in, which load when you load the data
####  If you are using a different datasets you'll need to run the commands
# hyper<-list()
## prior for the mortality rate is Gamma(alpha, beta)
# hyper$a.m <- 6  # replace the value with something sensible for your data
# hyper$b.m <-  10

## prior for average clutch size, lambda, is Gamma(alpha, beta)
# hyper$alpha.lambda <- 4
# hyper$beta.lambda  <- 1

## prior for p is beta(a, b)
hyper$a.p <-  1 ## 20   ## Hyper parameters for p's prior distribution
hyper$b.p <-  1##100

# prior for psi - assumed to be Gaussian
hyper$mu.psi <- 0  
hyper$sd.psi <- 1


plot.prior(hyper=hyper, show=TRUE, family="multbinom")
#plot.prior(hyper=hyper, show=FALSE, filename="Priors.pdf", family="multbinom") # save plot rather than printing it to screen

##########################################
### MCMC parameters

nbatch <- 1*10^3   ## number of MCMC samples to generate - usually we'll require at least 10^5 iterations 
thin = FALSE    
burnin <- 10^2 
thinby <- 3  


#######################################################################################################
# BINOMIAL MODEL
#######################################################################################################

## Define a start point for the MCMC chain - important to name the parameter and the column of the NM matrix
b.theta0 <-c(10, 0.1, 0.5) 
names(b.theta0) <-c("lambda", "p", "mort")

b.mcmc.out <- MCMCWithinGibbs( theta0=b.theta0,  data=GlegneriSecondary, hyper=hyper, nbatch=nbatch, family="binomial", keepNM=TRUE)


## Thin the chain or not
if(thin){
  b.mcmc.out.t  <- ThinChain(b.mcmc.out, thinby=thinby, burnin=burnin)
}
if(!thin){
  b.mcmc.out.t <- b.mcmc.out
}
rm(b.mcmc.out)  ## delete as not needed and takes a lot of memory


plot.posterior(chain=b.mcmc.out.t$chain, hyper=hyper, show=TRUE, family="binomial", theta.true=NULL)
plot.trace(chain=b.mcmc.out.t$chain, show=TRUE, family="binomial")




#######################################################################################################
# MULTIPLICATIVE BINOMIAL MODEL
#######################################################################################################

## Define the Metropolis-Hastings random walk step size
m.step.size<-c(0.3,0.2)   ## 0.3 and 0.2 aree 
names(m.step.size) <- c("p.logit", "psi")


m.theta0 <-c(10, 0.1, 0, 0.1) 
names(m.theta0) <-c("lambda", "p", "psi", "mort")

m.mcmc.out <- MCMCWithinGibbs( theta0=m.theta0,  data=GlegneriSecondary, hyper=hyper, nbatch=nbatch,  family="multbinom", step.size=m.step.size, keepNM=TRUE)

if(thin){
  m.mcmc.out.t  <- ThinChain(m.mcmc.out, thinby=thinby, burnin=burnin)
}
if(!thin){
  m.mcmc.out.t <- m.mcmc.out
}
rm(m.mcmc.out)   ## to save space


plot.posterior(chain=m.mcmc.out.t$chain, hyper=hyper, show=TRUE, family="multbinom", theta.true=NULL)
plot.trace(chain=m.mcmc.out.t$chain, show=TRUE,   family="multbinom")







#######################################################################################################
# DOUBLE BINOMIAL MODEL
#######################################################################################################

d.step.size<-c(0.3,0.2)   ## 0.3 and 0.2 aree 
names(d.step.size) <- c("p.logit", "psi")

d.theta0 <-c(10, 0.1, 0, 0.1) 
names(d.theta0) <-c("lambda", "p", "psi", "mort")

d.mcmc.out <- MCMCWithinGibbs( theta0=d.theta0,  data=GlegneriSecondary, hyper=hyper, nbatch=nbatch,  family="doublebinom", step.size=d.step.size, keepNM=TRUE)

if(thin){
  d.mcmc.out.t  <- ThinChain(d.mcmc.out, thinby=thinby, burnin=burnin)
}
if(!thin){
  d.mcmc.out.t <- d.mcmc.out
}
rm(d.mcmc.out)   ## to save space


plot.posterior(chain=d.mcmc.out.t$chain, hyper=hyper, show=TRUE,   family="doublebinom", theta.true=NULL)
plot.trace(chain=d.mcmc.out.t$chain, show=TRUE,   family="doublebinom")



############################## Estimate Model Evidences





b.log.evidence <- CalculateEvidence(mcmc.out=b.mcmc.out.t, data=GlegneriSecondary,  hyper=hyper, family="binomial", sd=FALSE)
m.log.evidence <- CalculateEvidence(mcmc.out=m.mcmc.out.t, data=GlegneriSecondary,  hyper=hyper, family="multbinom", sd=FALSE, nbatch=nbatch, step.size=m.step.size, thin=thin, thinby=thinby, burnin=burnin)
d.log.evidence <- CalculateEvidence(mcmc.out=d.mcmc.out.t, data=GlegneriSecondary,  hyper=hyper, family="doublebinom", sd=FALSE, nbatch=nbatch, step.size=d.step.size, thin=thin, thinby=thinby, burnin=burnin)

## For the  multiplicative and double models, we have to run an MCMC with theta fixed, so it will take longer

log.evidence <- c(b.log.evidence$log.evidence, m.log.evidence$log.evidence, d.log.evidence$log.evidence)


BF<-CalcBF(log.evidence)

chib.out <- list(BF=BF$BF, probH0 = BF$probH0 , 
                 ProbPosPsi = c("multbinom"=sum((m.mcmc.out.t$chain[,"psi"]>0))/length(m.mcmc.out.t$chain[,"psi"]), "doublebinom"=sum((d.mcmc.out.t$chain[,"psi"]>0))/length(d.mcmc.out.t$chain[,"psi"])),
                 log.BF=log(BF$BF), log.evidence=log.evidence, 
                 R= c("R"=meelis.out$R.av),
                 meelis = c("U"=meelis.out$U.av, "p"=meelis.out$p.av,  "conclusion"=ifelse(meelis.out$p.av<0.05, "RejectH0", "AcceptH0")),
                 james = c("U"=james.out$U, "p"=james.out$p.val,  "conclusion" = ifelse(james.out$p.val<0.05, "RejectH0", "AcceptH0") ) )   

print(chib.out)




