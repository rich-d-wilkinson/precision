\name{MCMCWithinGibbs}
\alias{MCMCWithinGibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrapper for the MCMC
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MCMCWithinGibbs(theta0, data, hyper, nbatch, family, step.size = NULL, keepNM = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta0}{Starting value for the Markov chain. Must be of length three for the binomial model, and length four for the multbinom or doublebinom models. The elements must also be named correctly. See example usage below.
%%     ~~Describe \code{theta0} here~~
}
  \item{data}{The data
%%     ~~Describe \code{data} here~~
}
  \item{hyper}{A list  containing the relevant hyperparameter information. hyper$a.m and hyper$b.m are the parameters in the Beta(a,b) prior on the mortality rate. \code{hyper$alpha.lambda} and \code{hyper$beta.lambda} are the parameters in a Gamma(alpha, beta) prior on lambda.  hyper$a.p and hyper$b.p are the parameters in the Beta(a,b) prior on the probability p. hyper$mu.psi and hyper$sd.psi are the mean and standard deviation parameters for the N(mu, sd) prior on psi.
%%     ~~Describe \code{hyper} here~~
}
  \item{nbatch}{Number of iterations in the MCMC
%%     ~~Describe \code{nbatch} here~~
}
  \item{family}{Either \code{binomial}, \code{multbinom}, or \code{doublebinom} depending on the model we wish to analyse.
%%     ~~Describe \code{family} here~~
}
  \item{step.size}{A vector of length 2 giving the step size for the Metropolis-Hastings algorithm. Only needed for the multbinom and doublebinom models. Needs to be name correctly (see example usage below)
%%     ~~Describe \code{step.size} here~~
}
  \item{keepNM}{logical; if TRUE, the Markov chains for N and M are returned.
%%     ~~Describe \code{keepNM} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{chain}{Output of the Markov chain}
\item{accep}{}
\item{acc.vec}{}
\item{N.chain}{The Markov chain output for the latent N values}
\item{M.chain}{The Markov chain output for the latent M values}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Richard Wilkinson
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#demo(MCMC, package="precision")

data(GlegneriSecondary)
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
##########################################
### MCMC parameters

nbatch <- 1*10^3   ## number of MCMC samples to generate - usually we'll require at least 10^5 iterations 
thin = FALSE    
burnin <- 10^4 
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MCMC}
