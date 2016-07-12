#########################################################################################
##
##
## Plotting tools
##
##
#########################################################################################


plot.prior <- function(hyper, show, filename=NULL, family="binomial"){
  stopifnot(!(!show &is.null(filename)))
  if(!show) {pdf(file=paste(filename, "priors.pdf", sep=""))}

  if(family=="binomial")  {
    par(mfrow=c(1,3))
    x<-seq(0,1, 0.01)
    y<-dbeta(x, shape1=hyper$a.p, shape2=hyper$b.p)
    plot(x,y, lty=2, col=2, type="l", main="p", sub=paste("a.p=", hyper$a.p, "b.p.=", hyper$b.p))


    x<-seq(0,1, 0.01)
    y<-dbeta(x, shape1=hyper$a.m, shape2=hyper$b.m)
    plot(x,y, lty=2, col=2, type="l",  main="mort",  sub=paste("a.m=", hyper$a.m, "b.m=", hyper$b.m))


    x<-seq(0,30, 0.1)
    y<-dgamma(x, shape=hyper$alpha.lambda, rate=hyper$beta.lambda)
    plot(x,y, lty=2, col=2, type="l", main="lambda", sub=paste("alpha=", hyper$alpha.lambda, "beta=", hyper$beta.lambda))
  }
  if(family=="multbinom" | family =="doublebinom"){
    
    par(mfrow=c(2,2))
    x<-seq(0,1, 0.01)
    y<-dbeta(x, shape1=hyper$a.p, shape2=hyper$b.p)
    plot(x,y, lty=2, col=2, type="l", main="p", sub=paste("a.p=", hyper$a.p, "b.p.=", hyper$b.p))
    
    
    x<-seq(0,1, 0.01)
    y<-dbeta(x, shape1=hyper$a.m, shape2=hyper$b.m)
    plot(x,y, lty=2, col=2, type="l",  main="mort",  sub=paste("a.m=", hyper$a.m, "b.m=", hyper$b.m))
    
    
    x<-seq(-3,3, 0.01)
    y<-dnorm(x, mean=hyper$mu.psi, sd=hyper$sd.psi)
    plot(x,y, lty=2, col=2, type="l",  main="psi", sub=paste("mu=", hyper$mu.psi, "sd=", hyper$sd.psi))
    
    x<-seq(0,30, 0.1)
    y<-dgamma(x, shape=hyper$alpha.lambda, rate=hyper$beta.lambda)
    plot(x,y, lty=2, col=2, type="l", main="lambda", sub=paste("alpha=", hyper$alpha.lambda, "beta=", hyper$beta.lambda))

  }
  if(!show) dev.off()
  if(show) {dev.copy2pdf(file=paste(filename, "priors.pdf", sep=""))}
}






#####
### if theta.true!=NULL then it plots the value on the graphs


plot.posterior <- function(chain, hyper, show,  filename=NULL, family="binomial", theta.true=NULL){
  stopifnot(!(!show &is.null(filename)))
  if(!show)  pdf(file=paste(filename, "posteriors.pdf", sep=""))

  if(family=="binomial")  par(mfrow=c(1,3))
  if(family=="multbinom" | family =="doublebinom") par(mfrow=c(2,2))

  plot(density(chain[,"lambda"]), main="Posterior lambda", lwd=2, sub=paste("alpha=", hyper$alpha.lambda, ", beta=",hyper$beta.lambda, sep=""))
  x<-seq(0,30, 0.1)
  y<-dgamma(x, shape=hyper$alpha.lambda, rate=hyper$beta.lambda)
  lines(x,y, lty=2, col=1)
  if(!is.null(theta.true)) lines(c(theta.true["lambda"],theta.true["lambda"]), c(0,10), lty=3, lwd=2, col=3)

  plot(density(chain[,"p"]),  main="Posterior p", lwd=2, sub=paste("a=", hyper$a.p, ", b=",hyper$b.p, sep=""))
  x<-seq(0,1, 0.01)
  y<-dbeta(x, shape1=hyper$a.p, shape2=hyper$b.p)
  lines(x,y, lty=2, col=1)
  if(!is.null(theta.true))  lines(c(theta.true["p"],theta.true["p"]), c(0,30), lty=3, lwd=2, col=3)


  if(family=="multbinom" | family =="doublebinom"){
    if(family=="multbinom")   plot(density(chain[,"psi"]),  main="Posterior psi - multbinom", lwd=2, , sub=paste("mu=", hyper$mu.psi, ", sd=",hyper$sd.psi, sep=""))
    if(family=="doublebinom")   plot(density(chain[,"psi"]),  main="Posterior psi - doublebinom", lwd=2, , sub=paste("mu=", hyper$mu.psi, ", sd=",hyper$sd.psi, sep=""))
    
    x<-seq(-3,3, 0.01)
    y<-dnorm(x, mean=hyper$mu.psi, sd=hyper$sd.psi)
    lines(x,y, lty=2, col=1)
    if(!is.null(theta.true)) lines(c(theta.true["psi"],theta.true["psi"]), c(0,30), lty=3, lwd=2, col=3)
  }



  plot(density(chain[,"mort"]),  main="Posterior mort", lwd=2, , sub=paste("a=", hyper$a.m, ", b=",hyper$b.m, sep=""))
  x<-seq(0,1, 0.01)
  y<-dbeta(x, shape1=hyper$a.m, shape2=hyper$b.m)
  lines(x,y, lty=2, col=1)
  if(!is.null(theta.true))  lines(c(theta.true["mort"],theta.true["mort"]), c(0,30), lty=3, lwd=2, col=3)

  if(!show) dev.off()

  if(show) {dev.copy2pdf(file=paste(filename, "posteriors.pdf",sep=""))}
}




plot.trace <- function(chain, show,  filename=NULL, family="binomial"){
  stopifnot(!(!show &is.null(filename)))
  if(!show)  pdf(file=paste(filename, "trace.pdf", sep=""))
  if(family=="binomial")  par(mfrow=c(1,3))
  if(family=="multbinom" | family =="doublebinom") par(mfrow=c(2,2))
  
  plot(chain[,"lambda"], type="l", main="lambda")
  plot(chain[,"p"], type="l", main="p")
  if(family=="multbinom" | family =="doublebinom")   plot(chain[,"psi"], type="l", main="psi")
  plot(chain[,"mort"], type="l", main="mort")
  
  if(!show) dev.off()
  if(show) {dev.copy2pdf(file=paste(filename,"trace.pdf", sep=""))}
  
}




