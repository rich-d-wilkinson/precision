####
####
####  Two functions for simulating mortality and multiple experiments for given parameters
####
#### 
####

Simulate.Data <- function(theta, exp.size, family="binomial"){
  
  clutch.sizes<-0
  while(!prod(clutch.sizes!=0)) {clutch.sizes <- rpois(exp.size, theta["lambda"]) }
  
  
  if(family=="binomial"){
    males <-sapply(clutch.sizes, function(x) rbinom(1, size=x, prob=theta["p"]))
    }
  if(family=="multbinom"){
    males <-sapply(clutch.sizes, function(x) rmultbinom(1, size=x, prob=theta["p"], psi=theta["psi"]))
    }
  if(family=="doublebinom"){
    males <-sapply(clutch.sizes, function(x) rdoublebinom(1, size=x, prob=theta["p"], psi=theta["psi"]))
  }
  experiment.prim <-cbind(clutch.sizes, males)
  colnames(experiment.prim) <-c("N", "M")
  experiment.sec <- Binom.Mort(experiment.prim, theta["mort"])
  
  data <- experiment.sec
  colnames(data)<-c("n", "m")
  
  return(data)
  
}

  
  

## Binomial mortality
Binom.Mort <- function(primary, mort){
	## Takes a primary data set, and returns it after being affected by binomial mortality.
	## primary of the form (clutch size, no. males)
	
	tmp<-c()
	
	secondary<-apply(primary, 1, 
		function(x) { 
			deaths <-rbinom(n=1, size=x[1], prob=mort)
			temp<-sample(1:x[1], size=deaths, replace=FALSE)
			m.deaths <- sum(temp <=x[2])
			tmp[1]<-x[1]-deaths
			tmp[2]<-x[2]-m.deaths
			return(tmp)
		}
	)			

return(t(secondary))
}



Experiments.MultBinom <- function(exp.size, mortality,  no.experiments, signif.level, param.multbin, av.clutch.size){

	####	Estimates the rejection rate for a given experiment size, mortality, and 
	####	significance level. no.experiments controls the number of MC replicates
	#### 	Assumes multiplicative binomial sex-ratios, and binomial mortality.

	p <-param.multbin[1]
	psi <-param.multbin[2]
	

	rejects<-0
	fail<-0	

	for(i in 1: no.experiments){
	
	clutch.sizes<-0
	while(!prod(clutch.sizes!=0)) {clutch.sizes <- rpois(exp.size, av.clutch.size) }
	tab <-table(clutch.sizes)	
	tab.mat <-cbind(as.numeric(names(tab)), as.vector(tab))

	males <-unlist(apply(tab.mat, 1, function(x) rmultbinom(x[2], size=x[1], prob=p, psi=psi)))
	experiment.prim <-cbind(sort(clutch.sizes), males)
	experiment.sec <- Binom.Mort(experiment.prim, mortality)
	
	M.out <- Meelis.test(experiment.sec)

	#print(table(experiment.sec[,1], experiment.sec[,2]))
	#print(M.out)

		if(is.na(M.out$p.av)) {fail <-fail+1}
		else{	if (M.out$p.av <= signif.level) rejects <- rejects+1}
	}
	print(paste("fail=", fail))
	return(list(rejects=rejects, prob.rej=rejects/(no.experiments-fail)))

}

