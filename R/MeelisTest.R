
## Returns the Meelis and R values for an experiment with different size clutches

Meelis.test <-function(experiment, TwoSided=FALSE){
	## Experiment should be of the form, (clutch size, males)

	## START by removing experiments of size 0
	experiment <- experiment[experiment[,1]!=0,]

	clutch.sizes<-experiment[,1]
	data <- as.data.frame(cbind(clutch.sizes, experiment[,2])) ##reorganise data into (clutch size, male) table

  
	tmp<-table(data[,1], data[,2])  # organise into a table 
	sizes <-as.numeric(row.names(tmp))

	if(sum(data[,2])==0) {
	  print("Warning: no males in this data set. R and Meelis test can not be applied")
	  list(vals=NA, R.av=NA, U.av=NA, p.av=NA, exp.table = tmp)
	  
	}
	
	out<-matrix(nrow=length(sizes), ncol=10)
	out[,1]<-sizes
	counter<-1

	for(i in sizes){
		n.clutches <- sum(clutch.sizes==i)  ## number of clutches of that size
		out[counter,2] <- n.clutches
		data.tmp<-data[clutch.sizes==i, ]
		s <-sum(data.tmp[,2])
		p.hat <- s/(i*n.clutches)
		out[counter,3]<-p.hat
		out[counter,4]<-i*p.hat*(1-p.hat)

		s_2 <- s*(s-1)
		s_3 <-s_2*(s-2)
		s_4 <-s_3*(s-3)
		N <- dim(experiment)[1]	
	
		M <- s*(s*(i-1)+i*(n.clutches-1))/(n.clutches*i-1)
	
		V<-s_4*(i-1)*(n.clutches*i*(i-1)-4*i+6)/((n.clutches*i-1)*(n.clutches*i-2)*(n.clutches*i-3)) + 4* s_3*(i-1)*(i-2)/((n.clutches*i-1)*(n.clutches*i-2)) + 2*s_2*(i-1)/(n.clutches*i-1) - s^2*(s-1)^2*(i-1)^2/(n.clutches*i-1)^2
		
		U <- (sum(data.tmp[,2]^2)-M)/sqrt(V)
		out[counter, 7] <-M
		out[counter, 8] <-V
		out[counter, 9] <-U
		out[counter, 10] <-pnorm(U)



		
		if(n.clutches >1){
			out[counter, 5] <- var(data.tmp[,2])
			out[counter, 6] <-out[counter,5]/out[counter,4]
			}
		else {
			out[counter, 5] <- 0
			}
		counter <-counter+1
	}

	colnames(out)<-c("clutch size", "no. clutches", "p hat", "binom var", "obs var", "R", "M", "V", "U", "p value")

	R.av <- sum(out[,2]*out[,5])/sum(out[,2]*out[,4])
	U.inc <-out[!is.na(out[,9]),9]
	U.av <- sum(U.inc)/sqrt(length(U.inc))
  s2 <- McCullagh(experiment)
  if(TwoSided) p.av <- 2*pnorm(-abs(U.av))
  else  p.av <- pnorm(U.av)



return(list(vals=out, R.av=R.av, s2=s2, U.av=U.av, p.av=p.av, exp.table = tmp))

}

### Can fail if all the clutch sizes end up with only 1 or 2 clutches with 1 male. If essentially there are no males, it will fail



McCullagh <- function(data){
  if(!isTRUE(all.equal(colnames(data), c('n', 'm')))){
    print('Warning - relabelling columns')
    colnames(data) <- c('n', 'm') 
  }
  phat <- sum(data[,'m'])/sum(data[,'n'])
  exp.size = dim(data)[1]
  num <- (data[,'m'] - phat*data[,'n'])^2
  denom <- data[,'n']*phat*(1-phat)
  sigma2 <- sum(num/denom)*1/(exp.size-1)
  return(sigma2)
}



#	exp.size<-50
#	mortality <-0.7
#	av.clutch.size<-10
#	p<-0.002
#	psi <- 0.445

#	clutch.sizes<-0
#	while(!prod(clutch.sizes!=0)) {clutch.sizes <- rpois(exp.size, av.clutch.size) }#
#	tab <-table(clutch.sizes)	
#	tab.mat <-cbind(as.numeric(names(tab)), as.vector(tab))

#	males <-unlist(apply(tab.mat, 1, function(x) rMultbinom(x[2], size=x[1], prob=p, psi=psi)))
#	experiment.prim <-cbind(sort(clutch.sizes), males)
#	experiment.sec <- Binom.Mort(experiment.prim, mortality)
	
#	M.out <- Meelis.t(experiment.sec)
#	M.out$U.av

















