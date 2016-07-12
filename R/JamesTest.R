
## Returns the Meelis and R values for an experiment with different size clutches

James.test <-function(experiment, TwoSided=FALSE){
	## Experiment should be of the form, (clutch size, males)
## Note - two-sided test - p-values would be smaller if one-sided test used
  
	## START by removing experiments of size 0
	experiment <- experiment[experiment[,1]!=0,]

	clutch.sizes<-experiment[,1]
	data <- as.data.frame(cbind(clutch.sizes, experiment[,2])) ##reorganise data into (clutch size, male) table

  male <- data[,2]
  female <- data[,1] -data[,2]

  tmp<-table(data[,1], data[,2])  # organise into a table 
	sizes <-as.numeric(row.names(tmp))

  phat <- sum(male)/sum(data[,1])   ## proportion of males
	qhat <- sum(female)/sum(data[,1]) ## proportion of females
	
  K <- c()
  K <- 1/2 *(female*(female-1)/qhat^2 +male*(male-1)/phat^2 - 2*female*male/(phat*qhat))
  
  I <- data[,1]*(data[,1]-1)/(2*phat^2*qhat^2)
  
  sumI <- sum(I)
  sumK <- sum(K)
  
  
  U <-  sumK/sqrt(sumI)
  
	
	if(TwoSided)   p.val <- 2*pnorm(-abs(U))
	else  p.val <- pnorm(U)
	

return(list(U=U, p.val = p.val , exp.table = tmp))

}

### Can fail if all the clutch sizes end up with only 1 or 2 clutches with 1 male. If essentially there are no males, it will fail



#U<-c()
#ps<-c()

#for(i in 1:100){
#exp.size<-50
#	mortality <-0.5
#	av.clutch.size<-10
#	p<-0.05
#	psi <- 0.445

#	clutch.sizes<-0
#	while(!prod(clutch.sizes!=0)) {clutch.sizes <- rpois(exp.size, av.clutch.size) }#
#	tab <-table(clutch.sizes)	
#	tab.mat <-cbind(as.numeric(names(tab)), as.vector(tab))
#
#	males <-unlist(apply(tab.mat, 1, function(x) rmultbinom(x[2], size=x[1], prob=p, psi=psi)))
#	experiment.prim <-cbind(sort(clutch.sizes), males)
#	experiment.sec <- Binom.Mort(experiment.prim, mortality)
#	
#	J.out <- James.Test(experiment.sec)
#	U[i] <- J.out$U
 # ps[i] <- J.out$p.val
#
#}
#
#hist(U)
#print(sum(U< -1.96 | U>1.96))

#print(sum(ps<0.05))










