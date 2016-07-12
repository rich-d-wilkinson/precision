
#include "random.h"
#include "factorial.h"


// should only be called by C functions. Not for use from R
void dDoubleBinomC_internal(int size, double p, double psi,  double *den){
  // returns entire density

  int i;
  double unnorm_logden[size+1];
  double unnorm_den[size+1];
	double sum=0;
	double log_norm_const, norm_const;

  // i==0
  unnorm_logden[0] = size*psi*log(size)+log(1-p)*size*(psi+1)-log(size)*psi*size;
  unnorm_den[0] = exp(unnorm_logden[0]);
  sum = unnorm_den[0];
  // i==size
  unnorm_logden[size] = size*psi*log(size)+size*(psi+1)*log(p)-log(size)*size*psi;
  unnorm_den[size] = exp(unnorm_logden[size]);
  sum = sum + unnorm_den[size];  
  
  
	for(i=1; i<size; i++){
	      unnorm_logden[i] = log(Binomial_Coefficient(size, i)) + size*psi*log(size) + i*(psi+1)*log(p) + log(1-p)*(size-i)*(psi+1) - log(i)*i*psi - log(size-i)*psi*(size-i);
	      unnorm_den[i] = exp(unnorm_logden[i]);
	      sum = sum + unnorm_den[i];	
	}
	
	norm_const = 1/sum;
	
	for(i=0; i<=size; i++){
	      den[i] = unnorm_den[i] * norm_const;
  }
// return(den);
}



/*
 Simulates *n random variables from the multiplicative binomial distribution.
*/

void rDoublebinomC(int *n, 
                 int *size_R, 
                 double *p_R, 
                 double *psi_R, 
                 int *rndm,  
                 int *seed, 
                 double *den){

      double cumden[*size_R+1];
      int i;	
      int j=0;
      double U;
      int seeds = seed[0];
      sdprand(seeds);  //seed the random number generator


	   dDoubleBinomC_internal(*size_R, *p_R, *psi_R, den);
	    cumden[0] =den[0];
	
	for(i=1; i<= (*size_R); i++){
	    cumden[i] = cumden[i-1]+den[i];
	}


	for(j=0; j<*n;j++){
      i=0;
      U = dprand();
	    while(U>cumden[i]) i++;
	    rndm[j]=i;
	}

}


int rDoublebinom_internal(int size, 
                 double p, 
                 double psi){

      double cumden[size+1];
      int i;  
      double U;
      int rndm;
      double den[size+1];

    	dDoubleBinomC_internal(size, p, psi, den);
	    cumden[0] =den[0];
	
    	for(i=1; i<= size; i++){
	        cumden[i] = cumden[i-1]+den[i];
	    }


      i=0;
      U = dprand();
	    while(U>cumden[i]) i++;
	    rndm=i;
	    
      return(rndm);
}


