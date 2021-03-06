
#include "random.h"
#include "factorial.h"


// should only be called by C functions. Not for use from R
void dMultBinomC_internal(int size, double p, double psi,  double *den){
  // returns entire density

  int i;
  double unnorm_logden[size+1];
	double unnorm_den[size+1];
	double sum=0;
	double log_norm_const, norm_const;

	for(i=0; i<=size; i++){
	      unnorm_logden[i] = log(Binomial_Coefficient(size, i)) + i*log(p) + log(1-p) * (size-i) +  psi * i * (size-i);
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

void rMultbinomC(int *n, 
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


	   dMultBinomC_internal(*size_R, *p_R, *psi_R, den);
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


int rMultbinom_internal(int size, 
                 double p, 
                 double psi){

      double cumden[size+1];
      int i;  
      double U;
      int rndm;
      double den[size+1];

    	dMultBinomC_internal(size, p, psi, den);
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


/*
 Simulates from the multiplicative binomial distribution conditional on the value being greater than or equal to lb. 

Possibly no longer needed - delete if so
*/
void rMultbinomC_lb(int *n,
                    int *size_R, 
                    double *p_R, 
                    double *psi_R, 
                    int *lb, 
                    int *rndm, 
                    int *seed, 
                    double *den){

      double cumden[*size_R+1];
      int i;  
      int j=0;
      double U;
      int seeds = seed[0];
      
      
      sdprand(seeds);  //seed the random number generator

	
      dMultBinomC_internal(*size_R, *p_R, *psi_R, den=den);
    	cumden[0] =den[0];

	
    	for(i=1; i<= (*size_R); i++){
    	    cumden[i] = cumden[i-1]+den[i];
    	}
  
  	  for(j=0; j<*n;j++){
        R_CheckUserInterrupt();  // allows R to interupt  with Ctrl-C.
      	rndm[j]= -1;
    	  while(rndm[j]<*lb)
      	 {
	          i=0;
	          U = dprand();
	         while(U > cumden[i]) i++;
	         rndm[j]=i;
      	 }
      }
}

/*
 Simulates from the multiplicative binomial distribution conditional on the value being less or equal to ub, and greater than or equal to lb. 


Possibly no longer needed - delete if so

*/
void rMultbinomC_lb_ub(int *n,
                       int *size_R, 
                       double *p_R, 
                       double *psi_R, 
                       int *lb, 
                       int *ub, 
                       int *rndm, 
                       int *seed, 
                       double *den){

      double cumden[*size_R+1];
      int i;  
      int j=0;
      double U;
      int seeds = seed[0];
      
      
      sdprand(seeds);  //seed the random number generator

	
  	  dMultBinomC_internal(*size_R, *p_R, *psi_R, den=den);
	    cumden[0] =den[0];
	
  	  for(i=1; i<= (*size_R); i++){
  	    cumden[i] = cumden[i-1]+den[i];
    	}

      for(j=0; j<*n;j++){
  	    R_CheckUserInterrupt();
    	  rndm[j] = -1;
	      while(rndm[j]<*lb || rndm[j]>*ub){
	        i=0;
	        U = dprand();
	        while(U > cumden[i]) i++;
	        rndm[j]=i;
        }
      }
}



