
#include "factorial.h"


void dDoublebinomC(int *size_R, double *p_R, double *psi_R, double *logden, double *den){
  // returns entire density

  int i;
	double p = *p_R;
	double psi = *psi_R;
	int size=*size_R;
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
	log_norm_const=log(norm_const);
	
	for(i=0; i<=size; i++){
	      den[i] = unnorm_den[i] * norm_const;
	      logden[i] = unnorm_logden[i] + log_norm_const;
	}
}
