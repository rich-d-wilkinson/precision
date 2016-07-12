
#include "factorial.h"


void dMultbinomC(int *size_R, double *p_R, double *psi_R, double *logden, double *den){
  // returns entire density

	int i;
	double p = *p_R;
	double psi = *psi_R;
	int size=*size_R;
	double unnorm_logden[size+1];
	double unnorm_den[size+1];
	double sum=0;
	double log_norm_const, norm_const;
	
	for(i=0; i<=size; i++){
	      unnorm_logden[i] = log(Binomial_Coefficient(size, i)) + i*log(p) + log(1-p) * (size-i) +  psi * i * (size-i);
	      unnorm_den[i] = exp(unnorm_logden[i]);
	      sum = sum + unnorm_den[i];	
	}
	//*psi_R= Binomial_Coefficient(305,15);//binomial(30,15);   //combination(size,1);
	
	norm_const = 1/sum;
	log_norm_const=log(norm_const);
	
	for(i=0; i<=size; i++){
	      den[i] = unnorm_den[i] * norm_const;
	      logden[i] = unnorm_logden[i] + log_norm_const;
	}
}

