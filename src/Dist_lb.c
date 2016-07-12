
/* C Code */

#include "factorial.h"
#include "random.h"

// checked seems to work fine
int *rpois_lb(int *n, double *lambda,  int *lb, int *rndm, int *seed){
  int smp;
  int i;
  int seeds = seed[0];
   
  sdprand(seeds);  //seed the random number generator
  
  for(i=1; i<=*n; i++){
    while((smp=Poisson(*lambda))<*lb){R_CheckUserInterrupt();};
    rndm[i-1]=smp;//smp;
   // printf("%d\n", rndm[i-1]);
  }
}



int *rbinom_lb(int *n, int *size, double *p,  int *lb, int *rndm, int *seed){
  int smp;
  int i;
  int seeds=seed[0];
  
  sdprand(seeds);
  
  for(i=1; i<=*n; i++){
    while((smp=Binomial(*size, *p))<*lb){R_CheckUserInterrupt();};
    rndm[i-1]=smp;
   // printf("%d\n", rndm[i-1]);
  }
}



int *rbinom_lb_ub(int *n, int *size, double *p,  int *lb, int *ub, int *rndm, int *seed){
  int smp;
  int i;
  int seeds=seed[0];
  
  sdprand(seeds);
  
  for(i=1; i<=*n; i++){
    while((smp=Binomial(*size, *p))<*lb || smp>*ub){R_CheckUserInterrupt();};
    rndm[i-1]=smp;
   // printf("%d\n", rndm[i-1]);
  }
}


/*
Returns n random variables (N,M) with
N>=size, M>=m, and M<=N-size+m 
where N ~ Po(lambda)
      M ~ Bin(size, p)
*/
int *rNM_given_nm_binomial(int *n, 
                           double *lambda,  
                           int *n_data, 
                           double *p, 
                           int *m, 
                           int *N, 
                           int *M, 
                           int *seed){
  
  int smp;
  int i;
  int seeds=seed[0];
  
  sdprand(seeds);
 
  for( i=0; i<*n; i++){
    N[i]=Poisson(*lambda);
    M[i]=Binomial(N[i], *p);
    while( (N[i] < *n_data || M[i] < *m ) || (M[i] > N[i] - *n_data + *m) ){
      R_CheckUserInterrupt();
      N[i]=Poisson(*lambda);
      M[i]=Binomial(N[i], *p);
    };
  }
}


/******************************************************************************


 Multiplicative binomial functions for when 



*******************************************************************************/


/*
Returns n random variables (N,M) with
N>=size, M>=m, and M<=N-size+m 
where N ~ Po(lambda)
      M ~ MultiBin(size, p)
*/
int *rNM_given_nm_multbinom(int *n, 
                           double *lambda,  
                           int *n_data, 
                           double *p,
                           double *psi,
                           int *m,
                           int *N, 
                           int *M,
                           int *seed){

      int i;	
      int j=0;
      int seeds = seed[0];
      sdprand(seeds);  //seed the random number generator
 

      for( i=0; i<*n; i++){
        N[i]=Poisson(*lambda);
        M[i]=rMultbinom_internal(N[i], *p, *psi);
        while( (N[i] < *n_data || M[i] < *m ) || (M[i] > N[i] - *n_data + *m) ){
          R_CheckUserInterrupt();
          N[i]=Poisson(*lambda);
          M[i]=rMultbinom_internal(N[i], *p, *psi);
        };
     }  	
}




/******************************************************************************


 Double binomial functions for when 



*******************************************************************************/


/*
Returns n random variables (N,M) with
N>=size, M>=m, and M<=N-size+m 
where N ~ Po(lambda)
      M ~ MultiBin(size, p)
*/
int *rNM_given_nm_doublebinom(int *n, 
                           double *lambda,  
                           int *n_data, 
                           double *p,
                           double *psi,
                           int *m,
                           int *N, 
                           int *M,
                           int *seed){

      int i;  
      int j=0;
      int seeds = seed[0];
      sdprand(seeds);  //seed the random number generator
 

      for( i=0; i<*n; i++){
        N[i]=Poisson(*lambda);
        M[i]=rDoublebinom_internal(N[i], *p, *psi);
        while( (N[i] < *n_data || M[i] < *m ) || (M[i] > N[i] - *n_data + *m) ){
          R_CheckUserInterrupt();
          N[i]=Poisson(*lambda);
          M[i]=rDoublebinom_internal(N[i], *p, *psi);
        };
     }  	
}




