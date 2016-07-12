/* Copyright (C) 1992  N.M. Maclaren
   Copyright (C) 1992  The University of Cambridge

    This software may be reproduced and used freely, provided that all
    users of it agree that the copyright holders are not liable for any
    damage or injury caused by use of this software and that this condition
    is passed onto all subsequent recipients of the software, whether
    modified or not.
*/


/*Modified by Rich Wilkinson*/


#include <stdlib.h>
#include <math.h>


#define xmod 1000009711.0
#define ymod 33554432.0




#if __STDC__
extern void sdprand (int seed);
extern double dprand (void);
#else
extern void sdprand ();
extern double dprand ();
#endif


/* The following is static only because some compilers do not optimise constant
double expressions. */

#ifdef __STDC__
const static double
#else
static double
#endif
    offset = 1.0/ymod, xmod2 = 2.0*xmod, xmod4 = 4.0*xmod, tiny = 1.0e-17;


#ifdef __STDC__
 double Exponential(double m);
 int Poisson(double m);
  int Bernoulli(double p);
 int Binomial(int n, double p);
 int rMultbinom_internal(int size, double p, double psi);
 #endif
