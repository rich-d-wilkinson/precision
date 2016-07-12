/* Copyright (C) 1992  N.M. Maclaren
   Copyright (C) 1992  The University of Cambridge

    This software may be reproduced and used freely, provided that all
    users of it agree that the copyright holders are not liable for any
    damage or injury caused by use of this software and that this condition
    is passed onto all subsequent recipients of the software, whether
    modified or not.
*/


/*  Modified by Rich Wilkinson

     I'm unsure of the reliability of this random number generator. I was given it by a professor in Cambridge who had          written it, but I've never tested it and don't know what type it is etc. No guarantees of reliability.





*/

#include "random.h"


/* The following is the workspace for the random number generator. */

static int initial = 1, index;
static double poly[101], other;

extern void sdprand (seed)
int seed;
{
    int ix, iy, iz, i;
    double x=0;

/* seed should be set to an integer between 0 and 9999 inclusive; a value of 0
will initialise the generator only if it has not already been done. */

    if (initial || seed != 0)
        initial = 0;
    else
        return;

/* index must be initialised to an integer between 1 and 101 inclusive,
poly[0...100] to integers between 0 and 1000009710 inclusive (not all 0), and
other to a non-negative proper fraction with denominator 33554432.  It uses the
Wichmann-Hill generator to do this. */

    ix = (seed >= 0 ? seed : -seed)%10000+1;
    iy = 2*ix+1;
    iz = 3*ix+1;
    for (i = -11; i < 101; ++i) {
        if (i >= 0) poly[i] = floor(xmod*x);
        ix = (171*ix)%30269;
        iy = (172*iy)%30307;
        iz = (170*iz)%30323;
        x = ((double)ix)/30269.0+((double)iy)/30307.0+((double)iz)/30323.0;
        x = x-floor(x);
    }
    other = floor(ymod*x)/ymod;
    index = 0;
}



extern double dprand ()
{
    int n;
    double x=0, y;

/* This returns a uniform (0,1) random number, with extremely good uniformity
properties.  It assumes that double precision provides at least 33 bits of
accuracy, and uses a power of two base. */

    if (initial) sdprand(0);

/* See [Knuth] for why this implements the algorithm described in the paper.
Note that this code is tuned for machines with fast double precision, but slow
multiply and divide; many, many other options are possible. */

    if ((n = index-64) < 0) n += 101;
    x = poly[index]+poly[index];
    x = xmod4-poly[n]-poly[n]-x-x-poly[index];
    if (x <= 0.0) {
        if (x < -xmod) x += xmod2;
        if (x < 0.0) x += xmod;
    } else {
        if (x >= xmod2) {
            x = x-xmod2;
            if (x >= xmod) x -= xmod;
        }
        if (x >= xmod) x -= xmod;
    }
    poly[index] = x;
    if (++index >= 101) index = 0;

/* Add in the second generator modulo 1, and force to be non-zero.  The
restricted ranges largely cancel themselves out. */

    do {
        y = 37.0*other+offset;
        other = y-floor(y);
    } while (other == 0.0);
    if ((x = x/xmod+other) >= 1.0) x -= 1.0;


    return x+tiny;
}


     double Exponential(double m)

{
  double aa;
  aa= dprand();
  if (aa==1.0) {return 10000000;}
  return (-m*log(1.0-aa));
}



   int Poisson(double m)
/* ================================================== 
 * Returns a Poisson distributed non-negative integer. 
 * NOTE: use m > 0
 * ==================================================
 */
{ 
  double t = 0.0;
  long   x = 0;

  while (t < m) {
    t += Exponential(1.0);
    x++;
  }
  return (x - 1);
}



/*************************Bernoulli(double p) ***********************/
   int Bernoulli(double p)

{
  double hh;
  int gg;
 hh=dprand();
gg= ((hh < (1.0 - p)) ? 0 : 1);
return gg;

}


/*************************Binomial(int n, double p) ***********************/
  // by RDW not original author 
  int Binomial(int n, double p)

{
  int ii;
  int gg=0;
  for(ii=1; ii<=n; ii++){
    gg = gg + Bernoulli(p);
  //  printf("%d\t", gg);
  }
return gg;
}







