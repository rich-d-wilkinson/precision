
#include "factorial.h"

////////////////////////////////////////////////////////////////////////////////
// double Factorial( int n )                                                  //
//                                                                            //
//  Description:                                                              //
//     This function computes n! for 0 <= n <= Factorial_Max_Arg().  If       //
//     n > Factorial_Max_Arg(), then DBL_MAX is returned and if n < 0, then   //
//     0 is returned.                                                         //
//                                                                            //
//  Arguments:                                                                //
//     int    n   Argument of the Factorial function.                         //
//                                                                            //
//  Return Values:                                                            //
//     If n is negative, then 0 is returned.  If n > Factorial_Max_Arg(),     //
//     then DBL_MAX is returned.  If 0 <= n <= Factorial_Max_Arg(), then      //
//     n! is returned.                                                        //
//                                                                            //
//  Example:                                                                  //
//     int n;                                                                 //
//     double x;                                                              //
//                                                                            //
//     x = Factorial( n );                                                    //
////////////////////////////////////////////////////////////////////////////////
double Factorial(int n) {

               // For a negative argument (n < 0) return 0.0 //

   if ( n < 0 ) return 0.0;


           // For a large postive argument (n >= N) return DBL_MAX //

   if ( n >= N ) return DBL_MAX;

                          // Otherwise return n!. //

   return (double) factorials[n];
}


////////////////////////////////////////////////////////////////////////////////
// long double xFactorial( int n )                                            //
//                                                                            //
//  Description:                                                              //
//     This function computes n! for 0 <= n <= Factorial_Max_Arg().  If       //
//     n > Factorial_Max_Arg(), then DBL_MAX is returned and if n < 0, then   //
//     0 is returned.                                                         //
//                                                                            //
//  Arguments:                                                                //
//     int    n   Argument of the Factorial function.                         //
//                                                                            //
//  Return Values:                                                            //
//     If n is negative, then 0 is returned.  If n > Factorial_Max_Arg(),     //
//     then DBL_MAX is returned.  If 0 <= n <= Factorial_Max_Arg(), then      //
//     n! is returned.                                                        //
//                                                                            //
//  Example:                                                                  //
//     int n;                                                                 //
//     long double x;                                                         //
//                                                                            //
//     x = xFactorial( n );                                                   //
////////////////////////////////////////////////////////////////////////////////
long double xFactorial(int n) {

               // For a negative argument (n < 0) return 0.0 //

   if ( n < 0 ) return 0.0L;


           // For a large postive argument (n >= N) return DBL_MAX //

   if ( n >= N ) return  (long double) DBL_MAX;

                          // Otherwise return n!. //

   return factorials[n];
}


////////////////////////////////////////////////////////////////////////////////
// int Factorial_Max_Arg( void )                                              //
//                                                                            //
//  Description:                                                              //
//     This function returns the maximum argument of the Factorial function   //
//     for which a number < DBL_MAX is returned.                              //
//                                                                            //
//  Arguments:                                                                //
//     none                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     N-1                                                                    //
//                                                                            //
//  Example:                                                                  //
//     int x;                                                                 //
//                                                                            //
//     x = Factorial_Max_Arg();                                               //
////////////////////////////////////////////////////////////////////////////////

int Factorial_Max_Arg( void ) { return N - 1; }



////////////////////////////////////////////////////////////////////////////////
// double Binomial_Coefficient( int n, int m )                                //
//                                                                            //
//  Description:                                                              //
//     This function returns C(n,m) for n >= 0, 0 <= m <= n and               //
//     C(n,m) < DBL_MAX.                                                      //
//                                                                            //
//     If C(n,m) >= DBL_MAX, then DBL_MAX is returned.  If n < 0 or           //
//     m < 0 or m > n, then 0 is returned.                                    //
//                                                                            //
//  Arguments:                                                                //
//     int    n   The number of objects.                                      //
//     int    m   Number of terms in the product.  m must be nonnegative.     //
//                                                                            //
//  Return Values:                                                            //
//     If n or m or both are negative,  then 0 is returned.  If m > n, then   //
//     0 is returned.  C(n,m) >= DBL_MAX, then DBL_MAX is returned.           //
//     Otherwise C(n,m) is returned.                                          //
//                                                                            //
//  Example:                                                                  //
//     int m, n;                                                              //
//     double x;                                                              //
//                                                                            //
//     x = Binomial_Coefficient( n, m );                                      //
////////////////////////////////////////////////////////////////////////////////
double Binomial_Coefficient(int n, int m) {
   return (double) xBinomial_Coefficient(n,m);
}


////////////////////////////////////////////////////////////////////////////////
// long double xBinomial_Coefficient( int n, int m )                          //
//                                                                            //
//  Description:                                                              //
//     This function returns C(n,m) for n >= 0, 0 <= m <= n and               //
//     C(n,m) < DBL_MAX.                                                      //
//                                                                            //
//     If C(n,m) >= DBL_MAX, then DBL_MAX is returned.  If n < 0 or           //
//     m < 0 or m > n, then 0 is returned.                                    //
//                                                                            //
//  Arguments:                                                                //
//     int    n   The number of objects.                                      //
//     int    m   Number of terms in the product.  m must be nonnegative.     //
//                                                                            //
//  Return Values:                                                            //
//     If n or m or both are negative,  then 0 is returned.  If m > n, then   //
//     0 is returned.  C(n,m) >= DBL_MAX, then DBL_MAX is returned.           //
//     Otherwise C(n,m) is returned.                                          //
//                                                                            //
//  Example:                                                                  //
//     int m, n;                                                              //
//     long double x;                                                         //
//                                                                            //
//     x = xBinomial_Coefficient( n, m );                                     //
////////////////////////////////////////////////////////////////////////////////
long double xBinomial_Coefficient(int n, int m) {

   long double ln_combination;
   long double combination;
   unsigned long c;

                 // For n < 0 or m < 0 or m > n return 0.0. //

   if ( (n < 0) || (m < 0) || (m > n) ) return 0.0L;

           // If n <= Factorial_Max_Arg(), simply return //
           // the quotient n! / m! (n - m)!.             //

   if ( n <= Factorial_Max_Arg() )
      return xFactorial(n) / ( xFactorial(m) * xFactorial(n-m) );

        // If ln(n!) - ln(m!) - ln((n-m)!) < ln(DBL_MAX) then return //
        // exp(ln(n!) - ln(m!) - ln((n-m)!) otherwise return DBL_MAX //
        // Further since C(n,m) is an integer, if possible correct   //
        // roundoff errors.                                          //

   ln_combination = xLn_Factorial(n) - xLn_Factorial(m) - xLn_Factorial(n-m);
   if ( ln_combination < ln_DBL_MAX ) {
      combination = expl(ln_combination);
      if (combination < (double) ULONG_MAX) {
         c = (unsigned long)(combination + 0.5L);
         combination = (long double) c;
      }
      return combination;
   }
   return (long double) DBL_MAX;
}

////////////////////////////////////////////////////////////////////////////////
// double Ln_Factorial( int n )                                               //
//                                                                            //
//  Description:                                                              //
//     This function computes ln(n!).  If n < 0, then -DBL_MAX is returned.   //
//                                                                            //
//  Arguments:                                                                //
//     int    n   Argument of the Factorial function.                         //
//                                                                            //
//  Return Values:                                                            //
//     If n is negative, then -DBL_MAX is returned otherwise ln(n!) is        //
//     returned.                                                              //
//                                                                            //
//  Example:                                                                  //
//     int n;                                                                 //
//     double x;                                                              //
//                                                                            //
//     x = Ln_Factorial( n );                                                 //
////////////////////////////////////////////////////////////////////////////////
double Ln_Factorial(int n) {

             // For a negative argument (n < 0) return -DBL_MAX //

   if ( n < 0 ) return -DBL_MAX;

                       // If n < N then return ln(n!). //

   if ( n < N2 ) return (double) ln_factorials[n];

            // Otherwise return asymptotic expansion of ln(n!). //

   return (double) xLn_Factorial_Asymptotic_Expansion( n );
}


////////////////////////////////////////////////////////////////////////////////
// long double xLn_Factorial( int n )                                         //
//                                                                            //
//  Description:                                                              //
//     This function computes ln(n!).  If n < 0, then -LDBL_MAX is returned.  //
//                                                                            //
//  Arguments:                                                                //
//     int    n   Argument of the Factorial function.                         //
//                                                                            //
//  Return Values:                                                            //
//     If n is negative, then -LDBL_MAX is returned otherwise ln(n!) is       //
//     returned.                                                              //
//                                                                            //
//  Example:                                                                  //
//     int n;                                                                 //
//     long double x;                                                         //
//                                                                            //
//     x = xLn_Factorial( n );                                                //
////////////////////////////////////////////////////////////////////////////////
long double xLn_Factorial(int n) {

            // For a negative argument (n < 0) return LDBL_MAX //

   if ( n < 0 ) return -LDBL_MAX;

                       // If n < N then return ln(n!). //

   if ( n < N2 ) return ln_factorials[n];

            // Otherwise return asymptotic expansion of ln(n!). //

   return xLn_Factorial_Asymptotic_Expansion( n );
}


////////////////////////////////////////////////////////////////////////////////
// static long double xLn_Factorial_Asymptotic_Expansion( int n )             //
//                                                                            //
//  Description:                                                              //
//     This function estimates log(n!) by evaluating the asymptotic           //
//     expression:                                                            //
//         ln(n!) ~ ln(2sqrt(2pi)) - (n+1) + (n + 1/2) ln (n+1)               //
//                    Sum B[2j] / [ 2j * (2j-1) * (n+1)^(2j-1) ], summed over //
//     j from 1 to 6 and where B[2j] is the 2j-th Bernoulli number.           //
//                                                                            //
//  Arguments:                                                                //
//     int n  Argument of the ln n!. The argument n must be non-negative.     //
//                                                                            //
//  Return Values:                                                            //
//     ln(n!) where n > N = 171.                                              //
//                                                                            //
//  Example:                                                                  //
//     int n;                                                                 //
//     long double g;                                                         //
//                                                                            //
//     g = xFactorial_Asymptotic_Expansion( n );                              //
////////////////////////////////////////////////////////////////////////////////

static long double xLn_Factorial_Asymptotic_Expansion( int n ) {
   const int  m = 6;
   long double term[6];
   long double sum;
   long double xx = (long double) ((n + 1) * (n + 1));
   long double xj = (long double) (n + 1);
   long double lnfactorial = log_sqrt_2pi - xj + (xj - 0.5L) * logl(xj);
   int i;

   for (i = 0; i < m; i++) { term[i] = B[i] / xj; xj *= xx; }
   sum = term[m-1];
   for (i = m - 2; i >= 0; i--)
      if (fabsl(sum) > fabsl(term[i])) sum = term[i];
      else break;
   for (; i >= 0; i--) sum += term[i]; 
   return lnfactorial + sum;
}





long int factorial(int n)
{
  if (n == 0)
    return 1;
  else
    return(n * factorial(n-1));
}

long int combination(int n, int r)
{
  return(factorial(n)/(factorial(r)*factorial(n-r)));
}

ULONG binomial(ULONG n, ULONG k)   // doesn't work for n>20 - s
{
	ULONG r = 1, d = n - k;
 
	/* choose the smaller of k and n - k */
	if (d > k) { k = d; d = n - k; }
 
	while (n > k) {
		if (r >= UINT_MAX / n) return 0; /* overflown */
		r *= n--;
 
		/* divide (n - k)! as soon as we can to delay overflows */
		while (d > 1 && !(r % d)) r /= d--;
	}
	return r;
}
 
