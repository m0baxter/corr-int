
#include <cmath>
#include <memory>
#include "newtoncotes.hpp"

/* 
Global Variables:
*/

double roots[4];   //roots used as abscissas in Gauss method.
double weights[4]; //weights used in the Gauss method.

void gauss_setup() {
   /*Calculates and stores the abscissas/weights for the 4th order Gauss integration method.*/

   //Calculate and store the roots:
   roots[0] =  sqrt( 525.0 - 70.0 * sqrt( 30.0 ) ) /35.0;
   roots[1] = -sqrt( 525.0 - 70.0 * sqrt( 30.0 ) ) /35.0;
   roots[2] =  sqrt( 525.0 + 70.0 * sqrt( 30.0 ) ) /35.0;
   roots[3] = -sqrt( 525.0 + 70.0 * sqrt( 30.0 ) ) /35.0;

   //calculate and store the weights:
   weights[0] = ( 18.0 + sqrt( 30.0 ) ) / 36.0;
   weights[1] = ( 18.0 + sqrt( 30.0 ) ) / 36.0;
   weights[2] = ( 18.0 - sqrt( 30.0 ) ) / 36.0;
   weights[3] = ( 18.0 - sqrt( 30.0 ) ) / 36.0;
}


double lagrangebasis( const double *x, const int k, const double p, const int n) {
   /*Return the value at p of the Lagrange basis element k on the points x of degree n.*/

   double result = 1.0;
   
   for(int i = 0; i < n; ++i) {
      if (i != k){
         result *= (p - x[i])/(x[k] - x[i]);
      }
   }
   
   return result;
}


double lagrangeinterp( const double *x, const double *y, const double p, const int n) {
   /*Returns the value of the lagrange interpolation polynomial of degree n on the data [x,y] evaluated at p.*/

   double result = 0.0;

   for (int i = 0; i < n; ++i) {
      result += y[i]*lagrangebasis(x, i, p, n);
   }

   return result;
}


std::unique_ptr<double[]> slice( const double *arr, const int first, const int last ) {
   /*Returns an array that consists of the elements arr[first] through arr[last-1]
      (just like python can do except you need to catch the pointer and delete when finished).*/

   std::unique_ptr<double[]> result( new double[last - first] );

   for (int i = first; i < last; ++i) {
      result[i - first] = arr[i]; 
   }

   return result;
}


double transform( const double a, const double b, const double x) {
   /*Transforms the interval [a,b] -> [-1,1] and returns the new value in [-1,1] at x in [a,b]. */

   return (b - a)*x/2.0 + (a + b)/2.0;
}


double gauss_interp(const double *x, const double *y, const int n) {
   /*Calculates the integral of the Lagrange interpolation polynomial on the data (x,y)
      using the 4th order gauss method. n is the number of elements in each array.*/

   double total = 0.0;

   //Prefactor of the transformed integral:
   double c = (x[n-1] - x[0])/2.0;

   for (int i = 0; i < 4; ++i) {
      total += lagrangeinterp(x, y, transform( x[0], x[n-1], roots[i]), n) * weights[i];
   }

   return total*c;
}


double newtoncotes( const double *x, const double *y, const int n) {
   /*Calculates the 7th order Newton Cotes approximation to the function (x,y).
      The lengths of the input arrays must be equal and the same as n.*/

   //Will hold the result:
   double total = 0.0;

   //the index of the last multiple of 7:
   int last = 7*((n - 1)/7);

   for (int i = 0; i < n - 7; i+=7) {

      std::unique_ptr<double[]> slx = slice(x, i, i + 7 + 1);
      std::unique_ptr<double[]> sly = slice(y, i, i + 7 + 1);

      total += gauss_interp(slx.get(), sly.get(), 8);
   }

   //add the tail:
   std::unique_ptr<double[]> slx = slice(x,last,n);
   std::unique_ptr<double[]> sly = slice(y,last,n);

   total += gauss_interp( slx.get(), sly.get(), n-last);

   return total;
}

