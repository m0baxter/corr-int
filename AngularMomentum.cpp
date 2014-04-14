
#include <algorithm>
#include <cmath>
#include "AngularMomentum.hpp"

const double PI = acos(-1.0);

int factorial( const int n ) {
   /*returns n!*/
   
   int fact = 1;
   
   for(int i = 2; i <= n; ++i) {
      fact *= i;
   }
   
   return fact;
}

bool TriangleBroken ( const int l1, const int l2, const int l3 ) {
   /*Tests whether the triangle condition |l1 - l2| <= l3 <= l1 - l2 is broken.
   If it does not it returns True, else False*/

   return ( abs( l1 - l2 ) > l3 ) or ( l1 + l2 < l3 );
}


bool IsEven( int i ) {
   /*Checks whether the integer i is even. Does not check all the conditions,
	most will be met from the contex.*/

   return (i % 2) == 0;
}


double _3j( const int l1, const int m1, const int l2, const int m2, const int l3, const int m3 ) {
/*Calculates the Wigner 3-j coefficient .*/
	

   double n1, n2, n3, d1, d2, d3;
   double norm, sum, phase;

   //calculate the normalization factor:
   n1 = factorial(l3 + l1 - l2) * factorial(l3 - l1 + l2);
   n2 = factorial(l1 + l2 - l3);
   n3 = factorial(l3 - m3) * factorial(l3 + m3);
   d1 = factorial(l1 + l2 + l3 + 1);
   d2 = factorial(l1 - m1) * factorial(l1 + m1);
   d3 = factorial(l2 - m2) * factorial(l2 + m2);
   norm = pow(-1.0, l1 - l2 - m3)  * sqrt( ( n1*n2*n3 ) / ( d1*d2*d3 ) );

   //determine the limits of the loop:
   int kmin = std::max(0, l2 - l1 - m3);
   int kmax = std::min(l3 - l1 + l2, l3 - m3);

   //The phase of the first term of the sum:
   phase = pow( -1.0 , kmin + l2 + m2 );

   sum = 0.0;

   for (int k = kmin; k <= kmax; ++k) {

      n1 = factorial(l2 + l3 + m1 - k);
      n2 = factorial(l1 - m1 + k);
      d1 = factorial(k) * factorial(l3 - l1 + l2 - k);
      d2 = factorial(l3 - m3 - k);
      d3 = factorial(k + l1 - l2 + m3);

      double temp = phase * ( n1*n2 ) / (  d1*d2*d3 ) ;

      phase = -phase;
      sum += norm*temp;
   }

   return sum;
}


double gaunt( const int l1, const int m1, const int l2, const int m2, const int l3, const int m3 ) {
   /*Calculates the Gaunt integral \int d\Omega Y^m1_l1 Y^m2_l2 Y^m3_l3. */

   double norm = sqrt( (2.0*l1 + 1.0)*(2.0*l2 + 1.0)*(2.0*l3 + 1.0)/( 4.0 * PI ) );

   return norm * _3j(l1,m1,l2,m2,l3,m3) * _3j(l1,0,l2,0,l3,0);

}
