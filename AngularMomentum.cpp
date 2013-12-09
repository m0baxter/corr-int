
#include <algorithm>
#include <math.h>
#include "AngularMomentum.hpp"

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982;

int factorial(int n) {
	/*returns n!*/

	if ( n <= 1 ) {

		return 1;
	}
	else {

		int fact = 1;

		for (int i = 1; i < n + 1; ++i) {
			fact = fact*i;
		}
		return fact;
	}
}

bool TriangleBroken ( int l1, int l2, int l3 ) {
	/*Tests whether the triangle condition |l1 - l2| <= l3 <= l1 - l2 is broken. If it does not it returns True, else False*/

   return ( abs( l1 - l2 ) > l3 ) or ( l1 + l2 < l3 );
}


bool IsEven( int i ) {
	/*Checks whether the integer i is even. Does not check all the conditions, most will be met from the contex.*/

	return (i % 2) == 0;
}


double _3j( int l1, int m1, int l2, int m2, int l3, int m3 ) {
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


double gaunt( int l1, int m1, int l2, int m2, int l3, int m3 ) {
	/*Calculates the Gaunt integral \int d\Omega Y^m1_l1 Y^m2_l2 Y^m3_l3. */

	double norm = sqrt( (2.0*l1 + 1.0)*(2.0*l2 + 1.0)*(2.0*l3 + 1.0)/( 4.0 * PI ) ); //4.0*atan(1.0)

	return norm * _3j(l1,m1,l2,m2,l3,m3) * _3j(l1,0,l2,0,l3,0);

}
