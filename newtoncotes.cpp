
#include "newtoncotes.hpp"

/* 
Global Variables:
*/

double roots[4];   //roots used as abscissas in Gauss method.
double weights[4]; //weights used in the Gauss method.

void gauss_setup() {
	/*Calculates and stores the abscissas/weights for the 4th order Gauss integration method.*/

	//Calculate and store the roots:
	roots[0] =   0.33998104358485626480266575910324468720057586977092; //sqrt( 525.0 - 70.0 * sqrt( 30.0 ) ) /35.0;
	roots[1] =  -0.33998104358485626480266575910324468720057586977092; //- sqrt( 525.0 - 70.0 * sqrt( 30.0 ) ) /35.0;
	roots[2] =   0.86113631159405257522394648889280950509572537962973; //sqrt( 525.0 + 70.0 * sqrt( 30.0 ) ) /35.0;
	roots[3] = - 0.86113631159405257522394648889280950509572537962973; //- sqrt( 525.0 + 70.0 * sqrt( 30.0 ) ) /35.0;

	//calculate and store the weights:
	weights[0] = ( 18.0 + sqrt( 30.0 ) ) / 36.0;
	weights[1] = ( 18.0 + sqrt( 30.0 ) ) / 36.0;
	weights[2] = ( 18.0 - sqrt( 30.0 ) ) / 36.0;
	weights[3] = ( 18.0 - sqrt( 30.0 ) ) / 36.0;
}


double lagrangebasis( double x[], int k, double p , int n) {
	/*Return the value at p of the Lagrange basis element k on the points x of degree n.*/

	double result = 1.0;
	
	for(int i = 0; i < n; ++i) {
		if (i != k){
			result *= (p - x[i])/(x[k] - x[i]);
		}
	}
	
	return result;
}


double lagrangeinterp( double x[],  double y[], double p , int n) {
	/*Returns the value of the lagrange interpolation polynomial of degree n on the data [x,y] evaluated at p.*/

	double result = 0.0;


	for (int i = 0; i < n; ++i) {
		result += y[i]*lagrangebasis(x, i, p, n);
	}

	return result;
}


double* slice( double arr[], int first, int last ) {
	/*Returns an array that consists of the elements arr[first] through arr[last-1]
      (just like python can do except you need to catch the pointer and delete when finished).*/

	double* result = new double[last - first];

	for (int i = first; i < last; ++i) {
		result[i - first] = arr[i]; 
	}

	return result;
}


double transform( double a, double b, double x) {
	/*Transforms the interval [a,b] -> [-1,1] and returns the new value in [-1,1] at x in [a,b]. */

	return (b - a)*x/2.0 + (a + b)/2.0;
}


double gauss_interp( double x[], double y[], int n) {
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


double newtoncotes( double x[], double y[], int n) {
	/*Calculates the 7th order Newton Cotes approximation to the function (x,y).
		The lengths of the input arrays must be equal and the same as n.*/

	//Will hold the result:
	double total = 0.0;

	//the index of the last multiple of 7:
	int last = 7*((n - 1)/7);

	for (int i = 0; i < n - 7; i+=7) {

		double* slx = slice(x, i, i + 7 + 1);
		double* sly = slice(y, i, i + 7 + 1);

		total += gauss_interp(slx, sly, 8);

		delete [] slx;
		delete [] sly;
	}

	//add the tail:
	double* slx = slice(x,last,n);
	double* sly = slice(y,last,n);

	total += gauss_interp( slx, sly, n-last);

	delete [] slx;
	delete [] sly;

	return total;
}

