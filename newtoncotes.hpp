
#include <math.h>

#ifndef NEWTONCOTES_HPP
#define NEWTONCOTES_HPP

	void gauss_setup();
	double lagrangebasis( double x[], int , double , int );
	double lagrangeinterp( double x[],  double y[], double , int );
	double* slice( double arr[], int , int );
	double transform( double , double , double );
	double gauss_interp( double x[], double y[], int );
	double newtoncotes( double x[], double y[], int );

#endif
