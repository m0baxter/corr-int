
#include <memory>

#ifndef NEWTONCOTES_HPP
#define NEWTONCOTES_HPP

   double lagrangebasis( const double *x, const int, const double, const int );
   double lagrangeinterp( const double *x, const double *y, const double, const int );
   std::unique_ptr<double[]> slice( const double *arr, const int, const int );
   double transform( const double, const double, const double );
   double gauss_interp( const double *x, const double *y, const int );
   double newtoncotes( const double *x, const double *y, const int );

#endif
