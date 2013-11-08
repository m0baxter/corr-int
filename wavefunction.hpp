//silverman
// LIMIT:     1601
//heopmxdyn
// LIMIT:     1534

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <complex>
#include <limits>
#include <stdlib.h>
#include "newtoncotes.hpp"
#include "AngularMomentum.hpp"


#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

class WaveFunction {

   public:
      WaveFunction( char* , int, wb );
      ~WaveFunction();
      float get_impact( int );
      double indi_electron( const char, int );
      double correlationintegral( const char, int );
      double correlationintegral_wb( const char, int );

   private:
      //Variables:
      static const int MAXn = 5;
      static const int LEN = 1601;
      static const int LIMIT = 1601; //1001;
      double r[LEN]; //Lattice points for the radial R_nl functions.
      double Rad[2][2][5][4][LEN];  //The R_nl function values on target.
      float b[30];  //holds the values of the impact parameter.
      std::complex<double> amps[2][30][5][4][7];  //Holds the amplitude values for the 30 impact parameters on the target.

      double int_table[2][5][4][5][4][5][4][5][4]; //Holds the integral table.
      double int_table_wb[2][5][4][5][4][5][4][5][4]; //Holds the integral table for the Wilken and Bauer version.

      int N_terms_T;  //The number of terms in the wave function on target.
      int (*terms_T)[2][3];  //Holds the quantum numbers of the terms of the wave function.
      double* factors_T;  //Holds the prefactors for each term of the wave function.
      int N_terms_P;  //The number of terms in the wave function on projectile.
      int (*terms_P)[2][3];  //Holds the quantum numbers of the terms of the wave function.
      double* factors_P;  //Holds the prefactors for each term of the wave function.      

      //Functions:
      int TP_toint( const char );
      int DG_toint( const char );
      void readlattice();
      void readradial( const char, const char );
      void readamplitudes( const char, int );
      void readinput( const char, char* );
      int r_index( double );
      double R( const char, const char, int , int , double );
      double Hep1( double );
      std::complex<double> a( const char, int , int, int, int );
      double table( const char, int, int, int, int, int, int, int, int );
      
            //Here
      
      //Modify these later:
      double* integrand( const char, int , int , int , int , int , int , int , int );
      void generate_integral_table( const, char );
      double* integrand_wb( const char, double , int , int , int , int , int , int , int , int );
      void generate_integral_table_wb( const char, double );
      
      // add TP integral functions

};

#endif
