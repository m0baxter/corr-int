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
#include <stdlib.h>
#include "newtoncotes.hpp"
#include "AngularMomentum.hpp"


#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

class WaveFunction {

   public:
      WaveFunction( char* , int );
      ~WaveFunction();
      float get_impact( int );
      double indi_electron( int );
      double correlationintegral( int );
      double correlationintegral_wb( int );

   private:
      //Variables:
      static const int MAXn = 5;
      static const int LEN = 1601;
      static const int LIMIT = 1601; //1001;
      double r[LEN]; //Lattice points for the radial R_nl functions.
      double Rad_T[5][4][LEN];  //The R_nl function values on target.
      double Rad_P[5][4][LEN];  //The R_nl function values on projectile.
      double Rad_gs_T[5][4][LEN];  //The ground state R_nl function values on target.
      double Rad_gs_P[5][4][LEN];  //The ground state R_nl function values on projectile.
      float b[30];  //holds the values of the impact parameter.
      std::complex<double> amps_T[30][5][4][7];  //Holds the amplitude values for the 30 impact parameters on the target.
      std::complex<double> amps_P[30][5][4][7];  //Holds the amplitude values for the 30 impact parameters.

      double int_table_T[5][4][5][4][5][4][5][4]; //Holds the integral table.
      double int_table_P[5][4][5][4][5][4][5][4]; //Holds the integral table.
      double int_table_wb_T[5][4][5][4][5][4][5][4]; //Holds the integral table for the Wilken and Bauer version.
      double int_table_wb_P[5][4][5][4][5][4][5][4]; //Holds the integral table for the Wilken and Bauer version.

      int N_terms_T;  //The number of terms in the wave function on target.
      int N_terms_P;  //The number of terms in the wave function on projectile.
      int (*terms_T)[2][3];  //Holds the quantum numbers of the terms of the wave function.
      int (*terms_P)[2][3];  //Holds the quantum numbers of the terms of the wave function.      
      double* factors_T;  //Holds the prefactors for each term of the wave function.
      double* factors_P;  //Holds the prefactors for each term of the wave function.      

      //Functions:
      void readlattice();
      void readradial_T();
      void readradial_gs_T();
      void readradial_P();
      void readradial_gs_P();
      void readamplitudes_T( int );
      void readamplitudes_P( int );
      void readinput_T( char* );
      void readinput_P( char* );      
      int r_index( double );
      double R_T( int , int , double );
      double R_T( int , int , double );
      double R_gs_P( int , int , double );
      double R_gs_P( int , int , double );
      double Hep1( double );
      //Modify these later:
      double* integrand(int , int , int , int , int , int , int , int );
      void generate_integral_table();
      double* integrand_wb( double , int , int , int , int , int , int , int , int );
      void generate_integral_table_wb( double );

};

#endif
