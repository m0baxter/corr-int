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
		double Rad[5][4][LEN];  //The R_nl function values.
		double Rad_gs[5][4][LEN];  //The ground state R_nl function values.
		float b[30];  //holds the values of the impact parameter.
		std::complex<double> amps[30][5][4][7];  //Holds the amplitude values for the 30 impact parameters.

		double int_table[5][4][5][4][5][4][5][4]; //Holds the integral table.
		double int_table_wb[5][4][5][4][5][4][5][4]; //Holds the integral table for the Wilken and Bauer version.

		int N_terms;  //The number of terms in the wave function.
		int (*terms)[2][3];  //Holds the quantum numbers of the terms of the wave function.
		double* factors;  //Holds the prefactors for each term of the wave function.

		//Functions:
		void readlattice();
		void readradial();
		void readradial_gs();
		void readamplitudes( int );
		void readinput( char* );
		int r_index( double );
		double R( int , int , double );
		double R_gs( int , int , double );
		double Hep1( double );
		double* integrand(int , int , int , int , int , int , int , int );
		void generate_integral_table();
		double* integrand_wb( double , int , int , int , int , int , int , int , int );
		void generate_integral_table_wb( double );

};

#endif
