//silverman
// LIMIT:     1601
//heopmxdyn
// LIMIT:     1534

#include <complex>
#include <memory>


#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

class WaveFunction {

   public:
      //Functions:
      WaveFunction( const int, const bool );
      WaveFunction()                  = delete;  //no default
      WaveFunction( WaveFunction & )  = delete;  //no copy
      WaveFunction( WaveFunction && ) = delete;  //no move
      ~WaveFunction();
      float get_impact( const int );
      double indi_electron( const char, const int );
      double correlationintegral( const char, const int );
      double correlationintegral_wb( const char, const int );

   private:
      //Variables:
      static const int MAXn = 5;
      static const int LEN = 1601;
      static const int LIMIT = 1543; //1001;
      double r[LEN];                                  //Lattice points for the radial R_nl functions.
      double Rad[2][2][5][4][LEN];                    //The R_nl functions.
      float b[30];                                    //holds the values of the impact parameter.
      std::complex<double> amps[2][30][5][4][7];      //Holds the amplitude values for the 30 impact parameters.
      double int_table[2][5][4][5][4][5][4][5][4];    //Holds the integral table.
      double int_table_wb[2][5][4][5][4][5][4][5][4]; //Holds the integral table for the Wilken and Bauer version.
      
      //Storage of ground state information:
      int N_terms_T;        //number of terms on target.
      int (*terms_T)[2][3]; //quantum numbers of terms on target.
      double *factors_T;    //prefactors for terms on target.
      int N_terms_P;        //number of terms on projectile.
      int (*terms_P)[2][3]; //quantum numbers of terms on projectile.
      double *factors_P;    //prefactors for terms on projectile.

      //Functions:
      int TP_toint( const char );
      int DG_toint( const char );
      int charge( const char );
      void readlattice();
      void readradial( const char, const char );
      void readamplitudes( const char, const int );
      void readinput( const char );
      double R( const char, const char, const int, const int, const int );
      double Hlike( int, int );
      std::complex<double> a( const char, const int, const int, const int, const int );
      std::unique_ptr<double[]> integrand( const char, const int, const int, const int, const int, const int, const int, const int, const int );
      void generate_integral_table( const char );
      std::unique_ptr<double[]> integrand_wb( const char, const double, const int, const int, const int, const int, const int, const int, const int, const int );
      void generate_integral_table_wb( const char, const double );
      double table( const char, const int, const int, const int, const int, const int, const int, const int, const int );
      double table_wb( const char, const int, const int, const int, const int, const int, const int, const int, const int );
      
};

#endif
