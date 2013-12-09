
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <sstream>
#include <algorithm>
#include <memory>
#include "StringManipulators.hpp"
#include "newtoncotes.hpp"
#include "AngularMomentum.hpp"
#include "wavefunction.hpp"


const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982;

WaveFunction::WaveFunction( const char* inputpath_T, const char* inputpath_P, int energy, bool wb ) {
   /*Constructor for objects of class WaveFunction.*/

   //initialize tables:
   readlattice();
   readradial('T', 'D');
   readradial('P', 'D');
   readradial('T', 'G');
   readradial('P', 'G');
   readamplitudes( 'T', energy );
   readamplitudes( 'P', energy );
   readinput( 'T', inputpath_T );
   readinput( 'P', inputpath_P );

   //Calculate table of integrals
   if (!wb){
      generate_integral_table('T');
      generate_integral_table('P');
   }
}


WaveFunction::~WaveFunction() {
   /*Destroys an object of class WaveFunction.*/

   delete [] terms_T;
   delete [] terms_P;
   delete [] factors_T;
   delete [] factors_P;
}


int WaveFunction::TP_toint( const char c ) {

   return ( 84 - (c + 0) )/4;
}


int WaveFunction::DG_toint( const char c ) {

   return ( (c + 0) - 68 )/3;
}


int WaveFunction::charge( const char c ) {

   return ( (c + 0) - 76 )/4;
}


float WaveFunction::get_impact( int i) {

   return b[i];
}


void WaveFunction::readlattice() {
   /*Reads/stores the file containing the lattice points.*/

   std::string line;
   const char* path;

   //Open file to be read:
     path = "/home/baxter/Documents/Observable/Input/radial/heopmxdyn/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/opm/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/silverman_1s1s/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/silverman full/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-2p/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s(all)/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d(all)/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4f/lattice.txt";

   std::ifstream readfile(path);

   r[0] = 0.0;

   for (int i = 1; i < LEN; ++i) {

      getline(readfile,line);
      r[i] = str_to<double>(line);
   }
   readfile.close();
}


void WaveFunction::readradial( const char centre, const char type ) {
   /*Reads/stores the R_{nl} function data.*/

   std::string path;
   std::string line;
   
   int c = TP_toint(centre);
   int t = DG_toint(type);

   for (int n = 1; n < 5; ++n) {
      for (int l = 0; l < n; ++l) {
      
         switch (centre + type) {
            case 'T' + 'D': //target dynamic
               //path = "/home/baxter/Documents/Observable/Input/radial/heopmxdyn/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/opm/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/silverman_1s1s/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/silverman full/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-2p/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s(all)/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d(all)/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4f/R_";
                 path = "/home/baxter/Documents/Observable/Input/radial/opmgrid/R_";
               break;
               
            case 'P' + 'D': //projectile dynamic
                 path = "/home/baxter/Documents/Observable/Input/radial/heopmxdyn/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/opm/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/silverman_1s1s/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/silverman full/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-2p/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s(all)/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d(all)/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4f/R_";
               break;
               
            case 'T' + 'G': //target ground state
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/silverman_1s1sp/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/silverman_full/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-2p/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4s/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4s(all)/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4d/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4d(all)/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4f/R_";
                 path = "/home/baxter/Documents/Observable/Input/radial/Ground state/opmgrid/R_";
               break;
               
            case 'P' + 'G': //projectile ground state
                 path = "/home/baxter/Documents/Observable/Input/radial/heopmxdyn/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/silverman_1s1sp/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/silverman_full/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-2p/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4s/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4s(all)/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4d/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4d(all)/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4f/R_";
               //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/hneg/R_";
               break;
               
            default:
               std::cout << "readradial error: Improper centre/type signifier" << std::endl;
               return;
         }

         std::stringstream ss;
         ss << n << l << ".txt";
         path += ss.str();

         //Open file:
         std::ifstream readfile(path.c_str());

         //Give R_nl(0) some value (doesn't really matter what it is as long as its not zero):
         Rad[c][t][n][l][0] = 1.0;

         for (int i = 1; i < LEN; ++i) {

            getline(readfile,line);
            Rad[c][t][n][l][i] = str_to<double>(line);
         }
         readfile.close();
      }
   }
}


void WaveFunction::readamplitudes( const char centre, int energy ) {
   /*Reads in the impact parameter and amplitude data and stores it in b and amp. The energy parameter must be 100 or 2000.*/

   std::string path, line;
   std::stringstream ss;
   
   int c = TP_toint(centre);
   
   //setup the path string:
   switch (centre) {
      case 'T':
         //path = "/home/baxter/Documents/Observable/Input/amplitudes/old/E";
         //path = "/home/baxter/Documents/Observable/Input/amplitudes/OPM/E";
         //path = "/home/baxter/Documents/Observable/Input/amplitudes/MCHF_1s4d/E";
         //path = "/home/baxter/Documents/Observable/Input/amplitudes/MCHF_1s4d(all)/E";
           path = "/home/baxter/Documents/Observable/Input/amplitudes/phe/target/E";
         break;

      case 'P':
           path = "/home/baxter/Desktop/E";   
         //path = "/home/baxter/Documents/Observable/Input/amplitudes/old/E";
         //std::string path = "/home/baxter/Documents/Observable/Input/amplitudes/OPM/E";
         //std::string path = "/home/baxter/Documents/Observable/Input/amplitudes/MCHF_1s4d/E";
         //std::string path = "/home/baxter/Documents/Observable/Input/amplitudes/MCHF_1s4d(all)/E";
           path = "/home/baxter/Documents/Observable/Input/amplitudes/phe/projectile/E";   
         break;
         
      default:
         std::cout << "readamplitudes error: Improper centre label" << std::endl;
         return;
   }

   ss << energy << ".txt";
   path += ss.str();

   //open file.
   std::ifstream readfile(path.c_str());

   for (int i = 0; i < 30; ++i) {

      //read impact parameter line.
      getline(readfile,line);
      std::vector<std::string> tokens = split(line,"=");
      b[i] = str_to<float>(tokens[1]);

      //store value in b.
      for (int n = 1; n < 5; ++n) {
         for (int l = 0; l < n; ++l) {
            for (int m = 0; m < 2*l + 1; ++m) {

               getline(readfile,line);
               tokens = split(line,"\t");
               double realpart = str_to<double>(tokens[0]);
               double impart   = str_to<double>(tokens[1]);
               amps[c][i][n][l][m] = std::complex<double>( realpart, impart );

            }
         }
      }
   }
   readfile.close();
}


void WaveFunction::readinput( const char centre, const char* inputpath ) {
   /*Reads and stores the wave function quantum number data from the file located at inputpath.*/

   std::ifstream readfile ( inputpath );
   std::string line;
   int N;

   //Read the number of terms:
   getline(readfile,line);
   
   switch (centre) {
      case 'T':
         N_terms_T = atoi(line.c_str());
         N = N_terms_T;
         break;
         
      case 'P':
         N_terms_P = atoi(line.c_str());
         N = N_terms_P;
         break;
         
      default:
         std::cout << "readinput error: Improper centre label" << std::endl;
         return;
   }

   //temperary varibles for holding factors and terms
   double *temp_f = new double[N];
   int (*temp_t)[2][3] = new int[N][2][3];

   //Read the other data:
   for (int i = 0; i < N; ++i) {

      //get next line:
      getline(readfile,line);
      std::vector<std::string> tokens = split(line," ");
      
      //Store prefactor:
      temp_f[i] = str_to<double>(tokens[0]);

      //store the quantum numbers
      for (int j = 0; j < 6; ++j) {
      
         //Determine which particle. store numbers:
         if (j < 3) {
            temp_t[i][0][ j % 3 ] = str_to<int>(tokens[j+1]);
         }
         else {
            temp_t[i][1][ j % 3 ] = str_to<int>(tokens[j+1]);
         }
      }
   }
   
   switch (centre) {
      case 'T':
         terms_T = temp_t;
         factors_T = temp_f;
         break;
      case 'P':
         terms_P = temp_t;
         factors_P = temp_f;
         break;
      default:
         std::cout << "readinput error: Improper centre label" << std::endl;
         return;
   }

   readfile.close();
}


int WaveFunction::r_index( double x) {
   /*Returns the index of the value x in  the lattice table. If x not in r not sure.*/

   for (int i = 0; i < LEN; ++i) {
      if (r[i] == x) {
         return i;
      }
   }

   //If the loop exits:
   std::cout << "Error: invalid grid point" << std::endl;
   return 0;
}


double WaveFunction::R( const char centre, const char type, int n, int l, double x) {
   /*Returns the value of the R_nl function evaluated at x. x must be in r (ie a lattice point).*/
   
   return Rad[TP_toint(centre)][DG_toint(type)][n][l][r_index(x)];
}


std::complex<double> WaveFunction::a( const char centre, int i, int n, int l, int m ) {
   /*Return the coefficient a_centre[n][l][m].*/
   
   return amps[TP_toint(centre)][i][n][l][m];
}


double WaveFunction::Hlike( int z, double x ) {
   /*Returns the of the He 1+ ground state at the point x.*/

   return sqrt(4.0*z*z*z) * exp(-z*x) ;
}


double WaveFunction::indi_electron( const char centre, int i ) {
   /*Calculates the indipendent electron model ionization probability (p_centre)
   for impact parameter b[i].*/

   std::complex<double> p = std::complex<double>(0.0,0.0);

   for (int n = 1; n < MAXn; ++n) {
      for (int l = 0; l < n; ++l) {
         for (int m = 0; m < 2*l + 1; ++m) {
                  p += a(centre,i,n,l,m) * conj( a(centre,i,n,l,m) );
         }
      }
   }
   return  real( p );
}


std::unique_ptr<double[]> WaveFunction::integrand(const char c, int n1, int l1, int n2, int l2 , int s_n_1 , int s_l_1, int s_n_2 , int s_l_2) {
   /*Generates the values of the integrand for the given values.*/

   std::unique_ptr<double[]> result(new double[LIMIT]);

   result[0] = 0.0;

   for (int i = 1; i < LIMIT; ++i) {

      result[i] = r[i]*r[i] * R(c,'G', s_n_1 , s_l_1, r[i]) * R(c,'G', s_n_2 , s_l_2, r[i]) * ( R(c,'D',n1,l1,r[i])/R(c,'D',1,0,r[i]) ) * ( R(c,'D',n2,l2,r[i])/R(c,'D',1,0,r[i]) );
   }

   return result;
}


std::unique_ptr<double[]> WaveFunction::integrand_wb( const char c, double N, int n1, int l1, int n2, int l2 , int s_n_1 , int s_l_1, int s_n_2 , int s_l_2) {
   /*Generates the values of the WB integrand for the given values.*/

   std::unique_ptr<double[]> result(new double[LIMIT]);
   int z = charge(c);

   result[0] = 0.0;

   for (int i = 1; i < LIMIT; ++i) {

      result[i] = r[i]*r[i] * R(c,'G', s_n_1 , s_l_1, r[i]) * R(c,'G', s_n_2 , s_l_2, r[i]) * ( R(c,'D',n1,l1,r[i])  * R(c,'D',n2,l2,r[i]) / ( (2 - N) * Hlike(z,r[i])*Hlike(z,r[i]) + 2.0 * (N - 1) * R(c,'D',1,0,r[i]) * R(c,'D',1,0,r[i]) ) );
   }

   return result;
}


void WaveFunction::generate_integral_table( const char centre ) {
   /*Fills the integral table with integrals.*/
   
   int c = TP_toint(centre);

   for (int n1 = 1; n1 < 5; ++n1) {
      for (int l1 = 0; l1 < n1; ++l1) {
         for (int n2 = n1; n2 < 5; ++n2) {
            for (int l2 = 0; l2 < n2; ++l2) {
               for (int n1_gs = 1; n1_gs < 5; ++n1_gs) {
                  for (int l1_gs = 0; l1_gs < n1_gs; ++l1_gs) {
                     for (int n2_gs = n1_gs; n2_gs < 5; ++n2_gs) {
                        for (int l2_gs = 0; l2_gs < n2_gs; ++l2_gs) {
                        
                           std::unique_ptr<double[]> temp = integrand( centre, n1, l1, n2, l2, n1_gs, l1_gs, n2_gs, l2_gs);

                           int_table[c][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs] = newtoncotes(r,temp.get(), LIMIT);
                           int_table[c][n2][l2][n1][l1][n1_gs][l1_gs][n2_gs][l2_gs] = int_table[c][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           int_table[c][n1][l1][n2][l2][n2_gs][l2_gs][n1_gs][l1_gs] = int_table[c][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           int_table[c][n2][l2][n1][l1][n2_gs][l2_gs][n1_gs][l1_gs] = int_table[c][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           
                        }
                     }
                  }
               }
            }
         }            
      }
   }
}


void WaveFunction::generate_integral_table_wb( const char centre, double N  ) {
   /*Fills the WB integral table with integrals.*/
   
   int c = TP_toint(centre);

   for (int n1 = 1; n1 < 5; ++n1) {
      for (int l1 = 0; l1 < n1; ++l1) {
         for (int n2 = n1; n2 < 5; ++n2) {
            for (int l2 = 0; l2 < n2; ++l2) {
               for (int n1_gs = 1; n1_gs < 5; ++n1_gs) {
                  for (int l1_gs = 0; l1_gs < n1_gs; ++l1_gs) {
                     for (int n2_gs = n1_gs; n2_gs < 5; ++n2_gs) {
                        for (int l2_gs = 0; l2_gs < n2_gs; ++l2_gs) {
                        
                           std::unique_ptr<double[]> temp = integrand( centre, n1, l1, n2, l2, n1_gs, l1_gs, n2_gs, l2_gs);

                           int_table_wb[c][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs] = newtoncotes(r,temp.get(), LIMIT);
                           int_table_wb[c][n2][l2][n1][l1][n1_gs][l1_gs][n2_gs][l2_gs] = int_table_wb[c][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           int_table_wb[c][n1][l1][n2][l2][n2_gs][l2_gs][n1_gs][l1_gs] = int_table_wb[c][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           int_table_wb[c][n2][l2][n1][l1][n2_gs][l2_gs][n1_gs][l1_gs] = int_table_wb[c][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           
                        }
                     }
                  }
               }
            }
         }            
      }
   }
}


double WaveFunction::table( const char centre, int n1, int l1, int n2, int l2 , int n1_gs , int l1_gs, int n2_gs , int l2_gs ) {
   /*Returns the integral for the given radial function identifiers on centre.*/
   
   return int_table[TP_toint(centre)][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
}


double WaveFunction::table_wb( const char centre, int n1, int l1, int n2, int l2 , int n1_gs , int l1_gs, int n2_gs , int l2_gs ) {
   /*Returns the integral for the given radial function identifiers on centre.*/
   
   return int_table_wb[TP_toint(centre)][n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
}


double  WaveFunction::correlationintegral( const char c, int k ) {
   /*Calculates the correlation integral for the impact parameter b[i].*/

   //Holds the result:
   std::complex<double> Ic = std::complex<double>(0.0,0.0);

   //holds the results of the angular integrals:
   double G1;
   double G2;

   //Used for temperary storage:
   double phase, phase1, phase2, G3, G4;
   
   //Holds the ground state configuration information:
   int N;
   double *t;
   int (*T)[2][3];
   
   //Pick the ground state configuration data on cnetre c:
   switch (c) {
      case 'T':
         N = N_terms_T;
         t = factors_T;
         T = terms_T;
         break;
         
      case 'P':
         N = N_terms_P;
         t = factors_P;
         T = terms_P;
         break;
         
      default:
         std::cout << "Correlation integral error: Improper centre signifier" << std::endl;
         return 0;
   }

   //Sweep through the terms of the wave function:
   for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {

         //perform the angular intergal of the first particle by sweeping through the possible Gaunt integrals:
         for (int L1 = abs( T[i][0][1] - T[j][0][1] ); L1 <= T[i][0][1] + T[j][0][1]; ++L1 ) {
            for (int M1 = -L1; M1 <= L1; ++M1) {

               if ( IsEven( T[i][0][1] + T[j][0][1] + L1) && ! (TriangleBroken( T[i][0][1], T[j][0][1], L1 ) ) && ( T[i][0][2] - T[j][0][2] + M1 == 0) ) {

                  G1 = gaunt ( T[i][0][1], T[i][0][2], T[j][0][1], -T[j][0][2], L1, M1 );

                  //perform the angular intergal of the second particle by sweeping through the possible Gaunt integrals:
                  for (int L2 = abs( T[i][1][1] - T[j][1][1] ); L2 <= T[i][1][1] + T[j][1][1]; ++L2 ) {
                     for (int M2 = -L2; M2 <= L2; ++M2) {

                        if ( IsEven( T[i][1][1] + T[j][1][1] + L2) && ! (TriangleBroken( T[i][1][1], T[j][1][1], L2 ) ) && ( T[i][1][2] - T[j][1][2] + M1 == 0) ) {

                           G2 =  gaunt ( T[i][1][1], T[i][1][2], T[j][1][1], -T[j][1][2], L2, M2);

                           //calculate the pahse we neglected untill now:
                           phase = pow(-1.0, T[j][0][2] + T[j][1][2] );

                           //Perform the two sums from the P_lm functions:
                           for (int n1 = 1; n1 < MAXn; ++n1) {
                              for (int l1 = 0; l1 < n1; ++l1) {
                                 for (int np1 = 1; np1 < MAXn; ++np1) {
                                    for (int lp1 = 0; lp1 < np1; ++lp1) {
                                       for (int n2 = 1; n2 < MAXn; ++n2) {
                                          for (int l2 = 0; l2 < n2; ++l2) {
                                             for (int np2 = 1; np2 < MAXn; ++np2) {
                                                for (int lp2 = 0; lp2 < np2; ++lp2) {
                                                   //Store the value of the products of the complete integrals:

                                                   double temp = phase * G1 * t[i] * table( c, n1, l1, np1, lp1, T[i][0][0], T[i][0][1], T[j][0][0] , T[j][0][1] ) * G2 * t[j] * table( c, n2, l2, np2, lp2, T[i][1][0], T[i][1][1], T[j][1][0] , T[j][1][1] );

                                                   for (int m1 = 0; m1 < 2*l1 + 1; ++m1) {
                                                      for (int mp1 = 0; mp1 < 2*lp1 + 1; ++mp1) {

                                                         if ( IsEven( l1 + lp1 + L1) && ! (TriangleBroken(l1, lp1, L1) ) && (m1 - l1 -(mp1 - lp1) - M1 == 0) ) { //here

                                                            G3 = gaunt( l1, m1 - l1, lp1, -(mp1 - lp1), L1, -M1 );

                                                            for (int m2 = 0; m2 < 2*l2 + 1; ++m2) {
                                                               for (int mp2 = 0; mp2 < 2*lp2 + 1; ++mp2) {

                                                                  if ( IsEven( l2 + lp2 + L2) && ! (TriangleBroken(l2, lp2, L2) ) && (m2 - l2 -(mp2 - lp2) - M2 == 0) ) { //here

                                                                     G4 = gaunt( l2, m2 - l2, lp2, -(mp2 - lp2), L2, -M2 );

                                                                     phase1 = pow(-1.0, mp1 - lp1 + M1);
                                                                     phase2 = pow(-1.0, mp2 - lp2 + M2);

                                                                     Ic += a(c,k,n1,l1,m1) * conj( a(c,k,np1,lp1,mp1) ) * a(c,k,n2,l2,m2) * conj( a(c,k,np2,lp2,mp2) ) * temp * phase1 * G3 * phase2 * G4;
                                                                  }
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return   32.0 * PI * PI * real( Ic ) ;
}


double  WaveFunction::correlationintegral_wb( const char c, int k ) {
   /*Calculates the correlation integral for the impact parameter b[i].*/
   
   //Holds the result:
   std::complex<double> Ic = std::complex<double>(0.0,0.0);

   //holds the results of the angular integrals:
   double G1;
   double G2;

   //Used for temperary storage:
   double phase, phase1, phase2, G3, G4;

   //Calculate the fractional ionization N(t):
   double N_e = 2.0 * indi_electron(c,k);
   
   //Holds the ground state configuration information:
   int N;
   double *t;
   int (*T)[2][3];
   
   //Pick the ground state configuration data on cnetre c:
   switch (c) {
      case 'T':
         N = N_terms_T;
         t = factors_T;
         T = terms_T;
         break;
         
      case 'P':
         N = N_terms_P;
         t = factors_P;
         T = terms_P;
         break;
         
      default:
         std::cout << "Correlation integral error: Improper centre signifier" << std::endl;
         return 0;
   }

   if ( N_e <= 1.0 ) {
      return 0.0;
   }
   else {

      //Generate the integral table for the given impact parameter/N value:
      generate_integral_table_wb(c,N_e);
      
      //Sweep through the terms of the wave function:
      for (int i = 0; i < N; ++i) {
         for (int j = 0; j < N; ++j) {

            //perform the angular intergal of the first particle by sweeping through the possible Gaunt integrals:
            for (int L1 = abs( T[i][0][1] - T[j][0][1] ); L1 <= T[i][0][1] + T[j][0][1]; ++L1 ) {
               for (int M1 = -L1; M1 <= L1; ++M1) {

                  if ( IsEven( T[i][0][1] + T[j][0][1] + L1) && ! (TriangleBroken( T[i][0][1], T[j][0][1], L1 ) ) && ( T[i][0][2] - T[j][0][2] + M1 == 0) ) {

                     G1 = gaunt ( T[i][0][1], T[i][0][2], T[j][0][1], -T[j][0][2], L1, M1 );

                     //perform the angular intergal of the second particle by sweeping through the possible Gaunt integrals:
                     for (int L2 = abs( T[i][1][1] - T[j][1][1] ); L2 <= T[i][1][1] + T[j][1][1]; ++L2 ) {
                        for (int M2 = -L2; M2 <= L2; ++M2) {

                           if ( IsEven( T[i][1][1] + T[j][1][1] + L2) && ! (TriangleBroken( T[i][1][1], T[j][1][1], L2 ) ) && ( T[i][1][2] - T[j][1][2] + M1 == 0) ) {

                              G2 =  gaunt ( T[i][1][1], T[i][1][2], T[j][1][1], -T[j][1][2], L2, M2);

                              //calculate the pahse we neglected untill now:
                              phase = pow(-1.0, T[j][0][2] + T[j][1][2] );

                              //Perform the two sums from the P_lm functions:
                              for (int n1 = 1; n1 < MAXn; ++n1) {
                                 for (int l1 = 0; l1 < n1; ++l1) {
                                    for (int np1 = 1; np1 < MAXn; ++np1) {
                                       for (int lp1 = 0; lp1 < np1; ++lp1) {
                                          for (int n2 = 1; n2 < MAXn; ++n2) {
                                             for (int l2 = 0; l2 < n2; ++l2) {
                                                for (int np2 = 1; np2 < MAXn; ++np2) {
                                                   for (int lp2 = 0; lp2 < np2; ++lp2) {

                                                      //Store the value of the products of the complete integrals:
                                                      double temp = phase * G1 * t[i] * table_wb( c, n1, l1, np1, lp1, T[i][0][0], T[i][0][1], T[j][0][0] , T[j][0][1] )  * G2 * t[j] * table_wb( c, n2, l2, np2, lp2, T[i][1][0], T[i][1][1], T[j][1][0] , T[j][1][1] );

                                                      for (int m1 = 0; m1 < 2*l1 + 1; ++m1) {
                                                         for (int mp1 = 0; mp1 < 2*lp1 + 1; ++mp1) {

                                                            if ( IsEven( l1 + lp1 + L1) && ! (TriangleBroken(l1, lp1, L1) ) && (m1 - l1 -(mp1 - lp1) - M1 == 0) ) { //here

                                                               G3 = gaunt( l1, m1 - l1, lp1, -(mp1 - lp1), L1, -M1 );

                                                               for (int m2 = 0; m2 < 2*l2 + 1; ++m2) {
                                                                  for (int mp2 = 0; mp2 < 2*lp2 + 1; ++mp2) {

                                                                     if ( IsEven( l2 + lp2 + L2) && ! (TriangleBroken(l2, lp2, L2) ) && (m2 - l2 -(mp2 - lp2) - M2 == 0) ) { //here

                                                                        G4 = gaunt( l2, m2 - l2, lp2, -(mp2 - lp2), L2, -M2 );

                                                                        phase1 = pow(-1.0, mp1 - lp1 + M1);
                                                                        phase2 = pow(-1.0, mp2 - lp2 + M2);

                                                                        Ic += a(c,k,n1,l1,m1) * conj( a(c,k,np1,lp1,mp1) ) * a(c,k,n2,l2,m2) * conj( a(c,k,np2,lp2,mp2) ) * temp * phase1 * G3 * phase2 * G4;
                                                                     }
                                                                  }
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      return   128.0 * ( N_e - 1 ) * PI * PI * real( Ic );
   }
}
