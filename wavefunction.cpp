
#include "wavefunction.hpp"

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982;

WaveFunction::WaveFunction( char* inputpath , int energy ) {
   /*Constructor for objects of class WaveFunction.*/

   //initialize tables:
   readlattice();
   readradial();
   readradial_gs();
   readamplitudes( energy );
   readinput( inputpath );

   //Calculate table of integrals
   generate_integral_table();

}


WaveFunction::~WaveFunction() {
   /*Destroys an object of class WaveFunction.*/

   delete [] terms;
   delete [] factors;
}


float WaveFunction::get_impact( int i) {

   return b[i];
}


void WaveFunction::readlattice() {
   /*Reads/stores the file containing the lattice points.*/

   std::string line;
   const char* path;

   //Open file to be read:
   //path = "/home/baxter/Documents/Observable/Input/radial/silverman_1s1s/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/silverman full/lattice.txt";
     path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-2p/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s(all)/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d(all)/lattice.txt";
   //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4f/lattice.txt";

   std::ifstream readfile(path);

   r[0] = 0.0;

   for (int i = 1; i < LEN; ++i) {

      //Read a line:
      getline(readfile,line);
      std::stringstream ss;
      ss << line;

      //store:
      ss >> r[i];

   }
   readfile.close();
}


void WaveFunction::readradial() {
   /*Reads/stores the R_{nl} function data.*/

   std::string path;
   std::string line;

   for (int n = 1; n < 5; ++n) {
      for (int l = 0; l < n; ++l) {
         //path = "/home/baxter/Documents/Observable/Input/radial/silverman_1s1s/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/silverman full/R_";
           path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-2p/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4s(all)/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4d(all)/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/MCHF_1s-4f/R_";

         std::stringstream ss;
         ss << n << l << ".txt";
         path += ss.str();

         //convert string to char*:
         char* pth = new char[path.size()+1];
         strcpy(pth,path.c_str());

         //Open file:
         std::ifstream readfile(pth);

         delete[] pth;

         //Give R_nl(0) some value (doesn't really matter what it is as long as its not zero:
         Rad[n][l][0] = 1.0;

         for (int i = 1; i < LEN; ++i) {
            //get a line
            getline(readfile,line);
            std::stringstream ssline;
            ssline << line;

            //Store:
            ssline >> Rad[n][l][i];
         }
         readfile.close();
      }
   }
}


void WaveFunction::readradial_gs() {
   /*Reads/stores the R_{nl} function data.*/

   std::string path;
   std::string line;

   for (int n = 1; n < 5; ++n) {
      for (int l = 0; l < n; ++l) {
         //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/silverman_1s1sp/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/silverman_full/R_";
           path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-2p/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4s/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4s(all)/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4d/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4d(all)/R_";
         //path = "/home/baxter/Documents/Observable/Input/radial/Ground state/MCHF_1s-4f/R_";

         std::stringstream ss;
         ss << n << l << ".txt";
         path += ss.str();

         //convert string to char*:
         char* pth = new char[path.size()+1];
         strcpy(pth,path.c_str());

         //Open file:
         std::ifstream readfile(pth);

         delete[] pth;

         //Give R_nl(0) some value (doesn't really matter what it is as long as its not zero:
         Rad_gs[n][l][0] = 1.0;

         for (int i = 1; i < LEN; ++i) {
            //get a line
            getline(readfile,line);
            std::stringstream ssline;
            ssline << line;

            //Store:
            ssline >> Rad_gs[n][l][i];

         }
         readfile.close();
      }
   }
}


void WaveFunction::readamplitudes( int energy ) {
   /*Reads in the impact parameter and amplitude data and stores it in b and amp. The energy parameter must be 100 or 2000.*/

   std::string line;
   std::stringstream ss;
   char* tok1;
   char* tok2;

   //setup the path string:
     std::string path = "/home/baxter/Documents/Observable/Input/amplitudes/old/E";
   //std::string path = "/home/baxter/Documents/Observable/Input/amplitudes/OPM/E";
   //std::string path = "/home/baxter/Documents/Observable/Input/amplitudes/MCHF_1s4d/E";
   //std::string path = "/home/baxter/Documents/Observable/Input/amplitudes/MCHF_1s4d(all)/E";
   ss << energy << ".txt";
   path += ss.str();

   //convert path to c_string.
   char* pth = new char[path.size()+1];
   strcpy(pth,path.c_str());

   //open file.
   std::ifstream readfile(pth);

   delete[] pth;

   for (int i = 0; i < 30; ++i) {

      //read impact parameter line.
      getline(readfile,line);

      //convert line to c_string.
      char* cstr_line = new char[line.size()+1];
      strcpy(cstr_line,line.c_str());

      //extract token containing impact parameter.
      strtok(cstr_line, "=");
      std::stringstream ss1;
      ss1 << strtok(NULL, " ");
      ss1 >> b[i];

      //store value in b.

      for (int n = 1; n < 5; ++n) {
         for (int l = 0; l < n; ++l) {
            for (int m = 0; m < 2*l + 1; ++m) {

               //get next line.
               getline(readfile,line);

               //convert line to c_string.
               cstr_line = new char[line.size()+1];
               strcpy(cstr_line,line.c_str());

               //get first token tok1.
               tok1 = strtok(cstr_line, "\t");

               //convert it to a double:
               std:: stringstream real;
               real << tok1;
               double realpart;
               real >> realpart;

               //get second token tok2.
               tok2 = strtok(NULL, "\t");

               //convert it to a double:
               std:: stringstream im;
               im << tok2;
               double impart;
               im >> impart;

               //convert to complex(tok1,tok2) and store in amps[i][n][l][m].
               amps[i][n][l][m] = std::complex<double>( realpart, impart );
            }
         }
      }
      delete[] cstr_line;
   }
   readfile.close();
}


void WaveFunction::readinput( char* inputpath ) {
   /*Reads and stores the wave function quantum number data from the file located at inputpath.*/

   std::ifstream readfile ( inputpath );
   std::string line;
   char* tok;

   //Read the number of terms:
   getline(readfile,line);
   N_terms = atoi(line.c_str());

   //temperary varibles for holding factors and terms
   double *temp_f = new double[N_terms];
   int (*temp_t)[2][3] = new int[N_terms][2][3];

   //Read the other data:
   for (int i = 0; i < N_terms; ++i) {

      //get next line:
      getline(readfile,line);

      //convert line to c_string.
      char* cstr_line = new char[line.size()+1];
      strcpy(cstr_line,line.c_str());

      //get first token tok1.
      tok = strtok(cstr_line, " ");

      //store prefactor:
      std::stringstream ss1;
      ss1 << tok;
      ss1 >> temp_f[i];

      //store the quantum numbers
      for (int j = 0; j < 6; ++j) {

         //get next token.
         tok = strtok(NULL, " \n");

         //Determine which particle. store numbers:
         if (j < 3) {
            temp_t[i][0][ j % 3 ] = atoi(tok);
         }
         else {
            temp_t[i][1][ j % 3 ] = atoi(tok);
         }
      }
      delete[] cstr_line;
   }
   
   terms = temp_t;
   factors = temp_f;

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


double WaveFunction::R(int n, int l, double x) {
   /*Returns the value of the R_nl function evaluated at x. x must be in r (ie a lattice point).*/

   return Rad[n][l][r_index(x)];
}


double WaveFunction::R_gs(int n, int l, double x) {
   /*Returns the value of the R_nl function evaluated at x. x must be in r (ie a lattice point).*/

   return Rad_gs[n][l][r_index(x)];
}


double WaveFunction::Hep1( double x ) {
   /*Returns the of the He 1+ ground state at the point x.*/

   return sqrt( 32.0 ) * exp( - 2.0 * x ) ;
}


double* WaveFunction::integrand(int n1, int l1, int n2, int l2 , int s_n_1 , int s_l_1, int s_n_2 , int s_l_2) {
   /*Generates the values of the integrand for the given values.*/

   static double result[LIMIT];

   result[0] = 0.0;

   for (int i = 1; i < LIMIT; ++i) {

      result[i] = r[i]*r[i] * R_gs( s_n_1 , s_l_1, r[i]) * R_gs( s_n_2 , s_l_2, r[i]) * ( R(n1,l1,r[i])/R(1,0,r[i]) ) * ( R(n2,l2,r[i])/R(1,0,r[i]) );
   }

   return result;
}


double* WaveFunction::integrand_wb( double N, int n1, int l1, int n2, int l2 , int s_n_1 , int s_l_1, int s_n_2 , int s_l_2) {
   /*Generates the values of the WB integrand for the given values.*/

   static double result[LIMIT];

   result[0] = 0.0;

   for (int i = 1; i < LIMIT; ++i) {

      result[i] = r[i]*r[i] * R_gs( s_n_1 , s_l_1, r[i]) * R_gs( s_n_2 , s_l_2, r[i]) * ( R(n1,l1,r[i])  * R(n2,l2,r[i]) / ( (2 - N) * Hep1(r[i])*Hep1(r[i]) + 2.0 * (N - 1) * R(1,0,r[i]) * R(1,0,r[i]) ) );
   }

   return result;
}


void WaveFunction::generate_integral_table() {
   /*Fills the integral table with integrals.*/

   for (int n1 = 1; n1 < 5; ++n1) {
      for (int l1 = 0; l1 < n1; ++l1) {
         for (int n2 = n1; n2 < 5; ++n2) {
            for (int l2 = 0; l2 < n2; ++l2) {
               for (int n1_gs = 1; n1_gs < 5; ++n1_gs) {
                  for (int l1_gs = 0; l1_gs < n1_gs; ++l1_gs) {
                     for (int n2_gs = n1_gs; n2_gs < 5; ++n2_gs) {
                        for (int l2_gs = 0; l2_gs < n2_gs; ++l2_gs) {

                           int_table[n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs] = newtoncotes(r,integrand( n1, l1, n2, l2, n1_gs, l1_gs, n2_gs, l2_gs), LIMIT);
                           int_table[n2][l2][n1][l1][n1_gs][l1_gs][n2_gs][l2_gs] = int_table[n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           int_table[n1][l1][n2][l2][n2_gs][l2_gs][n1_gs][l1_gs] = int_table[n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           int_table[n2][l2][n1][l1][n2_gs][l2_gs][n1_gs][l1_gs] = int_table[n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];

                        }
                     }
                  }
               }
            }
         }            
      }
   }
}


void WaveFunction::generate_integral_table_wb( double N  ) {
   /*Fills the WB integral table with integrals.*/

   for (int n1 = 1; n1 < 5; ++n1) {
      for (int l1 = 0; l1 < n1; ++l1) {
         for (int n2 = n1; n2 < 5; ++n2) {
            for (int l2 = 0; l2 < n2; ++l2) {
               for (int n1_gs = 1; n1_gs < 5; ++n1_gs) {
                  for (int l1_gs = 0; l1_gs < n1_gs; ++l1_gs) {
                     for (int n2_gs = n1_gs; n2_gs < 5; ++n2_gs) {
                        for (int l2_gs = 0; l2_gs < n2_gs; ++l2_gs) {

                           int_table_wb[n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs] = newtoncotes(r,integrand_wb( N, n1, l1, n2, l2, n1_gs, l1_gs, n2_gs, l2_gs), LIMIT);
                           int_table_wb[n2][l2][n1][l1][n1_gs][l1_gs][n2_gs][l2_gs] = int_table_wb[n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           int_table_wb[n1][l1][n2][l2][n2_gs][l2_gs][n1_gs][l1_gs] = int_table_wb[n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                           int_table_wb[n2][l2][n1][l1][n2_gs][l2_gs][n1_gs][l1_gs] = int_table_wb[n1][l1][n2][l2][n1_gs][l1_gs][n2_gs][l2_gs];
                        }
                     }
                  }
               }
            }
         }            
      }
   }
}


double WaveFunction::indi_electron( int i ) {
   /*Calculates the indipendent electron model ionization probability for impact parameter b[i].*/

   std::complex<double> p = std::complex<double>(0.0,0.0);

   for (int n = 1; n < MAXn; ++n) {
      for (int l = 0; l < n; ++l) {
         for (int m = 0; m < 2*l + 1; ++m) {
            p += amps[i][n][l][m] * conj( amps[i][n][l][m] );
         }
      }
   }

   return  real( 1.0 - p );
}


double  WaveFunction::correlationintegral( int a ) {
   /*Calculates the correlation integral for the impact parameter b[i].*/

   //Holds the result:
   std::complex<double> Ic = std::complex<double>(0.0,0.0);

   //holds the results of the angular integrals:
   double G1;
   double G2;

   //Used for temperary storage:
   double phase, phase1, phase2, G3, G4;

   //Sweep through the terms of the wave function:
   for (int i = 0; i < N_terms; ++i) {
      for (int j = 0; j < N_terms; ++j) {

         //perform the angular intergal of the first particle by sweeping through the possible Gaunt integrals:
         for (int L1 = abs( terms[i][0][1] - terms[j][0][1] ); L1 <= terms[i][0][1] + terms[j][0][1]; ++L1 ) {
            for (int M1 = -L1; M1 <= L1; ++M1) {

               if ( IsEven( terms[i][0][1] + terms[j][0][1] + L1) && ! (TriangleBroken( terms[i][0][1], terms[j][0][1], L1 ) ) && (terms[i][0][2] -terms[j][0][2] + M1 == 0) ) { //here

                  G1 = gaunt (terms[i][0][1], terms[i][0][2], terms[j][0][1], -terms[j][0][2], L1, M1);

                  //perform the angular intergal of the second particle by sweeping through the possible Gaunt integrals:
                  for (int L2 = abs( terms[i][1][1] - terms[j][1][1] ); L2 <= terms[i][1][1] + terms[j][1][1]; ++L2 ) {
                     for (int M2 = -L2; M2 <= L2; ++M2) {

                        if ( IsEven( terms[i][1][1] + terms[j][1][1] + L2) && ! (TriangleBroken( terms[i][1][1], terms[j][1][1], L2 ) ) && (terms[i][1][2] -terms[j][1][2] + M1 == 0) ) { //here

                           G2 =  gaunt (terms[i][1][1], terms[i][1][2], terms[j][1][1], -terms[j][1][2], L2, M2);

                           //calculate the pahse we neglected untill now:
                           phase = pow(-1.0, terms[j][0][2] + terms[j][1][2] );

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

                                                   double temp = phase * G1 * factors[i] * int_table[n1][l1][np1][lp1][ terms[i][0][0] ][ terms[i][0][1] ][ terms[j][0][0] ][ terms[j][0][1] ] * G2 * factors[j] * int_table[n2][l2][np2][lp2][ terms[i][1][0] ][ terms[i][1][1] ][ terms[j][1][0] ][ terms[j][1][1] ];

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

                                                                     Ic += amps[a][n1][l1][m1] * conj( amps[a][np1][lp1][mp1] ) * amps[a][n2][l2][m2] * conj( amps[a][np2][lp2][mp2] ) * temp * phase1 * G3 * phase2 * G4;
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


double  WaveFunction::correlationintegral_wb( int a ) {
   /*Calculates the correlation integral for the impact parameter b[i].*/
   
   //Holds the result:
   std::complex<double> Ic = std::complex<double>(0.0,0.0);

   //holds the results of the angular integrals:
   double G1;
   double G2;

   //Used for temperary storage:
   double phase, phase1, phase2, G3, G4;

   //Calculate the fractional ionization N(t):
   double N = 2.0 * ( 1.0 - indi_electron(a) );

   if ( N <= 1.0    ) {
      return 0.0;
   }
   else {

      //Generate the integral table for the given impact parameter/N value:
      generate_integral_table_wb(N);

      //Sweep through the terms of the wave function:
      for (int i = 0; i < N_terms; ++i) {
         for (int j = 0; j < N_terms; ++j) {

            //perform the angular intergal of the first particle by sweeping through the possible Gaunt integrals:
            for (int L1 = abs( terms[i][0][1] - terms[j][0][1] ); L1 <= terms[i][0][1] + terms[j][0][1]; ++L1 ) {
               for (int M1 = -L1; M1 <= L1; ++M1) {

                  if ( IsEven( terms[i][0][1] + terms[j][0][1] + L1) && ! (TriangleBroken( terms[i][0][1], terms[j][0][1], L1 ) ) && (terms[i][0][2] -terms[j][0][2] + M1 == 0) ) { //here

                     G1 = gaunt (terms[i][0][1], terms[i][0][2], terms[j][0][1], -terms[j][0][2], L1, M1);

                     //perform the angular intergal of the second particle by sweeping through the possible Gaunt integrals:
                     for (int L2 = abs( terms[i][1][1] - terms[j][1][1] ); L2 <= terms[i][1][1] + terms[j][1][1]; ++L2 ) {
                        for (int M2 = -L2; M2 <= L2; ++M2) {

                           if ( IsEven( terms[i][1][1] + terms[j][1][1] + L2) && ! (TriangleBroken( terms[i][1][1], terms[j][1][1], L2 ) ) && (terms[i][1][2] -terms[j][1][2] + M1 == 0) ) { //here

                              G2 =  gaunt (terms[i][1][1], terms[i][1][2], terms[j][1][1], -terms[j][1][2], L2, M2);

                              //calculate the pahse we neglected untill now:
                              phase = pow(-1.0, terms[j][0][2] + terms[j][1][2] );

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

                                                      double temp = phase * G1 * factors[i] * int_table_wb[n1][l1][np1][lp1][ terms[i][0][0] ][ terms[i][0][1] ][ terms[j][0][0] ][ terms[j][0][1] ] * G2 * factors[j] * int_table_wb[n2][l2][np2][lp2][ terms[i][1][0] ][ terms[i][1][1] ][ terms[j][1][0] ][ terms[j][1][1] ];

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

                                                                        Ic += amps[a][n1][l1][m1] * conj( amps[a][np1][lp1][mp1] ) * amps[a][n2][l2][m2] * conj( amps[a][np2][lp2][mp2] ) * temp * phase1 * G3 * phase2 * G4;
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

      return   128.0 * ( N - 1 ) * PI * PI * real( Ic );
   }
}
