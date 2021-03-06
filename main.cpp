
#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include "newtoncotes.hpp"
#include "wavefunction.hpp"

using namespace std;

int main( int argc, char* argv[] ) {

   //Create an array to hold the impact energies:
   int E = stoi(argv[1]);
   int z = stoi(argv[2]);
   bool WB = true;

   cout << "\n\n";

   cout.precision(10);

   time_t begin, setup, t1, t2, end;

   time(&begin);

   std::stringstream ss;
   ss << "./output/MCHF-un/" << "z" << z << "/He2pHe_MCHF_resp_E" << E << "_z" << z << ".txt";

   //Open file to be written to:
   ofstream writefile  ( ss.str().c_str() );
   writefile << "b p_T p_P Ic_TT Ic_PP Ic_TP ie_pTT ie_pTI ie_pII ie_pTP ie_pIP ie_pPP p_TT p_TI p_II p_TP p_PI p_PP\n";

   cout << "Starting E = " << E << " keV\n";

   //Create Wave function:
   WaveFunction wfn( E, z, WB );

   time(&setup);

   cout << "setup done (s): " << difftime(setup,begin) << "\n";

   for (int j = 0; j < wfn.get_Nb(); ++j) {

      time(&t1);
      double Ic_TT, Ic_PP, Ic_TP, p_TT, p_TI, p_II, p_TP, p_PI, p_PP;

      //Calculate IEM p:
      double p_T = wfn.indi_electron( 'T', j );
      double p_P = wfn.indi_electron( 'P', j );

      //Calculate IEM probabilities:
      double ie_pTT = p_T * p_T; 
      double ie_pTI = 2.0 * p_T * ( 1.0 - p_T - p_P );
      double ie_pII = ( 1.0 - p_T - p_P ) * ( 1.0 - p_T - p_P );
      double ie_pTP = 2.0 * p_T * p_P;
      double ie_pIP = 2.0 * p_P * ( 1.0 - p_T - p_P );
      double ie_pPP = p_P * p_P;

      //Calculate the Ic's:
      if ( WB ) {
         Ic_TT = wfn.correlationintegral_wb( 'T', j );
         Ic_PP = wfn.correlationintegral_wb( 'P', j );
      }
      else {
         Ic_TT = wfn.correlationintegral( 'T', j );
         Ic_PP = wfn.correlationintegral( 'P', j );
      }
      
      Ic_TP = corrint_TP( p_T, p_P, Ic_TT, Ic_PP );

      //"Exact" probabilities:
      p_TT = 0.5 * Ic_TT;
      p_TI = 2.0 * p_T - Ic_TT - Ic_TP;
      p_II = 1.0 - 2.0 * p_T - 2.0 * p_P + 0.5 * Ic_TT + Ic_TP + 0.5 * Ic_PP;
      p_TP = Ic_TP;
      p_PI = 2.0 * p_P - Ic_PP - Ic_TP;
      p_PP = 0.5 * Ic_PP;
      
      //Write to file:
      writefile << wfn.get_impact(j) << " " << p_T    << " " << p_P    << " " << Ic_TT  << " "
                << Ic_PP             << " " << Ic_TP  << " " << ie_pTT << " " << ie_pTI << " "
                << ie_pII            << " " << ie_pTP << " " << ie_pIP << " " << ie_pPP << " "
                << p_TT              << " " << p_TI   << " " << p_II   << " " << p_TP   << " "
                << p_PI              << " " << p_PP   << "\n";

      time(&t2);

      cout << "done #" << j + 1 << " (s): " << difftime(t2,t1) << "\n";
   }

   writefile.close();

   time(&end);

   cout << "\ntotal time (s): " << difftime(end,begin) << "\n";
   cout << "\n" << endl;

   return 0;
}
