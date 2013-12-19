
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include "newtoncotes.hpp"
#include "wavefunction.hpp"

using namespace std;

int main() {

   //Create an array to hold the impact energies:
   int E = 100;
   bool WB = true;

   cout << "\n" << endl;

   cout.precision(10);

   time_t begin, setup, t1, t2, end;

   time(&begin);

   std::stringstream ss;

   std::string Z   = "/home/baxter/Documents/TPZ_E";
   std::string EQ4 = "/home/baxter/Documents/TPEQ4_E";
   std::string HL  = "/home/baxter/Documents/TPHL_E";

	ss << E << ".txt";
   EQ4 += ss.str();
   Z += ss.str();
   HL += ss.str();

   //Open file to be written to:
   ofstream writeZ  ( Z.c_str()   );
   ofstream writeEQ4( EQ4.c_str() );
   ofstream writeHL ( HL.c_str()  );

   string header = "b p_T p_P Ic_TT Ic_PP Ic_TP ie_pTT ie_pTI ie_pII ie_pTP ie_pIP ie_pPP p_TT p_TI p_II p_TP p_PI p_PP";
   writeZ   << header << endl;
   writeEQ4 << header << endl;
   writeHL  << header << endl;

   cout << "Starting E = " << E << " keV" << endl;

   //Create Wave function:
   WaveFunction wfn( E , WB );

   time(&setup);

   cout << "setup done (s): " << difftime(setup,begin) << endl;

   for (int j = 0; j < 30; ++j) {

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

      //Probabilities that do not change:
      p_TT = Ic_TT/2.0;
      p_PP = Ic_PP/2.0;

      //I_TP = 0:
      Ic_TP = 0;

      //"Exact" probabilities:
      p_TI = 2.0 * p_T * ( 1 - p_P ) - Ic_TT;
      p_II = 1.0 - 2.0 * p_T - 2.0 * p_P + 2.0 * p_T * p_P + Ic_TT/2.0 + Ic_PP/2.0;
      p_TP = 2.0 * p_P * p_T;
      p_PI = 2.0 * p_P * ( 1 - p_T) - Ic_PP;

      //Write to file:
      writeZ << wfn.get_impact(j) << " " << p_T    << " " << p_P    << " " << Ic_TT  << " "
             << Ic_PP             << " " << Ic_TP  << " " << ie_pTT << " " << ie_pTI << " "
             << ie_pII            << " " << ie_pTP << " " << ie_pIP << " " << ie_pPP << " "
             << p_TT              << " " << p_TI   << " " << p_II   << " " << p_TP   << " "
             << p_PI              << " " << p_PP   << endl;

      //EQ (4):
      Ic_TP = 2 * p_P * ( 1 - p_T ) - Ic_PP;

      //"Exact" probabilities:
      p_TI = 2.0 * ( p_T - p_P ) + Ic_PP- Ic_TT;
      p_II = 1.0 - 2.0 * p_T + Ic_TT/2.0 - Ic_PP/2.0;
      p_TP = 2.0 * p_P - Ic_PP;
      p_PI = 0.0;

      //Write to file:
      writeEQ4 << wfn.get_impact(j) << " " << p_T    << " " << p_P    << " " << Ic_TT  << " "
               << Ic_PP             << " " << Ic_TP  << " " << ie_pTT << " " << ie_pTI << " "
               << ie_pII            << " " << ie_pTP << " " << ie_pIP << " " << ie_pPP << " "
               << p_TT              << " " << p_TI   << " " << p_II   << " " << p_TP   << " "
               << p_PI              << " " << p_PP   << endl;
		
      //Heitler-London TP integral:
      if ( p_T + p_P > .5 ) {
         Ic_TP = (2*p_P*p_T)/( p_T + p_P - 0.5);
      }
      else {
         Ic_TP = 0;
      }

      //Calculate "exact" probabilities:
      p_TI = 2.0 * p_T - Ic_TT - Ic_TP;
      p_II = 1.0 - 2.0 * p_T - 2.0 * p_P + Ic_TT/2.0 + Ic_TP + Ic_PP/2.0;
      p_TP = Ic_TP;
      p_PI = 2.0 * p_P - Ic_PP - Ic_TP;

      //Write to file:
      writeHL  << wfn.get_impact(j) << " " << p_T    << " " << p_P    << " " << Ic_TT  << " "
               << Ic_PP             << " " << Ic_TP  << " " << ie_pTT << " " << ie_pTI << " "
               << ie_pII            << " " << ie_pTP << " " << ie_pIP << " " << ie_pPP << " "
               << p_TT              << " " << p_TI   << " " << p_II   << " " << p_TP   << " "
               << p_PI              << " " << p_PP   << endl;

      time(&t2);

      cout << "done #" << j + 1 << " (s): " << difftime(t2,t1) << endl;
   }

   writeZ.close();
   writeEQ4.close();
   writeZ.close();writeHL.close();

   time(&end);

   cout << "\ntotal time (s): " << difftime(end,begin) << endl;
   cout << "\n" << endl;

   return 0;
}
