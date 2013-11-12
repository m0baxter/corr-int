
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include "wavefunction.hpp"
#include <complex>

using namespace std;

int main() {

	//Create an array to hold the impact energies:
	int E;

	cout << "\n" << endl;

	cout.precision(10);

	//setup gauss method:
	gauss_setup();

	time_t begin, setup, t1, t2, end;

	cout << "Enter impact energy:  " << endl;
	cin >> E;

	time(&begin);

	std::stringstream ss;

	std::string path = "/home/baxter/Documents/test_std_E";
	ss << E << ".txt";
	path += ss.str();

	//Open file to be written to:
	ofstream writefile( path.c_str() );
	writefile << "b p_T p_P Ic_TT Ic_PP Ic_TP ie_pTT ie_pTI ie_pII ie_pTP ie_pIP ie_pPP p_TT p_TI p_II p_TP p_PI p_PP\n";

	cout << "Starting E = " << E << " keV" << endl;
	
	//Input paths:
	string input_T = "/home/baxter/Documents/Observable/Input/wfn/MCHF_1s-4d(all).txt";
	string input_P = "/home/baxter/Documents/Observable/Input/wfn/MCHF_1s-4d(all).txt"; 

	//Create Wave function:
	WaveFunction MCHF_1s4d ( input_T.c_str(), input_P.c_str() , E , true );

	time(&setup);

	cout << "setup done (s): " << difftime(setup,begin) << endl;

	for (int j = 0; j < 1; ++j) {

		time(&t1);

		//Calculate the Ic's:
		double Ic_TT = MCHF_1s4d.correlationintegral_wb( 'T', j );
		double Ic_PP = MCHF_1s4d.correlationintegral_wb( 'P', j );
		double Ic_TP = 0;

		//Calculate IEM p:
		double p_T = MCHF_1s4d.indi_electron( 'T', j );
		double p_P = MCHF_1s4d.indi_electron( 'P', j );

		//Calculate IEM probabilities:
		double ie_pTT = p_T * p_T; 
		double ie_pTI = 2.0 * p_T * ( 1.0 - p_T - p_P );
		double ie_pII = ( 1.0 - p_T - p_P ) * ( 1.0 - p_T - p_P );
		double ie_pTP = 2.0 * p_T * p_P;
		double ie_pIP = 2.0 * p_P * ( 1.0 - p_T - p_P );
      double ie_pPP = p_P * p_P;
		
		//Calculate "exact" probabilities:
		double p_TT = Ic_TT/2.0;
		double p_TI = 2.0 * p_T - Ic_TT - Ic_TP;
		double p_II = 1.0 - 2.0 * p_T - 2.0 * p_P + Ic_TT/2.0 + Ic_TP + Ic_PP/2.0;
		double p_TP = Ic_TP;
		double p_PI = 2.0 * p_P - Ic_PP - Ic_TP;
		double p_PP = Ic_PP/2.0;

		writefile << MCHF_1s4d.get_impact(j) << " " << p_T    << " " << p_P    << " " << Ic_TT  << " "
		          << Ic_PP                   << " " << Ic_TP  << " " << ie_pTT << " " << ie_pTI << " "
		          << ie_pII                  << " " << ie_pTP << " " << ie_pIP << " " << ie_pPP << " "
		          << p_TT                    << " " << p_TI   << " " << p_II   << " " << p_TP   << " "
		          << p_PI                    << " " << p_PP <<"\n";

		time(&t2);

		cout << "done #" << j + 1 << " (s): " << difftime(t2,t1) << endl;
		/*cout << MCHF_1s4d.get_impact(j) << "     " << p << "     " <<  Ic << "     " << ie_p0 <<
		 "     " << ie_p1 << "     " << ie_p2 << "     " << p0 << "     " << p1 << "     " << p2 <<
		 "\n\n"; */
	}

	writefile.close();

	time(&end);

	cout << "\ntotal time (s): " << difftime(end,begin) << endl;
	cout << "\n" << endl;

	return 0;
}
