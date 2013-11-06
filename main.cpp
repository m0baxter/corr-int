
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

	std::string path = "/home/baxter/Documents/Observable/results/test_std_E";
	ss << E << ".txt";
	path += ss.str();

	//convert path to c_string.
	char* pth = new char[path.size()+1];
	strcpy(pth,path.c_str());

	//Open file to be written to:
	ofstream writefile( pth );
	writefile << "b      p       Ic       ie_p0       ie_p1       ie_p2       p0       p1       p2\n";

	cout << "Starting E = " << E << " keV" << endl;

	//Create Wave function:
	WaveFunction MCHF_1s4d ( (char*) "/home/baxter/Documents/Observable/Input/wfn/MCHF_1s-4d(all).txt", E ); //MCHF_1s-4d(all)

	time(&setup);

	cout << "setup done (s): " << difftime(setup,begin) << endl;

	for (int j = 0; j < 1; ++j) {

		time(&t1);

		//Calculate Ic:
		double Ic = MCHF_1s4d.correlationintegral_wb( j );

		//Calculate IEM p:
		double p = MCHF_1s4d.indi_electron( j );

		//Calculate IEM probabilities:
		double ie_p0 = ( 1.0 - p ) * ( 1.0 - p ); 
		double ie_p1 = 2.0 * p * ( 1.0 - p );
		double ie_p2 = p * p;

		//Calculate "exact" probabilities:
		double p0 = Ic/2.0;
		double p1 = 2.0*(1.0-p) - Ic;
		double p2 = 2.0*p - 1.0 + Ic/2.0;

		writefile << MCHF_1s4d.get_impact(j) << "     " << p << "     " <<  Ic << "     " << ie_p0 << "     " << ie_p1 << "     " << ie_p2 << "     " << p0 << "     " << p1 << "     " << p2 << "\n";

		time(&t2);

		cout << "done #" << j + 1 << " (s): " << difftime(t2,t1) << endl;
		cout << MCHF_1s4d.get_impact(j) << "     " << p << "     " <<  Ic << "     " << ie_p0 << "     " << ie_p1 << "     " << ie_p2 << "     " << p0 << "     " << p1 << "     " << p2 << "\n\n";
	}

	writefile.close();

	time(&end);

	cout << "\ntotal time (s): " << difftime(end,begin) << endl;
	cout << "\n" << endl;

	return 0;
}
