#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include "monte_carlo.h"
#include "mellin_pdf.h"
#include "parameters.h"
#include "qq_process.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <string.h>
#include <sstream>
#include "inout.h"
using namespace std;


string to_string_round(double Q){
	ostringstream q_to_str;
	q_to_str << Q;
	return q_to_str.str();
}


string to_string2(string Q){
	ostringstream q_to_str;
	q_to_str << Q;
	return q_to_str.str();
}

double closest(vector<double> const& vec, double value) {
    auto it = lower_bound(vec.begin(), vec.end(), value);
		auto itm = prev(it);

		double itlow = abs(*itm - value);
		double ithigh = abs(*it - value);
		if (itlow < ithigh){ return *itm;}
    return *it;
}


int main(int argc, char* argv[]){
	//////////////////////////////////////////////
	/// predefinition of everything, setting it up
	//////////////////////////////////////////////
	configure(argc,argv, to_string2("input.cfg"), false);
	results out_result;
	double z =0.1;
	lumni_params params = {z, Q, 2*Q/S, 0, 0, 0,0,0};
    
    muF = 235.;
	muR = 235.;
	update_defaults();
	vector<double> muF_values_pdf;
	for (std::unordered_map<double, std::vector<std::vector<double>>>::iterator it=fitcoeff.begin(); it!=fitcoeff.end(); ++it)
	muF_values_pdf.push_back((double) it->first);
	sort(muF_values_pdf.begin(),muF_values_pdf.end());

	//////////////////////////////////////////////
	/// creating the output file
	//////////////////////////////////////////////
	ofstream output;
	string homedir = "4top_29042021";
	string q_str = "results/"+homedir+"/output_4top_res_CMP_"+to_string_round(CMP)+"_phiMP_"+to_string_round(phiMP);
	if(LO){q_str = "results/"+homedir+"/output_4top_LO_CMP_"+to_string_round(CMP)+"_phiMP_"+to_string_round(phiMP);}
	//////////////////////////////////////////////
	/// code for total xsec
	//////////////////////////////////////////////
	cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	Q2 = pow(4.*mt,2); Q = sqrt(Q2);
	muF = Q;
	muR = Q;
	tau = Q2/S2;
	muF = closest(muF_values_pdf, muF);
	muR = closest(muF_values_pdf, muR);
	update_defaults();
	
	cout << " muF = " << muF << ", muR = " << muR  << endl;
	output << " muF = " << muF << ", muR = " << muR << endl;
	if(LO){
		//out_result = call_vegas(init_vegas_4top("resum"),params, true, true);
		out_result = call_vegas(init_vegas_4top("LO"),params, true, true);
		cout << "Res: " << out_result.res << " " << out_result.err << endl;
		output << "Res: " << out_result.res << " " << out_result.err << endl;
	}
		/*	else{
				diagsoft = true;
				higgs1 = call_vegas(init_vegas_ttH("tot_N1"),params, true, true);
				cout << "N1(abs)_diag, " << higgs1.res << " " << higgs1.err << endl;
				output << "N1(abs)_diag, " << higgs1.res << " " << higgs1.err << endl;
				diagsoft = false;
				higgs1 = call_vegas(init_vegas_ttH("tot_N1"),params, true, true);
				cout << "N1(abs), " << higgs1.res << " " << higgs1.err << endl;
				output << "N1(abs), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N2"),params, true, true);
				cout << "N2(stt), " << higgs1.res << " " << higgs1.err << endl;
				output << "N2(stt), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N3"),params, true, true);
				cout << "N3(Q), " << higgs1.res << " " << higgs1.err << endl;
				output << "N3(Q), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N5"),params, true, true);
				cout << "N5(sttpT), " << higgs1.res << " " << higgs1.err << endl;
				output << "N5(sttpT), " << higgs1.res << " " << higgs1.err << endl;


			}*/
	return 0;
}
	

