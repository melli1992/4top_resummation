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
	configure(argc,argv, to_string2("4top.cfg"), false);
	// use a different file with argument -c
	results out_result;
	double z =0.1;
	lumni_params params = {z, Q, 2*Q/S, 0, 0, 0,0,0};
    
    muF = 235.;
	muR = 235.;
	mt  = 172.5;
	mt2 = mt*mt;
	update_defaults();
	vector<double> muF_values_pdf;
	for (std::unordered_map<double, std::vector<std::vector<double>>>::iterator it=fitcoeff.begin(); it!=fitcoeff.end(); ++it)
	muF_values_pdf.push_back((double) it->first);
	sort(muF_values_pdf.begin(),muF_values_pdf.end());

	//////////////////////////////////////////////
	/// creating the output file
	//////////////////////////////////////////////
	ofstream output;
	string homedir = "4top_07052021";
	string q_str = "results/"+homedir+"/output_4top_res_";
	if((include_gg == true) and (include_qqbar==false))  q_str += "gg_";
	if((include_gg == false) and (include_qqbar==true))  q_str += "qq_";
	if(expansion) q_str += "expanded_";
	q_str += "CMP_"+to_string_round(CMP)+"_phiMP_"+to_string_round(phiMP);
	output.open(q_str.c_str()); 
	//////////////////////////////////////////////
	/// code for total xsec
	//////////////////////////////////////////////
	cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	
	vector<double> top_quark_values = {100,125,150,175,200,225,250,275,300,350,400,450,500,550,600,650,700,750};
	for(int k = 0; k < 18; k++){
		mt = top_quark_values[k];
		mt2 = pow(mt,2);
		Q2 = pow(4.*mt,2); Q = sqrt(Q2);
		muF = 4.*mt;
		muR = 4.*mt;
		tau = Q2/S2;
		muF = closest(muF_values_pdf, muF);
		muR = closest(muF_values_pdf, muR);
		update_defaults();
		output << "mt = " << mt << ", muF = " << muF << ", muR = " << muR << endl;
		output << " qqbar_included = " << include_qqbar << " gg_included = " << include_gg << " expanded = " << expansion << endl;
		out_result = call_vegas(init_vegas_4top("resum"),params, true, true);
		cout << "Res: " << out_result.res << " " << out_result.err << endl;
		output << "Res: " << out_result.res << " " << out_result.err << endl;
	}
		
	output.close();
	return 0;
}
	

