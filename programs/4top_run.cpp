#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include "monte_carlo.h"
#include "mellin_pdf.h"
#include "parameters.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <string.h>
#include <sstream>
#include "inout.h"
#include "4top_vegas.h"
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
	//muR = 235.;
	//mt  = 172.5;
	//mt2 = mt*mt;
	update_defaults();
	vector<double> muF_values_pdf;
	for (std::unordered_map<double, std::vector<std::vector<double>>>::iterator it=fitcoeff.begin(); it!=fitcoeff.end(); ++it)
	muF_values_pdf.push_back((double) it->first);
	sort(muF_values_pdf.begin(),muF_values_pdf.end());

	//////////////////////////////////////////////
	/// creating the output file
	//////////////////////////////////////////////
	//ISLL  = 1, ISNLL = 0, include_S1 = false, include_C1 = false;
	//string homedir = "4top_26052022";
	ofstream output;
	//include_S1 = false; include_C1 = false;
	string homedir = "4top_mass_differences";//"4top_additional_scale_piece_13TeV"; //"4top_17062022";
	//string q_str = "results/"+homedir+"/output_4top_"+method+"_LL_";
	string q_str = "results/"+homedir+"/output_4top_"+method+"_";
	if((include_gg == true) and (include_qqbar==false))  q_str += "gg_";
	if((include_gg == false) and (include_qqbar==true))  q_str += "qq_";
	if(expansion) q_str += "expanded_";
	if(INCEULER) q_str += "nbaryes_";
	else q_str += "nbarno_";
	if(include_S1) q_str += "S1yes_";
	if(include_C1) q_str += "C1yes_";
	if(fitPDF) q_str+= "fitpdf_";
	if(!fitPDF) q_str+= "realpdf_";
	if(ONLY_SF==0) q_str+= "ONLY_SF_";
	q_str += "CMP_"+to_string_round(CMP)+"_phiMP_"+to_string_round(phiMP)+"_muR_"+to_string_round(muR)+"_mt_"+to_string_round(mt);
	output.open(q_str.c_str()); 
	cout << "Opening the file " << q_str << endl;
	//////////////////////////////////////////////
	/// code for total xsec
	//////////////////////////////////////////////
	//cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	//output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	
	//vector<double> top_quark_values = {100,125,150,172.5, 173, 175,200,225,250,275,300,350,400,450,500,550,625,750};
	//for(int k = 0; k < 19; k++){
		mt2 = pow(mt,2);
		Q2 = pow(4.*mt,2); Q = sqrt(Q2);
		tau = Q2/S2;
		//vector<double> scales = {muR};
		vector<double> scales = {mt, 2*mt, 4*mt};
		/*muF = scales[0];
		muR = scales[0];
		if(fitPDF){
			muF = closest(muF_values_pdf, muF);
			muR = closest(muF_values_pdf, muR);
		}*/
		
		update_defaults();
		/*                 
		double s = 8134873.7525206236 ;
		double s12 = 823711.00265290437 ;
		double s34 = 495636.87293792958;
		double thetaCM = 1.9619335642393936;
		double phiCM = 0.0;
		double theta12 = 2.4125718167750776;
		double phi12 = 2.1718290474314728 ;
		double theta34 = 0.69835376399455884 ;
		double phi34 = 1.1490681137878767;
		complex<double> N = 1.7682287539440571+I*0.23177121566281500;
		
		test_res(s, s12, s34, thetaCM, theta12, theta34, phiCM, phi12, phi34, N);*/
		// exit(0);
	//	vector<double> CMPval = {1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
	//	vector<double> phival = {2./3.*M_PI, 0.70833333333*M_PI, 3./4.*M_PI};
	//	for(int k = 0; k < 3; k++){
	//		phiMP = phival[k];
	//		for(int j = 0; j < 11; j++){
	//			if((method == "LO") && (k > 0 || j > 0)) continue;
	//			CMP = CMPval[j];
	
					for(int scalF = 0; scalF < 3; scalF++){
					muF = scales[scalF];
					//muR = scales[scal];
					if(fitPDF){
						muF = closest(muF_values_pdf, muF);
						//muR = closest(muF_values_pdf, muR);
					}
					if(muR == muF){
						continue;
					}
					update_defaults();
					out_result = call_vegas(init_vegas_4top(method),params, true, true);
					cout << "Res: " << out_result.res << " " << out_result.err << endl;
					output << "CMP = " << CMP << " phiMP = " << phiMP << " mt = " << mt << ", muF = " << muF << ", muR = " << muR << " qqbar_included = " << include_qqbar << " gg_included = " << include_gg << " expanded = " << expansion << " includeS1 = " << include_S1 << " includeC1 = " << include_C1 << " Res = " << out_result.res << " " << out_result.err << endl;
					}
				
		
	output.close();
	return 0;
}
	

