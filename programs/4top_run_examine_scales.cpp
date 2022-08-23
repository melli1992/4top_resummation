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
	full_sad      = true;
	NLL_truncated = true;
	expansion     = true;
	ONLY_SF = 1;
	
	string homedir = "4top_230822_nondiagonalpieces";//"4top_each_scale_piece_individually";//"4top_additional_scale_piece_13TeV"; //"4top_17062022";
	//string q_str = "results/"+homedir+"/output_4top_"+method+"_LL_";
	string q_str = "results/"+homedir+"/output_4top_"+method+"_";
	if((include_gg == true) and (include_qqbar==false))  q_str += "gg_";
	if((include_gg == false) and (include_qqbar==true))  q_str += "qq_";
	if(expansion && NLL_truncated) q_str += "NLLtruncated_";
	if(INCEULER) q_str += "nbaryes_";
	else q_str += "nbarno_";
	if(include_S1) q_str += "S1yes_";
	if(include_C1) q_str += "C1yes_";
	if(fitPDF) q_str+= "fitpdf_";
	if(!fitPDF) q_str+= "realpdf_";
	if(ONLY_SF==0) q_str+= "ONLY_SF_";
	if     ((cusp_piece_LO == 1) && (cusp_piece_NLO == 0) && (wide_soft_piece== 0) && (s1_piece == 0) && (c1_piece == 0) && (b0_piece == 0)) q_str+= "cLO_";
	else if((cusp_piece_LO == 0) && (cusp_piece_NLO == 1) && (wide_soft_piece== 0) && (s1_piece == 0) && (c1_piece == 0) && (b0_piece == 0)) q_str+= "cNLO_";
	else if((cusp_piece_LO == 0) && (cusp_piece_NLO == 0) && (wide_soft_piece== 1) && (s1_piece == 0) && (c1_piece == 0) && (b0_piece == 0)) q_str+= "ws_";
	else if((cusp_piece_LO == 0) && (cusp_piece_NLO == 0) && (wide_soft_piece== 0) && (s1_piece == 1) && (c1_piece == 0) && (b0_piece == 0)) q_str+= "s1_";
	else if((cusp_piece_LO == 0) && (cusp_piece_NLO == 0) && (wide_soft_piece== 0) && (s1_piece == 0) && (c1_piece == 1) && (b0_piece == 0)) q_str+= "c1_";
	else if((cusp_piece_LO == 0) && (cusp_piece_NLO == 0) && (wide_soft_piece== 0) && (s1_piece == 0) && (c1_piece == 0) && (b0_piece == 1)) q_str+= "b0_";
	else if((cusp_piece_LO == 1) && (cusp_piece_NLO == 1) && (wide_soft_piece== 1) && (s1_piece == 1) && (c1_piece == 1) && (b0_piece == 1)) q_str+= "";
	else{
		cout << "Cannot be! Look at turning on all the pieces" << endl;
		exit(0);
	}
	
	q_str += "CMP_"+to_string_round(CMP)+"_phiMP_"+to_string_round(phiMP)+"_muR_"+to_string_round(muR)+"_mt_"+to_string_round(mt);
	output.open(q_str.c_str()); 
	cout << "Opening the file " << q_str << endl;
	//////////////////////////////////////////////
	/// code for total xsec
	//////////////////////////////////////////////
		mt2 = pow(mt,2);
		Q2 = pow(4.*mt,2); Q = sqrt(Q2);
		tau = Q2/S2;
		vector<double> scales = {mt, 2*mt, 4*mt};
		
		update_defaults();
		
		for(int scalF = 0; scalF < 3; scalF++){
			muF = scales[scalF];
			update_defaults();
			out_result = call_vegas(init_vegas_4top(method),params, true, true);
			cout << "Res: " << out_result.res << " " << out_result.err << endl;
			output << "CMP = " << CMP << " phiMP = " << phiMP << " mt = " << mt << ", muF = " << muF << ", muR = " << muR << " qqbar_included = " << include_qqbar << " gg_included = " << include_gg << " expanded = " << expansion << " includeS1 = " << include_S1 << " includeC1 = " << include_C1 << " Res = " << out_result.res << " " << out_result.err << endl;
		}
				
		
	output.close();
	return 0;
}
	

