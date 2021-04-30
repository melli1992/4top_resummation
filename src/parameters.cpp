#include <bits/stdc++.h>
#include <iostream> //for stirng
#include <fstream>
#include <sstream>
#include <unistd.h> //for getting options
#include <gsl/gsl_math.h>
#include <complex>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdlib.h> // for exit()
#include <getopt.h> //for getting options
#include <vector>
#include "LHAPDF/LHAPDF.h"
#include "fit_coefficients.h"
#include "qq_process.h"

using namespace std;

//////////////////////////////////////////////////////////////////
///
/// parameters.cpp contains the functions to read in the parameters
/// that are used in the code and unchanged
/// contains stuff like invariant mass, com energy, tau, color
/// factors and the electromagnetic coupling constant
///
//////////////////////////////////////////////////////////////////


complex<double> I(0,1);
double CMP(2.5), phiMP(3./4.*M_PI);
string result_map = "/home/mbeekveld/4top/results";

//scales
double M2;
double Q;
double Q2(pow(Q,2.));
double S(13000.);
double S2(pow(S,2.));
double tau(Q2/S2);
double muF(Q);
double muR(Q);
double muF2(muF*muF);
double muR2(muR*muR);
//(6.62607015*10^(-34)*299792458/(2*pi))^2/(1.602176634*10^(−19))^2
// 3.8937937e
double pbunits(0.38937966*pow(10.,9.));
double fbunits(0.38937966*pow(10.,12.));
double mt(173.0);//mt(173.0);
double mt2(pow(mt,2));
double M_gammaE(0.57721566490153286060651209008240243104215933593992);
double zeta2(1.6449340668482264);
double zeta3(1.2020569031595942);

double CF(4./3.);
double CA(3.);
double TF(1./2.);
double alphas_muF(0);
double alphas_muR(0);
double alphas_Q(0);
double alphas_mt(0);
double LambdaQCD(0.208364837218848);
double nF(5); 

double b0((11*CA/3.-2*nF/3.)/(4.*M_PI));
double b1((17.*pow(CA,2)-10.*CA*TF*nF-6.*CF*TF*nF)/(24.*pow(M_PI,2)));

//switches

double ISLL(1);
double ISNLL(1);
bool fitPDF=true, realPDF=false, LO=true;
double INCEULER(1.);

// anomalous dimensions for dQCD
double A1q(CF); // 1405.4827 eq. 12
double A1g(CA);// 1405.4827 eq. 12
double A2q(CF/2.*(CA*(67./18.-pow(M_PI,2)/6.)-10./9.*TF*nF));// 1405.4827 eq. 12
double A2g(CA/2.*(CA*(67./18.-pow(M_PI,2)/6.)-10./9.*TF*nF));// 1405.4827 eq. 12

double B1q(-3./4*CF);
double B1g(-M_PI*b0);

//PDF declarations
string setname("PDF4LHC15_nnlo_100");
std::vector<LHAPDF::PDF*> pdfs;
double xmin_pdfs(0.), xmax_pdfs(0.); //min x, max x and alphas
int use_member(0);
std::unordered_map<double, std::vector<std::vector<double>>> fitcoeff;

//process declarations
qq_process qqhard;
// Read param_card and set parameters
// qqhard.initProc("param_card.dat");

string observable;
///////////////////////////////////////
/// update defaults of the programm
///////////////////////////////////////
void update_defaults(bool printout , bool pdfset){
    Q2 = pow(Q,2);
    muF2 = muF*muF;
    muR2 = muR*muR;
    mt2 = mt*mt;
    tau = Q2/S2;
	if(pdfset){
		LHAPDF::PDFSet setk(setname);
		int nmem(0.); //number of members
		vector<int> pids; //number of flavors, span from -5 to 5 with 0 = 21 gluon
		nmem = setk.size()-1;
		if(use_member>nmem){cout << "PDF member not there, using default value of 0" << endl; use_member = 0;}
		pdfs = setk.mkPDFs();
		pids = pdfs[use_member]->flavors();
		xmin_pdfs = pdfs[use_member]->xMin();
		xmax_pdfs = pdfs[use_member]->xMax();
	}
	
		
	realPDF = false;
	fitPDF = true;
	if(fitPDF){
    cout << "Using fitted PDFs" << endl;
		if (setname=="PDF4LHC15_nnlo_100"){fitcoeff = fitcoeff_PDF4LHC15_nnlo_100;}
		else {cout << "PDFset " << setname << " not implemented, using default" << endl; fitcoeff = fitcoeff_PDF4LHC15_nnlo_100;}
		if (fitcoeff.find(muF) == fitcoeff.end())
		{  cout << "Scale " << muF << " not present in the fitted coefficients!" << endl;
				 cout << "Exiting program!" << endl;
				exit(0);}
	}
	alphas_muF = pdfs[use_member]->alphasQ2(muF*muF);
	alphas_Q = pdfs[use_member]->alphasQ2(Q*Q);
	alphas_muR = pdfs[use_member]->alphasQ2(muR*muR);
    if(printout){
		cout << "=========================================" << endl;
		cout << "PDFset: 				" << setname << endl
			<< "PDFmember: 				" << use_member << endl
			<< "Center of Mass energy [TeV]: 		" << S/1000. << endl
			<< "Momentum [GeV]: 			" << Q << endl
			<< "Renormalization scale [GeV]: 		" << muR << endl
			<< "Factorization scale [GeV]:		" << muF << endl;
	   cout << "alphas_muR: 		 		" << alphas_muR << endl;
		}
        if(fitPDF){cout << "Using the fitted PDFs" << endl; realPDF = false;}
		else{ cout << "Using the real PDFs" << endl; realPDF = true;}
		cout << "INCEULER (resumming gammaE) = " << INCEULER << endl;
		cout << "=========================================" << endl;
	qqhard.initProc("param_card.dat");
}
