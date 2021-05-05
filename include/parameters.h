#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdlib.h> // for exit()
#include <getopt.h> //for getting options
#include <vector>
#include "LHAPDF/LHAPDF.h"
#include "qq_process.h"
#include "gg_process.h"


#ifndef PARAM_H 
#define PARAM_H

extern std::string result_map;

extern std::complex<double> I;
extern double CMP, phiMP;

extern double tau;
extern double S;
extern double S2;
extern double Q;
extern double M2;
extern double Q2;
extern double muF;
extern double muF2;
extern double muR;
extern double muR2;

extern double pbunits;
extern double fbunits;

extern double mt;
extern double mt2;

extern double M_gammaE;
extern double zeta2;
extern double zeta3;

extern double CF;
extern double CA;
extern double TF;
extern double alphas_muF;
extern double alphas_muR;
extern double alphas_Q;
extern double LambdaQCD;
extern double nF;

extern double b0;
extern double b1;

//pQCD parameters

extern double A1q;
extern double A1g;
extern double A2q;
extern double A2g;

//PDF parameters
extern std::string setname;
extern std::vector<LHAPDF::PDF*> pdfs; //pdf vector
extern double xmin_pdfs, xmax_pdfs; //min x, max x
extern int use_member; //the member that one needs to use
struct lumni_params {double z; double pT; double xT; double epeta; double emeta; int power; int flavor; int coefficient;};

//settings
extern bool realPDF, fitPDF, LO;
extern double INCEULER;
extern double ISLL;
extern double ISNLL;

//process
extern qq_process qqhard;
extern gg_process gghard;

void update_defaults(bool printout = true , bool pdfset = true);

extern std::unordered_map<double, std::vector<std::vector<double>>> fitcoeff;
#endif
