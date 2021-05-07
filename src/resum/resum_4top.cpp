//#include "4top_softanom.h" //not needed for now
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameters.h"
#include "resum_4top.h"
#include "qq_process.h"
#include "gg_process.h"

using namespace std;

// LP LL function h0 (or g1) hep-ph/0306211 eqn 39 (note that gammaE is not part there of lambda and the factor 2 (as they have 2*h0 = g1!) or 1905.11771 eqn 6
complex<double> g1_M2(double A1,complex<double>lambda){
	return A1/(2.*M_PI*pow(b0,2))*(2.*lambda+(1.-2.*lambda)*log(1.-2.*lambda));
}

/////////////////////////////////////////////////////////////
// LP NLL function h1 (or g2) hep-ph/0306211 (note factor 2 and gammaE) eqn 40 or 1905.11771 eqn 61
// (checked with mathematica)
complex<double> g2_M2(double A1,double A2,complex<double>lambda){
	double INCeuler = 0.;
	if(INCEULER == 0) {INCeuler = 1.;}
	return 1./(2.*M_PI*b0)*(-A2/(M_PI*b0)+A1*log(M2/muR2))*(2.*lambda+log(1.-2.*lambda))
	+ A1*b1/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	- A1/(M_PI*b0)*lambda*log(M2/muF2)
	- INCeuler*2.*M_gammaE*log(1.-2.*lambda)*A1/(2.*M_PI*b0);
}

//expanded version of the cusp
complex<double> delidelj_exp(complex<double> N, double A1){
	//if(!INCEULER)
	return alphas_muR*2.*A1/M_PI*log(N)*(log(N)-ISNLL*log(M2/muF2));//+ISNLL*INCEULER*M_gammaE);
}

complex<double> qq_res_abs(complex<double> N, vector<double*> mom){
	complex<double> lambda = b0*alphas_muR*log(N);
	qqhard.setMomenta(mom);
	// Evaluate matrix element
	qqhard.sigmaKin();
	const double* matrix_elements = qqhard.getMatrixElements();
	// Do the resummation
	complex<double> wide_soft = exp(ISNLL*-1.*CA*log(1.-2.*lambda)/(2.*M_PI*b0));
	complex<double> cusp_factor = exp(2.*(1./alphas_muR*ISLL*g1_M2(A1q,lambda)+ISNLL*g2_M2(A1q,A2q,lambda)));
	complex<double> sumqq = 0;
	for(int i=0; i<qqhard.ncol;i++)
	{
		if(i < 2) sumqq+=matrix_elements[i];
		else sumqq+=matrix_elements[i]*wide_soft;
	}
	complex<double> result = sumqq*cusp_factor;
	//compute pure correction on top of NLO
	if(expansion){
		 	complex<double> wide_soft_exp = ISNLL*CA*alphas_muR*log(N)/(M_PI);
			complex<double> cusp_expanded = 1.+delidelj_exp(N,A1q);
			complex<double> sumqq_exp = 0;
			for(int i=0; i<qqhard.ncol;i++)
			{
				if(i < 2) sumqq_exp+=matrix_elements[i]*cusp_expanded;
				else sumqq_exp+=matrix_elements[i]*(cusp_expanded+wide_soft_exp);
			}
			return result - sumqq_exp;
	}
	return result;
}

complex<double> gg_res_abs(complex<double> N, vector<double*> mom){
	complex<double> lambda = b0*alphas_muR*log(N);
	
	// Evaluate matrix element
	gghard.setMomenta(mom);
	gghard.sigmaKin();
	const double* matrix_elements = gghard.getMatrixElements();
	
	// Do the resummation
	complex<double> res_factor = -1.*ISNLL*log(1.-2.*lambda)/(2.*M_PI*b0);
	// written in the diagbasis, not in the orthogonal basis
	vector<complex<double>> wide_soft = {exp(8.*res_factor),
											exp(6.*res_factor),
											exp(6.*res_factor),
											exp(4.*res_factor),
											exp(3.*res_factor),
											exp(3.*res_factor),
											exp(3.*res_factor),
											exp(3.*res_factor),
											exp(3.*res_factor),
											exp(3.*res_factor),
											exp(3.*res_factor),
											exp(3.*res_factor),
											1.,
											1.};
	complex<double> cusp_factor = exp(2.*(1./alphas_muR*ISLL*g1_M2(A1g,lambda)+ISNLL*g2_M2(A1g,A2g,lambda)));
	complex<double> sumgg = 0;
	for(int i=0; i<gghard.ncol;i++)
	{
		sumgg+=matrix_elements[i]*wide_soft[i];
	}
	complex<double> result = sumgg*cusp_factor;
	//compute pure correction on top of NLO
	if(expansion){
		 	complex<double> res_factor_exp = ISNLL*alphas_muR*log(N)/(M_PI);
		 	vector<complex<double>> wide_soft_exp = {8.*res_factor_exp,
											6.*res_factor_exp,
											6.*res_factor_exp,
											4.*res_factor_exp,
											3.*res_factor_exp,
											3.*res_factor_exp,
											3.*res_factor_exp,
											3.*res_factor_exp,
											3.*res_factor_exp,
											3.*res_factor_exp,
											3.*res_factor_exp,
											3.*res_factor_exp,
											0.,
											0.};
			complex<double> cusp_expanded = 1.+delidelj_exp(N,A1g);
			complex<double> sumgg_exp = 0;
			for(int i=0; i<gghard.ncol;i++)
			{
				sumgg_exp+=matrix_elements[i]*(cusp_expanded+wide_soft_exp[i]);
			}
			return result - sumgg_exp;
	}
	return result;
}

