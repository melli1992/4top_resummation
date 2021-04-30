//#include "4top_softanom.h" //not needed for now
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameters.h"
#include "resum_4top.h"
#include "qq_process.h"

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

//expanded version of the wide-angle
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
	/*if(expansion){
		 	complex<double> wide_soft2 = 1.-ISNLL*CA*(-1.)*alphas_muR*log(N)/(M_PI);
			complex<double> exponent = delidelj_exp(N,A1q);
			//cout << "difference in result qq " << result << " " << (H22qq*S22qq*wide_soft2+H22qq*S22qq*delidelj_exp(N,A1q)) << endl;
			return result - (H22qq*S22qq*(wide_soft2+exponent));
	}*/
	return result;
}
/*
complex<double> gg_res_abs(complex<double> N, double s, double t13, double t14, double t23, double t24)
{
	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> wide_soft = ISNLL*CA*(-1.)*log(1.-2.*lambda)/(2.*M_PI*b0);
	complex<double> H22gg = H22_gg_c(s, t13, t14, t23, t24);
	complex<double> H33gg = H33_gg_c(s, t13, t14, t23, t24);
	double S11gg = S11_gg();
	double S22gg = S22_gg();
	double S33gg = S33_gg();

	complex<double> result = (H22gg*(S22gg*exp(wide_soft)+S11gg/(pow(CA,2)))+ exp(wide_soft)*H33gg*S33gg)
														*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1g,lambda)+ISNLL*g2_M2(A1g,A2g,lambda)));
	if(expansion){
		complex<double> wide_soft2 = 1.-ISNLL*CA*(-1.)*alphas_muR*log(N)/(M_PI);
		//cout << "difference in result gg " << result << " " << ((H22gg*(S22gg*wide_soft2+S11gg/(pow(CA,2)))+ wide_soft2*H33gg*S33gg)
		//									+ (H22gg*(S22gg+S11gg/(pow(CA,2)))+ H33gg*S33gg)*delidelj_exp(N,A1g)) << endl;
		complex<double> exponent = delidelj_exp(N,A1g);
		return result - ((H22gg*(S22gg*(exponent+wide_soft2)+S11gg/(pow(CA,2))*(1.+exponent)))+ (exponent+wide_soft2)*H33gg*S33gg);
					}
	return result;
}*/

