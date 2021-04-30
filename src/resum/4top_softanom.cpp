#include "4top_softanom.h"
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameters.h"

using namespace std;
//FEED t wiggle to all the functions!
double beta34(double s34)
{
	return sqrt(1.-4.*mt2/s34);
}//check
complex<double> log_beta34(double s34)
{
	double betatt = beta34(s34);
	double result = (1.+pow(betatt,2))/(2.*betatt)*(log((1.-betatt)/(1.+betatt)));
	if(isnan(result) || isinf(result)){return -1.+I*M_PI;}
	return result +I*M_PI;
}//check

complex<double> Tij(complex<double> tij, complex<double> s)
{
	return log(-tij/(mt*sqrt(s)))+(I*M_PI-1.)/2.;
} //check

complex<double> Omega3(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return (Tij(t13, s)+Tij(t24,s)-Tij(t14,s)-Tij(t23,s))/2.;
}//LET OP DIT MOET t13 = (p1+p3)^2-mt^2 zijn!! Zorg dat dat klopt

complex<double> Lambda3(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return (Tij(t13, s)+Tij(t24,s)+Tij(t14,s)+Tij(t23,s))/2.;
}//LET OP DIT MOET t13 = (p1+p3)^2-mt^2 zijn!! Zorg dat dat klopt

complex<double> lambda_qq_11(double s34)
{
	//cout << log_beta34(s34) << endl;
	return -alphas_muR/(M_PI)*CF*(log_beta34(s34)+1.);//+alphas_muR/(M_PI)*CF*(1.-I*M_PI); //klopt dit?
}

complex<double> lambda_qq_12(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(M_PI)*(CF)/CA*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_qq_21(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(M_PI)*2.*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_qq_22(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, double s34, complex<double> s)
{
	return alphas_muR/(2.*M_PI)*((CA-2.*CF)*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s)+(8.*CF-3.*CA)*Omega3(t13,t24,t14,t23,s));//+alphas_muR/(M_PI)*CF*(1.-I*M_PI);
}




complex<double> lambda_gg_11(double s34)
{
	return -alphas_muR/(M_PI)*CF*(log_beta34(s34)+1.);//+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
}

complex<double> lambda_gg_12()
{
	return 0;
}

complex<double> lambda_gg_13(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(M_PI)*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_gg_21()
{
	return 0;
}

complex<double> lambda_gg_22(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s, double s34)
{
	return alphas_muR/(2.*M_PI)*(1./CA*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s));//+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
}

complex<double> lambda_gg_23(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(2.*M_PI)*CA*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_gg_31(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(M_PI)*2.*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_gg_32(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR*(pow(CA,2)-4.)/(2.*CA*M_PI)*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_gg_33(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s, double s34)
{
	return alphas_muR/(2.*M_PI)*(1./CA*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s));//+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
}
