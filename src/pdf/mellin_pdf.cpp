#include <bits/stdc++.h>
#include <stdlib.h>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gamma.h>
#include "LHAPDF/LHAPDF.h"
#include "mellin_pdf.h"
#include "parameters.h"
using namespace std;

//////////////////////////////////////////////////////////
/// contains the N space pdfs 
///
//////////////////////////////////////////////////////////


double fit_single_div_sophis(int i, double x)
{
	double epsilon = 1E-5*x;
	if((x-2.*epsilon ) < 0){return 0;}
	else if((x+2.*epsilon) > 1){return 0;}
	else {return (-real(xfit_pdfs(i,x+2.*epsilon)) + 8*real(xfit_pdfs(i,x+epsilon))-8.*real(xfit_pdfs(i,x-epsilon))+real(xfit_pdfs(i,x-2.*epsilon)))/(12.*epsilon);}
}
double fit_double_div_sophis(int i, double x)
{
	double epsilon = 1E-5*x;
	if((x-2.*epsilon ) < 0){return 0;}
	else if((x+2.*epsilon) > 1){return 0;}
	else {return (-fit_single_div_sophis(i, x + 2.*epsilon) + 8*fit_single_div_sophis(i, x + epsilon)-8.*fit_single_div_sophis(i, x - epsilon)+fit_single_div_sophis(i, x - 2.*epsilon))/(12.*epsilon);}
}

double fit_single_div(int i, double x)
{
	double eps = 1E-5*x;
	if((x-eps/2. ) < 0){return 0;}
	else if((x+eps/2.) > 1){return 0;}
	else {return (real(xfit_pdfs(i, x+eps/2.))-real(xfit_pdfs(i, x-eps/2.)))/(eps);}
}
double fit_double_div(int i, double x)
{
	double epsilon = 1E-5*x;
	if((x-epsilon/2. ) < 0){return 0;}
	else if((x+epsilon/2.) > 1){return 0;}
	else {
		double fitplus = (real(xfit_pdfs(i, x+epsilon))-real(xfit_pdfs(i, x)))/(epsilon);
		double fitmin = (real(xfit_pdfs(i, x))-real(xfit_pdfs(i, x-epsilon)))/(epsilon);
		return (fitplus-fitmin)/(epsilon);}

}

double fit_gluon_d0PDF(double x1, double x2)
{
	//no div
	return (real(xfit_pdfs(5,x1))*real(xfit_pdfs(5,x2)))/(x1*x2);
}
double fit_quark_d0PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)//it says so in the .info file!
			{
				// we can have qqbar or qbarq, so add those together!
				//no div
				factor = factor+(real(xfit_pdfs(i+5,x1))*real(xfit_pdfs(5-i,x2)) + real(xfit_pdfs(i+5,x2))*real(xfit_pdfs(5-i,x1)))/(x1*x2);
			}
	return factor;

}
double fit_gluon_d1PDF(double x1, double x2)
{
	//one div
	return fit_single_div(5, x1)*fit_single_div(5, x2);
}
double fit_gluon_d2PDF(double x1, double x2)
{
	//double div
	return (fit_single_div(5, x1)+x1*fit_double_div(5,x1))*(fit_single_div(5, x2)+x2*fit_double_div(5,x2));

}
double fit_quark_d1PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)
			{
				//one div
				factor = (fit_single_div(5+i, x1)*fit_single_div(5-i, x2)+fit_single_div(5+i, x2)*fit_single_div(5-i, x1)) + factor;
			}
	return factor;

}
double fit_quark_d2PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)
			{
				//double div
				factor = (fit_single_div(5+i, x1)+x1*fit_double_div(5+i,x1))*(fit_single_div(5-i, x2)+x2*fit_double_div(5-i,x2)) +  (fit_single_div(5-i, x1)+x1*fit_double_div(5-i,x1))*(fit_single_div(5+i, x2)+x2*fit_double_div(5+i,x2)) + factor;
			}
	return factor;

}

complex<double> fit_mellin_pdf_sum_qqbar(complex<double> Nint){
	complex<double> sum_pdf(0,0);
	for(int i = 1; i <=5; i++){
		sum_pdf+= 2.*(xfit_Nspace_pdfs(5-i,Nint)*xfit_Nspace_pdfs(5+i,Nint));
	}
	return sum_pdf;
}

// gluon gluon convolution
complex<double> fit_mellin_pdf_sum_gg(complex<double> Nint){
	complex<double> sum_pdf(0,0);
	sum_pdf+= xfit_Nspace_pdfs(5,Nint)*xfit_Nspace_pdfs(5,Nint);
	return sum_pdf;
}

/// ---------------------------------------------------------------------------------------------


/////////////////////////////////////////////////////
/// and the functions that read from this structure
/////////////////////////////////////////////////////
// this is xfx(x)
complex<double> xfit_pdfs(int i, complex<double> x){
	complex<double> y = 1.-2.*pow(x,0.5);
	return fitcoeff[muF][i][0]*pow(1.-x,fitcoeff[muF][i][1])*pow(x,fitcoeff[muF][i][2])*(1.+fitcoeff[muF][i][3]*y+fitcoeff[muF][i][4]*(2.*pow(y,2)-1.))+fitcoeff[muF][i][5]*pow(1.-x,fitcoeff[muF][i][6])*pow(x,fitcoeff[muF][i][7])*(1.+ fitcoeff[muF][i][8]*pow(x,fitcoeff[muF][i][9]));
}
//this is the derivative
complex<double> Dxfit_pdfs(int i, complex<double> x){
	double A = fitcoeff[muF][i][0];
	double x3 = fitcoeff[muF][i][1];
	double x4 = fitcoeff[muF][i][2];
	double x5 = fitcoeff[muF][i][3];
	double x6 = fitcoeff[muF][i][4];
	double B = fitcoeff[muF][i][5];
	double x7 = fitcoeff[muF][i][6];
	double x8 = fitcoeff[muF][i][7];
	double C = fitcoeff[muF][i][8];
	double x9 = fitcoeff[muF][i][9];
	return A*pow(1. - x,x3)*pow(x,x4)*(-x5*pow(x,-0.5) + 8.*x6 - (4.*x6)*pow(x,-0.5))
			- A*pow((1. - x),(-1. + x3))*pow(x,x4)*x3*(1. + x5 - 2.*pow(x,0.5)*x5 + x6 - 8.*pow(x,0.5)*x6 + 8.*x*x6)
			+ A*pow((1. - x),x3)*pow(x,(-1. + x4))*x4*(1. + x5 - 2.*pow(x,0.5)*x5 + x6 - 8.*pow(x,0.5)*x6 + 8.*x*x6)
			- B*pow((1. - x),(-1. + x7))*pow(x,x8)*(1. + C*pow(x,x9))*x7
			+ B*pow((1. - x),x7)*pow(x,(-1. + x8))*(1. + C*pow(x,x9))*x8
			+ B*C*pow((1. - x),x7)*pow(x,(-1. + x8 + x9))*x9;
}

//this is the double derivative
complex<double> DDxfit_pdfs(int i, complex<double> x){
	double A = fitcoeff[muF][i][0];
	double x3 = fitcoeff[muF][i][1];
	double x4 = fitcoeff[muF][i][2];
	double x5 = fitcoeff[muF][i][3];
	double x6 = fitcoeff[muF][i][4];
	double B = fitcoeff[muF][i][5];
	double x7 = fitcoeff[muF][i][6];
	double x8 = fitcoeff[muF][i][7];
	double C = fitcoeff[muF][i][8];
	double x9 = fitcoeff[muF][i][9];
	return (A*pow(1. - x,x3)*pow(x,-1.5 + x4)*(x5 + 4.*x6))/2.
				+ 2.*A*pow(1. - x,-1. + x3)*pow(x,-0.5 + x4)*x3*(x5 + (4. - 8.*pow(x,0.5))*x6)
				- 2.*A*pow(1. - x,x3)*pow(x,-1.5 + x4)*x4*(x5 + (4. - 8.*pow(x,0.5))*x6)
				+ A*pow(1. - x,-2. + x3)*pow(x,x4)*(-1. + x3)*x3*(1. + x5 - 2.*pow(x,0.5)*x5 + x6 - 8.*pow(x,0.5)*x6 + 8.*x*x6)
				- 2.*A*pow(1. - x,-1. + x3)*pow(x,-1. + x4)*x3*x4*(1. + x5 - 2.*pow(x,0.5)*x5 + x6 - 8.*pow(x,0.5)*x6 + 8.*x*x6)
				+ A*pow(1. - x,x3)*pow(x,-2. + x4)*(-1. + x4)*x4*(1. + x5 - 2.*pow(x,0.5)*x5 + x6 - 8.*pow(x,0.5)*x6 + 8.*x*x6)
				+ B*pow(1. - x,-2. + x7)*pow(x,x8)*(1. + C*pow(x,x9))*(-1. + x7)*x7
				- 2.*B*pow(1. - x,-1. + x7)*pow(x,-1. + x8)*(1. + C*pow(x,x9))*x7*x8
				+ B*pow(1. - x,x7)*pow(x,-2. + x8)*(1. + C*pow(x,x9))*(-1. + x8)*x8
				- 2.*B*C*pow(1. - x,-1. + x7)*pow(x,-1. + x8 + x9)*x7*x9
				+ B*C*pow(1. - x,x7)*pow(x,-2. + x8 + x9)*x8*x9
				+ B*C*pow(1. - x,x7)*pow(x,-2. + x8 + x9)*x9*(-1. + x8 + x9);
}

// this is fx(x)
complex<double> fit_pdfs(int i, complex<double> x){
	complex<double> y = 1.-2.*pow(x,0.5);
	return 1./x*(fitcoeff[muF][i][0]*pow(1.-x,fitcoeff[muF][i][1])*pow(x,fitcoeff[muF][i][2])*(1.+fitcoeff[muF][i][3]*y+fitcoeff[muF][i][4]*(2.*pow(y,2)-1.))+fitcoeff[muF][i][5]*pow(1.-x,fitcoeff[muF][i][6])*pow(x,fitcoeff[muF][i][7])*(1.+ fitcoeff[muF][i][8]*pow(x,fitcoeff[muF][i][9])));
}
//this is the derivative
complex<double> Dfit_pdfs(int i, complex<double> x){
	double A = fitcoeff[muF][i][0];
	double x3 = fitcoeff[muF][i][1];
	double x4 = fitcoeff[muF][i][2];
	double x5 = fitcoeff[muF][i][3];
	double x6 = fitcoeff[muF][i][4];
	double B = fitcoeff[muF][i][5];
	double x7 = fitcoeff[muF][i][6];
	double x8 = fitcoeff[muF][i][7];
	double C = fitcoeff[muF][i][8];
	double x9 = fitcoeff[muF][i][9];
	return -((1./pow(x,2))*(A*pow((1. - x),(-1. + x3))*pow(x,x4)*(8.*pow(x,2)*(x3 + x4)*x6 - (-1. + x4)*(1. + x5 + x6)
																	+ pow(x,0.5)*(-1. + 2.*x4)*(x5 + 4.*x6)
																	- pow(x,(3./2.))*(-1. + 2.*x3 + 2.*x4)*(x5 + 4.*x6)
																	+ x*(-1. + x4 - x5 + x4*x5 - x6 - 7.*x4*x6 +  x3*(1. + x5 + x6)))
						+ B*pow((1. - x),(-1. + x7))*pow(x,x8)*(1. - x8 + x*(-1. + x7 + x8) - C*pow(x,x9)*(-1. + x8 + x9) + C*pow(x,(1. + x9))*(-1. + x7 + x8 + x9))));
}

// this is xfx(N) (mellin transform) checked that it gives the same answer when it is tranformed back to x-space
// note that if fx(N) is needed, then N-> N-1 works
complex<double> xfit_Nspace_pdfs(int i, complex<double> N){
	double A = fitcoeff[muF][i][0];
	double x3 = fitcoeff[muF][i][1];
	double x4 = fitcoeff[muF][i][2];
	double x5 = fitcoeff[muF][i][3];
	double x6 = fitcoeff[muF][i][4];
	double B = fitcoeff[muF][i][5];
	double x7 = fitcoeff[muF][i][6];
	double x8 = fitcoeff[muF][i][7];
	double C = fitcoeff[muF][i][8];
	double x9 = fitcoeff[muF][i][9];

	// these are optional, it looks like the bended contour now doesnt perform so well anymore
	// although it gets better for higher CMP (as expected)
	// checked that we can take these out, then it still works (it is quicker that way, as it doesn't have to evaluate the if statements)
	// also then we can bend the contour without making an error, as the function still exists there (it can be continued without introducing an error so it seems)
	//if(real(x4+N)<=0){return 0;}
	//if(real(x3)<=-1){return 0;}
	//if(real(x8+N)<=0){return 0;}
	//if(real(x7)<=-1){return 0;}
	return (A*Gamma(x3 + 1.)*(Gamma(N + x4)*Gamma(N + x3 + x4 + 3./2.)*((x5 + 1.)*Gamma(N + x3 + x4 + 2.) + x6*(9.*N + x3 + 9.*x4 + 1.)*Gamma(N + x3 + x4 + 1.)) - 2.*(x5 + 4.*x6)*Gamma(N + x4 + 1./2.)*Gamma(N + x3 + x4 + 1.)*Gamma(N + x3 + x4 + 2.)))/(Gamma(N + x3 + x4 + 1.)*Gamma(N + x3 + x4 + 3./2.)*Gamma(N + x3 + x4 + 2.))
	+ B*Gamma(x7 + 1.)*((C*Gamma(N + x8 + x9))/Gamma(N + x7 + x8 + x9 + 1.) + Gamma(N + x8)/Gamma(N + x7 + x8 + 1.));
}

//needed for the functions above
complex<double> Gamma(complex<double> x){
	gsl_sf_result lnr;
	gsl_sf_result arg;
	complex<double> I(0,1);
	gsl_sf_lngamma_complex_e(real(x),imag(x),&lnr, &arg);
	return exp(lnr.val+I*arg.val);
}

