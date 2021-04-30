#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "LHAPDF/LHAPDF.h"
#include "deriv_pdf.h"
#include "parameters.h"
using namespace std;

////////////////////////////////////////////////////////////////////////
///
/// this file contains all pdf related stuff, like sums of qqbar
/// and derivatives of pdfs
/// note that the xfxQ function of LHAPDF returns the x*f(x) value
///
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/// derivatives for 4top
//////////////////////////////////////////////////////////////

double single_div_sophis(int i, double x)
{
	double epsilon = 1E-5*x;
	if((x-2.*epsilon ) < 0){return 0;}
	else if((x+2.*epsilon) > 1){return 0;}
	else {return (-pdfs[use_member]->xfxQ(i,x+2.*epsilon,muF) + 8*pdfs[use_member]->xfxQ(i,x+epsilon,muF)-8.*pdfs[use_member]->xfxQ(i,x-epsilon,muF)+pdfs[use_member]->xfxQ(i,x-2.*epsilon,muF))/(12.*epsilon);}
}
double double_div_sophis(int i, double x)
{
	double epsilon = 1E-5*x;
	if((x-2.*epsilon ) < 0){return 0;}
	else if((x+2.*epsilon) > 1){return 0;}
	else {return (-single_div_sophis(i, x + 2.*epsilon) + 8*single_div_sophis(i, x + epsilon)-8.*single_div_sophis(i, x - epsilon)+single_div_sophis(i, x - 2.*epsilon))/(12.*epsilon);}
}
double single_div(int i, double x)
{
	double epsilon = 1E-3*x;
	if((x-epsilon/2. ) < 0){return 0;}
	else if((x+epsilon/2.) > 1){return 0;}
	else {return (pdfs[use_member]->xfxQ(i,x+epsilon/2.,muF) - pdfs[use_member]->xfxQ(i,x-epsilon/2.,muF))/(epsilon);}
}
double double_div(int i, double x)
{
	double epsilon = 1E-6*x;
	if((x-epsilon ) < 0){return 0;}
	else if((x+epsilon) > 1){return 0;}
	else {return (pdfs[use_member]->xfxQ(i,x+epsilon,muF) -2.*pdfs[use_member]->xfxQ(i,x,muF)+pdfs[use_member]->xfxQ(i,x-epsilon,muF))/(epsilon*epsilon);}

}
double gluon_d0PDF(double x1, double x2)
{
	//no div
	return (pdfs[use_member]->xfxQ(0,x1,muF))*(pdfs[use_member]->xfxQ(0,x2,muF))/(x1*x2);
}
double gluon_d1PDF(double x1, double x2)
{
	//one div
	return single_div(0, x1)*single_div(0, x2);
}
double gluon_d2PDF(double x1, double x2)
{
	//double div
	return (single_div(0, x1)+x1*double_div(0,x1))*(single_div(0, x2)+x2*double_div(0,x2));

}
double quark_d0PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)
			{
				// we can have qqbar or qbarq, so add those together!
				factor = factor+((pdfs[use_member]->xfxQ(i,x1,muF))*(pdfs[use_member]->xfxQ(-i,x2,muF)) + (pdfs[use_member]->xfxQ(-i,x1,muF))*(pdfs[use_member]->xfxQ(i,x2,muF)))/(x1*x2);
			}
	return factor;

}
double quark_d1PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)
			{
				//one div
				factor = (single_div(i, x1)*single_div(-i, x2)+single_div(i, x2)*single_div(-i, x1)) + factor;
			}
	return factor;

}
double quark_d2PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)
			{
				//double div
				factor = (single_div(i, x1)+x1*double_div(i,x1))*(single_div(-i, x2)+x2*double_div(-i,x2)) +  (single_div(-i, x1)+x1*double_div(-i,x1))*(single_div(i, x2)+x2*double_div(i,x2)) + factor;
			}
	return factor;

}
