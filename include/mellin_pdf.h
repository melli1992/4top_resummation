#include <stdlib.h>
#include <complex>
#include <gsl/gsl_math.h>
#include "LHAPDF/LHAPDF.h"

#ifndef MELLINPDF_H
#define MELLINPDF_H

//
double fit_single_div_sophis(int i, double x);
double fit_double_div_sophis(int i, double x);
double fit_single_div(int i, double x);
double fit_double_div(int i, double x);
double fit_gluon_d0PDF(double x1, double x2);
double fit_quark_d0PDF(double x1, double x2);
double fit_gluon_d1PDF(double x1, double x2);
double fit_gluon_d2PDF(double x1, double x2);
double fit_quark_d1PDF(double x1, double x2);
double fit_quark_d2PDF(double x1, double x2);

//Mellin space
std::complex<double> fit_mellin_pdf_sum_gg(std::complex<double> Nint);
std::complex<double> fit_mellin_pdf_sum_qqbar(std::complex<double> Nint);

//additionals
std::complex<double> xfit_pdfs(int i, std::complex<double> x);
std::complex<double> Dxfit_pdfs(int i, std::complex<double> x);
std::complex<double> DDxfit_pdfs(int i, std::complex<double> x);
std::complex<double> fit_pdfs(int i, std::complex<double> x);
std::complex<double> Dfit_pdfs(int i, std::complex<double> x);
std::complex<double> xfit_Nspace_pdfs(int i, std::complex<double> N);
std::complex<double> Gamma(std::complex<double> x);
#endif
