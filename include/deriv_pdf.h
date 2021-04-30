#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "LHAPDF/LHAPDF.h"

#ifndef DERIVPDF_H
#define DERIVPDF_H

double single_div_sophis(int i, double x);
double double_div_sophis(int i, double x);
double single_div(int i, double x);
double double_div(int i, double x);
double gluon_d0PDF(double x1, double x2);
double gluon_d1PDF(double x1, double x2);
double gluon_d2PDF(double x1, double x2);
double quark_d0PDF(double x1, double x2);
double quark_d1PDF(double x1, double x2);
double quark_d2PDF(double x1, double x2);

#endif
