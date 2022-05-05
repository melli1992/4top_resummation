#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#ifndef VEGAS_4top_H
#define VEGAS_4top_H

double vegas_4top_LO(double *k, size_t dim, void *params);
double vegas_4top_LO_N(double *k, size_t dim, void *params);
double vegas_4top_LO_1st(double *k, size_t dim, void *params);
double vegas_4top_LO_2nd(double *k, size_t dim, void *params);
double vegas_4top_res(double *k, size_t dim, void *params);
double vegas_4top_res_1st(double *k, size_t dim, void *params);
double vegas_4top_res_2nd(double *k, size_t dim, void *params);
void test_res(double s,  double s12, double s34, double thetaCM, double phiCM, double theta12, double phi12, double theta34, double phi34, std::complex<double> N);

double kallen(double x, double y, double z);
double kallen12(double x, double y, double z);
std::vector<double> xsec_LO(double rho,  double s12, double s34, double thetaCM, double phiCM, double theta12, double phi12, double theta34, double phi34);
std::vector<std::complex<double>> xsec_res(std::complex<double> N, double rho,  double s12, double s34, double thetaCM, double phiCM, double theta12, double phi12, double theta34, double phi34);
std::vector<double> check_res(double s,  double s12, double s34, double thetaCM, double phiCM, double theta12, double phi12, double theta34, double phi34);

#endif
