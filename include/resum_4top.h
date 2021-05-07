#include <complex>
#include <vector>
#ifndef res_4top_H
#define res_4top_H

std::complex<double> g1_M2(double A1,std::complex<double>lambda);
std::complex<double> g2_M2(double A1, double A2, std::complex<double>lambda);
std::complex<double> delidelj_exp(std::complex<double> N, double A1);

std::complex<double> qq_res_abs(std::complex<double> N, vector<double*> mom);
std::complex<double> gg_res_abs(std::complex<double> N, vector<double*> mom);


#endif
