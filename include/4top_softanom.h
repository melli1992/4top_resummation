#include <complex>
#include <vector>
#ifndef softanom_4top_H
#define softanom_4top_H

double beta34(double s34);
std::complex<double> log_beta34(double s34);
std::complex<double> Tij(std::complex<double> tij, std::complex<double> s);
std::complex<double> Omega3(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> Lambda3(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);

std::complex<double> lambda_qq_11(double s34);
std::complex<double> lambda_qq_12(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_qq_21(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_qq_22(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, double s34, std::complex<double> s);
std::complex<double> lambda_gg_11(double s34);
std::complex<double> lambda_gg_12();
std::complex<double> lambda_gg_13(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_21();
std::complex<double> lambda_gg_22(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s, double s34);
std::complex<double> lambda_gg_23(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_31(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_32(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_33(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s, double s34);

#endif
