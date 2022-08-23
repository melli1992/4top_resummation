#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameters.h"
#include "resum_4top.h"
#include "qq_process.h"
#include "gg_process.h"
#include "4top_softanom.h"

using namespace std;

// make sure this is not reinitialised every time
Eigen::ComplexEigenSolver<Eigen::Matrix<complex<double>,6,6>> ces;
    

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
	M2 = (4.*mt)*(4.*mt);
	
	return 1./(2.*M_PI*b0)*(-A2/(M_PI*b0)+A1*log(M2/muR2))*(2.*lambda+log(1.-2.*lambda))
	+ A1*b1/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	- A1/(M_PI*b0)*lambda*log(M2/muF2)
	- INCeuler*2.*M_gammaE*log(1.-2.*lambda)*A1/(2.*M_PI*b0);
}

//expanded version of the cusp
complex<double> delidelj_exp(complex<double> N, double A1){
	double INCeuler = 0.;
	if(INCEULER == 0) {INCeuler = 1.;}
	return alphas_muR*2.*A1/M_PI*log(N)*(log(N)+2.*INCeuler*M_gammaE-ISNLL*log(M2/muF2));
}

// C1 coefficent for collinear radiation
// needs to be multiplied with H0
double C1_coefficient(double Ci, double gammai){
	double INCeuler = 0.;
	if(INCEULER == 0) {INCeuler = 1.;}
	return alphas_muR*2*Ci/M_PI*(1./2.*zeta2 + 1./4.*pow(log(muR2/Q2),2) + 
					(1./2.*gammai/Ci-INCeuler*M_gammaE)*log(muR2/muF2) +
					INCeuler*M_gammaE*log(muR2/Q2) + pow(INCeuler*M_gammaE, 2)
					);
}

// scale coefficient
// needs to be multiplied with H0
double scale_coefficient(double Ci, double gammai){
	// double INCeuler = 0.;
	// if(INCEULER == 0) {INCeuler = 1.;}
	return -4.*b0*alphas_muR*log(muF2/muR2); // - alphas_muR*2*Ci/M_PI*((1./2.*gammai/Ci-INCeuler*M_gammaE)*log(muR2/muF2));
}

complex<double> qq_res_abs(complex<double> N, vector<double*> mom){

	double INCeuler = 0.;
	if(INCEULER == 0) {INCeuler = 1.;}

	complex<double> lambda = b0*alphas_muR*log(N);
	qqhard.setMomenta(mom);
	// Evaluate matrix element
	qqhard.sigmaKin();
	const double* matrix_elements = qqhard.getMatrixElements();
	// Do the resummation
	// cusp anomalous dimension is (0,0, -CA/2, -CA/2, -CA/2, -CA/2)
	complex<double> wide_soft = exp(ISNLL*-1.*CA*log(1.-2.*lambda)/(2.*M_PI*b0));
	complex<double> cusp_factor = exp(2.*(1./alphas_muR*ISLL*g1_M2(A1q,lambda)+ISNLL*g2_M2(A1q,A2q,lambda)));
	complex<double> sumqq = 0;
	double C1 = 0, SF = 0;
	int size_qqhard = qqhard.ncol;
	vector<double> S1(size_qqhard, 0.0);
	if(include_C1){
		C1 = ISNLL*C1_coefficient(CF, gamma_q)*ONLY_SF;
		SF = ISNLL*scale_coefficient(CF, gamma_q);
	}
	// handle the colour part
	for(int i=0; i<size_qqhard;i++)
	{
		if(include_S1){
			if(i < 2) S1[i] = 0;
			else      S1[i] = ISNLL*alphas_muR/M_PI*CA*(1 + INCeuler*M_gammaE + 1./2.*log(muR2/Q2));
		}
		if(i < 2) sumqq+=matrix_elements[i]*((1.+C1)*ONLY_SF+SF);
		else sumqq+=matrix_elements[i]*((1+S1[i]+C1)*ONLY_SF+SF)*wide_soft;
	}
	complex<double> result = sumqq*cusp_factor;
	//compute pure correction on top of NLO
	if(expansion){
		 	complex<double> wide_soft_exp = ISNLL*CA*alphas_muR*log(N)/(M_PI);
			complex<double> cusp_expanded = 1.*cusp_piece_LO+delidelj_exp(N,A1q)*cusp_piece_NLO;
			complex<double> sumqq_exp = 0;
			for(int i=0; i<size_qqhard;i++)
			{
				//if(i < 2) sumqq_exp+=matrix_elements[i]*(((1.+C1)*ONLY_SF+SF)*cusp_factor-cusp_expanded*ONLY_SF-C1*ONLY_SF-SF);
				//else      sumqq_exp+=matrix_elements[i]*(((1.+S1[i]+C1)*ONLY_SF+SF)*cusp_factor*wide_soft-cusp_expanded*ONLY_SF-wide_soft_exp*ONLY_SF-C1*ONLY_SF-S1[i]*ONLY_SF-SF);
				if(i < 2) sumqq_exp+=matrix_elements[i]*(cusp_expanded*ONLY_SF+C1*ONLY_SF*c1_piece+SF*b0_piece);
				else      sumqq_exp+=matrix_elements[i]*(cusp_expanded*ONLY_SF+wide_soft_exp*ONLY_SF*wide_soft_piece+C1*ONLY_SF*c1_piece+S1[i]*ONLY_SF*s1_piece+SF*b0_piece);
			}
			if (NLL_truncated) return sumqq_exp;
			return result - sumqq_exp;
	}
	return result;
}


complex<double> qq_res_full_sad(complex<double> N, vector<double*> mom){

	double INCeuler = 0.;
	// if(INCEULER == 0) {INCeuler = 1.;}

	complex<double> lambda = b0*alphas_muR*log(N);
	qqhard.setMomenta(mom);
	// Evaluate matrix element
	qqhard.sigmaKin();
	const double* matrix_elements = qqhard.getMatrixElements();
	// Do the resummation
	complex<double> wide_soft   = ISNLL*log(1.-2.*lambda)/(2.*M_PI*b0);
	complex<double> cusp_factor = exp(2.*(1./alphas_muR*ISLL*g1_M2(A1q,lambda)+ISNLL*g2_M2(A1q,A2q,lambda)));
	complex<double> sumqq = 0;
	double C1 = 0, SF = 0;
	int size_qqhard = qqhard.ncol;
	if(include_C1){
		C1 = ISNLL*C1_coefficient(CF, gamma_q)*ONLY_SF;
		SF = ISNLL*scale_coefficient(CF, gamma_q);
	}
	kinematic_invariants KI(mom);
	Eigen::Matrix<complex<double>, 6, 6> sad     = get_Gamma_R_eigen(KI);
	ces.compute(sad);
	// we have
	// GammaR = Rinv.Gamma.R
	// ces.eigenvalues() = ces.eigenvectors().inverse() * sad * ces.eigenvectors()
	// HR = Rinv* H * (Rinv)^dagger
	// (dagger = conjugate transpose = adjoint)
	// SR = R^dagger * S * R
	// need Tr[HR*SR]
	Eigen::DiagonalMatrix<complex<double>, 6> Hard(matrix_elements[0], matrix_elements[1], matrix_elements[2],
												   matrix_elements[3], matrix_elements[4], matrix_elements[5] );
	Eigen::DiagonalMatrix<complex<double>, 6> Stilde_LO(1,1,1,1,1,1);
	Eigen::Matrix<complex<double>, 6, 6> Stilde_NLO = get_S1tilde(KI);

	Eigen::Matrix<complex<double>, 6, 6> HR     = (ces.eigenvectors().inverse()) * Hard.diagonal().asDiagonal()      * (ces.eigenvectors().inverse()).adjoint();
	Eigen::Matrix<complex<double>, 6, 6> SR_LO  = (ces.eigenvectors().adjoint()) * Stilde_LO.diagonal().asDiagonal() * (ces.eigenvectors());
	Eigen::Matrix<complex<double>, 6, 6> SR_NLO = (ces.eigenvectors().adjoint()) * Stilde_NLO                        * (ces.eigenvectors());
	
	// checks
	Eigen::Matrix<complex<double>, 6, 6> HR_SR = HR*(SR_LO + SR_NLO);
	// cout << "========" << endl;
    // cout << "The eigenvalues of sad are:" << endl << ces.eigenvalues() << endl;
    // cout << "The eigenvalues of sad are:" << endl << ces.eigenvalues()[0] << endl;
    // cout << "sad = V * D * V^(-1) = " << endl
    //     << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << endl;
	// cout << "H " << Hard.diagonal() << endl;
	// cout << "S " << Stilde.diagonal() << endl;
	// cout << "HR " << HR << endl;
	// cout << "SR " << SR << endl;
	// cout << "HR_SR " << HR_SR << endl;
	// complex<double> sum_HR_SR;
	// complex<double> sum_H_S;
	// for(int i = 0; i < size_qqhard; i++){
	// 	//cout << "HR_SR " << HR_SR.diagonal()[i] << endl;
	// 	sum_HR_SR +=  HR_SR.diagonal()[i];
	// 	sum_H_S   +=  Hard.diagonal()[i];
	// }
	// // indeed the same up to complex part!
	// std::cout << "sum_HR_SR " << sum_HR_SR << " sum_H_S " << sum_H_S << endl;
	
	// handle the colour part
	for(int i=0; i<size_qqhard;i++){
	  sumqq+=HR_SR(i,i)*((1+C1)*ONLY_SF+SF)*exp(wide_soft*(2.*real(ces.eigenvalues()[i])));
	}
	complex<double> result = sumqq*cusp_factor;
	//compute pure correction on top of NLO
	if(expansion){
			// redefine HR_SR 
			Eigen::Matrix<complex<double>, 6, 6> HR_SR_LO  = HR*SR_LO;
			Eigen::Matrix<complex<double>, 6, 6> HR_SR_NLO = HR_SR - HR_SR_LO;
		 	complex<double> wide_soft_exp = -ISNLL*alphas_muR*log(N)/(M_PI);
			complex<double> cusp_expanded = 1.*cusp_piece_LO+delidelj_exp(N,A1q)*cusp_piece_NLO;
			complex<double> sumqq_exp = 0;
			for(int i=0; i<size_qqhard;i++)
			{
				sumqq_exp+=(HR_SR_LO(i,i)    * cusp_expanded*ONLY_SF
				            + HR_SR_LO(i,i)  * wide_soft_exp*(2.*real(ces.eigenvalues()[i]))*ONLY_SF*wide_soft_piece
							+ HR_SR_LO(i,i)  * C1*ONLY_SF*c1_piece
							+ HR_SR_NLO(i,i) * ONLY_SF*s1_piece
							+ HR_SR_LO(i,i)  * SF*b0_piece);
			}
			if (NLL_truncated) return sumqq_exp;
			return result - sumqq_exp;
	}
	return result;
}


complex<double> gg_res_abs(complex<double> N, vector<double*> mom){

	double INCeuler = 0.;
	if(INCEULER == false) {INCeuler = 1.;}
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
	//cout << "cusp_factor " << cusp_factor << " wo exponent " << log(cusp_factor) << endl;
	complex<double> sumgg = 0;

	// hard collinear correction
	double C1 = 0, SF = 0;
	int size_gghard = gghard.ncol;
	if(include_C1){
		C1 = ISNLL*C1_coefficient(CA, gamma_g)*ONLY_SF;
		SF = ISNLL*scale_coefficient(CA, gamma_g);
	}

	// soft correction
	vector<double> S1(size_gghard, 0.0);
	if(include_S1){
		double mult_factor_s1 = ISNLL*alphas_muR/(M_PI)*(1 + INCeuler*M_gammaE + 1./2.*log(muR2/Q2));
		S1 = {8.*mult_factor_s1,
				6*mult_factor_s1,
				6*mult_factor_s1,
				4*mult_factor_s1,
				3*mult_factor_s1,
				3*mult_factor_s1,
				3*mult_factor_s1,
				3*mult_factor_s1,
				3*mult_factor_s1,
				3*mult_factor_s1,
				3*mult_factor_s1,
				3*mult_factor_s1,
				0,
				0
		};
	}
	
	// handle the colour part
	//cout << "C1=" << C1 << " C1/(alphas/pi) " << C1/(alphas_muR/M_PI) << endl;
	//cout << "INCeuler*M_gammaE=" << INCeuler*M_gammaE << " log(muR2/Q2) " << log(muR2/Q2) << endl;
	for(int i=0; i<size_gghard;i++)
	{	
		//cout << "matrix_elements[" << i << "]=" << matrix_elements[i] << " wide_soft[" << i << "]" << wide_soft[i] << "S1 " << S1[i] << " S1/(alphas_muR/pi) " << S1[i]/(alphas_muR/M_PI) << endl;
		sumgg+=matrix_elements[i]*((1+S1[i]+C1)*ONLY_SF+SF)*wide_soft[i];
	}
	//cout << "Tr[H*Ubar*S*U]= " << sumgg << endl;
	//cout << "cusp_factor=" << cusp_factor << endl;
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
			complex<double> cusp_expanded = 1.*cusp_piece_LO+delidelj_exp(N,A1g)*cusp_piece_NLO;
			complex<double> sumgg_exp = 0;
			for(int i=0; i<gghard.ncol;i++)
			{
				sumgg_exp+=matrix_elements[i]*((cusp_expanded+wide_soft_exp[i]*wide_soft_piece+S1[i]*s1_piece+C1*c1_piece)*ONLY_SF + SF*b0_piece);
			}
			if (NLL_truncated) return sumgg_exp;
			return result - sumgg_exp;
	}
	return result;
}

