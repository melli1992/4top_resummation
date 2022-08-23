#include "4top_vegas.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "parameters.h"
#include "qq_process.h"
#include "resum_4top.h"
#include <limits>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

double vegas_4top_LO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double s   = k[0]*k[1]*S2;
	M2 = 16.*mt2;
	//cut-off on s
	if(s  < M2) {return 0;}
	//invariants s12 and s34 for the top-quark pairs
	double s34 = 4.*mt2+(pow(sqrt(s) - 2.*mt,2)-4*mt2)*k[2];
	double s12 = 4.*mt2+(pow(sqrt(s) - sqrt(s34),2)-4*mt2)*k[3];
	//jacobian for the transformation
	double Jac = (pow(sqrt(s) - 2.*mt,2)-4*mt2)*(pow(sqrt(s) - sqrt(s34),2)-4*mt2);
	//computation of the partonic cross section
	// partonic[1] = gg channel
	// partonic[0] = qqbar channel
   	vector<double> partonic = xsec_LO(M2/s, s12, s34, k[4], k[5], k[6], k[7], k[8], k[9]) ;
	double result = 0;
	if(!fitPDF) result = Jac*(gluon_d0PDF(k[0],  k[1])*partonic[1] + quark_d0PDF(k[0],  k[1])*partonic[0]);
	if(fitPDF) result = Jac*(fit_gluon_d0PDF(k[0],  k[1])*partonic[1] + fit_quark_d0PDF(k[0],  k[1])*partonic[0]);
	if (isnan(result)){return 0;}
	else{ return result;}
}


double vegas_4top_LO_N(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 ,	 rho, s12   s34, 	angles, 
	//k[0]	 ,	k[1], k[2] k[3],	k[4]-k[8]
	if(k[0] == 1){return 0;}
	M2 = 16.*mt2;
	tau = M2/S2;
	//jacobian enforces cut-off on rho
	double rho = M2/S2+(1.-M2/S2)*k[1];
	double s   = M2/rho;
	//invariants s12 and s34 for the top-quark pairs
	double s34 = 4.*mt2+(pow(sqrt(s) - 2.*mt,2)-4*mt2)*k[2];
	double s12 = 4.*mt2+(pow(sqrt(s) - sqrt(s34),2)-4*mt2)*k[3];
	//jacobian for the transformation
	double Jac = (1.-M2/S2)*(pow(sqrt(s) - 2.*mt,2)-4*mt2)*(pow(sqrt(s) - sqrt(s34),2)-4*mt2);
	//N for inverse Mellin transform
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
	double result = 0;
	vector<double> partonic = xsec_LO(rho, s12, s34, k[4], k[5], k[6], k[7], k[8], k[9]);
	result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
			              (fit_mellin_pdf_sum_gg(N)*partonic[1] 
			                + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
			               *pow(rho,N-1.)*pow(tau, -N));
	if (isnan(result)){return 0;}
	return result;
}


double vegas_4top_LO_1st(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 ,	 rho, s12   s34, 	angles, 
	//k[0]	 ,	k[1], k[2] k[3],	k[4]-k[8]
	if(k[0] == 1){return 0;}
	M2 = 16.*mt2;
	tau = M2/S2;
	//jacobian enforces cut-off on rho
	double rho = M2/S2+(1.-M2/S2)*k[1];
	double s   = M2/rho;
	if(rho == 0){return 0;}
	//to be sure
	if((M2/rho > k[10]*k[11]*S2)){return 0;}
	//if((M2/rho <= k[10]*k[11]*S2)){return 0;} // for the boundary term
	
	//invariants s12 and s34 for the top-quark pairs
	double s34 = 4.*mt2+(pow(sqrt(s) - 2.*mt,2)-4*mt2)*k[2];
	double s12 = 4.*mt2+(pow(sqrt(s) - sqrt(s34),2)-4*mt2)*k[3];
	//jacobian for the transformation
	double Jac = (1.-M2/S2)*(pow(sqrt(s) - 2.*mt,2)-4*mt2)*(pow(sqrt(s) - sqrt(s34),2)-4*mt2);
	//N for inverse Mellin transform
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
	double result = 0;
	vector<double> partonic = xsec_LO(rho, s12, s34, k[4], k[5], k[6], k[7], k[8], k[9]);
	if(!fitPDF){
		result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				    (gluon_d1PDF(k[10],  k[11])*partonic[1] + quark_d1PDF(k[10],  k[11])*partonic[0])
					*pow(rho, N - 1.)*pow(N,-2.)*pow(k[10]*k[11], N)*pow(tau, -N));
	}
	if(fitPDF){ 
		result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				   (fit_gluon_d1PDF(k[10],  k[11])*partonic[1] + fit_quark_d1PDF(k[10],  k[11])*partonic[0])
					*pow(rho, N - 1.)*pow(N,-2.)*pow(k[10]*k[11], N)*pow(tau, -N));
	}
	if (isnan(result)){return 0;}
	return result;
}

double vegas_4top_LO_2nd(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 ,	 rho, s12   s34, 	angles, 
	//k[0]	 ,	k[1], k[2] k[3],	k[4]-k[8]
	if(k[0] == 1){return 0;}
	M2 = 16.*mt2;
	tau = M2/S2;
	//jacobian enforces cut-off on rho
	double rho = M2/S2+(1.-M2/S2)*k[1];
	double s   = M2/rho;
	if(rho == 0){return 0;}
	//to be sure
	if((M2/rho > k[10]*k[11]*S2)){return 0;}
	//if((M2/rho <= k[10]*k[11]*S2)){return 0;} // for the boundary term
	
	//invariants s12 and s34 for the top-quark pairs
	double s34 = 4.*mt2+(pow(sqrt(s) - 2.*mt,2)-4*mt2)*k[2];
	double s12 = 4.*mt2+(pow(sqrt(s) - sqrt(s34),2)-4*mt2)*k[3];
	//jacobian for the transformation
	double Jac = (1.-M2/S2)*(pow(sqrt(s) - 2.*mt,2)-4*mt2)*(pow(sqrt(s) - sqrt(s34),2)-4*mt2);
	//N for inverse Mellin transform
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
	double result = 0;
	vector<double> partonic = xsec_LO(rho, s12, s34, k[4], k[5], k[6], k[7], k[8], k[9]);
	if(!fitPDF){
		result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				    (gluon_d2PDF(k[10],  k[11])*partonic[1] + quark_d2PDF(k[10],  k[11])*partonic[0])
					*pow(rho, N - 1.)*pow(N,-4.)*pow(k[10]*k[11], N)*pow(tau, -N));
	}
	if(fitPDF){ 
		result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				   (fit_gluon_d2PDF(k[10],  k[11])*partonic[1] + fit_quark_d2PDF(k[10],  k[11])*partonic[0])
					*pow(rho, N - 1.)*pow(N,-4.)*pow(k[10]*k[11], N)*pow(tau, -N));
	}
	if (isnan(result)){return 0;}
	return result;
}


double vegas_4top_res(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 ,	 rho, s12   s34, 	angles, 
	//k[0]	 ,	k[1], k[2] k[3],	k[4]-k[8]
	if(k[0] == 1){return 0;}
	M2 = 16.*mt2;
	tau = M2/S2;
	//jacobian enforces cut-off on rho
	double rho = M2/S2+(1.-M2/S2)*k[1];
	double s   = M2/rho;
	//invariants s12 and s34 for the top-quark pairs
	double s34 = 4.*mt2+(pow(sqrt(s) - 2.*mt,2)-4*mt2)*k[2];
	double s12 = 4.*mt2+(pow(sqrt(s) - sqrt(s34),2)-4*mt2)*k[3];
	//jacobian for the transformation
	double Jac = (1.-M2/S2)*(pow(sqrt(s) - 2.*mt,2)-4*mt2)*(pow(sqrt(s) - sqrt(s34),2)-4*mt2);
	//N for inverse Mellin transform
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
	double result = 0;
	vector<complex<double>> partonic = xsec_res(N*exp(INCEULER*M_gammaE), rho, s12, s34, k[4], k[5], k[6], k[7], k[8], k[9]);
	result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
			              (fit_mellin_pdf_sum_gg(N)*partonic[1] 
			                + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
			               *pow(rho,N-1.)*pow(tau, -N));
	if (isnan(result)){return 0;}
	return result;
}


double vegas_4top_res_2nd(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 ,	 rho, s12   s34, 	angles, 
	//k[0]	 ,	k[1], k[2] k[3],	k[4]-k[8]
	if(k[0] == 1){return 0;}
	M2 = 16.*mt2;
	tau = M2/S2;
	//jacobian enforces cut-off on rho
	double rho = M2/S2+(1.-M2/S2)*k[1];
	double s   = M2/rho;
	if(rho == 0){return 0;}
	//to be sure
	if((M2/rho > k[10]*k[11]*S2)){return 0;}
	//if((M2/rho <= k[10]*k[11]*S2)){return 0;} // for the boundary term
	
	//invariants s12 and s34 for the top-quark pairs
	double s34 = 4.*mt2+(pow(sqrt(s) - 2.*mt,2)-4*mt2)*k[2];
	double s12 = 4.*mt2+(pow(sqrt(s) - sqrt(s34),2)-4*mt2)*k[3];
	//jacobian for the transformation
	double Jac = (1.-M2/S2)*(pow(sqrt(s) - 2.*mt,2)-4*mt2)*(pow(sqrt(s) - sqrt(s34),2)-4*mt2);
	//N for inverse Mellin transform
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
	double result = 0;
	vector<complex<double>> partonic = xsec_res(N*exp(INCEULER*M_gammaE), rho, s12, s34, k[4], k[5], k[6], k[7], k[8], k[9]);
	if(!fitPDF){result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				    (gluon_d2PDF(k[10],  k[11])*partonic[1] + quark_d2PDF(k[10],  k[11])*partonic[0])
					*pow(rho, N - 1.)*pow(N,-4.)*pow(k[10]*k[11], N)*pow(tau, -N));
	}
	if(fitPDF){result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				    (fit_gluon_d2PDF(k[10],  k[11])*partonic[1] + fit_quark_d2PDF(k[10],  k[11])*partonic[0])
					*pow(rho, N - 1.)*pow(N,-4.)*pow(k[10]*k[11], N)*pow(tau, -N));
	}
	if (isnan(result)){return 0;}
	return result;
}



double vegas_4top_res_1st(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 ,	 rho, s12   s34, 	angles, 
	//k[0]	 ,	k[1], k[2] k[3],	k[4]-k[8]
	if(k[0] == 1){return 0;}
	M2 = 16.*mt2;
	tau = M2/S2;
	//jacobian enforces cut-off on rho
	double rho = M2/S2+(1.-M2/S2)*k[1];
	double s   = M2/rho;
	if(rho == 0){return 0;}
	//to be sure
	if((M2/rho > k[10]*k[11]*S2)){return 0;}
	//if((M2/rho <= k[10]*k[11]*S2)){return 0;} // for the boundary term
	
	//invariants s12 and s34 for the top-quark pairs
	double s34 = 4.*mt2+(pow(sqrt(s) - 2.*mt,2)-4*mt2)*k[2];
	double s12 = 4.*mt2+(pow(sqrt(s) - sqrt(s34),2)-4*mt2)*k[3];
	//jacobian for the transformation
	double Jac = (1.-M2/S2)*(pow(sqrt(s) - 2.*mt,2)-4*mt2)*(pow(sqrt(s) - sqrt(s34),2)-4*mt2);
	//N for inverse Mellin transform
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
	double result = 0;
	vector<complex<double>> partonic = xsec_res(N*exp(INCEULER*M_gammaE), rho, s12, s34, k[4], k[5], k[6], k[7], k[8], k[9]);
	if(!fitPDF){result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				    (gluon_d1PDF(k[10],  k[11])*partonic[1] + quark_d1PDF(k[10],  k[11])*partonic[0])
					*pow(rho, N - 1.)*pow(N,-2.)*pow(k[10]*k[11], N)*pow(tau, -N));
	}
	if(fitPDF){result = 2.*Jac*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				    (fit_gluon_d1PDF(k[10],  k[11])*partonic[1] + fit_quark_d1PDF(k[10],  k[11])*partonic[0])
					*pow(rho, N - 1.)*pow(N,-2.)*pow(k[10]*k[11], N)*pow(tau, -N));
	}
	if (isnan(result)){return 0;}
	return result;
}
//kallen function - note that the sqrt still has to be taken
double kallen(double x, double y, double z){
return (pow(x,2)+pow(y,2)+pow(z,2)-2*x*y-2*x*z-2*y*z);
}
//kallen function - do not take the sqrt
double kallen12(double x, double y, double z){
return sqrt(pow(x,2)+pow(y,2)+pow(z,2)-2*x*y-2*x*z-2*y*z);
}
// returns momenta in the partonic CM frame
vector<double*> get_momenta(double s, double s12, double s34, double thetaCM, double phiCM, double theta12, double phi12, double theta34, double phi34){
  // these are the initial-state momenta
  vector<double*> p_new(1, new double[4]);
  p_new[0][0] = sqrt(s)/2.;
  p_new[0][1] = 0;
  p_new[0][2] = 0;
  p_new[0][3] = sqrt(s)/2.;
  p_new.push_back(new double[4]);
  p_new[1][0] = sqrt(s)/2.;
  p_new[1][1] = 0;
  p_new[1][2] = 0;
  p_new[1][3] = -sqrt(s)/2.;
  
  //now we have the s12 frame
  double E12 = sqrt(s12)/2.;
  double p12 = sqrt(s12-4.*mt2)/2.;
  //these need to be boosted to the CM frame
  double E12CM = (s+s12-s34)/(2.*sqrt(s));
  double p12CM = sqrt(kallen(s,s12,s34))/(2.*sqrt(s));
  //setting up the boost matrix
  double gamma = E12CM/sqrt(s12);
  double gammaom = (E12CM-sqrt(s12))/sqrt(s12);
  double betax = -p12CM*sin(thetaCM)*cos(phiCM)/E12CM;
  double betay = -p12CM*sin(thetaCM)*sin(phiCM)/E12CM;
  double betaz = -p12CM*cos(thetaCM)/E12CM;
  double beta2 = pow(betax,2)+pow(betay,2)+pow(betaz,2);
  //performing the boost
  double p12x = p12*cos(phi12)*sin(theta12);
  double p12y = p12*sin(phi12)*sin(theta12);
  double p12z = p12*cos(theta12);
  
  double Et1_CM = gamma*E12 - gamma*betax*p12x - gamma*betay*p12y - gamma*betaz*p12z;
  double ptx1_CM = -gamma*betax*E12 + (gammaom*betax*betax/beta2+1.)*p12x + gammaom*betay*betax/beta2*p12y      + gammaom*betax*betaz/beta2*p12z;
  double pty1_CM = -gamma*betay*E12 + gammaom*betax*betay/beta2*p12x      + (gammaom*betay*betay/beta2+1.)*p12y + gammaom*betay*betaz/beta2*p12z;
  double ptz1_CM = -gamma*betaz*E12 + gammaom*betax*betaz/beta2*p12x      + gammaom*betay*betaz/beta2*p12y      + (gammaom*betaz*betaz/beta2+1.)*p12z;
  double Etbar2_CM = gamma*E12 + gamma*betax*p12x + gamma*betay*p12y + gamma*betaz*p12z;
  double ptbarx2_CM = -gamma*betax*E12 - (gammaom*betax*betax/beta2+1.)*p12x - gammaom*betay*betax/beta2*p12y - gammaom*betax*betaz/beta2*p12z;
  double ptbary2_CM = -gamma*betay*E12 - gammaom*betay*betax/beta2*p12x - (gammaom*betay*betay/beta2+1.)*p12y - gammaom*betay*betaz/beta2*p12z;
  double ptbarz2_CM = -gamma*betaz*E12 - gammaom*betax*betaz/beta2*p12x - gammaom*betay*betaz/beta2*p12y - (gammaom*betaz*betaz/beta2+1.)*p12z;
  
  p_new.push_back(new double[4]);
  p_new[2][0] = Et1_CM;
  p_new[2][1] = ptx1_CM;
  p_new[2][2] = pty1_CM;
  p_new[2][3] = ptz1_CM;
  p_new.push_back(new double[4]);
  p_new[3][0] = Etbar2_CM;
  p_new[3][1] = ptbarx2_CM;
  p_new[3][2] = ptbary2_CM;
  p_new[3][3] = ptbarz2_CM;
  
  // and now the same for the s34 frame
  double E34 = sqrt(s34)/2.;
  double p34 = sqrt(s34-4.*mt2)/2.;
  //these need to be boosted to the CM frame
  double E34CM = (s+s34-s12)/(2.*sqrt(s));
  //setting up the boost matrix
  gamma = E34CM/sqrt(s34);
  gammaom = (E34CM-sqrt(s34))/sqrt(s34);
  betax = p12CM*sin(thetaCM)*cos(phiCM)/E34CM;
  betay = p12CM*sin(thetaCM)*sin(phiCM)/E34CM;
  betaz = p12CM*cos(thetaCM)/E34CM;
  beta2 = pow(betax,2)+pow(betay,2)+pow(betaz,2);
  //performing the boost
  double p34x = p34*cos(phi34)*sin(theta34);
  double p34y = p34*sin(phi34)*sin(theta34);
  double p34z = p34*cos(theta34);
  
  double Et3_CM = gamma*E34 - gamma*betax*p34x - gamma*betay*p34y - gamma*betaz*p34z;
  double ptx3_CM = -gamma*betax*E34 + (gammaom*betax*betax/beta2+1.)*p34x + gammaom*betay*betax/beta2*p34y      + gammaom*betax*betaz/beta2*p34z;
  double pty3_CM = -gamma*betay*E34 + gammaom*betax*betay/beta2*p34x      + (gammaom*betay*betay/beta2+1.)*p34y + gammaom*betay*betaz/beta2*p34z;
  double ptz3_CM = -gamma*betaz*E34 + gammaom*betax*betaz/beta2*p34x      + gammaom*betay*betaz/beta2*p34y      + (gammaom*betaz*betaz/beta2+1.)*p34z;
  double Etbar4_CM = gamma*E34 + gamma*betax*p34x + gamma*betay*p34y + gamma*betaz*p34z;
  double ptbarx4_CM = -gamma*betax*E34 - (gammaom*betax*betax/beta2+1.)*p34x - gammaom*betay*betax/beta2*p34y - gammaom*betax*betaz/beta2*p34z;
  double ptbary4_CM = -gamma*betay*E34 - gammaom*betay*betax/beta2*p34x - (gammaom*betay*betay/beta2+1.)*p34y - gammaom*betay*betaz/beta2*p34z;
  double ptbarz4_CM = -gamma*betaz*E34 - gammaom*betax*betaz/beta2*p34x - gammaom*betay*betaz/beta2*p34y - (gammaom*betaz*betaz/beta2+1.)*p34z;
  
  
  p_new.push_back(new double[4]);
  p_new[4][0] = Et3_CM;
  p_new[4][1] = ptx3_CM;
  p_new[4][2] = pty3_CM;
  p_new[4][3] = ptz3_CM;
  p_new.push_back(new double[4]);
  p_new[5][0] = Etbar4_CM;
  p_new[5][1] = ptbarx4_CM;
  p_new[5][2] = ptbary4_CM;
  p_new[5][3] = ptbarz4_CM;
  
  /*for(int i = 0; i < 6; i++){
		cout << i << " " << p_new[i][0] << "," << p_new[i][1] << "," << p_new[i][2] << "," << p_new[i][3] << " " << p_new[i][0]*p_new[i][0]-p_new[i][1]*p_new[i][1]-p_new[i][2]*p_new[i][2]-p_new[i][3]*p_new[i][3] << endl;
	}
  cout << "Check mom conservation " << endl;
  cout << "Energy: " << p_new[0][0]+p_new[1][0]-(p_new[2][0]+p_new[3][0]+p_new[4][0]+p_new[5][0]) << endl;
  cout << "px: " << p_new[0][1]+p_new[1][1]-(p_new[2][1]+p_new[3][1]+p_new[4][1]+p_new[5][1]) << endl;
  cout << "py: " << p_new[0][2]+p_new[1][2]-(p_new[2][2]+p_new[3][2]+p_new[4][2]+p_new[5][2]) << endl;
  cout << "pz: " << p_new[0][3]+p_new[1][3]-(p_new[2][3]+p_new[3][3]+p_new[4][3]+p_new[5][3]) << endl;
  cout << "Check massive particles " << endl;
  cout << "t1+t2: " << p_new[2][0]+p_new[3][0] << "," << p_new[2][1]+p_new[3][1] << "," << p_new[2][2]+p_new[3][2] << "," << p_new[2][3]+p_new[3][3] << endl;
  cout << "p12: " << E12CM << "," << p12CM*sin(thetaCM)*cos(phiCM) << "," << p12CM*sin(thetaCM)*sin(phiCM) << "," << p12CM*cos(thetaCM) << endl; 
  cout << "t3+t4: " << p_new[4][0]+p_new[5][0] << "," << p_new[4][1]+p_new[5][1] << "," << p_new[4][2]+p_new[5][2] << "," << p_new[4][3]+p_new[5][3] << endl;
  cout << "p34: " << E34CM << "," << -p12CM*sin(thetaCM)*cos(phiCM) << "," << -p12CM*sin(thetaCM)*sin(phiCM) << "," << -p12CM*cos(thetaCM) << endl; 
  cout << "p12+p34: "<< E12CM+E34CM << ", p1+p2: " << p_new[0][0]+p_new[1][0] << endl;
  */
  return p_new;
}
	
std::vector<double*> get_momenta_wo_phiCM(double s, double s12, double s34, double thetaCM, double theta12, double phi12, double theta34, double phi34){
  
  double sqrts = sqrt(s);
  double sqrts12 = sqrt(s12);
  double sqrts34 = sqrt(s34);
  
  double p12s = s12*sqrt(s12-4.*mt2);
  double p34s = s34*sqrt(s34-4.*mt2);
  double pCMs = kallen12(s,s12,s34);
  
  double cosCM = cos(thetaCM);
  double cos12 = cos(theta12);
  double cos34 = cos(theta34);
  double sinCM = sin(thetaCM);
  double sin12 = sin(theta12);
  double sin34 = sin(theta34);
  
  double sinphi12 = sin(phi12);
  double cosphi12 = cos(phi12);
  double sinphi34 = sin(phi34);
  double cosphi34 = cos(phi34);

  
  vector<double*> p_new(1, new double[4]);
    
  p_new[0][0] = sqrts/2.;
  p_new[0][1] = 0;
  p_new[0][2] = 0;
  p_new[0][3] = sqrts/2.;
  p_new.push_back(new double[4]);
  p_new[1][0] = sqrt(s)/2.;
  p_new[1][1] = 0;
  p_new[1][2] = 0;
  p_new[1][3] = -sqrts/2.;
  
  
  p_new.push_back(new double[4]);
  p_new[2][0] = (1./(4.*sqrts*s12))*((mt2-mt2+s12)*(s+s12-s34) + p12s*pCMs*cos12);
  p_new[2][1] = p12s*sin12*cosphi12/(2.*sqrts12);
  p_new[2][2] = (1./(4.*sqrts*s12))*(((mt2-mt2+s12)*pCMs + p12s*(s+s12-s34)*cos12)*sinCM + 2.*sqrts*sqrts12*p12s*cosCM*sinphi12*sin12);
  p_new[2][3] = (1./(4.*sqrts*s12))*(((mt2-mt2+s12)*pCMs + p12s*(s+s12-s34)*cos12)*cosCM - 2.*sqrts*sqrts12*p12s*sinCM*sinphi12*sin12);


  p_new.push_back(new double[4]);
  p_new[3][0] = (1./(4.*sqrts*s12))*(-(mt2-mt2-s12)*(s+s12-s34) - p12s*pCMs*cos12);
  p_new[3][1] = -p12s*sin12*cosphi12/(2.*sqrts12);
  p_new[3][2] = (1./(4.*sqrts*s12))*(((-mt2+mt2+s12)*pCMs - p12s*(s+s12-s34)*cos12)*sinCM - 2.*sqrts*sqrts12*p12s*cosCM*sinphi12*sin12);
  p_new[3][3] = (1./(4.*sqrts*s12))*(((-mt2+mt2+s12)*pCMs - p12s*(s+s12-s34)*cos12)*cosCM + 2.*sqrts*sqrts12*p12s*sinCM*sinphi12*sin12);


  p_new.push_back(new double[4]);
  p_new[4][0] = (1./(4.*sqrts*s34))*((mt2-mt2+s34)*(s-s12+s34) + p34s*pCMs*cos34);
  p_new[4][1] = p34s*sin34*cosphi34/(2.*sqrts34);
  p_new[4][2] = (1./(4.*sqrts*s34))*(-((mt2-mt2+s34)*pCMs + p34s*(s-s12+s34)*cos34)*sinCM - 2.*sqrts*sqrts34*p34s*cosCM*sinphi34*sin34);
  p_new[4][3] = (1./(4.*sqrts*s34))*(-((mt2-mt2+s34)*pCMs + p34s*(s-s12+s34)*cos34)*cosCM + 2.*sqrts*sqrts34*p34s*sinCM*sinphi34*sin34);


  p_new.push_back(new double[4]);
  p_new[5][0] = (1./(4.*sqrts*s34))*(-(mt2-mt2-s34)*(s-s12+s34) - p34s*pCMs*cos34);
  p_new[5][1] = -p34s*sin34*cosphi34/(2.*sqrts34);
  p_new[5][2] =(1./(4.*sqrts*s34))*((-(-mt2+mt2+s34)*pCMs + p34s*(s-s12+s34)*cos34)*sinCM + 2.*sqrts*sqrts34*p34s*cosCM*sinphi34*sin34);
  p_new[5][3] = (1./(4.*sqrts*s34))*((-(-mt2+mt2+s34)*pCMs + p34s*(s-s12+s34)*cos34)*cosCM - 2.*sqrts*sqrts34*p34s*sinCM*sinphi34*sin34);
	
	
  return p_new;
	
}	


std::vector<double> xsec_LO(double rho, double s12, double s34, double thetaCM, double theta12, double theta34, double phiCM, double phi12, double phi34){

	double s = M2/rho;
     
	std::vector<double> result = {0,0};
	double s12kal = pow(1.-4.*mt2/s12,0.5);
	double s34kal = pow(1.-4.*mt2/s34,0.5);
	double skal = pow(kallen(s, s12, s34),0.5)/(s);
	double flux = 1./(2.*s);
	double jac_angles = sin(thetaCM)*sin(theta12)*sin(theta34);

	double constants = fbunits/512./pow(2*M_PI,8)*flux*s12kal*s34kal*skal*jac_angles;
	double Mqqbar2 = 0.;
	double Mgg2 = 0.;
	
	vector<double*> mom = get_momenta(s, s12, s34, thetaCM, phiCM, theta12, phi12, theta34, phi34);
	//vector<double*> mom = get_momenta_wo_phiCM(s, s12, s34, thetaCM, theta12, phi12, theta34, phi34);
	
	//for(int i =0; i < mom.size(); i++){
	//	cout << mom[i][0] << " " << mom[i][1] << " " << mom[i][2] << " " << mom[i][3] <<  endl;}
	
	// Evaluate matrix element qqbar
	if(include_qqbar){
		qqhard.setMomenta(mom);
		qqhard.sigmaKin();
		const double* matrix_elements_qq = qqhard.getMatrixElements();
		for(int i=0; i<qqhard.ncol;i++) Mqqbar2+=matrix_elements_qq[i];
	}
	
	// Evaluate matrix element gg
	if(include_gg){
		gghard.setMomenta(mom);
		gghard.sigmaKin();
		const double* matrix_elements_gg = gghard.getMatrixElements();
		for(int i=0; i<gghard.ncol;i++) Mgg2+=matrix_elements_gg[i];
	}
	
	result[0] = constants*Mqqbar2;
	result[1] = constants*Mgg2;
	  
	return result;


}


std::vector<complex<double>> xsec_res(complex<double> N, double rho, double s12, double s34, double thetaCM, double theta12, double theta34, double phiCM, double phi12, double phi34){
    double s = M2/rho;
	
	//cout << "s = " << s << ", s12 = " << s12 << ", s34 = "<<  s34 << ", thetaCM = " << thetaCM << ", theta12 = " << theta12 << ", theta34 = " << theta34 << ", phiCM = " << phiCM << ", phi12 = "<< phi12 << "  
	std::vector<complex<double>> result = {0,0};
	double s12kal = pow(1.-4.*mt2/s12,0.5);
	double s34kal = pow(1.-4.*mt2/s34,0.5);
	double skal = pow(kallen(s, s12, s34),0.5)/(s);
	double flux = 1./(2.*s);
	double jac_angles = sin(thetaCM)*sin(theta12)*sin(theta34);

	double constants = fbunits/512./pow(2*M_PI,8)*flux*s12kal*s34kal*skal*jac_angles;
	
	vector<double*> mom = get_momenta(s, s12, s34, thetaCM, phiCM, theta12, phi12, theta34, phi34);
	complex<double> Mqqbar2 = 0.;
    complex<double> Mgg2 = 0.;
	if(include_qqbar){
		if(full_sad) Mqqbar2 = qq_res_full_sad(N+1., mom);
		else         Mqqbar2 = qq_res_abs(N+1., mom);
	}
	if(include_gg) Mgg2 = gg_res_abs(N+1.,mom);

	result[0] = constants*Mqqbar2;
	result[1] = constants*Mgg2;
	  
	return result;
}


void test_res(double s, double s12, double s34, double thetaCM, double theta12, double theta34, double phiCM, double phi12, double phi34, complex<double> N){
	  /*vector<double*> p_new(1, new double[4]);
	  p_new[0][0] = 0.6920000E+03 ;
	  p_new[0][1] = 0;
	  p_new[0][2] = 0;
	  p_new[0][3] = 0.6920000E+03;
	  p_new.push_back(new double[4]);
	  p_new[1][0] = 0.6920000E+03;
	  p_new[1][1] = 0;
	  p_new[1][2] = 0;
	  p_new[1][3] = -0.6920000E+03;
	  p_new.push_back(new double[4]);
	  p_new[2][0] = 0.2005889E+03;
	  p_new[2][1] = -0.2533817E+02;
	  p_new[2][2] = 0.4595164E+02;
	  p_new[2][3] = -0.8691000E+02;
	  p_new.push_back(new double[4]);
	  p_new[3][0] = 0.4142768E+03;
	  p_new[3][1] = -0.1190623E+03;
	  p_new[3][2] = -0.3461634E+03;
	  p_new[3][3] = 0.8770050E+02;
	  p_new.push_back(new double[4]);
	  p_new[4][0] = 0.2458474E+03;
	  p_new[4][1] = -0.1213912E+03;
	  p_new[4][2] = -0.1120229E+03;
	  p_new[4][3] = 0.5680662E+02;
	  p_new.push_back(new double[4]);
	  p_new[5][0] = 0.5232868E+03;
	  p_new[5][1] = 0.2657917E+03;
	  p_new[5][2] = 0.4122346E+03;
	  p_new[5][3] = -0.5759711E+02;*/
	N*=exp(INCEULER*M_gammaE); 
	complex<double> Mqqbar2 = 0.;
    complex<double> Mgg2 = 0.;
	
	std::vector<complex<double>> result = {0,0};
	double s12kal = pow(1.-4.*mt2/s12,0.5);
	double s34kal = pow(1.-4.*mt2/s34,0.5);
	double skal = pow(kallen(s, s12, s34),0.5)/(s);
	double flux = 1./(2.*s);
	double jac_angles = sin(thetaCM)*sin(theta12)*sin(theta34);

	double constants = fbunits/512./pow(2*M_PI,8)*flux*s12kal*s34kal*skal*jac_angles;
	
	vector<double*> mom = get_momenta(s, s12, s34, thetaCM, phiCM, theta12, phi12, theta34, phi34);
	for(int i =0; i < mom.size(); i++){
		cout << mom[i][0] << " " << mom[i][1] << " " << mom[i][2] << " " << mom[i][3] <<  endl;}
	
	if(include_qqbar) Mqqbar2 = qq_res_abs(N+1., mom);
	if(include_gg) Mgg2 = gg_res_abs(N+1.,mom);
	cout << "cusp*soft*hard= " << Mgg2 << endl;
	
	
	cout << "including constants " << constants*Mgg2 << endl;
	cout << "only constants " << constants << endl;

	
}



std::vector<double> check_res(double s, double s12, double s34, double thetaCM, double theta12, double theta34, double phiCM, double phi12, double phi34){

	cout << s << "," << s12 <<"," << s34 <<"," << thetaCM <<"," << phiCM <<"," << theta12 <<"," << phi12 << "," << theta34 << "," << phi34 << endl;

	std::vector<double> result = {0,0};
	double s12kal = pow(1.-4.*mt2/s12,0.5);
	double s34kal = pow(1.-4.*mt2/s34,0.5);
	double skal = pow(kallen(s, s12, s34),0.5)/(s);
	double flux = 1./(2.*s);
	double jac_angles = sin(thetaCM)*sin(theta12)*sin(theta34);

	double constants = fbunits/512./pow(2*M_PI,8)*flux*s12kal*s34kal*skal*jac_angles;
	double Mqqbar2 = 0.;
	double Mgg2 = 0.;

	vector<double*> mom = get_momenta(s, s12, s34, thetaCM, phiCM, theta12, phi12, theta34, phi34);
	for(int i =0; i < mom.size(); i++){
		cout << mom[i][0] << " " << mom[i][1] << " " << mom[i][2] << " " << mom[i][3] <<  endl;}
	
	// Evaluate matrix element qqbar
	if(include_qqbar){
		qqhard.setMomenta(mom);
		qqhard.sigmaKin();
		const double* matrix_elements_qq = qqhard.getMatrixElements();
		for(int i=0; i<qqhard.ncol;i++) Mqqbar2+=matrix_elements_qq[i];
	}
	
	// Evaluate matrix element gg
	if(include_gg){
		gghard.setMomenta(mom);
		gghard.sigmaKin();
		const double* matrix_elements_gg = gghard.getMatrixElements();
		for(int i=0; i<gghard.ncol;i++) Mgg2+=matrix_elements_gg[i];
	}
	
	cout << Mqqbar2 << " " << Mgg2 << endl;

	result[0] = constants*Mqqbar2;
	result[1] = constants*Mgg2;
	  
	return result;


}
