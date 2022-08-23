#ifndef softanom_4top_H
#define softanom_4top_H

#include <complex>
#include <vector>
#include <Eigen/Dense>
#include "gsl/gsl_sf_dilog.h"
#include <iostream>

using namespace std;
template<typename T>
T pow2(T x){
    return x*x;
}


double beta34(double s34);
std::complex<double> log_beta34(double s34);
std::complex<double> Tij(std::complex<double> tij, std::complex<double> s);
std::complex<double> Omega3(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> Lambda3(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
void check();

class kinematic_invariants {
    public:

    // kinematic functions
    double DLc13, DLc14, DLc15, DLc16, DLc23, DLc24, DLc25, DLc26;
    double DL34, DL35, DL36, DL45, DL46, DL56;
    double L2beta34, L2beta35, L2beta36, L2beta45, L2beta46, L2beta56;    
    double Lbeta3, Lbeta4, Lbeta5, Lbeta6;
    double T13, T14, T15, T16, T23, T24, T25, T26;
    double L34, L35, L36, L45, L46, L56;
    
    //ctor 
    kinematic_invariants(vector<double*> mom){
        vector<double> p1 = {mom[0][0], mom[0][1], mom[0][2], mom[0][3]};
        vector<double> p2 = {mom[1][0], mom[1][1], mom[1][2], mom[1][3]};
        vector<double> p3 = {mom[2][0], mom[2][1], mom[2][2], mom[2][3]};
        vector<double> p4 = {mom[3][0], mom[3][1], mom[3][2], mom[3][3]};
        vector<double> p5 = {mom[4][0], mom[4][1], mom[4][2], mom[4][3]};
        vector<double> p6 = {mom[5][0], mom[5][1], mom[5][2], mom[5][3]};
        fill_beta_and_angles_XVij(p1, p2, p3, p4, p5, p6);
        fill_Tij();
        fill_Lij();
        fill_DLcij();
        fill_Lbetai();
        fill_L2betaij();
        fill_DLij();
    }

    double dot(vector<double> pi, vector<double> pj){
        return pi[0]*pj[0] - pi[1]*pj[1] - pi[2]*pj[2] - pi[3]*pj[3]; 
    }
    //private: 

    // all kinematic objects needed
    // angles between particles
    double c13, c14, c15, c16, c23, c24, c25, c26, c34, c35, c36, c45, c46, c56;
    // beta's
    double beta3, beta4, beta5, beta6;
    // 1 - betai * betaj * cij
    double v34, v35, v36, v45, v46, v56;
    // = Xij func
    double X34, X35, X36, X45, X46, X56;
    // vijT
    double v34T, v43T, v35T, v53T, v36T, v63T, v45T, v54T, v46T, v64T, v56T, v65T;

    double Xij(double & betai, double & betaj, double & cij){
        return pow2(1 - betai*betaj*cij) - (1 - pow2(betai)) * (1 - pow2(betaj));
    }


    void fill_beta_and_angles_XVij(vector<double> & p1, vector<double> & p2, vector<double> & p3, vector<double> & p4, vector<double> & p5, vector<double> & p6){
        //std::cout << p1[0] << std::endl;
        beta3 = sqrt(1 - dot(p3, p3)/pow2(p3[0]));
        beta4 = sqrt(1 - dot(p4, p4)/pow2(p4[0]));
        beta5 = sqrt(1 - dot(p5, p5)/pow2(p5[0]));
        beta6 = sqrt(1 - dot(p6, p6)/pow2(p6[0]));     
        c13   = 1/beta3 * (1 - dot(p1, p3)/(p1[0]*p3[0]));
        c14   = 1/beta4 * (1 - dot(p1, p4)/(p1[0]*p4[0]));
        c15   = 1/beta5 * (1 - dot(p1, p5)/(p1[0]*p5[0]));
        c16   = 1/beta6 * (1 - dot(p1, p6)/(p1[0]*p6[0])); 
        c23   = 1/beta3 * (1 - dot(p2, p3)/(p2[0]*p3[0]));
        c24   = 1/beta4 * (1 - dot(p2, p4)/(p2[0]*p4[0]));
        c25   = 1/beta5 * (1 - dot(p2, p5)/(p2[0]*p5[0]));
        c26   = 1/beta6 * (1 - dot(p2, p6)/(p2[0]*p6[0])); 
        v34   = dot(p3, p4)/(p3[0]*p4[0]);
        v35   = dot(p3, p5)/(p3[0]*p5[0]);
        v36   = dot(p3, p6)/(p3[0]*p6[0]);
        v45   = dot(p4, p5)/(p4[0]*p5[0]);
        v46   = dot(p4, p6)/(p4[0]*p6[0]);
        v56   = dot(p5, p6)/(p5[0]*p6[0]);
        c34   = 1/(beta3 * beta4) * (1 - v34);
        c35   = 1/(beta3 * beta5) * (1 - v35);
        c36   = 1/(beta3 * beta6) * (1 - v36);
        c45   = 1/(beta4 * beta5) * (1 - v45);
        c46   = 1/(beta4 * beta6) * (1 - v46);
        c56   = 1/(beta5 * beta6) * (1 - v56);
        X34   = pow2(v34) - (1 - pow2(beta3)) * (1 - pow2(beta4));
        X35   = pow2(v35) - (1 - pow2(beta3)) * (1 - pow2(beta5));
        X36   = pow2(v36) - (1 - pow2(beta3)) * (1 - pow2(beta6));
        X45   = pow2(v45) - (1 - pow2(beta4)) * (1 - pow2(beta5));
        X46   = pow2(v46) - (1 - pow2(beta4)) * (1 - pow2(beta6));
        X56   = pow2(v56) - (1 - pow2(beta5)) * (1 - pow2(beta6));
        vector<double> vi_jT = vijT(beta3, beta4, c34, X34);
        v34T  = vi_jT[0];
        v43T  = vi_jT[1];
        vi_jT = vijT(beta3, beta5, c35, X35);
        v35T  = vi_jT[0];
        v53T  = vi_jT[1];
        vi_jT = vijT(beta3, beta6, c36, X36);
        v36T  = vi_jT[0];
        v63T  = vi_jT[1];
        vi_jT = vijT(beta4, beta6, c46, X46);
        v46T  = vi_jT[0];
        v64T  = vi_jT[1];
        vi_jT = vijT(beta4, beta5, c45, X45);
        v45T  = vi_jT[0];
        v54T  = vi_jT[1];
        vi_jT = vijT(beta5, beta6, c56, X56);
        v56T  = vi_jT[0];
        v65T  = vi_jT[1];
        // std::cout << "beta3 = " << beta3 << "\n"
        //           << "beta4 = " << beta4 << "\n"
        //           << "beta5 = " << beta5 << "\n"
        //           << "beta6 = " << beta6 << "\n"
        //           << "c13   = " << c13   << "\n"
        //           << "c14   = " << c14   << "\n"
        //           << "c15   = " << c15   << "\n"
        //           << "c16   = " << c16   << "\n"
        //           << "c23   = " << c23   << "\n"
        //           << "c24   = " << c24   << "\n"
        //           << "c25   = " << c25   << "\n"
        //           << "c26   = " << c26   << "\n"
        //           << "c34   = " << c34   << "\n"
        //           << "c35   = " << c35   << "\n"
        //           << "c36   = " << c36   << "\n"
        //           << "c45   = " << c45   << "\n"
        //           << "c46   = " << c46   << "\n"
        //           << "c56   = " << c56   << "\n"
        //           << "X34   = " << X34   << "\n"
        //           << "X35   = " << X35   << "\n"
        //           << "X36   = " << X36   << "\n"
        //           << "X45   = " << X45   << "\n"
        //           << "X46   = " << X46   << "\n"
        //           << "X56   = " << X56   << "\n"
        //           << "v34   = " << v34   << "\n"
        //           << "v35   = " << v35   << "\n"
        //           << "v36   = " << v36   << "\n"
        //           << "v45   = " << v45   << "\n"
        //           << "v46   = " << v46   << "\n"
        //           << "v56   = " << v56   << "\n" << std::endl;
    }
    vector<double> vijT(double & betai, double & betaj, double & cij, double & Xij){
        double sqrt_Xij = sqrt(Xij); 
        double viT = ((1 - pow2(betai))*(1 - pow2(betaj) + sqrt_Xij) - (1 - betai * betaj * cij) * (1 - betai * betaj * cij + sqrt_Xij)) / ( - pow2(betai) - pow2(betaj) + 2 * betai * betaj * cij);
        double vjT = ((1 - pow2(betaj))*(1 - pow2(betai) - sqrt_Xij) - (1 - betai * betaj * cij) * (1 - betai * betaj * cij - sqrt_Xij)) / ( - pow2(betai) - pow2(betaj) + 2 * betai * betaj * cij);
        return (vector<double>){viT, vjT};
    }
    // kinematic functions
    void fill_Lij(){
        L34 = Lij_real(beta3, beta4, X34, v34, v34T, v43T);
        L35 = Lij_real(beta3, beta5, X35, v35, v35T, v53T);
        L36 = Lij_real(beta3, beta6, X36, v36, v36T, v63T);
        L45 = Lij_real(beta4, beta5, X45, v45, v45T, v54T);
        L46 = Lij_real(beta4, beta6, X46, v46, v46T, v64T);
        L56 = Lij_real(beta5, beta6, X56, v56, v56T, v65T);
        // std::cout << "L34   = " << L34 << "\n"
        //           << "L35   = " << L35 << "\n"
        //           << "L36   = " << L36 << "\n"
        //           << "L45   = " << L45 << "\n"
        //           << "L46   = " << L46 << "\n"
        //           << "L56   = " << L56 << "\n" << std::endl;
    }
    double Lij_real(double & betai, double & betaj, double & Xij, double & vij, double & viT, double & vjT){
        return vij/(2.*sqrt(Xij))*(log((1-pow2(betai))/pow2(viT)) - log((1-pow2(betaj))/pow2(vjT)));
    }
    //
    void fill_Tij(){
        T13 = Tij_real(c13, beta3);
        T14 = Tij_real(c14, beta4);
        T15 = Tij_real(c15, beta5);
        T16 = Tij_real(c16, beta6);
        T23 = Tij_real(c23, beta3);
        T24 = Tij_real(c24, beta4);
        T25 = Tij_real(c25, beta5);
        T26 = Tij_real(c26, beta6);
        // std::cout << "T13   = " << T13   << "\n"
        //           << "T14   = " << T14   << "\n"
        //           << "T15   = " << T15   << "\n"
        //           << "T16   = " << T16   << "\n"
        //           << "T23   = " << T23   << "\n"
        //           << "T24   = " << T24   << "\n"
        //           << "T25   = " << T25   << "\n"
        //           << "T26   = " << T26   << "\n" << std::endl;
    }
    double Tij_real(double & cij, double & betaj){
        return -1./2.*log((1 - pow2(betaj))/pow2(1 - betaj * cij)) - 1./2.;
    }
    //
    void fill_DLcij(){
        DLc13 = DLcij(c13, beta3);
        DLc14 = DLcij(c14, beta4);
        DLc15 = DLcij(c15, beta5);
        DLc16 = DLcij(c16, beta6);
        DLc23 = DLcij(c23, beta3);
        DLc24 = DLcij(c24, beta4);
        DLc25 = DLcij(c25, beta5);
        DLc26 = DLcij(c26, beta6);
        // std::cout << "DLc13   = " << DLc13   << "\n"
        //           << "DLc14   = " << DLc14   << "\n"
        //           << "DLc15   = " << DLc15   << "\n"
        //           << "DLc16   = " << DLc16   << "\n"
        //           << "DLc23   = " << DLc23   << "\n"
        //           << "DLc24   = " << DLc24   << "\n"
        //           << "DLc25   = " << DLc25   << "\n"
        //           << "DLc26   = " << DLc26   << "\n" << std::endl;
    }
    double DLcij(double & cij, double & betaj){
        return gsl_sf_dilog(- betaj*(1-cij)/(1-betaj)) + gsl_sf_dilog(betaj*(1+cij)/(1+betaj)); 
    }
    //
    void fill_Lbetai(){
        Lbeta3 = Lbetai(beta3);
        Lbeta4 = Lbetai(beta4);
        Lbeta5 = Lbetai(beta5);
        Lbeta6 = Lbetai(beta6);
    }    
    double Lbetai(double & betaj){
        return 1./betaj * log((1 + betaj) / (1 - betaj));
    }
    //
    void fill_L2betaij(){
        L2beta34 = L2betaij(beta3, beta4, X34, v34, v34T, v43T);
        L2beta35 = L2betaij(beta3, beta5, X35, v35, v35T, v53T);
        L2beta36 = L2betaij(beta3, beta6, X36, v36, v36T, v63T);
        L2beta45 = L2betaij(beta4, beta5, X45, v45, v45T, v54T);
        L2beta46 = L2betaij(beta4, beta6, X46, v46, v46T, v64T);
        L2beta56 = L2betaij(beta5, beta6, X56, v56, v56T, v65T);
    }
    double L2betaij(double & betai, double & betaj, double & Xij, double & vij, double & viT, double & vjT){
        return (vij/sqrt(Xij))*(pow2(log((1-pow2(betai))/pow2(viT))) - pow2(log((1-pow2(betaj))/pow2(vjT))));
    }
    //
    void fill_DLij(){
        DL34 = DLij(beta3, beta4, X34, v34, v34T, v43T);
        DL35 = DLij(beta3, beta5, X35, v35, v35T, v53T);
        DL36 = DLij(beta3, beta6, X36, v36, v36T, v63T);
        DL45 = DLij(beta4, beta5, X45, v45, v45T, v54T);
        DL46 = DLij(beta4, beta6, X46, v46, v46T, v64T);
        DL56 = DLij(beta5, beta6, X56, v56, v56T, v65T);
    }
    double DLij(double & betai, double & betaj, double & Xij, double & vij, double & viT, double & vjT){
        return (vij/sqrt(Xij))*(  (gsl_sf_dilog((1-betai-viT)/(1-betai)) + gsl_sf_dilog((1+betai-viT)/(1+betai)))
                                - (gsl_sf_dilog((1-betaj-vjT)/(1-betaj)) + gsl_sf_dilog((1+betaj-vjT)/(1+betaj))));
    }
};

Eigen::Matrix<complex<double>, 6, 6> get_Gamma_R_eigen(kinematic_invariants & KI);
Eigen::Matrix<complex<double>, 6, 6> get_S1tilde(kinematic_invariants & KI);


#endif
