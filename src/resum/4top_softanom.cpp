#include "4top_softanom.h"
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameters.h"

using namespace std;

complex<double> i(0,1);
// void set_threshold(){
//     double beta = 1;
//     Lbeta34 = -1. + i*M_PI/(2*beta);
//     Lbeta35 = Lbeta34;
//     Lbeta36 = Lbeta34;
//     Lbeta45 = Lbeta34;
//     Lbeta46 = Lbeta34;
//     Lbeta56 = Lbeta34;
//     T13 = (i*M_PI - 1.)/2.;
//     T14 = T13;
//     T15 = T13;
//     T16 = T13;
//     T23 = T13;
//     T24 = T13;
//     T25 = T13;
//     T26 = T13;
// }

double dot(vector<double> pi, vector<double> pj){
    return pi[0]*pj[0] - pi[1]*pj[1] - pi[2]*pj[2] - pi[3]*pj[3]; 
}
//double mass(vector<double> pi){
//    return sqrt(pi[0]*pi[0] - pi[1]*pi[1] - pi[2]*pi[2] - pi[3]*pi[3]); 
//}

vector<double> betaij_kappaij_fcn(vector<double> pi, vector<double> pj, double mi, double mj){
    double pi_plus_pj_sq = pow(mi,2) + pow(mj, 2) + 2.*dot(pi, pj);

    double betaij  =  sqrt(1 - pow(mi + mj, 2)/pi_plus_pj_sq);
    double kappaij =  sqrt(1 - pow(mi - mj, 2)/pi_plus_pj_sq);
    return (vector<double>){betaij, kappaij};
}

complex<double> Lbetaij(vector<double> pi, vector<double> pj, double mi = mt, double mj = mt){
    vector<double> betaij_kappaij = betaij_kappaij_fcn(pi, pj, mi, mj);
    double betaij  = betaij_kappaij[0];
    double kappaij = betaij_kappaij[1];
    return (pow(kappaij,2) + pow(betaij,2))/(2.*kappaij*betaij)*(log((kappaij-betaij)/(kappaij+betaij)) + i*M_PI);
}

complex<double> Tij(vector<double> pi, vector<double> pj, double sqrt_s, double mj){
    return log(2.*dot(pi,pj)/(sqrt_s*mj)) + (i*M_PI - 1.)/2.;
}

//vector<Eigen::Matrix<complex<double>, 6, 6>, Eigen::Matrix<complex<double>, 6, 6>, Eigen::Matrix<complex<double>, 6, 6>> get_Gamma_R_eigen(vector<double*> mom){
// Eigen::Matrix<complex<double>, 6, 6> get_Gamma_R_eigen(vector<double*> mom){
    
//     vector<double> p1 = {mom[0][0], mom[0][1], mom[0][2], mom[0][3]};
//     vector<double> p2 = {mom[1][0], mom[1][1], mom[1][2], mom[1][3]};
//     vector<double> p3 = {mom[2][0], mom[2][1], mom[2][2], mom[2][3]};
//     vector<double> p4 = {mom[3][0], mom[3][1], mom[3][2], mom[3][3]};
//     vector<double> p5 = {mom[4][0], mom[4][1], mom[4][2], mom[4][3]};
//     vector<double> p6 = {mom[5][0], mom[5][1], mom[5][2], mom[5][3]};
    
//     complex<double> Lbeta34, Lbeta35, Lbeta36, Lbeta45, Lbeta46, Lbeta56;
//     complex<double> KI.T13, KI.T14, KI.T15, KI.T16, KI.T23, KI.T24,KI.T25, KI.T26;

//     Lbeta34 = Lbetaij(p3, p4);
//     Lbeta35 = Lbetaij(p3, p5);
//     Lbeta36 = Lbetaij(p3, p6);
//     Lbeta45 = Lbetaij(p4, p5);
//     Lbeta46 = Lbetaij(p4, p6);
//     Lbeta56 = Lbetaij(p5, p6);
//     KI.T13     = KI.Tij(p1, p3, Q, mt);
//     KI.T14     = KI.Tij(p1, p4, Q, mt);
//     KI.T15     = KI.Tij(p1, p5, Q, mt);
//     KI.T16     = KI.Tij(p1, p6, Q, mt);
//     KI.T23     = KI.Tij(p2, p3, Q, mt);
//     KI.T24     = KI.Tij(p2, p4, Q, mt);
//     KI.T25     = KI.Tij(p2, p5, Q, mt);
//     KI.T26     = KI.Tij(p2, p6, Q, mt);

//     complex<double> G11 ,G12 ,G13 ,G14, G15, G16,
//     G21 ,G22 ,G23 ,G24, G25, G26,
//     G31 ,G32 ,G33 ,G34, G35, G36,
//     G41 ,G42 ,G43 ,G44, G45, G46,
//     G51 ,G52 ,G53 ,G54, G55, G56,
//     G61 ,G62 ,G63 ,G64, G65, G66;
//     //set_threshold();
//     G11 = -CF*(Lbeta34 + Lbeta56 + 2.);
//     G12 = sqrt(pow(CA,2) - 1)/(2*CA)*(Lbeta35+Lbeta46-Lbeta36-Lbeta45);
//     G21 = G12;
//     G13 = sqrt(pow(CA,2) - 1)/(2*CA)*(KI.T15 - KI.T16 - KI.T25 + KI.T26);
//     G31 = G13;
//     G14 = sqrt(pow(CA,2) - 1)/(2*CA)*(KI.T13 - KI.T14 - KI.T23 + KI.T24);
//     G41 = G14;
//     G15 = 0;
//     G51 = 0;
//     G16 = 0;
//     G61 = 0;
//     G22 = -2*CF + 1/(2*CA)*(Lbeta34+Lbeta56) - 1/CA*(Lbeta35+Lbeta46) - (pow(CA,2)-2)/(2*CA)*(Lbeta36 + Lbeta45);
//     G23 = 1/(2*CA)*(KI.T13 - KI.T14 - KI.T23 + KI.T24);
//     G32 = G23;
//     G24 = 1/(2*CA)*(KI.T15 - KI.T16 - KI.T25 + KI.T26);
//     G42 = G24;
//     G25 = sqrt(pow(CA,2) - 4)/(2*sqrt(2)*CA)*(KI.T13 - KI.T14 +KI.T15 - KI.T16 - KI.T23 + KI.T24 - KI.T25 + KI.T26);
//     G52 = G25;
//     G26 = 1/(2*sqrt(2))*(-KI.T13 - KI.T14 +KI.T15 + KI.T16 + KI.T23 + KI.T24 - KI.T25 - KI.T26);
//     G62 = G26;
//     G33 = -(pow(CA,2) - 2)/(2*CA) - CF*Lbeta34 + 1/(2*CA)*Lbeta56 + (pow(CA,2)-2)/(2*CA)*(KI.T15+KI.T26) + 1/CA*(KI.T25+KI.T16);
//     G34 = G12/sqrt(pow(CA,2) - 1);
//     G43 = G34;
//     G35 = sqrt(pow(CA,2)-4)/(2*sqrt(2)*CA)*(Lbeta35-Lbeta36-Lbeta45+Lbeta46+KI.T13 - KI.T14 - KI.T23 + KI.T24);
//     G53 = G35;
//     G36 = 1/(2*sqrt(2))*(-Lbeta35-Lbeta36+Lbeta45+Lbeta46 - KI.T13 + KI.T14 - KI.T23 + KI.T24);
//     G63 = G36;
//     G44 = -CF + 1/(2*CA) - CF*Lbeta56 + 1/(2*CA)*Lbeta34 + (pow(CA,2)-2)/(2*CA)*(KI.T13+KI.T24) + 1/CA*(KI.T14+KI.T23);
//     G45 = sqrt(pow(CA,2)-4)/(2*sqrt(2)*CA)*(Lbeta35-Lbeta36-Lbeta45+Lbeta46+KI.T15 - KI.T16 - KI.T25 + KI.T26);
//     G54 = G45;
//     G46 = 1/(2*sqrt(2))*(Lbeta35-Lbeta36+Lbeta45-Lbeta46 + KI.T15 - KI.T16 + KI.T25 - KI.T26);
//     G64 = G46;
//     G55 = -CF + 1/(2*CA) +1/(2*CA)*(Lbeta34 - 3.*Lbeta35 - 3.*Lbeta46 + Lbeta56 - (pow(CA,2) - 6)/2*(Lbeta36 + Lbeta45))
//               + 3/(2*CA)*(KI.T14+KI.T16+KI.T23+KI.T25) + (pow(CA,2)-6)/(4*CA)*(KI.T13+KI.T15+KI.T24+KI.T26);
//     G56 = sqrt(pow(CA,2)-4)/(4)*(-Lbeta36+Lbeta45-KI.T13 + KI.T15 + KI.T24 - KI.T26);
//     G65 = G56;
//     G66 = -CF + 1/(2*CA) +1/(2*CA)*(Lbeta34 - Lbeta35 - Lbeta46 + Lbeta56 - (pow(CA,2) - 2)/(2)*(Lbeta36 + Lbeta45))
//               + 1/(2*CA)*(KI.T14+KI.T16+KI.T23+KI.T25) + (pow(CA,2)-2)/(4*CA)*(KI.T13+KI.T15+KI.T24+KI.T26);

//     Eigen::Matrix<complex<double>, 6, 6> sad;
//     sad <<  G11 ,G12 ,G13 ,G14, G15, G16,
//             G21 ,G22 ,G23 ,G24, G25, G26,
//             G31 ,G32 ,G33 ,G34, G35, G36,
//             G41 ,G42 ,G43 ,G44, G45, G46,
//             G51 ,G52 ,G53 ,G54, G55, G56,
//             G61 ,G62 ,G63 ,G64, G65, G66;
    
//     return sad;

// }

Eigen::Matrix<complex<double>, 6, 6> get_Gamma_R_eigen(kinematic_invariants & KI){
    
    complex<double> G11 ,G12 ,G13 ,G14, G15, G16,
    G21 ,G22 ,G23 ,G24, G25, G26,
    G31 ,G32 ,G33 ,G34, G35, G36,
    G41 ,G42 ,G43 ,G44, G45, G46,
    G51 ,G52 ,G53 ,G54, G55, G56,
    G61 ,G62 ,G63 ,G64, G65, G66;
    //set_threshold();
    G11 = -CF*(KI.L34 + KI.L56 + 2.);
    G12 = sqrt(pow(CA,2) - 1)/(2*CA)*(KI.L35+KI.L46-KI.L36-KI.L45);
    G21 = G12;
    G13 = sqrt(pow(CA,2) - 1)/(2*CA)*(KI.T15 - KI.T16 - KI.T25 + KI.T26);
    G31 = G13;
    G14 = sqrt(pow(CA,2) - 1)/(2*CA)*(KI.T13 - KI.T14 - KI.T23 + KI.T24);
    G41 = G14;
    G15 = 0;
    G51 = 0;
    G16 = 0;
    G61 = 0;
    G22 = -2*CF + 1/(2*CA)*(KI.L34+KI.L56) - 1/CA*(KI.L35+KI.L46) - (pow(CA,2)-2)/(2*CA)*(KI.L36 + KI.L45);
    G23 = 1/(2*CA)*(KI.T13 - KI.T14 - KI.T23 + KI.T24);
    G32 = G23;
    G24 = 1/(2*CA)*(KI.T15 - KI.T16 - KI.T25 + KI.T26);
    G42 = G24;
    G25 = sqrt(pow(CA,2) - 4)/(2*sqrt(2)*CA)*(KI.T13 - KI.T14 +KI.T15 - KI.T16 - KI.T23 + KI.T24 - KI.T25 + KI.T26);
    G52 = G25;
    G26 = 1/(2*sqrt(2))*(-KI.T13 - KI.T14 +KI.T15 + KI.T16 + KI.T23 + KI.T24 - KI.T25 - KI.T26);
    G62 = G26;
    G33 = -(pow(CA,2) - 2)/(2*CA) - CF*KI.L34 + 1/(2*CA)*KI.L56 + (pow(CA,2)-2)/(2*CA)*(KI.T15+KI.T26) + 1/CA*(KI.T25+KI.T16);
    G34 = G12/sqrt(pow(CA,2) - 1);
    G43 = G34;
    G35 = sqrt(pow(CA,2)-4)/(2*sqrt(2)*CA)*(KI.L35-KI.L36-KI.L45+KI.L46+KI.T13 - KI.T14 - KI.T23 + KI.T24);
    G53 = G35;
    G36 = 1/(2*sqrt(2))*(-KI.L35-KI.L36+KI.L45+KI.L46 - KI.T13 + KI.T14 - KI.T23 + KI.T24);
    G63 = G36;
    G44 = -CF + 1/(2*CA) - CF*KI.L56 + 1/(2*CA)*KI.L34 + (pow(CA,2)-2)/(2*CA)*(KI.T13+KI.T24) + 1/CA*(KI.T14+KI.T23);
    G45 = sqrt(pow(CA,2)-4)/(2*sqrt(2)*CA)*(KI.L35-KI.L36-KI.L45+KI.L46+KI.T15 - KI.T16 - KI.T25 + KI.T26);
    G54 = G45;
    G46 = 1/(2*sqrt(2))*(KI.L35-KI.L36+KI.L45-KI.L46 + KI.T15 - KI.T16 + KI.T25 - KI.T26);
    G64 = G46;
    G55 = -CF + 1/(2*CA) +1/(2*CA)*(KI.L34 - 3.*KI.L35 - 3.*KI.L46 + KI.L56 - (pow(CA,2) - 6)/2*(KI.L36 + KI.L45))
              + 3/(2*CA)*(KI.T14+KI.T16+KI.T23+KI.T25) + (pow(CA,2)-6)/(4*CA)*(KI.T13+KI.T15+KI.T24+KI.T26);
    G56 = sqrt(pow(CA,2)-4)/(4)*(-KI.L36+KI.L45-KI.T13 + KI.T15 + KI.T24 - KI.T26);
    G65 = G56;
    G66 = -CF + 1/(2*CA) +1/(2*CA)*(KI.L34 - KI.L35 - KI.L46 + KI.L56 - (pow(CA,2) - 2)/(2)*(KI.L36 + KI.L45))
              + 1/(2*CA)*(KI.T14+KI.T16+KI.T23+KI.T25) + (pow(CA,2)-2)/(4*CA)*(KI.T13+KI.T15+KI.T24+KI.T26);

    Eigen::Matrix<complex<double>, 6, 6> sad;
    sad <<  G11 ,G12 ,G13 ,G14, G15, G16,
            G21 ,G22 ,G23 ,G24, G25, G26,
            G31 ,G32 ,G33 ,G34, G35, G36,
            G41 ,G42 ,G43 ,G44, G45, G46,
            G51 ,G52 ,G53 ,G54, G55, G56,
            G61 ,G62 ,G63 ,G64, G65, G66;
    
    return sad;

}

Eigen::Matrix<complex<double>, 6, 6> get_S1tilde(kinematic_invariants & KI){
    assert(INCEULER == 1 && "Cannot handle the euler factor now!");
    complex<double> S1tilde11 ,S1tilde12 ,S1tilde13 ,S1tilde14, S1tilde15, S1tilde16,
                    S1tilde21 ,S1tilde22 ,S1tilde23 ,S1tilde24, S1tilde25, S1tilde26,
                    S1tilde31 ,S1tilde32 ,S1tilde33 ,S1tilde34, S1tilde35, S1tilde36,
                    S1tilde41 ,S1tilde42 ,S1tilde43 ,S1tilde44, S1tilde45, S1tilde46,
                    S1tilde51 ,S1tilde52 ,S1tilde53 ,S1tilde54, S1tilde55, S1tilde56,
                    S1tilde61 ,S1tilde62 ,S1tilde63 ,S1tilde64, S1tilde65, S1tilde66;
    //set_threshold();
    S1tilde11 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*((-1 + pow(CA,2))*((4*KI.DL34 + 4*KI.DL56 + KI.L2beta34 + KI.L2beta56 + 2*KI.Lbeta3 + 2*KI.Lbeta4 + 2*KI.Lbeta5 + 2*KI.Lbeta6)/8. + ((2 + KI.L34 + KI.L56)*log(muR2/M2))/2.))/CA;
    S1tilde12 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*(-((4*KI.DL35 - 4*KI.DL36 - 4*KI.DL45 + 4*KI.DL46 + KI.L2beta35 - KI.L2beta36 - KI.L2beta45 + KI.L2beta46)*sqrt(-1 + pow(CA,2)))/8. - ((KI.L35 - KI.L36 - KI.L45 + KI.L46)*sqrt(-1 + pow(CA,2))*log(muR2/M2))/2.)/CA;
    S1tilde13 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*((sqrt(-1 + pow(CA,2))*(KI.DLc15 - KI.DLc16 - KI.DLc25 + KI.DLc26 + KI.T15 + pow(KI.T15,2) - KI.T16 - pow(KI.T16,2) - KI.T25 - pow(KI.T25,2) + KI.T26 + pow(KI.T26,2)))/2. + (sqrt(-1 + pow(CA,2))*(-KI.T15 + KI.T16 + KI.T25 - KI.T26)*log(muR2/M2))/2.)/CA;
    S1tilde14 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*((sqrt(-1 + pow(CA,2))*(KI.DLc13 - KI.DLc14 - KI.DLc23 + KI.DLc24 + KI.T13 + pow(KI.T13,2) - KI.T14 - pow(KI.T14,2) - KI.T23 - pow(KI.T23,2) + KI.T24 + pow(KI.T24,2)))/2. + (sqrt(-1 + pow(CA,2))*(-KI.T13 + KI.T14 + KI.T23 - KI.T24)*log(muR2/M2))/2.)/CA;
    S1tilde15 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*0;
    S1tilde16 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*0;
    S1tilde22 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*CA*((4*KI.DL36 + 4*KI.DL45 + KI.L2beta36 + KI.L2beta45 + 2*KI.Lbeta3 + 2*KI.Lbeta4 + 2*KI.Lbeta5 + 2*KI.Lbeta6)/8. - ((-8 - 4*KI.L36 - 4*KI.L45)*log(muR2/M2))/8.) + ((-4*KI.DL34 + 8*KI.DL35 - 8*KI.DL36 - 8*KI.DL45 + 8*KI.DL46 - 4*KI.DL56 - KI.L2beta34 + 2*KI.L2beta35 - 2*KI.L2beta36 - 2*KI.L2beta45 + 2*KI.L2beta46 - KI.L2beta56 - 2*KI.Lbeta3 - 2*KI.Lbeta4 - 2*KI.Lbeta5 - 2*KI.Lbeta6)/8. - ((8 + 4*KI.L34 - 8*KI.L35 + 8*KI.L36 + 8*KI.L45 - 8*KI.L46 + 4*KI.L56)*log(muR2/M2))/8.)/CA;
    S1tilde23 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*((KI.DLc13 - KI.DLc14 - KI.DLc23 + KI.DLc24 + KI.T13 + pow(KI.T13,2) - KI.T14 - pow(KI.T14,2) - KI.T23 - pow(KI.T23,2) + KI.T24 + pow(KI.T24,2))/2. + ((-KI.T13 + KI.T14 + KI.T23 - KI.T24)*log(muR2/M2))/2.)/CA;
    S1tilde24 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*((KI.DLc15 - KI.DLc16 - KI.DLc25 + KI.DLc26 + KI.T15 + pow(KI.T15,2) - KI.T16 - pow(KI.T16,2) - KI.T25 - pow(KI.T25,2) + KI.T26 + pow(KI.T26,2))/2. + ((-KI.T15 + KI.T16 + KI.T25 - KI.T26)*log(muR2/M2))/2.)/CA;
    S1tilde25 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*((sqrt(-4 + pow(CA,2))*(KI.DLc13 - KI.DLc14 + KI.DLc15 - KI.DLc16 - KI.DLc23 + KI.DLc24 - KI.DLc25 + KI.DLc26 + KI.T13 + pow(KI.T13,2) - KI.T14 - pow(KI.T14,2) + KI.T15 + pow(KI.T15,2) - KI.T16 - pow(KI.T16,2) - KI.T23 - pow(KI.T23,2) + KI.T24 + pow(KI.T24,2) - KI.T25 - pow(KI.T25,2) + KI.T26 + pow(KI.T26,2)))/(2.*sqrt(2)) + (sqrt(-4 + pow(CA,2))*(-KI.T13 + KI.T14 - KI.T15 + KI.T16 + KI.T23 - KI.T24 + KI.T25 - KI.T26)*log(muR2/M2))/(2.*sqrt(2)))/CA;
    S1tilde26 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*(-KI.DLc13 - KI.DLc14 + KI.DLc15 + KI.DLc16 + KI.DLc23 + KI.DLc24 - KI.DLc25 - KI.DLc26 - KI.T13 - pow(KI.T13,2) - KI.T14 - pow(KI.T14,2) + KI.T15 + pow(KI.T15,2) + KI.T16 + pow(KI.T16,2) + KI.T23 + pow(KI.T23,2) + KI.T24 + pow(KI.T24,2) - KI.T25 - pow(KI.T25,2) - KI.T26 - pow(KI.T26,2))/(2.*sqrt(2)) + ((KI.T13 + KI.T14 - KI.T15 - KI.T16 - KI.T23 - KI.T24 + KI.T25 + KI.T26)*log(muR2/M2))/(2.*sqrt(2));
    S1tilde33 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*CA*((2 + 4*KI.DL34 + 4*KI.DLc15 + 4*KI.DLc26 + KI.L2beta34 + 2*KI.Lbeta3 + 2*KI.Lbeta4 + 2*KI.Lbeta5 + 2*KI.Lbeta6 + 4*KI.T15 + 4*pow(KI.T15,2) + 4*KI.T26 + 4*pow(KI.T26,2))/8. + ((4 + 4*KI.L34 - 4*KI.T15 - 4*KI.T26)*log(muR2/M2))/8.) + ((-4*KI.DL34 - 4*KI.DL56 - 8*KI.DLc15 + 8*KI.DLc16 + 8*KI.DLc25 - 8*KI.DLc26 - KI.L2beta34 - KI.L2beta56 - 2*KI.Lbeta3 - 2*KI.Lbeta4 - 2*KI.Lbeta5 - 2*KI.Lbeta6 - 8*KI.T15 - 8*pow(KI.T15,2) + 8*KI.T16 + 8*pow(KI.T16,2) + 8*KI.T25 + 8*pow(KI.T25,2) - 8*KI.T26 - 8*pow(KI.T26,2))/8. + ((-8 - 4*KI.L34 - 4*KI.L56 + 8*KI.T15 - 8*KI.T16 - 8*KI.T25 + 8*KI.T26)*log(muR2/M2))/8.)/CA;
    S1tilde34 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*((-4*KI.DL35 + 4*KI.DL36 + 4*KI.DL45 - 4*KI.DL46 - KI.L2beta35 + KI.L2beta36 + KI.L2beta45 - KI.L2beta46)/8. - ((KI.L35 - KI.L36 - KI.L45 + KI.L46)*log(muR2/M2))/2.)/CA;
    S1tilde35 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*(-(sqrt(-4 + pow(CA,2))*(4*KI.DL35 - 4*KI.DL36 - 4*KI.DL45 + 4*KI.DL46 - 4*KI.DLc13 + 4*KI.DLc14 + 4*KI.DLc23 - 4*KI.DLc24 + KI.L2beta35 - KI.L2beta36 - KI.L2beta45 + KI.L2beta46 - 4*KI.T13 - 4*pow(KI.T13,2) + 4*KI.T14 + 4*pow(KI.T14,2) + 4*KI.T23 + 4*pow(KI.T23,2) - 4*KI.T24 - 4*pow(KI.T24,2)))/(8.*sqrt(2)) - (sqrt(-4 + pow(CA,2))*(KI.L35 - KI.L36 - KI.L45 + KI.L46 + KI.T13 - KI.T14 - KI.T23 + KI.T24)*log(muR2/M2))/(2.*sqrt(2)))/CA;
    S1tilde36 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*(4*KI.DL35 + 4*KI.DL36 - 4*KI.DL45 - 4*KI.DL46 - 4*KI.DLc13 + 4*KI.DLc14 - 4*KI.DLc23 + 4*KI.DLc24 + KI.L2beta35 + KI.L2beta36 - KI.L2beta45 - KI.L2beta46 - 4*KI.T13 - 4*pow(KI.T13,2) + 4*KI.T14 + 4*pow(KI.T14,2) - 4*KI.T23 - 4*pow(KI.T23,2) + 4*KI.T24 + 4*pow(KI.T24,2))/(8.*sqrt(2)) + ((KI.L35 + KI.L36 - KI.L45 - KI.L46 + KI.T13 - KI.T14 + KI.T23 - KI.T24)*log(muR2/M2))/(2.*sqrt(2));
    S1tilde44 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*CA*((2 + 4*KI.DL56 + 4*KI.DLc13 + 4*KI.DLc24 + KI.L2beta56 + 2*KI.Lbeta3 + 2*KI.Lbeta4 + 2*KI.Lbeta5 + 2*KI.Lbeta6 + 4*KI.T13 + 4*pow(KI.T13,2) + 4*KI.T24 + 4*pow(KI.T24,2))/8. + ((4 + 4*KI.L56 - 4*KI.T13 - 4*KI.T24)*log(muR2/M2))/8.) + ((-4*KI.DL34 - 4*KI.DL56 - 8*KI.DLc13 + 8*KI.DLc14 + 8*KI.DLc23 - 8*KI.DLc24 - KI.L2beta34 - KI.L2beta56 - 2*KI.Lbeta3 - 2*KI.Lbeta4 - 2*KI.Lbeta5 - 2*KI.Lbeta6 - 8*KI.T13 - 8*pow(KI.T13,2) + 8*KI.T14 + 8*pow(KI.T14,2) + 8*KI.T23 + 8*pow(KI.T23,2) - 8*KI.T24 - 8*pow(KI.T24,2))/8. + ((-8 - 4*KI.L34 - 4*KI.L56 + 8*KI.T13 - 8*KI.T14 - 8*KI.T23 + 8*KI.T24)*log(muR2/M2))/8.)/CA;
    S1tilde45 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*(-(sqrt(-4 + pow(CA,2))*(4*KI.DL35 - 4*KI.DL36 - 4*KI.DL45 + 4*KI.DL46 - 4*KI.DLc15 + 4*KI.DLc16 + 4*KI.DLc25 - 4*KI.DLc26 + KI.L2beta35 - KI.L2beta36 - KI.L2beta45 + KI.L2beta46 - 4*KI.T15 - 4*pow(KI.T15,2) + 4*KI.T16 + 4*pow(KI.T16,2) + 4*KI.T25 + 4*pow(KI.T25,2) - 4*KI.T26 - 4*pow(KI.T26,2)))/(8.*sqrt(2)) - (sqrt(-4 + pow(CA,2))*(KI.L35 - KI.L36 - KI.L45 + KI.L46 + KI.T15 - KI.T16 - KI.T25 + KI.T26)*log(muR2/M2))/(2.*sqrt(2)))/CA;
    S1tilde46 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*(-4*KI.DL35 + 4*KI.DL36 - 4*KI.DL45 + 4*KI.DL46 + 4*KI.DLc15 - 4*KI.DLc16 + 4*KI.DLc25 - 4*KI.DLc26 - KI.L2beta35 + KI.L2beta36 - KI.L2beta45 + KI.L2beta46 + 4*KI.T15 + 4*pow(KI.T15,2) - 4*KI.T16 - 4*pow(KI.T16,2) + 4*KI.T25 + 4*pow(KI.T25,2) - 4*KI.T26 - 4*pow(KI.T26,2))/(8.*sqrt(2)) - ((KI.L35 - KI.L36 + KI.L45 - KI.L46 + KI.T15 - KI.T16 + KI.T25 - KI.T26)*log(muR2/M2))/(2.*sqrt(2));
    S1tilde55 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*CA*((4 + 4*KI.DL36 + 4*KI.DL45 + 4*KI.DLc13 + 4*KI.DLc15 + 4*KI.DLc24 + 4*KI.DLc26 + KI.L2beta36 + KI.L2beta45 + 4*KI.Lbeta3 + 4*KI.Lbeta4 + 4*KI.Lbeta5 + 4*KI.Lbeta6 + 4*KI.T13 + 4*pow(KI.T13,2) + 4*KI.T15 + 4*pow(KI.T15,2) + 4*KI.T24 + 4*pow(KI.T24,2) + 4*KI.T26 + 4*pow(KI.T26,2))/16. + ((8 + 4*KI.L36 + 4*KI.L45 - 4*KI.T13 - 4*KI.T15 - 4*KI.T24 - 4*KI.T26)*log(muR2/M2))/16.) + ((-8*KI.DL34 + 24*KI.DL35 - 24*KI.DL36 - 24*KI.DL45 + 24*KI.DL46 - 8*KI.DL56 - 24*KI.DLc13 + 24*KI.DLc14 - 24*KI.DLc15 + 24*KI.DLc16 + 24*KI.DLc23 - 24*KI.DLc24 + 24*KI.DLc25 - 24*KI.DLc26 - 2*KI.L2beta34 + 6*KI.L2beta35 - 6*KI.L2beta36 - 6*KI.L2beta45 + 6*KI.L2beta46 - 2*KI.L2beta56 - 4*KI.Lbeta3 - 4*KI.Lbeta4 - 4*KI.Lbeta5 - 4*KI.Lbeta6 - 24*KI.T13 - 24*pow(KI.T13,2) + 24*KI.T14 + 24*pow(KI.T14,2) - 24*KI.T15 - 24*pow(KI.T15,2) + 24*KI.T16 + 24*pow(KI.T16,2) + 24*KI.T23 + 24*pow(KI.T23,2) - 24*KI.T24 - 24*pow(KI.T24,2) + 24*KI.T25 + 24*pow(KI.T25,2) - 24*KI.T26 - 24*pow(KI.T26,2))/16. + ((-16 - 8*KI.L34 + 24*KI.L35 - 24*KI.L36 - 24*KI.L45 + 24*KI.L46 - 8*KI.L56 + 24*KI.T13 - 24*KI.T14 + 24*KI.T15 - 24*KI.T16 - 24*KI.T23 + 24*KI.T24 - 24*KI.T25 + 24*KI.T26)*log(muR2/M2))/16.)/CA;
    S1tilde56 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*(sqrt(-4 + pow(CA,2))*(4*KI.DL36 - 4*KI.DL45 - 4*KI.DLc13 + 4*KI.DLc15 + 4*KI.DLc24 - 4*KI.DLc26 + KI.L2beta36 - KI.L2beta45 - 4*KI.T13 - 4*pow(KI.T13,2) + 4*KI.T15 + 4*pow(KI.T15,2) + 4*KI.T24 + 4*pow(KI.T24,2) - 4*KI.T26 - 4*pow(KI.T26,2)))/16. + (sqrt(-4 + pow(CA,2))*(KI.L36 - KI.L45 + KI.T13 - KI.T15 - KI.T24 + KI.T26)*log(muR2/M2))/4.;
    S1tilde66 = include_S1*ISNLL*alphas_muR/M_PI*include_S1*CA*((4 + 4*KI.DL36 + 4*KI.DL45 + 4*KI.DLc13 + 4*KI.DLc15 + 4*KI.DLc24 + 4*KI.DLc26 + KI.L2beta36 + KI.L2beta45 + 4*KI.Lbeta3 + 4*KI.Lbeta4 + 4*KI.Lbeta5 + 4*KI.Lbeta6 + 4*KI.T13 + 4*pow(KI.T13,2) + 4*KI.T15 + 4*pow(KI.T15,2) + 4*KI.T24 + 4*pow(KI.T24,2) + 4*KI.T26 + 4*pow(KI.T26,2))/16. + ((8 + 4*KI.L36 + 4*KI.L45 - 4*KI.T13 - 4*KI.T15 - 4*KI.T24 - 4*KI.T26)*log(muR2/M2))/16.) + ((-8*KI.DL34 + 8*KI.DL35 - 8*KI.DL36 - 8*KI.DL45 + 8*KI.DL46 - 8*KI.DL56 - 8*KI.DLc13 + 8*KI.DLc14 - 8*KI.DLc15 + 8*KI.DLc16 + 8*KI.DLc23 - 8*KI.DLc24 + 8*KI.DLc25 - 8*KI.DLc26 - 2*KI.L2beta34 + 2*KI.L2beta35 - 2*KI.L2beta36 - 2*KI.L2beta45 + 2*KI.L2beta46 - 2*KI.L2beta56 - 4*KI.Lbeta3 - 4*KI.Lbeta4 - 4*KI.Lbeta5 - 4*KI.Lbeta6 - 8*KI.T13 - 8*pow(KI.T13,2) + 8*KI.T14 + 8*pow(KI.T14,2) - 8*KI.T15 - 8*pow(KI.T15,2) + 8*KI.T16 + 8*pow(KI.T16,2) + 8*KI.T23 + 8*pow(KI.T23,2) - 8*KI.T24 - 8*pow(KI.T24,2) + 8*KI.T25 + 8*pow(KI.T25,2) - 8*KI.T26 - 8*pow(KI.T26,2))/16. + ((-16 - 8*KI.L34 + 8*KI.L35 - 8*KI.L36 - 8*KI.L45 + 8*KI.L46 - 8*KI.L56 + 8*KI.T13 - 8*KI.T14 + 8*KI.T15 - 8*KI.T16 - 8*KI.T23 + 8*KI.T24 - 8*KI.T25 + 8*KI.T26)*log(muR2/M2))/16.)/CA;
    S1tilde21 = S1tilde12;
    S1tilde31 = S1tilde13;
    S1tilde41 = S1tilde14;
    S1tilde51 = 0;
    S1tilde61 = 0;
    S1tilde32 = S1tilde23;
    S1tilde42 = S1tilde24;
    S1tilde52 = S1tilde25;
    S1tilde62 = S1tilde26;
    S1tilde43 = S1tilde34;
    S1tilde53 = S1tilde35;
    S1tilde63 = S1tilde36;
    S1tilde54 = S1tilde45;
    S1tilde64 = S1tilde46;
    S1tilde65 = S1tilde56;

    Eigen::Matrix<complex<double>, 6, 6> S1tilde_fin;
    S1tilde_fin <<  S1tilde11 ,S1tilde12 ,S1tilde13 ,S1tilde14, S1tilde15, S1tilde16,
                    S1tilde21 ,S1tilde22 ,S1tilde23 ,S1tilde24, S1tilde25, S1tilde26,
                    S1tilde31 ,S1tilde32 ,S1tilde33 ,S1tilde34, S1tilde35, S1tilde36,
                    S1tilde41 ,S1tilde42 ,S1tilde43 ,S1tilde44, S1tilde45, S1tilde46,
                    S1tilde51 ,S1tilde52 ,S1tilde53 ,S1tilde54, S1tilde55, S1tilde56,
                    S1tilde61 ,S1tilde62 ,S1tilde63 ,S1tilde64, S1tilde65, S1tilde66;
    
    return S1tilde_fin;

}


void check(){
    
    complex<double> Lbeta34, Lbeta35, Lbeta36, Lbeta45, Lbeta46, Lbeta56;
    complex<double> T13, T14, T15, T16, T23, T24, T25,  T26;
    double sqrts = sqrt(1000000.);
    //double mt    = 172.5;
    vector<double> p1 = {sqrts/2,0,0,sqrts/2};
    vector<double> p2 = {sqrts/2,0,0,-sqrts/2};
    vector<double> p3 = {sqrt(pow(mt,2) + pow(158.24,2) + pow(159.74,2)  + pow(200.5,2)), 158.24, 159.74, 200.5};
    vector<double> p4 = {sqrt(pow(mt,2) + pow(120.,2)   + pow(89.74,2)   + pow(100.5,2)), 120., 89.74, 100.5};
    vector<double> p5 = {sqrt(pow(mt,2) + pow(158.24,2) + pow(159.74,2)  + pow(100.5,2)), -158.24, 159.74, 100.5};
    vector<double> p6 = {sqrt(pow(mt,2) + pow(120.,2)   + pow(409.22,2)  + pow(401.5,2)), -120., -409.22, -401.5};
    Lbeta34 = Lbetaij(p3, p4);
    Lbeta35 = Lbetaij(p3, p5);
    Lbeta36 = Lbetaij(p3, p6);
    Lbeta45 = Lbetaij(p4, p5);
    Lbeta46 = Lbetaij(p4, p6);
    Lbeta56 = Lbetaij(p5, p6);
    T13     = Tij(p1, p3, sqrts, mt);
    T14     = Tij(p1, p4, sqrts, mt);
    T15     = Tij(p1, p5, sqrts, mt);
    T16     = Tij(p1, p6, sqrts, mt);
    T23     = Tij(p2, p3, sqrts, mt);
    T24     = Tij(p2, p4, sqrts, mt);
    T25     = Tij(p2, p5, sqrts, mt);
    T26     = Tij(p2, p6, sqrts, mt);

    complex<double> G11 ,G12 ,G13 ,G14, G15, G16,
    G21 ,G22 ,G23 ,G24, G25, G26,
    G31 ,G32 ,G33 ,G34, G35, G36,
    G41 ,G42 ,G43 ,G44, G45, G46,
    G51 ,G52 ,G53 ,G54, G55, G56,
    G61 ,G62 ,G63 ,G64, G65, G66;
    //set_threshold();
    G11 = -CF*(Lbeta34 + Lbeta56 + 2.);
    G12 = sqrt(pow(CA,2) - 1)/(2*CA)*(Lbeta35+Lbeta46-Lbeta36-Lbeta45);
    G21 = G12;
    G13 = sqrt(pow(CA,2) - 1)/(2*CA)*(T15 - T16 - T25 + T26);
    G31 = G13;
    G14 = sqrt(pow(CA,2) - 1)/(2*CA)*(T13 - T14 - T23 + T24);
    G41 = G14;
    G15 = 0;
    G51 = 0;
    G16 = 0;
    G61 = 0;
    G22 = -2*CF + 1/(2*CA)*(Lbeta34+Lbeta56) - 1/CA*(Lbeta35+Lbeta46) - (pow(CA,2)-2)/(2*CA)*(Lbeta36 + Lbeta45);
    G23 = 1/(2*CA)*(T13 - T14 - T23 + T24);
    G32 = G23;
    G24 = 1/(2*CA)*(T15 - T16 - T25 + T26);
    G42 = G24;
    G25 = sqrt(pow(CA,2) - 4)/(2*sqrt(2)*CA)*(T13 - T14 +T15 - T16 - T23 + T24 - T25 + T26);
    G52 = G25;
    G26 = 1/(2*sqrt(2))*(-T13 - T14 +T15 + T16 + T23 + T24 - T25 - T26);
    G62 = G26;
    G33 = -(pow(CA,2) - 2)/(2*CA) - CF*Lbeta34 + 1/(2*CA)*Lbeta56 + (pow(CA,2)-2)/(2*CA)*(T15+T26) + 1/CA*(T25+T16);
    G34 = G12/sqrt(pow(CA,2) - 1);
    G43 = G34;
    G35 = sqrt(pow(CA,2)-4)/(2*sqrt(2)*CA)*(Lbeta35-Lbeta36-Lbeta45+Lbeta46+T13 - T14 - T23 + T24);
    G53 = G35;
    G36 = 1/(2*sqrt(2))*(-Lbeta35-Lbeta36+Lbeta45+Lbeta46 - T13 + T14 - T23 + T24);
    G63 = G36;
    G44 = -CF + 1/(2*CA) - CF*Lbeta56 + 1/(2*CA)*Lbeta34 + (pow(CA,2)-2)/(2*CA)*(T13+T24) + 1/CA*(T14+T23);
    G45 = sqrt(pow(CA,2)-4)/(2*sqrt(2)*CA)*(Lbeta35-Lbeta36-Lbeta45+Lbeta46+T15 - T16 - T25 + T26);
    G54 = G45;
    G46 = 1/(2*sqrt(2))*(Lbeta35-Lbeta36+Lbeta45-Lbeta46 + T15 - T16 + T25 - T26);
    G64 = G46;
    G55 = -CF + 1/(2*CA) +1/(2*CA)*(Lbeta34 - 3.*Lbeta35 - 3.*Lbeta46 + Lbeta56 - (pow(CA,2) - 6)/2*(Lbeta36 + Lbeta45))
              + 3/(2*CA)*(T14+T16+T23+T25) + (pow(CA,2)-6)/(4*CA)*(T13+T15+T24+T26);
    G56 = sqrt(pow(CA,2)-4)/(4)*(-Lbeta36+Lbeta45-T13 + T15 + T24 - T26);
    G65 = G56;
    G66 = -CF + 1/(2*CA) +1/(2*CA)*(Lbeta34 - Lbeta35 - Lbeta46 + Lbeta56 - (pow(CA,2) - 2)/(2)*(Lbeta36 + Lbeta45))
              + 1/(2*CA)*(T14+T16+T23+T25) + (pow(CA,2)-2)/(4*CA)*(T13+T15+T24+T26);

    Eigen::Matrix<complex<double>, 6, 6> sad;
    sad <<  G11 ,G12 ,G13 ,G14, G15, G16,
            G21 ,G22 ,G23 ,G24, G25, G26,
            G31 ,G32 ,G33 ,G34, G35, G36,
            G41 ,G42 ,G43 ,G44, G45, G46,
            G51 ,G52 ,G53 ,G54, G55, G56,
            G61 ,G62 ,G63 ,G64, G65, G66;
    cout << "SAD " << endl << sad << endl << endl;


    Eigen::ComplexEigenSolver<Eigen::Matrix<complex<double>,6,6>> ces;
    ces.compute(sad);
    cout << "KI.The eigenvalues of sad are:" << endl << ces.eigenvalues() << endl;
    cout << "sad = V * D * V^(-1) = " << endl
        << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << endl;
    cout << "And    , D = V^(-1) * sad * V^(-1) = " << endl << ces.eigenvectors().inverse() * sad * ces.eigenvectors() << endl;
    cout << "DiffereCAe = " << endl << ces.eigenvectors().inverse() * sad * ces.eigenvectors() << endl;

 }

