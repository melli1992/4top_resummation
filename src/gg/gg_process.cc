//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 3.1.0, 2021-03-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "gg_process.h"
#include "gg_helamps.h"
#include "Parameters_sm.h"

using namespace gg_MG5_sm;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ t t~ WEIGHTED<=4 @1

//--------------------------------------------------------------------------
// Initialize process.

void gg_process::initProc()
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm::getInstance();
  pars->setIndependentParameters();
  pars->setIndependentCouplings();
  pars->printIndependentParameters();
  pars->printIndependentCouplings();
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO);
  mME.push_back(pars->ZERO);
  mME.push_back(pars->mdl_MT);
  mME.push_back(pars->mdl_MT);
  mME.push_back(pars->mdl_MT);
  mME.push_back(pars->mdl_MT);
  jamp2[0] = new double[ncol];

  matrix_element[0] = new double[ncol];
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void gg_process::sigmaKin()
{
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
  static bool firsttime = true;
  if (firsttime)
  {
    pars->printDependentParameters();
    pars->printDependentCouplings();
    firsttime = false;
  }

  // Reset color flows
  for(int i = 0; i < 12; i++ )
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 64;
  static bool goodhel[ncomb] = {ncomb * false};
  static int ntry = 0, sum_hel = 0, ngood = 0;
  static int igood[ncomb];
  static int jhel;
  std::complex<double> * * wfs;
  double t[nprocesses];
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {1024};

  ntry = ntry + 1;

  // Reset the matrix elements
  // for(int i = 0; i < nprocesses; i++){
  // matrix_element[i] = 0.;
  // }
  for(int i = 0; i < ncol; i++ )
  {
    matrix_element[0][i] = 0.;
  }

  // Define permutation
  int perm[nexternal];
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i;
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]);
        t[0] = matrix_1_gg_ttxttx();

        for(int icol = 0; icol < ncol; icol++ )
        {
          matrix_element[0][icol] += jamp2[0][icol];
        }
        double tsum = 0;
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          // matrix_element[iproc]+=t[iproc];
          tsum += t[iproc];
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true;
          ngood++;
          igood[ngood] = ihel;
        }
      }
    }
    jhel = 0;
    sum_hel = min(sum_hel, ngood);
  }
  // else
  // {
  // Only use the "good" helicities
  // for(int j=0; j < sum_hel; j++){
  // jhel++;
  // if (jhel >= ngood) jhel=0;
  // double hwgt = double(ngood)/double(sum_hel);
  // int ihel = igood[jhel];
  // calculate_wavefunctions(perm, helicities[ihel]);
  // t[0]=matrix_1_gg_ttxttx();
  //
  // for(int iproc = 0;iproc < nprocesses; iproc++){
  // matrix_element[iproc]+=t[iproc]*hwgt;
  // }
  // }
  // }

  // for (int i=0;i < nprocesses; i++)
  // matrix_element[i] /= denominators[i];
  for (int i = 0; i < ncol; i++ )
    matrix_element[0][i] /= denominators[0];
    

}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void gg_process::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j;

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]);
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]);
  VVV1P0_1(w[0], w[1], pars->GC_10, pars->ZERO, pars->ZERO, w[6]);
  FFV1P0_3(w[3], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[7]);
  FFV1_1(w[4], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[8]);
  FFV1_2(w[5], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[9]);
  FFV1P0_3(w[5], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[10]);
  FFV1P0_3(w[5], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[11]);
  FFV1_2(w[3], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[12]);
  FFV1P0_3(w[3], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[13]);
  FFV1_1(w[2], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[14]);
  FFV1_1(w[2], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[15]);
  FFV1_2(w[3], w[1], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[16]);
  FFV1P0_3(w[5], w[15], pars->GC_11, pars->ZERO, pars->ZERO, w[17]);
  FFV1_1(w[4], w[1], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[18]);
  FFV1P0_3(w[3], w[15], pars->GC_11, pars->ZERO, pars->ZERO, w[19]);
  FFV1_2(w[5], w[1], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[20]);
  FFV1_1(w[15], w[1], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[21]);
  FFV1_2(w[3], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[22]);
  FFV1_1(w[2], w[1], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[23]);
  FFV1P0_3(w[22], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[24]);
  FFV1P0_3(w[22], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[25]);
  FFV1_2(w[22], w[1], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[26]);
  FFV1_1(w[4], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[27]);
  FFV1P0_3(w[3], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[28]);
  FFV1P0_3(w[5], w[27], pars->GC_11, pars->ZERO, pars->ZERO, w[29]);
  FFV1_1(w[27], w[1], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[30]);
  FFV1_2(w[5], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[31]);
  FFV1P0_3(w[31], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[32]);
  FFV1P0_3(w[31], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[33]);
  FFV1_2(w[31], w[1], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[34]);
  FFV1_1(w[23], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[35]);
  VVV1P0_1(w[0], w[13], pars->GC_10, pars->ZERO, pars->ZERO, w[36]);
  VVV1P0_1(w[0], w[10], pars->GC_10, pars->ZERO, pars->ZERO, w[37]);
  FFV1_2(w[16], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[38]);
  VVV1P0_1(w[0], w[11], pars->GC_10, pars->ZERO, pars->ZERO, w[39]);
  FFV1_1(w[18], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[40]);
  VVV1P0_1(w[0], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[41]);
  FFV1_2(w[20], w[0], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[42]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[5], w[8], w[7], pars->GC_11, amp[0]);
  FFV1_0(w[9], w[4], w[7], pars->GC_11, amp[1]);
  VVV1_0(w[6], w[7], w[10], pars->GC_10, amp[2]);
  FFV1_0(w[12], w[4], w[11], pars->GC_11, amp[3]);
  FFV1_0(w[3], w[8], w[11], pars->GC_11, amp[4]);
  VVV1_0(w[6], w[11], w[13], pars->GC_10, amp[5]);
  FFV1_0(w[5], w[14], w[13], pars->GC_11, amp[6]);
  FFV1_0(w[9], w[2], w[13], pars->GC_11, amp[7]);
  FFV1_0(w[3], w[14], w[10], pars->GC_11, amp[8]);
  FFV1_0(w[12], w[2], w[10], pars->GC_11, amp[9]);
  FFV1_0(w[16], w[4], w[17], pars->GC_11, amp[10]);
  FFV1_0(w[16], w[15], w[10], pars->GC_11, amp[11]);
  FFV1_0(w[5], w[18], w[19], pars->GC_11, amp[12]);
  FFV1_0(w[3], w[18], w[17], pars->GC_11, amp[13]);
  FFV1_0(w[20], w[4], w[19], pars->GC_11, amp[14]);
  FFV1_0(w[20], w[15], w[13], pars->GC_11, amp[15]);
  FFV1_0(w[5], w[21], w[13], pars->GC_11, amp[16]);
  VVV1_0(w[1], w[13], w[17], pars->GC_10, amp[17]);
  FFV1_0(w[3], w[21], w[10], pars->GC_11, amp[18]);
  VVV1_0(w[1], w[10], w[19], pars->GC_10, amp[19]);
  FFV1_0(w[5], w[23], w[24], pars->GC_11, amp[20]);
  FFV1_0(w[22], w[23], w[10], pars->GC_11, amp[21]);
  FFV1_0(w[5], w[18], w[25], pars->GC_11, amp[22]);
  FFV1_0(w[22], w[18], w[11], pars->GC_11, amp[23]);
  FFV1_0(w[20], w[4], w[25], pars->GC_11, amp[24]);
  FFV1_0(w[20], w[2], w[24], pars->GC_11, amp[25]);
  FFV1_0(w[26], w[4], w[11], pars->GC_11, amp[26]);
  VVV1_0(w[1], w[11], w[24], pars->GC_10, amp[27]);
  FFV1_0(w[26], w[2], w[10], pars->GC_11, amp[28]);
  VVV1_0(w[1], w[10], w[25], pars->GC_10, amp[29]);
  FFV1_0(w[5], w[23], w[28], pars->GC_11, amp[30]);
  FFV1_0(w[3], w[23], w[29], pars->GC_11, amp[31]);
  FFV1_0(w[16], w[2], w[29], pars->GC_11, amp[32]);
  FFV1_0(w[16], w[27], w[11], pars->GC_11, amp[33]);
  FFV1_0(w[20], w[2], w[28], pars->GC_11, amp[34]);
  FFV1_0(w[20], w[27], w[7], pars->GC_11, amp[35]);
  FFV1_0(w[5], w[30], w[7], pars->GC_11, amp[36]);
  VVV1_0(w[1], w[7], w[29], pars->GC_10, amp[37]);
  FFV1_0(w[3], w[30], w[11], pars->GC_11, amp[38]);
  VVV1_0(w[1], w[11], w[28], pars->GC_10, amp[39]);
  FFV1_0(w[3], w[23], w[32], pars->GC_11, amp[40]);
  FFV1_0(w[31], w[23], w[13], pars->GC_11, amp[41]);
  FFV1_0(w[16], w[4], w[33], pars->GC_11, amp[42]);
  FFV1_0(w[16], w[2], w[32], pars->GC_11, amp[43]);
  FFV1_0(w[3], w[18], w[33], pars->GC_11, amp[44]);
  FFV1_0(w[31], w[18], w[7], pars->GC_11, amp[45]);
  FFV1_0(w[34], w[4], w[7], pars->GC_11, amp[46]);
  VVV1_0(w[1], w[7], w[32], pars->GC_10, amp[47]);
  FFV1_0(w[34], w[2], w[13], pars->GC_11, amp[48]);
  VVV1_0(w[1], w[13], w[33], pars->GC_10, amp[49]);
  FFV1_0(w[5], w[35], w[13], pars->GC_11, amp[50]);
  FFV1_0(w[5], w[23], w[36], pars->GC_11, amp[51]);
  FFV1_0(w[3], w[35], w[10], pars->GC_11, amp[52]);
  FFV1_0(w[3], w[23], w[37], pars->GC_11, amp[53]);
  FFV1_0(w[38], w[4], w[11], pars->GC_11, amp[54]);
  FFV1_0(w[16], w[4], w[39], pars->GC_11, amp[55]);
  FFV1_0(w[38], w[2], w[10], pars->GC_11, amp[56]);
  FFV1_0(w[16], w[2], w[37], pars->GC_11, amp[57]);
  FFV1_0(w[5], w[40], w[7], pars->GC_11, amp[58]);
  FFV1_0(w[5], w[18], w[41], pars->GC_11, amp[59]);
  FFV1_0(w[3], w[40], w[11], pars->GC_11, amp[60]);
  FFV1_0(w[3], w[18], w[39], pars->GC_11, amp[61]);
  FFV1_0(w[42], w[4], w[7], pars->GC_11, amp[62]);
  FFV1_0(w[20], w[4], w[41], pars->GC_11, amp[63]);
  FFV1_0(w[42], w[2], w[13], pars->GC_11, amp[64]);
  FFV1_0(w[20], w[2], w[36], pars->GC_11, amp[65]);
  VVVV1_0(w[0], w[1], w[7], w[10], pars->GC_12, amp[66]);
  VVVV3_0(w[0], w[1], w[7], w[10], pars->GC_12, amp[67]);
  VVVV4_0(w[0], w[1], w[7], w[10], pars->GC_12, amp[68]);
  VVV1_0(w[1], w[10], w[41], pars->GC_10, amp[69]);
  VVV1_0(w[1], w[7], w[37], pars->GC_10, amp[70]);
  VVVV1_0(w[0], w[1], w[11], w[13], pars->GC_12, amp[71]);
  VVVV3_0(w[0], w[1], w[11], w[13], pars->GC_12, amp[72]);
  VVVV4_0(w[0], w[1], w[11], w[13], pars->GC_12, amp[73]);
  VVV1_0(w[1], w[13], w[39], pars->GC_10, amp[74]);
  VVV1_0(w[1], w[11], w[36], pars->GC_10, amp[75]);

}
double gg_process::matrix_1_gg_ttxttx()
{
  int i, j;
  // Local variables
  const int ngraphs = 76;
  const int ncolor = 12;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncol];
  // The color matrix;
  /*static const double denom[ncolor] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  static const double cf[ncolor][ncolor] = {{48, 16, 16, 6, 0, 16, -2, 0, -6,
      -2, -2, 6}, {16, 48, 6, 16, 16, 0, 0, -2, -2, -6, 6, -2}, {16, 6, 48, 16,
      -2, 0, 0, 16, -2, 6, -6, -2}, {6, 16, 16, 48, 0, -2, 16, 0, 6, -2, -2,
      -6}, {0, 16, -2, 0, 48, 16, 16, 6, 0, -2, 16, 0}, {16, 0, 0, -2, 16, 48,
      6, 16, -2, 0, 0, 16}, {-2, 0, 0, 16, 16, 6, 48, 16, 16, 0, 0, -2}, {0,
      -2, 16, 0, 6, 16, 16, 48, 0, 16, -2, 0}, {-6, -2, -2, 6, 0, -2, 16, 0,
      48, 16, 16, 6}, {-2, -6, 6, -2, -2, 0, 0, 16, 16, 48, 6, 16}, {-2, 6, -6,
      -2, 16, 0, 0, -2, 16, 6, 48, 16}, {6, -2, -2, -6, 0, 16, -2, 0, 6, 16,
      16, 48}};
  /*static const double fromTrtoOrtho[ncol][ncol] = {{0, (4 * 1.)/3, -1./6, 0, 4
      * 1., (4 * 1.)/3, (4 * 1.)/3, 1./2, 0, -1./6, (4 * 1.)/3, 0, 0, 0},
      {sqrt(2.) * 1., (sqrt(2.) * 1)/3, (sqrt(2.) * 1.)/3, sqrt(2.) * 1., 0,
      (sqrt(2.) * 1.)/3, (sqrt(2.) * 1.)/3, 0, sqrt(2.) * 1., (sqrt(2.) *
      1.)/3, (sqrt(2.) * 1.)/3, sqrt(2.) * 1., 0, 0}, {0, (sqrt(5.) * 1)/3,
      (sqrt(5.) * 1.)/3, sqrt(5.) * 1., 0, (sqrt(5.) * 1.)/3, (sqrt(5.) *
      1.)/3, 0, 0, (sqrt(5.) * 1.)/3, (sqrt(5.) * 1.)/3, sqrt(5.) * 1., 0, 0},
      {0, 1., 1., 3 * 1., 0, -1., 1., 0, 0, -1., -1., -3 * 1., 0, 0}, {sqrt(5.)
      * 1., (sqrt(5.) * 1)/3, (sqrt(5.) * 1.)/3, 0, 0, (sqrt(5.) * 1.)/3,
      (sqrt(5.) * 1.)/3, 0, sqrt(5.) * 1., (sqrt(5.) * 1.)/3, (sqrt(5.) *
      1.)/3, 0, 0, 0}, {0, (5 * 1)/(3 * sqrt(2.)), (-2 * sqrt(2.) * 1.)/3, 0,
      0, (5 * 1.)/(3 * sqrt(2.)), (5 * 1.)/(3 * sqrt(2.)), sqrt(2.) * 1., 0,
      (-2 * sqrt(2.) * 1.)/3, (5 * 1.)/(3 * sqrt(2.)), 0, 0, 0}, {0,
      sqrt(5./2.) * 1, 0, 0, 0, -(sqrt(5./2.) * 1.), sqrt(5./2.) * 1., 0, 0, 0,
      -(sqrt(5./2.) * 1.), 0, 0, 0}, {3., 1., 1., 0, 0, 1., -1., 0, -3 * 1.,
      -1., -1., 0, 0, 0}, {0, sqrt(5./2.) * 1, 0, 0, 0, sqrt(5./2.) * 1.,
      -(sqrt(5./2.) * 1.), 0, 0, 0, -(sqrt(5./2.) * 1.), 0, 0, 0}, {0, (3 *
      1)/sqrt(2.), 0, 0, 0, (-3 * 1.)/sqrt(2.), (-3 * 1.)/sqrt(2.), -(sqrt(2.)
      * 1.), 0, 0, (3 * 1.)/sqrt(2.), 0, 0, 0}, {0, 0, -(sqrt(5./2.) * 1.), 0,
      0, 0, 0, -(sqrt(5./2.) * 1.), 0, sqrt(5./2.) * 1., 0, 0, 0, 0}, {0, 0,
      sqrt(5./2.) * 1., 0, 0, 0, 0, -(sqrt(5./2.) * 1.), 0, -(sqrt(5./2.) *
      1.), 0, 0, 0, 0}, {0, 0, (3 * sqrt(3.) * 1.)/2, 0, 0, 0, 0, (3 * sqrt(3.)
      * 1.)/2, 0, (3 * sqrt(3.) * 1.)/2, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0}};*/

  static const double fromTrtoDiag[ncol][ncol] = {{0, 0, 0, 0, (3 *
      sqrt(3.))/2., (3 * sqrt(3.))/2., (3 * sqrt(3.))/2., (3 * sqrt(3.))/2., 0,
      0, 0, 0, 0, 0}, {0, 0, 0, 0, -sqrt(5.), 0, 0, sqrt(5.), 0, 0, 0, 0, 0,
      0}, {0, 0, 0, 0, 0, sqrt(5.), -sqrt(5.), 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -3/sqrt(2.), 3/sqrt(2.), 0,
      -sqrt(2.), 0, 0, sqrt(2.), 0, 3/sqrt(2.), -3/sqrt(2.), 0, 0, 0}, {0, (-13
      * sqrt(2./11.))/3., (8 * sqrt(2./11.))/3., 0, -4 * sqrt(2./11.), (5 *
      sqrt(2./11.))/3., (5 * sqrt(2./11.))/3., -sqrt(2./11.), 0, (-7 *
      sqrt(2./11.))/3., (-13 * sqrt(2./11.))/3., 0, 0, 0}, {0, (-13 *
      sqrt(2./11.))/3., (-7 * sqrt(2./11.))/3., 0, -4 * sqrt(2./11.), (5 *
      sqrt(2./11.))/3., (5 * sqrt(2./11.))/3., -sqrt(2./11.), 0, (8 *
      sqrt(2./11.))/3., (-13 * sqrt(2./11.))/3., 0, 0, 0}, {0, sqrt(5.), 0, 0,
      0, 0, 0, 0, 0, 0, -sqrt(5.), 0, 0, 0}, {3., 1., 1., 0, 0, 1., -1., 0,
      -3., -1., -1., 0, 0, 0}, {sqrt(5.), sqrt(5.)/3., sqrt(5.)/3., 0, 0,
      sqrt(5.)/3., sqrt(5.)/3., 0, sqrt(5.), sqrt(5.)/3., sqrt(5.)/3., 0, 0,
      0}, {0, 1., 1., 3., 0, -1., 1., 0, 0, -1., -1., -3., 0, 0}, {0,
      sqrt(5.)/3., sqrt(5.)/3., sqrt(5.), 0, sqrt(5.)/3., sqrt(5.)/3., 0, 0,
      sqrt(5.)/3., sqrt(5.)/3., sqrt(5.), 0, 0}, {0, 4/3., 4/3., 0, 1./2.,
      -1./6., -1./6., 1./2., 0, 4/3., 4/3., 0, 0, 0}, {sqrt(2.), sqrt(2.)/3.,
      sqrt(2.)/3., sqrt(2.), 0, sqrt(2.)/3., sqrt(2.)/3., 0, sqrt(2.),
      sqrt(2.)/3., sqrt(2.)/3., sqrt(2.), 0, 0}};

  static const double multfact[ncol][ncol] = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0,
      0, 0, 0, 34./25., (-3. * sqrt(11))/25., (-3. * sqrt(11.))/25., 0, 0, 0,
      0, 0, 0, 0}, {0, 0, 0, 0, (-3. * sqrt(11.))/25., 77./50., -33./50., 0, 0,
      0, 0, 0, 0, 0}, {0, 0, 0, 0, (-3 * sqrt(11.))/25., -33./50., 77./50., 0,
      0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
      0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}};

  // Calculate color flows
  jamp[0] = +1./2. * (+std::complex<double> (0, 1) * amp[3] - amp[5] +
      std::complex<double> (0, 1) * amp[6] + 1./3. * std::complex<double> (0,
      1) * amp[8] + 1./3. * std::complex<double> (0, 1) * amp[9] - amp[10] -
      1./3. * amp[11] - amp[16] + std::complex<double> (0, 1) * amp[17] - 1./3.
      * amp[18] - amp[54] - std::complex<double> (0, 1) * amp[55] - 1./3. *
      amp[56] - amp[71] - amp[72] - amp[74]);
  jamp[1] = +1./2. * (-std::complex<double> (0, 1) * amp[1] + amp[2] - 1./3. *
      std::complex<double> (0, 1) * amp[6] - 1./3. * std::complex<double> (0,
      1) * amp[7] - std::complex<double> (0, 1) * amp[8] + amp[14] + 1./3. *
      amp[15] + 1./3. * amp[16] + amp[18] - std::complex<double> (0, 1) *
      amp[19] + amp[62] + std::complex<double> (0, 1) * amp[63] + 1./3. *
      amp[64] + amp[66] + amp[67] + amp[69]);
  jamp[2] = +1./2. * (-std::complex<double> (0, 1) * amp[0] - amp[2] - 1./3. *
      std::complex<double> (0, 1) * amp[3] - 1./3. * std::complex<double> (0,
      1) * amp[4] - std::complex<double> (0, 1) * amp[9] + amp[32] + 1./3. *
      amp[33] + amp[36] - std::complex<double> (0, 1) * amp[37] + 1./3. *
      amp[38] + 1./3. * amp[54] + amp[56] + std::complex<double> (0, 1) *
      amp[57] - amp[66] + amp[68] + amp[70]);
  jamp[3] = +1./2. * (+1./3. * std::complex<double> (0, 1) * amp[0] + 1./3. *
      std::complex<double> (0, 1) * amp[1] + std::complex<double> (0, 1) *
      amp[4] + amp[5] + std::complex<double> (0, 1) * amp[7] - amp[34] - 1./3.
      * amp[35] - 1./3. * amp[36] - amp[38] + std::complex<double> (0, 1) *
      amp[39] - 1./3. * amp[62] - amp[64] - std::complex<double> (0, 1) *
      amp[65] + amp[71] - amp[73] - amp[75]);
  jamp[4] = +1./2. * (-1./3. * amp[12] - amp[13] - 1./3. * amp[14] - amp[15] -
      std::complex<double> (0, 1) * amp[17] - 1./3. * amp[22] - amp[23] - 1./3.
      * amp[24] - amp[25] + std::complex<double> (0, 1) * amp[27] -
      std::complex<double> (0, 1) * amp[61] + std::complex<double> (0, 1) *
      amp[65] + amp[72] + amp[73] + amp[74] + amp[75]);
  jamp[5] = +1./2. * (+1./3. * amp[10] + amp[11] + amp[12] + 1./3. * amp[13] +
      std::complex<double> (0, 1) * amp[19] + 1./3. * amp[42] + amp[43] + 1./3.
      * amp[44] + amp[45] - std::complex<double> (0, 1) * amp[47] -
      std::complex<double> (0, 1) * amp[57] + std::complex<double> (0, 1) *
      amp[59] - amp[67] - amp[68] - amp[69] - amp[70]);
  jamp[6] = +1./2. * (+1./3. * amp[20] + amp[21] + amp[24] + 1./3. * amp[25] -
      std::complex<double> (0, 1) * amp[29] + 1./3. * amp[30] + amp[31] + 1./3.
      * amp[34] + amp[35] + std::complex<double> (0, 1) * amp[37] +
      std::complex<double> (0, 1) * amp[53] - std::complex<double> (0, 1) *
      amp[63] - amp[67] - amp[68] - amp[69] - amp[70]);
  jamp[7] = +1./2. * (-amp[30] - 1./3. * amp[31] - 1./3. * amp[32] - amp[33] -
      std::complex<double> (0, 1) * amp[39] - 1./3. * amp[40] - amp[41] -
      amp[42] - 1./3. * amp[43] + std::complex<double> (0, 1) * amp[49] -
      std::complex<double> (0, 1) * amp[51] + std::complex<double> (0, 1) *
      amp[55] + amp[72] + amp[73] + amp[74] + amp[75]);
  jamp[8] = +1./2. * (-std::complex<double> (0, 1) * amp[3] + amp[5] -
      std::complex<double> (0, 1) * amp[6] - 1./3. * std::complex<double> (0,
      1) * amp[8] - 1./3. * std::complex<double> (0, 1) * amp[9] - amp[20] -
      1./3. * amp[21] - amp[26] - std::complex<double> (0, 1) * amp[27] - 1./3.
      * amp[28] - amp[50] + std::complex<double> (0, 1) * amp[51] - 1./3. *
      amp[52] + amp[71] - amp[73] - amp[75]);
  jamp[9] = +1./2. * (+std::complex<double> (0, 1) * amp[1] - amp[2] + 1./3. *
      std::complex<double> (0, 1) * amp[6] + 1./3. * std::complex<double> (0,
      1) * amp[7] + std::complex<double> (0, 1) * amp[8] + amp[40] + 1./3. *
      amp[41] + amp[46] + std::complex<double> (0, 1) * amp[47] + 1./3. *
      amp[48] + 1./3. * amp[50] + amp[52] - std::complex<double> (0, 1) *
      amp[53] - amp[66] + amp[68] + amp[70]);
  jamp[10] = +1./2. * (+std::complex<double> (0, 1) * amp[0] + amp[2] + 1./3. *
      std::complex<double> (0, 1) * amp[3] + 1./3. * std::complex<double> (0,
      1) * amp[4] + std::complex<double> (0, 1) * amp[9] + amp[22] + 1./3. *
      amp[23] + 1./3. * amp[26] + amp[28] + std::complex<double> (0, 1) *
      amp[29] + amp[58] - std::complex<double> (0, 1) * amp[59] + 1./3. *
      amp[60] + amp[66] + amp[67] + amp[69]);
  jamp[11] = +1./2. * (-1./3. * std::complex<double> (0, 1) * amp[0] - 1./3. *
      std::complex<double> (0, 1) * amp[1] - std::complex<double> (0, 1) *
      amp[4] - amp[5] - std::complex<double> (0, 1) * amp[7] - amp[44] - 1./3.
      * amp[45] - 1./3. * amp[46] - amp[48] - std::complex<double> (0, 1) *
      amp[49] - 1./3. * amp[58] - amp[60] + std::complex<double> (0, 1) *
      amp[61] - amp[71] - amp[72] - amp[74]);

  jamp[12] = 0.;
  jamp[13] = 0.;
  // Sum and square the color flows to get the matrix element
  /*double matrix_orig = 0;
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.;
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j];
    matrix_orig = matrix_orig + real(ztemp * conj(jamp[i]))/denom[i];
  }*/

  // Store the leading color flows for choice of color
  // for(i=0;i < ncolor; i++)
  // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

  double matrix = 0;
  std::complex<double> new_basis[ncol];

  for(i = 0; i < ncol; i++ )
  {
    new_basis[i] = 0.;
    for(j = 0; j < ncol; j++ )
      new_basis[i] = new_basis[i] + fromTrtoDiag[i][j] * jamp[j];
  }

  matrix = 0.;
  for(i = 0; i < ncol; i++ )
  {
    ztemp = 0.;
    for(j = 0; j < ncol; j++ )
    {
      ztemp = ztemp + (multfact[i][j]) * new_basis[j];
    }
    jamp2[0][i] = real(ztemp * conj(new_basis[i]));
    matrix = matrix + real(ztemp * conj(new_basis[i]));
  }

  return matrix;
}
