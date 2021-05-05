//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 3.1.0, 2021-03-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "qq_process.h"
#include "qq_helamps.h"
#include "Parameters_sm.h"

using namespace qq_MG5_sm;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: u u~ > t t~ t t~ WEIGHTED<=4 @1

//--------------------------------------------------------------------------
// Initialize process.

void qq_process::initProc()
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

void qq_process::sigmaKin()
{
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
  /*bool firsttime = true;
  if (firsttime)
  {
    pars->printDependentParameters();
    pars->printDependentCouplings();
    firsttime = false;
  }*/

  // Reset color flows
  for(int i = 0; i < 6; i++ )
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 64;
  bool goodhel[ncomb] = {ncomb && false};
  int ntry = 0, sum_hel = 0, ngood = 0;
  int igood[ncomb];
  int jhel;
  std::complex<double> * * wfs;
  double t[nprocesses];
  // Helicities for the process
  const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
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
  const int denominators[nprocesses] = {144};

  ntry = ntry + 1;

  // Reset the matrix elements
  // for(int i = 0; i < nprocesses; i++){
  // matrix_element[i] = 0.;
  // }

  for(int i = 0; i < 6; i++ )
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
        t[0] = matrix_1_uux_ttxttx();

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
  // t[0]=matrix_1_uux_ttxttx();
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

void qq_process::calculate_wavefunctions(const int perm[], const int hel[])
{
  
  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]);
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]);
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]);
  FFV1P0_3(w[0], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[6]);
  FFV1P0_3(w[3], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[7]);
  FFV1_1(w[4], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[8]);
  FFV1_2(w[5], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[9]);
  FFV1P0_3(w[5], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[10]);
  FFV1P0_3(w[5], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[11]);
  FFV1_2(w[3], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[12]);
  FFV1P0_3(w[3], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[13]);
  FFV1_1(w[2], w[6], pars->GC_11, pars->mdl_MT, pars->mdl_WT, w[14]);
  FFV1_2(w[0], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[15]);
  FFV1_2(w[0], w[10], pars->GC_11, pars->ZERO, pars->ZERO, w[16]);
  FFV1_2(w[0], w[11], pars->GC_11, pars->ZERO, pars->ZERO, w[17]);
  FFV1_2(w[0], w[13], pars->GC_11, pars->ZERO, pars->ZERO, w[18]);

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
  FFV1_0(w[15], w[1], w[10], pars->GC_11, amp[10]);
  FFV1_0(w[16], w[1], w[7], pars->GC_11, amp[11]);
  FFV1_0(w[17], w[1], w[13], pars->GC_11, amp[12]);
  FFV1_0(w[18], w[1], w[11], pars->GC_11, amp[13]);

}
double qq_process::matrix_1_uux_ttxttx()
{
  int i, j;
  // Local variables
  const int ncolor = 6;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  const double denom[ncolor] = {1, 1, 1, 1, 1, 1};
  const double cf_new[ncolor][ncolor] = {{3. * sqrt(3.), sqrt(3.),
      sqrt(3.), 1./sqrt(3.), 1./sqrt(3.), sqrt(3.)}, {0., 2. * sqrt(6.), 0, 2.
      * sqrt(2./3), 2. * sqrt(2./3), 0}, {0., 0, 0, 2. * sqrt(2./3), 2. *
      sqrt(2./3), 2. * sqrt(6.)}, {0., 0, 2. * sqrt(6.), 2. * sqrt(2./3), 2. *
      sqrt(2./3.), 0.}, {0., 0, 0., 2. * sqrt(5./3), 2. * sqrt(5./3.), 0.},
      {0., 0, 0., 2. * sqrt(3.), -2. * sqrt(3.), 0.}};
  // Calculate color flows
  jamp[0] = +1./4. * (+1./9. * amp[0] + 1./9. * amp[1] + 1./3. * amp[3] + 1./3.
      * amp[4] + 1./3. * amp[6] + 1./3. * amp[7] + 1./9. * amp[8] + 1./9. *
      amp[9] + 1./9. * amp[10] + 1./9. * amp[11]);
  jamp[1] = +1./4. * (-1./3. * amp[0] - 1./3. * amp[1] - 1./9. * amp[3] - 1./9.
      * amp[4] - 1./9. * amp[6] - 1./9. * amp[7] - 1./3. * amp[8] - 1./3. *
      amp[9] - 1./9. * amp[12] - 1./9. * amp[13]);
  jamp[2] = +1./4. * (-amp[3] - std::complex<double> (0, 1) * amp[5] - amp[6] -
      1./3. * amp[8] - 1./3. * amp[9] - 1./3. * amp[10] - 1./3. * amp[11] -
      amp[12]);
  jamp[3] = +1./4. * (+amp[0] - std::complex<double> (0, 1) * amp[2] + 1./3. *
      amp[3] + 1./3. * amp[4] + amp[9] + amp[11] + 1./3. * amp[12] + 1./3. *
      amp[13]);
  jamp[4] = +1./4. * (+amp[1] + std::complex<double> (0, 1) * amp[2] + 1./3. *
      amp[6] + 1./3. * amp[7] + amp[8] + amp[10] + 1./3. * amp[12] + 1./3. *
      amp[13]);
  jamp[5] = +1./4. * (-1./3. * amp[0] - 1./3. * amp[1] - amp[4] +
      std::complex<double> (0, 1) * amp[5] - amp[7] - 1./3. * amp[10] - 1./3. *
      amp[11] - amp[13]);

  // Sum and square the color flows to get the matrix element
  // double matrix = 0;
  // for(i=0;i < ncolor; i++){
  // ztemp = 0.;
  // for(j = 0; j < ncolor; j++)
  // ztemp = ztemp + cf[i][j]*jamp[j];
  // matrix = matrix+real(ztemp*conj(jamp[i]))/denom[i];
  // }

  // Store the leading color flows for choice of color
  // for(i=0;i < ncolor; i++)
  // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

  double matrix = 0;
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.;
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf_new[i][j] * jamp[j];
    jamp2[0][i] = real(ztemp * conj(ztemp))/denom[i];
    matrix = matrix + real(ztemp * conj(ztemp))/denom[i];
  }

  return matrix;
}
