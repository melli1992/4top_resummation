//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 3.1.0, 2021-03-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "Parameters_sm.h"
#include "parameters.h"

// Initialize static instance
Parameters_sm * Parameters_sm::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_sm * Parameters_sm::getInstance()
{
  if (instance == 0)
    instance = new Parameters_sm(); 

  return instance; 
}

void Parameters_sm::setIndependentParameters()
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  mdl_WH = 6.382339e-03; 
  mdl_WW = 2.047600e+00; 
  mdl_WZ = 2.441404e+00; 
  mdl_WT = 1.491500e+00; 
  mdl_ymtau = 1.777000e+00; 
  mdl_ymt = 1.730000e+02; 
  mdl_ymb = 4.700000e+00; 
  aS = alphas_muR; 
  mdl_Gf = 1.166390e-05; 
  aEWM1 = 1.325070e+02; 
  mdl_MH = 1.250000e+02; 
  mdl_MZ = 9.118800e+01; 
  mdl_MTA = mdl_ymtau; 
  mdl_MT = mt; 
  mdl_MB = mdl_ymb; 
  mdl_CKM3x3 = 1.; 
  mdl_conjg__CKM1x1 = 1.; 
  mdl_conjg__CKM3x3 = 1.; 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_MZ__exp__2 = ((mdl_MZ) * (mdl_MZ)); 
  mdl_MZ__exp__4 = ((mdl_MZ) * (mdl_MZ) * (mdl_MZ) * (mdl_MZ)); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_MH__exp__2 = ((mdl_MH) * (mdl_MH)); 
  mdl_aEW = 1./aEWM1; 
  mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = sqrt(mdl_aEW); 
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
  mdl_MW__exp__2 = ((mdl_MW) * (mdl_MW)); 
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
  mdl_cw = sqrt(1. - mdl_sw2); 
  mdl_sqrt__sw2 = sqrt(mdl_sw2); 
  mdl_sw = mdl_sqrt__sw2; 
  mdl_g1 = mdl_ee/mdl_cw; 
  mdl_gw = mdl_ee/mdl_sw; 
  mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee; 
  mdl_vev__exp__2 = ((mdl_vev) * (mdl_vev)); 
  mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2); 
  mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_vev; 
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev; 
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev; 
  mdl_muH = sqrt(mdl_lam * mdl_vev__exp__2); 
  mdl_I1x33 = mdl_yb * mdl_conjg__CKM3x3; 
  mdl_I2x33 = mdl_yt * mdl_conjg__CKM3x3; 
  mdl_I3x33 = mdl_CKM3x3 * mdl_yt; 
  mdl_I4x33 = mdl_CKM3x3 * mdl_yb; 
  mdl_ee__exp__2 = ((mdl_ee) * (mdl_ee)); 
  mdl_sw__exp__2 = ((mdl_sw) * (mdl_sw)); 
  mdl_cw__exp__2 = ((mdl_cw) * (mdl_cw)); 
}
void Parameters_sm::setIndependentCouplings()
{

}
void Parameters_sm::setDependentParameters()
{
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = ((G) * (G)); 
}
void Parameters_sm::setDependentCouplings()
{
  GC_10 = -G; 
  GC_11 = mdl_complexi * G; 
  GC_12 = mdl_complexi * mdl_G__exp__2; 
}

// Routines for printing out parameters
void Parameters_sm::printIndependentParameters()
{
  cout <<  "sm model parameters independent of event kinematics:" << endl; 
  cout << setw(20) <<  "mdl_WH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WH << endl;
  cout << setw(20) <<  "mdl_WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WW << endl;
  cout << setw(20) <<  "mdl_WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WZ << endl;
  cout << setw(20) <<  "mdl_WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WT << endl;
  cout << setw(20) <<  "mdl_ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymtau << endl;
  cout << setw(20) <<  "mdl_ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymt << endl;
  cout << setw(20) <<  "mdl_ymb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymb << endl;
  cout << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  cout << setw(20) <<  "mdl_Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gf << endl;
  cout << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  cout << setw(20) <<  "mdl_MH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MH << endl;
  cout << setw(20) <<  "mdl_MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MZ << endl;
  cout << setw(20) <<  "mdl_MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MTA << endl;
  cout << setw(20) <<  "mdl_MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MT << endl;
  cout << setw(20) <<  "mdl_MB " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MB << endl;
  cout << setw(20) <<  "mdl_CKM3x3 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM3x3 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM1x1 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM1x1 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM3x3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM3x3 << endl;
  cout << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi << endl;
  cout << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2 << endl;
  cout << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4 << endl;
  cout << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2 << endl;
  cout << setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__2 << endl;
  cout << setw(20) <<  "mdl_aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_aEW << endl;
  cout << setw(20) <<  "mdl_MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MW << endl;
  cout << setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__aEW << endl;
  cout << setw(20) <<  "mdl_ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ee << endl;
  cout << setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw2 << endl;
  cout << setw(20) <<  "mdl_cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cw << endl;
  cout << setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__sw2 << endl;
  cout << setw(20) <<  "mdl_sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw << endl;
  cout << setw(20) <<  "mdl_g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_g1 << endl;
  cout << setw(20) <<  "mdl_gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gw << endl;
  cout << setw(20) <<  "mdl_vev " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_vev << endl;
  cout << setw(20) <<  "mdl_vev__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_vev__exp__2 << endl;
  cout << setw(20) <<  "mdl_lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_lam << endl;
  cout << setw(20) <<  "mdl_yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yb << endl;
  cout << setw(20) <<  "mdl_yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yt << endl;
  cout << setw(20) <<  "mdl_ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ytau << endl;
  cout << setw(20) <<  "mdl_muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_muH << endl;
  cout << setw(20) <<  "mdl_I1x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1x33 << endl;
  cout << setw(20) <<  "mdl_I2x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2x33 << endl;
  cout << setw(20) <<  "mdl_I3x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3x33 << endl;
  cout << setw(20) <<  "mdl_I4x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4x33 << endl;
  cout << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2 << endl;
  cout << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2 << endl;
}
void Parameters_sm::printIndependentCouplings()
{
  cout <<  "sm model couplings independent of event kinematics:" << endl; 

}
void Parameters_sm::printDependentParameters()
{
  cout <<  "sm model parameters dependent on event kinematics:" << endl; 
  cout << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS << endl;
  cout << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  cout << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2 << endl;
}
void Parameters_sm::printDependentCouplings()
{
  cout <<  "sm model couplings dependent on event kinematics:" << endl; 
  cout << setw(20) <<  "GC_10 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_10 << endl;
  cout << setw(20) <<  "GC_11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_11 << endl;
}


