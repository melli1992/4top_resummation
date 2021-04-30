#include <iostream>
#include <iomanip>

#include "qq_process.h"

int main(int argc, char** argv){

  // Create a process object
  qq_process process;

  // Read param_card and set parameters
  process.initProc("Cards/param_card.dat");
  
  //double energy = 1500;
  //double weight;

  // Get phase space point
  //vector<double*> p = get_momenta(process.ninitial, energy,
//				 process.getMasses(), weight);
  //for(int k = 0; k < 1; k++)
  //{
    //p = get_momenta(process.ninitial, energy,process.getMasses(), weight);
    // Set momenta for this event
    //process.setMomenta(p);

    // Evaluate matrix element
    //process.sigmaKin();


    /*cout << "Momenta:" << endl;
    for(int i=0;i < process.nexternal; i++)
      cout << setw(4) << i+1
  	 << setiosflags(ios::scientific) << setw(14) << p[i][0]
  	 << setiosflags(ios::scientific) << setw(14) << p[i][1]
  	 << setiosflags(ios::scientific) << setw(14) << p[i][2]
  	 << setiosflags(ios::scientific) << setw(14) << p[i][3] << endl;
    cout << " -----------------------------------------------------------------------------" << endl;

    // Display matrix elements
    double sum = 0.;
    for(int i=0; i<process.ncol;i++)
    {sum+=matrix_elements[i];
      cout << " Matrix element = "
  //	 << setiosflags(ios::fixed) << setprecision(17)
  	 << matrix_elements[i]
  	 << " GeV^" << -(2*process.nexternal-8) << endl;}
    cout << "SUM " << sum << endl;
    cout << " -----------------------------------------------------------------------------" << endl;
  }*/

  vector<double*> p_new(1, new double[4]);
  p_new[0][0] = 7.500000E2;
  p_new[0][1] = 0;
  p_new[0][2] = 0;
  p_new[0][3] = 7.500000E2;
  p_new.push_back(new double[4]);
  p_new[1][0] = 7.500000E2;
  p_new[1][1] = 0;
  p_new[1][2] = 0;
  p_new[1][3] = -7.500000E2;
  p_new.push_back(new double[4]);
  p_new[2][0] = 2.311780E2;
  p_new[2][1] = -7.905776E1;
  p_new[2][2] = 9.276470E1;
  p_new[2][3] = 9.305303E1;
  p_new.push_back(new double[4]);
  p_new[3][0] = 4.886960E2;
  p_new[3][1] = -4.139796E1;
  p_new[3][2] = 1.093149E2;
  p_new[3][3] = 4.418498E2;
  p_new.push_back(new double[4]);
  p_new[4][0] = 4.867211E2;
  p_new[4][1] = 2.396366E2;
  p_new[4][2] = -6.710502E1;
  p_new[4][3] = -3.808407E2;
  p_new.push_back(new double[4]);
  p_new[5][0] = 2.934048E2;
  p_new[5][1] = -1.191809E2;
  p_new[5][2] = -1.349746E2;
  p_new[5][3] = -1.540621E2;
  process.setMomenta(p_new);
  // Evaluate matrix element
  process.sigmaKin();
  const double* matrix_elements = process.getMatrixElements();
  double sum = 0.;
  for(int i=0; i<process.ncol;i++)
  {sum+=matrix_elements[i];
    cout << " Matrix element = "
//	 << setiosflags(ios::fixed) << setprecision(17)
   << matrix_elements[i]
   << " GeV^" << -(2*process.nexternal-8) << endl;}
  cout << "SUM " << sum << endl;
  cout << " -----------------------------------------------------------------------------" << endl;
/*  7.500000e+02  0.000000e+00  0.000000e+00  7.500000e+02
  2  7.500000e+02  0.000000e+00  0.000000e+00 -7.500000e+02
  3  2.311780e+02 -7.905776e+01  9.276470e+01  9.305303e+01
  4  4.886960e+02 -4.139796e+01  1.093149e+02  4.418498e+02
  5  4.867211e+02  2.396366e+02 -6.710502e+01 -3.808407e+02
  6  2.934048e+02 -1.191809e+02 -1.349746e+02 -1.540621e+02*/


}
