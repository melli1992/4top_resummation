#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "parameters.h"
#include "monte_carlo.h"
#include "4top_vegas.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"

using namespace std;

/////////////////////////////////////////////////////////////////////
/// this contains a call to vegas and a functionint to write the results
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
/// display results of the MC integration
/////////////////////////////////////////
void display_results(string title, double &result, double &error, double chisq){
	cout << title << " result = " << result <<  " sigma = " << error  << ", chisq = " << chisq << endl;
}

functionint init_vegas_4top(std::string process){
	functionint integrand;

	if(process == "LO"){
			integrand.G.f =&vegas_4top_LO;
			integrand.G.dim = 10;
			integrand.xl = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
			integrand.xu = {1.,1.,1.,1.,M_PI, M_PI,M_PI, 2.*M_PI, 2.*M_PI, 2.*M_PI};
			//order is x1 x2 s12 s34 thetaCM theta12 theta34 phiCM phi12 phi34
		}
	else if(process == "LO_N"){
		integrand.G.f =&vegas_4top_LO_N;
		if(fitPDF){
			integrand.G.dim = 10;
			integrand.xl = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
			integrand.xu = {1.,1.,1.,1.,M_PI, M_PI,M_PI, 2.*M_PI, 2.*M_PI, 2.*M_PI};
			//order is N rho s12 s34 thetaCM theta12 theta34 phiCM phi12 phi34
		}
		else{cout << "This option for integration does not exist" << endl;
			 exit(0);}
		//else{
			//integrand.G.dim = 11;
			//integrand.xl = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
			//integrand.xu = {1.,1.,1.,1.,1.,1.,M_PI, M_PI,M_PI, 2.*M_PI, 2.*M_PI};
			//order is N x1, x2, rho s12 s34 thetaCM theta12 theta34 phi12 phi34
		//	}
		}
	else if(process == "resum"){
		integrand.G.f =&vegas_4top_res;
		if(fitPDF){
			integrand.G.dim = 10;
			integrand.xl = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
			integrand.xu = {1.,1.,1.,1.,M_PI, M_PI,M_PI, 2.*M_PI, 2.*M_PI, 2.*M_PI};
			//order is N rho s12 s34 thetaCM theta12 theta34 phiCM phi12 phi34
		}
		else{cout << "This option for integration does not exist" << endl;
			 exit(0);}
		}
	else{cout << "this option does not exist" << endl;
			exit(0);}
		return integrand;

}

//////////////////////////////////////////////////////////
/// call vegas, store the results in a result struc
/// handles all the declaration
/// input is the functionint to integrate
/// the parameter structure (may change)
/// the number of dimensions
/// the upper and lower bounds of the integrand
//////////////////////////////////////////////////////////
results call_vegas(functionint integrand, lumni_params params, bool verbal, bool high_prec){
	  double res(0),err(0);
	  int MAX_ITER = 20;
	  const gsl_rng_type *T;
	  gsl_rng *r; //the random number generator
	  gsl_rng_env_setup ();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);


	  integrand.G.params = &params;
	  size_t calls = 100000;

	  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (integrand.G.dim);
	  gsl_monte_vegas_init(s); //whatever, just do it
		gsl_monte_vegas_params *params_run = (gsl_monte_vegas_params*)malloc( sizeof(gsl_monte_vegas_params) );
		gsl_monte_vegas_params_get(s, params_run);
	  //params_run->alpha = 2.0;
	  gsl_monte_vegas_params_set(s, params_run);

		gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim,  10000, r, s, &res, &err); //integrate g over dim-dim region defined by lower and upper limits in arrays xl and xu, each of size dim. r is the random number generator. The result is given to res and err. This function is to prepare/warm up the grid (1000 function calls). After this the main run is called with calls/5 function calls.
		display_results("Initialization - " ,res, err);

	  int n_iter = 0, k = 0;
	  double chisq;
	  double minerr = 1.E-2;
		//params_run->alpha = 1.8;
		gsl_monte_vegas_params_set(s, params_run);
  	do
       {
       	gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls, r, s, &res, &err);
		chisq = gsl_monte_vegas_chisq (s); //gsl_monte_vegas_chisq (s) provides chisq per d.o.f.
		//gsl_monte_vegas_runval(s, &res, &sigma); //returns raw(unaveraged) values of integral result and error sigma from the most recent iteration of algorithm
		display_results("Intermediate step - " ,res, err, chisq);
		n_iter++;
		if(n_iter > MAX_ITER)
			{minerr = minerr*10;
			 n_iter = 0;
			 k+=1;}
		if(k>3){break;}

       }
       while (((isnan(res))||(fabs (chisq - 1.0) > 0.1)||(fabs (err/res) > minerr)));
      display_results("Final result - ", res, err, chisq);


	  /*
	  verbal = true;
	  if(high_prec==true)
	  {gsl_monte_vegas_params *params_run = (gsl_monte_vegas_params*)malloc( sizeof(gsl_monte_vegas_params) );
	  gsl_monte_vegas_params_get(s, params_run);
	  params_run->iterations = MAX_ITER;
	  params_run->alpha = 2.0;
	  params_run->stage = 0;
	  gsl_monte_vegas_params_set(s, params_run);

	  gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, 10000, r, s, &res, &err);
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);
	  gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls/10, r, s, &res, &err);
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);
	  if(verbal){display_results ("vegas warm-up", res, err);}
	  do{
			n_iter+=1;
			params_run->stage = 2;
	        params_run->alpha = 1.5;
	        gsl_monte_vegas_params_set(s, params_run);
			gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls, r, s, &res, &err);
			if(verbal){ cout << "round " << n_iter << ", result = " << res << " error = " << err << " chisq/dof = " << gsl_monte_vegas_chisq (s) << " err/res = " << err/res << endl;}
			if(n_iter >= MAX_ITER){break;}
		}
	  while ((isnan(res))||(fabs (err/res) > 1E-4));
	  params_run->stage = 1;
	  gsl_monte_vegas_params_set(s, params_run);

	  gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls*10, r, s, &res, &err);

	  params_run->stage = 3;
	  gsl_monte_vegas_params_set(s, params_run);
	  n_iter = 0;
	  do{
			n_iter+=1;
			gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls/5, r, s, &res, &err);
			if(verbal){ cout << "round " << n_iter << ", result = " << res << " error = " << err << " chisq/dof = " << gsl_monte_vegas_chisq (s) << " err/res = " << err/res << endl;}
			if(n_iter > 100){break;}
		}
	  while ((isnan(res))||(fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.05)||(fabs (err/res) > 1E-4));

	  if(verbal){display_results ("vegas final", res, err);}
	  }
	 if(high_prec==false)
	 {
		 do{
			n_iter+=1;
			gsl_monte_vegas_integrate (&integrand.G, &integrand.xl[0], &integrand.xu[0], integrand.G.dim, calls/5, r, s, &res, &err);
			if(verbal){ cout << "round " << n_iter << ", result = " << res << " error = " << err << " chisq/dof = " << gsl_monte_vegas_chisq (s) << " err/res = " << err/res << endl;}
			if(n_iter > 100){break;}
	  }
	  while ((isnan(res))||(fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.05)||(fabs (err/res) > 1E-4));
      }
      */

	  gsl_monte_vegas_free (s);

	  gsl_rng_free (r);
	  results result = {res,err};
	  return result;

}
