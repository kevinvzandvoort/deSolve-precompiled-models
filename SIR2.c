#include <R.h>
//#include "Rcpp.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

//extern "C" void SIRcpp_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip);
//extern "C" void initmod(void (* odeparms)(int *, double *));

//static double parms[2];
//#define beta parms[0]
//#define gamma parms[1]

//std::vector<double> parms;
//Rcpp::List parms;

//double test;

/* initializer */

static SEXP *cparms;
void initmod(void (* odeparms)(int *, double *))
{
  DL_FUNC get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
	SEXP cparms = get_deSolve_gparms();
}

/*void initmod(void (* odeparms)(int *, double *))
{
  int N=2;
  odeparms(&N, parms);

  test = 2;
}*/

void SIRcpp_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {

  //Rcpp::Rcout << "neq: " << neq << "; t: " << t << "; y[0]: " << y[0] << "; ydot[0]: " << ydot[0] << "; yout[0]" << yout[0] << "; ip[0]: " << ip[0] << std::endl;

  //initialize variables
  double S;
  double I;
  double R;

  //double beta = parms[0];

  //double beta = parms["beta"];
  //Rcpp::List parameters = parms;
  //int size = parameters.size();
  //double gamma = parms["gamma"];

  //Rcpp::Rcout << "Size of list: " << size << "; Size of parms: " << parms.size() << std::endl;

  //Rcpp::Rcout << "Size of parms: " << parms.size() << std::endl;

  S = y[0];
  I = y[1];
  R = y[2];

  double dS;
  double dI;
  double dR;

  if (ip[0] <1) error("nout should be at least 1");

  //ODEs
  //dS = -beta*S*I;
  //dI = beta*S*I -gamma*I;
  //dR = gamma*I + test + parms[1];

  double beta = 0.0005;

  dS = -beta*S*I;
  dI = beta*S*I -0.25*I;
  dR = 0.25*I;

  //Return
  ydot[0] = dS;//dS;
  ydot[1] = dI;//dI;
  ydot[2] = dR;//dR;

  yout[0] = S+I+R;
}
