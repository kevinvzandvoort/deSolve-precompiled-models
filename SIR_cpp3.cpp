#include <R.h>
#include "Rcpp.h"
#include <Rinternals.h>

SEXP attribute_hidden get_deSolve_gparms_Rcpp() {
  static SEXP(*fun)() = NULL;
  if (fun == NULL)
      fun = (SEXP(*)()) R_GetCCallable("deSolve","get_deSolve_gparms");
  return fun();
}

extern "C" {
  void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip);
  void initmod(void (* odeparms)(int *, double *));
}

//parms as passed to deSolve will be written to this variable as SEXP
Rcpp::List parameters;

//initialize states
double Sus, Inf, Rec;

//initialize differentials
double dS, dI, dR;

//initialize parameters as global vars
double contact, recovery, mrate;

void initmod(void (* odeparms)(int *, double *)) {
  //get parms argument passed to deSolve as SEXP object
  SEXP parms = get_deSolve_gparms_Rcpp();

  try {
    //parse parameters passed to deSolve as Rcpp::List
    parameters = Rcpp::as<Rcpp::List>(parms);

    //initialize parameters
    contact = parameters["beta"];
    recovery = parameters["gamma"];
    mrate = parameters["mrate"];

  } catch(std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}

void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
  if (ip[0] <1) error("nout should be at least 1");

  //initialize state values
  Sus = y[0];
  Inf = y[1];
  Rec = y[2];

  dS = -contact*Sus*Inf +mrate*(Inf+Rec);
  dI = contact*Sus*Inf -recovery*Inf -mrate*Inf;
  dR = recovery*Inf -mrate*Rec;

  //Rcpp::Rcout << "S: " << S << "; I: " << I << "; R: " << R << "; beta: " << beta << "; gamma: " << gamma << "; dS: " << dS << "; dI: " << dI << "; dR: " << dR << std::endl;

  //Return
  ydot[0] = dS;
  ydot[1] = dI;
  ydot[2] = dR;

  yout[0] = dS + dI + dR;
}
