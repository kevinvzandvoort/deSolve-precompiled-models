#include "SIR3.h"

void SIRcpp_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
  if (ip[0] <1) error("nout should be at least 1");

  //declare variables
  double beta;
  double gamma;

  //initialize state values
  double S = y[0];
  double I = y[1];
  double R = y[2];

  //initialize differentials
  double dS;
  double dI;
  double dR;

  try {
    //parse parameters passed to deSolve as Rcpp::List
    Rcpp::List parameters(parms);

    //elements passed as parameters can now be selected by their name
    //assign values to respective variables
    beta = parameters["beta"];
    gamma = parameters["gamma"];
  } catch(std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  dS = -beta*S*I;
  dI = beta*S*I -gamma*I;
  dR = gamma*I;

  //Return
  ydot[0] = dS;
  ydot[1] = dI;
  ydot[2] = dR;

  yout[0] = dS + dI + dR;
}
