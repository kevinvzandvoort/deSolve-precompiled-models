#include "Rcpp.h"
#include <R.h>

// [[Rcpp::export]]
Rcpp::List SIRcpp(int t, Rcpp::NumericVector state, Rcpp::NumericVector parameters) {

  double beta = parameters["beta"];
  double gamma = parameters["gamma"];

  double S = state["S"];
  double I = state["I"];
  double R = state["R"];

  Rcpp::List result = Rcpp::List::create(
    Rcpp::NumericVector::create(
      Rcpp::Named("S") = S,
      Rcpp::Named("I") = I,
      Rcpp::Named("R") = R
    )
  );
  return result;
}
