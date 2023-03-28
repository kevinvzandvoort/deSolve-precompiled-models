#include <R.h>

#ifdef __cplusplus
  #include "Rcpp.h"
  // the C++ compiler will see this causing all declarations in the block to have C linkage,
  // the C compiler won't see it (skipped by the preprocessor)
  extern "C" {
#endif

    void SIRcpp_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip);
    void initmod(void (* odeparms)(int *, double *));

#ifdef __cplusplus
  }
#endif

#include <Rinternals.h>

//parms as passed to deSolve will be written to this variable as SEXP
SEXP parms;
