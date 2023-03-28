#include "SIR3.h"

void initmod(void (* odeparms)(int *, double *))
{
  DL_FUNC get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");

  //this returns parms as SEXP and writes the SEXP to pars defined in the header file
  parms = get_deSolve_gparms();
}
