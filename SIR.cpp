#include <R.h>
static double parms[2];

#define beta parms[0]
#define gamma parms[1]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
  int N=2;
  odeparms(&N, parms);
}

void SIRcpp_derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {

  //initialize variables
  double S;
  double I;
  double R;

  S = y[0];
  I = y[1];
  R = y[2];

  double dS;
  double dI;
  double dR;

  if (ip[0] <1) error("nout should be at least 1");

  //ODEs
  dS = -beta*S*I;
  dI = beta*S*I -gamma*I;
  dR = gamma*I;

  //Return
  ydot[0] = dS;
  ydot[1] = dI;
  ydot[2] = dR;

  yout[0] = S+I+R;
}
