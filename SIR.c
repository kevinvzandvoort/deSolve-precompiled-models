
#include <R.h>
static double parms[2];

#define beta parms[0]
#define gamma parms[1]

/* initializer */
void initmod(void (* odeparms)(int *, double *)){
  int N=2;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
  if(ip[0] < 1)
    error("nout should be at least 1");
    //dS
    ydot[0] = -beta*y[0]*y[1];
    //dI
    ydot[1] = beta*y[0]*y[1] -gamma*y[1];
    //dR
    ydot[2] = gamma*y[1];
    
    yout[0] = y[0]+y[1]+y[2];
}
