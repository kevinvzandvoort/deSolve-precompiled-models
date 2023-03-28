
#include <R.h>
static double parms[3];

//assign variables by index
#define beta parms[0]
#define gamma parms[1]
#define mrate parms[2]

double S, I, R;
double dS, dI, dR;

/* initializer */
void initmod(void (* odeparms)(int *, double *)){
  int N=3;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
  if(ip[0] < 1)
    error("nout should be at least 1");
    
    S = y[0];
    I = y[1];
    R = y[2];
    
    //model is more readable
    dS = -beta*S*I +mrate*(I+R);
    dI = beta*S*I -gamma*I -mrate*I;
    dR = gamma*I -mrate*R;
    
    //return values
    ydot[0] = dS;
    ydot[1] = dI;
    ydot[2] = dR;
    
    yout[0] = S+I+R;
}
