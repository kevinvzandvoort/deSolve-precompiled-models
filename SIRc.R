library("deSolve")
library("Rcpp")

setwd("~/workspace/espicc_model")
dir.create("build")

#' Use a very simple SIR model
N = 100000

dur_infectiousness = 10 #days
r0 = 2.5

times = seq(1, 10 * 365, 0.1)

parameters = list(
  "beta" = r0/( N/dur_infectiousness ),
  "gamma" = 1/dur_infectiousness
)

state = c(
  "S" = N-1,
  "I" = 1,
  "R" = 0
)

SIR_r = function(t, state, pars){
  with( as.list(c(state, pars)), {
    
    dS = -beta*S*I
    dI = beta*S*I -gamma*I
    dR = gamma*I
      
    return(list(c(dS, dI, dR)))
  })
}
SIR_r_out = deSolve::ode(state, times, SIR_r, parameters)

SIR_r2 = function(t, state, pars){
    
  dS = -pars[["beta"]]*state["S"]*state["I"]
  dI = pars[["beta"]]*state["S"]*state["I"] -pars[["gamma"]]*state["I"]
  dR = pars[["gamma"]]*state["I"]
    
  return(list(c(dS, dI, dR)))
}
SIR_r2_out = deSolve::ode(state, times, SIR_r2, parameters)

SIR_r3 = function(t, state, pars){
  beta = pars[["beta"]]
  gamma = pars[["gamma"]]
  
  S = state["S"]
  I = state["I"]
  R = state["R"]
  
  dS = -beta*S*I
  dI = beta*S*I -gamma*I
  dR = gamma*I
  
  return(list(c(dS, dI, dR)))
}
SIR_r3_out = deSolve::ode(state, times, SIR_r3, parameters)

SIR_r4 = function(t, state, pars){
  beta = pars["beta"]
  gamma = pars["gamma"]
  
  S = state["S"]
  I = state["I"]
  R = state["R"]
  
  dS = -beta*S*I
  dI = beta*S*I -gamma*I
  dR = gamma*I
  
  return(list(c(dS, dI, dR)))
}
parameters_vec = unlist(parameters)
SIR_r4_out = deSolve::ode(state, times, SIR_r4, parameters_vec)

SIR_c = '
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
'
#write to file
tmpfile <- file("SIR_c.c")
cat(SIR_c, file = tmpfile);
close(tmpfile)

#compile model
system("R CMD SHLIB SIR_c.c")

#load in R
dyn.load("SIR_c.so")

#check that functions are available
is.loaded("initmod", "SIR_c")
is.loaded("derivs", "SIR_c")

#run model
SIR_c_out = deSolve::ode(state, times, func = "derivs", initfunc = "initmod", dllname = "SIR_c", nout = 1, outnames = "Sum", parms = parameters_vec)

SIR_c2 = '
#include <R.h>
static double parms[2];

//assign variables by index
#define beta parms[0]
#define gamma parms[1]

double S, I, R;
double dS, dI, dR;

/* initializer */
void initmod(void (* odeparms)(int *, double *)){
  int N=2;
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
    dS = -beta*S*I;
    dI = beta*S*I -gamma*I;
    dR = gamma*I;
    
    //return values
    ydot[0] = dS;
    ydot[1] = dI;
    ydot[2] = dR;
    
    yout[0] = S+I+R;
}
'
#write to file
tmpfile <- file("SIR_c2.c")
cat(SIR_c2, file = tmpfile);
close(tmpfile)

#compile model
system("R CMD SHLIB SIR_c2.c")

#load in R
dyn.load("SIR_c2.so")

#check that functions are available
is.loaded("initmod", "SIR_c2")
is.loaded("derivs", "SIR_c2")

#run model
SIR_c2_out = deSolve::ode(state, times, func = "derivs", initfunc = "initmod", dllname = "SIR_c2", nout = 1, outnames = "Sum", parms = parameters_vec)

SIR_cpp = '
Rcpp::List SIR_cpp(double t, Rcpp::NumericVector state, Rcpp::List pars){
  
  double beta = pars["beta"];
  double gamma = pars["gamma"];
  
  double S = state["S"];
  double I = state["I"];
  double R = state["R"];
  
  double dS, dI, dR;
  
  dS = -beta*S*I;
  dI = beta*S*I -gamma*I;
  dR = gamma*I;
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::NumericVector::create(
      Rcpp::Named("dS", dS),
      Rcpp::Named("dI", dI),
      Rcpp::Named("dR", dR)
    )
  );
  
  return result;
}
'
SIR_cpp_func = Rcpp::cppFunction(SIR_cpp, rebuild = T);
SIR_cpp_out = deSolve::ode(state, times, SIR_cpp_func, parameters)

#' Header file
SIR_h = '
#include <R.h>
#ifdef __cplusplus
  #include "Rcpp.h"
  // the C++ compiler will see this causing all declarations in the block to have C linkage,
  //  (skipped by the preprocessor for C compiler)
  extern "C" {
#endif
  void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip);
  void initmod(void (* odeparms)(int *, double *));
#ifdef __cplusplus
}
#endif
#include <Rinternals.h>
//parms as passed to deSolve will be written to this variable as SEXP
SEXP parms;
'
#write to file
tmpfile <- file("SIR_cpp.h")
cat(SIR_h, file = tmpfile);
close(tmpfile)

#' C file to read parameters
SIR_cpp_c = '
#include "SIR_cpp.h"
void initmod(void (* odeparms)(int *, double *))
{
  DL_FUNC get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");

  //this returns parms as SEXP and writes the SEXP to pars defined in the header file
  parms = get_deSolve_gparms();
}
'
#write to file
tmpfile <- file("SIR_cpp.c")
cat(SIR_cpp_c, file = tmpfile);
close(tmpfile)

#' Actual model
SIR_cpp2 = '
#include "SIR_cpp.h"
void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
  if (ip[0] <1) error("nout should be at least 1");
  
  //declare variables
  double beta, gamma;
  
  //initialize state values
  double S = y[0];
  double I = y[1];
  double R = y[2];
  
  //initialize differentials
  double dS, dI, dR;
  
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
'
#write to file
tmpfile <- file("SIR_cpp2.cpp")
cat(SIR_cpp2, file = tmpfile);
close(tmpfile)

#compile
system('gcc -c -o SIR_cpp_readpar.o SIR_cpp.c -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -fPIC');
system('g++ -std=gnu++11 -c -o SIR_cpp2.o SIR_cpp2.cpp -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -fPIC')
#link
system('g++ -std=gnu++11 -o SIR_cpp2.so SIR_cpp_readpar.o SIR_cpp2.o -shared -L/usr/lib/R/lib -lR')

dyn.load("SIR_cpp2.so")
is.loaded("derivs", "SIR_cpp2")
is.loaded("initmod", "SIR_cpp2")

SIR_cpp2_out = ode(state, times, func = "derivs", parms = parameters, dllname = "SIR_cpp2", nout = 1, initfunc = "initmod")


(model_bmark = rbenchmark::benchmark(
  deSolve::ode(state, times, SIR_r, parameters),
  deSolve::ode(state, times, SIR_r2, parameters),
  deSolve::ode(state, times, SIR_r3, parameters),
  deSolve::ode(state, times, SIR_r4, parameters_vec),
  deSolve::ode(state, times, func = "derivs", initfunc = "initmod", dllname = "SIR_c", nout = 1, outnames = "Sum", parms = parameters_vec),
  deSolve::ode(state, times, func = "derivs", initfunc = "initmod", dllname = "SIR_c2", nout = 1, outnames = "Sum", parms = parameters_vec),
  deSolve::ode(state, times, SIR_cpp_func, parameters),
  ode(state, times, func = "derivs", parms = parameters, dllname = "SIR_cpp2", nout = 1, initfunc = "initmod"),
  replications=rep(20, 1),
  columns=c('test', 'elapsed', 'relative', 'replications')
))

system('gcc -c -o SIR3_init.o SIR3.c -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -fPIC')
system('g++ -c -o SIR3_derivs.o SIR3.cpp -L/usr/lib/R/lib -lR -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions,-z,relro -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -fPIC')

# link them together, using the C++ library libmycpplib:
system('g++ -oSIR3.so SIR3_init.o SIR3_derivs.o -L/usr/lib/R/lib -lR -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions,-z,relro -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -fPIC')
#g++ -oprogram main.o module.o -lmycpplib

dyn.load("SIR3.so")
is.loaded("SIRcpp_derivs", "SIR3")
is.loaded("initmod", "SIR3")

(out <- ode(state, times, func = "SIRcpp_derivs", parms = parameters, dllname = "SIR3", nout = 1, initfunc = "initmod"))
