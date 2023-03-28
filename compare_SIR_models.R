library("deSolve")
library("Rcpp")

setwd("~/workspace/espicc_model")
#dir.create("build")

#' Use a very simple SIR model
N = 100000

dur_infectiousness = 10 #days
r0 = 3.5

times = seq(1, 100 * 365, 1)

parameters = list(
  "beta" = r0/( N/(1/dur_infectiousness )),
  "gamma" = 1/dur_infectiousness,
  "mrate" = 1/(75*365)
)

state = c(
  "S" = N-1,
  "I" = 1,
  "R" = 0
)

SIR_r = function(t, state, pars){
  with( as.list(c(state, pars)), {
    
    dS = -beta*S*I +mrate*(I+R)
    dI = beta*S*I -gamma*I -mrate*I
    dR = gamma*I -mrate*R
      
    return(list(c(dS, dI, dR)))
  })
}
SIR_r_out = deSolve::lsoda(state, times, SIR_r, parameters)

SIR_r2 = function(t, state, pars){
    
  dS = -pars[["beta"]]*state["S"]*state["I"] +pars[["mrate"]]*(state["I"]+state["R"])
  dI = pars[["beta"]]*state["S"]*state["I"] -pars[["gamma"]]*state["I"] -pars[["mrate"]]*state["I"]
  dR = pars[["gamma"]]*state["I"] -pars[["mrate"]]*state["R"]
    
  return(list(c(dS, dI, dR)))
}
SIR_r2_out = deSolve::lsoda(state, times, SIR_r2, parameters)

SIR_r3 = function(t, state, pars){
  beta = pars[["beta"]]
  gamma = pars[["gamma"]]
  mrate = pars[["mrate"]]
  
  S = state["S"]
  I = state["I"]
  R = state["R"]
  
  dS = -beta*S*I +mrate*(I+R)
  dI = beta*S*I -gamma*I -mrate*I
  dR = gamma*I -mrate*R
  
  return(list(c(dS, dI, dR)))
}
SIR_r3_out = deSolve::lsoda(state, times, SIR_r3, parameters)

SIR_r4 = function(t, state, pars){
  beta = pars["beta"]
  gamma = pars["gamma"]
  mrate = pars["mrate"]
  
  S = state["S"]
  I = state["I"]
  R = state["R"]
  
  dS = -beta*S*I +mrate*(I+R)
  dI = beta*S*I -gamma*I -mrate*I
  dR = gamma*I -mrate*R
  
  return(list(c(dS, dI, dR)))
}
parameters_vec = unlist(parameters)
SIR_r4_out = deSolve::lsoda(state, times, SIR_r4, parameters_vec)

SIR_c = '
#include <R.h>
static double parms[3];

#define beta parms[0]
#define gamma parms[1]
#define mrate parms[2]

/* initializer */
void initmod(void (* odeparms)(int *, double *)){
  int N=3;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
  if(ip[0] < 1)
    error("nout should be at least 1");
    //dS
    ydot[0] = -beta*y[0]*y[1] +mrate*(y[1]+y[2]);
    //dI
    ydot[1] = beta*y[0]*y[1] -gamma*y[1] -mrate*y[1];
    //dR
    ydot[2] = gamma*y[1] -mrate*y[2];
    
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
SIR_c_out = deSolve::lsoda(state, times, func = "derivs", initfunc = "initmod", dllname = "SIR_c", nout = 1, outnames = "Sum", parms = parameters_vec)

SIR_c2 = '
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
SIR_c2_out = deSolve::lsoda(state, times, func = "derivs", initfunc = "initmod", dllname = "SIR_c2", nout = 1, outnames = "Sum", parms = parameters_vec)

SIR_cpp = '
Rcpp::List SIR_cpp(double t, Rcpp::NumericVector state, Rcpp::List pars){
  
  double beta = pars["beta"];
  double gamma = pars["gamma"];
  double mrate = pars["mrate"];
  
  double S = state["S"];
  double I = state["I"];
  double R = state["R"];
  
  double dS, dI, dR;
  
  dS = -beta*S*I +mrate*(I+R);
  dI = beta*S*I -gamma*I -mrate*I;
  dR = gamma*I -mrate*R;
  
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
SIR_cpp_out = deSolve::lsoda(state, times, SIR_cpp_func, parameters)

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
  double beta, gamma, mrate;
  
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
    mrate = parameters["mrate"];
  } catch(std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  dS = -beta*S*I +mrate*(I+R);
  dI = beta*S*I -gamma*I -mrate*I;
  dR = gamma*I -mrate*R;

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

SIR_cpp2_out = lsoda(state, times, func = "derivs", parms = parameters, dllname = "SIR_cpp2", nout = 1, initfunc = "initmod")

system('g++ -c -o SIR_cpp3.o SIR_cpp3.cpp -L/usr/lib/R/lib -lR -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions,-z,relro -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -fPIC')
system('g++ -o SIR_cpp3.so SIR_cpp3.o -L/usr/lib/R/lib -lR -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions,-z,relro -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -fPIC')

dyn.load("SIR_cpp3.so")
is.loaded("derivs", "SIR_cpp3")
is.loaded("initmod", "SIR_cpp3")

SIR_cpp3_out = lsoda(state, times, func = "derivs", parms = parameters, dllname = "SIR_cpp3", nout = 1, initfunc = "initmod")

library(data.table)
library(ggplot2)

compare_output = rbindlist(list(
  as.data.table(SIR_r_out)[, name := "R - basic"],
  as.data.table(SIR_r2_out)[, name := "R - optimised 1"],
  as.data.table(SIR_r3_out)[, name := "R - optimised 2"],
  as.data.table(SIR_r4_out)[, name := "R - optimised 3"],
  as.data.table(SIR_c_out)[, name := "C - basic"],
  as.data.table(SIR_c2_out)[, name := "C - clarified"],
  as.data.table(SIR_cpp_out)[, name := "C++ (as R-function) - basic"],
  as.data.table(SIR_cpp2_out)[, name := "C++ (as compiled-function) - basic"],
  as.data.table(SIR_cpp3_out)[, name := "C++ (as compiled-function) - optimised"]
), fill=T)
compare_output = melt(compare_output[, -c("Sum", c("4")), with=F], measure.vars=c("S","I","R"))

ggplot(
  compare_output,
  aes(x=time, colour=name, y=value, group=paste(name,variable))
)+geom_line()+facet_grid(variable~factor(
  name,
  c(
    "R - basic", "R - optimised 1", "R - optimised 2", "R - optimised 3",
    "C - basic", "C - clarified",
    "C++ (as R-function) - basic", "C++ (as compiled-function) - basic", "C++ (as compiled-function) - optimised"
  ),
  c(
    "R\nbasic", "R\noptimised 1", "R\noptimised 2", "R\noptimised 3",
    "C\nbasic", "C\nclarified",
    "C++ (R)\nbasic", "C++ (compiled)\nbasic", "C++ (compiled)\noptimised"
  )
), scales="free")+theme_classic()+theme(legend.position="bottom")+scale_x_continuous(
  breaks = seq(0,100,25) * 365, 
  labels = seq(0,100,25)
)+labs(x = "year")

ggplot(
  compare_output,
  aes(x=time, colour=name, y=value, group=paste(name,variable))
)+geom_line()+facet_grid(.~variable, scales="free")+theme_classic()+theme(legend.position="bottom")+scale_x_continuous(
  breaks = seq(0,100,25) * 365, 
  labels = seq(0,100,25)
)+labs(x = "year")

(model_bmark = rbenchmark::benchmark(
  "R - basic" = deSolve::lsoda(state, times, SIR_r, parameters),
  "R - optimised 1" = deSolve::lsoda(state, times, SIR_r2, parameters),
  "R - optimised 2" = deSolve::lsoda(state, times, SIR_r3, parameters),
  "R - optimised 3" = deSolve::lsoda(state, times, SIR_r4, parameters_vec),
  "C - basic" = deSolve::lsoda(state, times, func = "derivs", initfunc = "initmod", dllname = "SIR_c", nout = 1, outnames = "Sum", parms = parameters_vec),
  "C - clarified" = deSolve::lsoda(state, times, func = "derivs", initfunc = "initmod", dllname = "SIR_c2", nout = 1, outnames = "Sum", parms = parameters_vec),
  "C++ (as R-function) - basic" = deSolve::lsoda(state, times, SIR_cpp_func, parameters),
  "C++ (as compiled-function) - basic" = lsoda(state, times, func = "derivs", parms = parameters, dllname = "SIR_cpp2", nout = 1, initfunc = "initmod"),
  "C++ (as compiled-function) - optimised" = lsoda(state, times, func = "derivs", parms = parameters, dllname = "SIR_cpp3", nout = 1, initfunc = "initmod"),
  replications=rep(100, 1),
  columns=c('test', 'elapsed', 'relative', 'replications')
))

as.data.table(model_bmark)[order(relative)]
