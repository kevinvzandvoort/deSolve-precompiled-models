library("Rcpp")
library("deSolve")

setwd("~/workspace/espicc_model")
dir.create("build")

Rcpp::sourceCpp(
  "espicc_model.cpp",
  cacheDir = "build", rebuild = T, verbose = T
)

(dlls=getLoadedDLLs())
is.loaded("sourceCpp_1_SIRcpp", "sourceCpp_9")

parameters = c(
  "beta" = 0.00005,
  "gamma" = 1/3
)

times = seq(1,10,0.1)

state = c(
  "S" = 100000,
  "I" = 1,
  "R" = 0
)

deSolve::lsoda(y = state, times = times, func = "sourceCpp_1_SIRcpp", parms = parameters, dllname = "sourceCpp_9", nout=1)




sourceCpp(
  "espicc_model.cpp",
  cacheDir = "build", rebuild = T, verbose = T
)

dyn.unload('/home/lsh1604011/workspace/espicc_model/build/sourceCpp-x86_64-pc-linux-gnu-1.0.6/sourcecpp_89f53fd7bc8b/sourceCpp_2.so')
`.sourceCpp_1_DLLInfo` <- dyn.load('/home/lsh1604011/workspace/espicc_model/build/sourceCpp-x86_64-pc-linux-gnu-1.0.6/sourcecpp_89f53fd7bc8b/sourceCpp_2.so')
SIRcpp2 <- Rcpp:::sourceCppFunction(function(t, state, parameters) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_SIRcpp')

(dlls=getLoadedDLLs())
getDLLRegisteredRoutines(dlls$sourceCpp_5)
is.loaded("SIRcpp")
is.loaded("sourceCpp_1_SIRcpp", "sourceCpp_1_SIRcpp")

system('g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG   -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -I"/home/lsh1604011/workspace/espicc_model"    -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-8T8CYO/r-base-4.0.3=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c espicc_model.cpp -o espicc_model.o')
system('g++ -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o espicc_model.so espicc_model.o -L/usr/lib/R/lib -lR')
x = dyn.load('espicc_model.so')
is.loaded("SIRcpp")
is.loaded("espicc_model_SIRcpp")
is.loaded("sourceCpp_1_SIRcpp")


system("R CMD SHLIB mymod.c")
dyn.load("mymod.so")

system("R CMD SHLIB mymod2.cpp")
dyn.load("mymod2.so")

(dlls = getLoadedDLLs())
#getDLLRegisteredRoutines(dlls$espicc_model)
getDLLRegisteredRoutines(dlls$mymod)

(x=getLoadedDLLs())
is.loaded("derivs", "mymod")
is.loaded("derivs", "mymod2")
is.loaded("derivs", "espicc_model")
is.loaded("SIRcpp")

parms <- c(k1 = 0.04, k2 = 1e4, k3=3e7)
Y <- c(y1 = 1.0, y2 = 0.0, y3 = 0.0)
times <- c(0, 0.4*10^(0:11) )

out <- ode(Y, times, func = "derivs", parms = parms,
           jacfunc = "jac", dllname = "mymod",
           initfunc = "initmod", nout = 1, outnames = "Sum")

out <- ode(Y, times, func = "derivs", parms = parms, dllname = "mymod", nout = 1)

SIR = function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dS = -beta*S*I
    dI = beta*S*I -gamma*I
    dR = gamma*I
    
    return(list(c(dS, dI, dR)))
  })
}

deSolve::lsoda(state, times, SIR, parameters)



system('R CMD SHLIB -o "espicc_model.so" "espicc_model.cpp" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -I"/home/lsh1604011/workspace/espicc_model" -L/usr/lib/R/lib -lR')



deSolve::ode(y = state, times = times, func = "SIRcpp", parms = parameters, dllname = "espicc_model", nout=1)

( SIR(1, state, parameters) )
( SIRcpp(1, state, parameters) )

is.loaded("sourceCpp_1_SIRcpp", "sourceCpp_8")
deSolve::ode(y = state, times = times, func = "sourceCpp_1_SIRcpp", parms = parameters, dllname = "sourceCpp_8", nout=1)
