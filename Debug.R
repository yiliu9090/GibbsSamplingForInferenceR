library(DirichletReg)
library(jsonlite)
library(matrixStats)
library(Rcpp)
sourceCpp("gsl.cpp")


Aj = cbind(rep(0.3,1000), rep(0.1,1000),rep(0.2,1000),rep(0.1,1000),rep(0.1,1000),rep(0.2,1000))
print(Aj)
D = gsl_mmm(Aj)
print(D)
print(colMeans(D))