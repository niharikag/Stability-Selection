#to check GIC for correlated variables
library(glmnet)
library(MASS)
library(lars)
library(hdi)
require(stabs)
library(c060)

#set.seed(10102016)
#q = .001
n=500
p=5000
err.sigma = 1
s = 20 #number of active features
true.beta = c(rep(3,s),rep(0,p-s))
S = c(1:s)

########generate orthonormal data
cov.mat = diag(1,p)
generate.data.orth.old = function(p, n, true.beta, err.sigma)
{
  x.train= mvrnorm(n=n, mu = rep(0,p), Sigma = cov.mat)
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train)
  return(result)
}

generate.data.orth = function(p, n, true.beta, err.sigma)
{
  x.train= matrix(0,n,p)
  for(i in 1:n)
     for(j in 1:p)
     {
         x.train[i,j] = rnorm(1)
     }
	
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train)
  return(result)
}


#with different number if observations
runAll = function(n)
{
  p=5000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=10000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=15000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=20000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=25000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=30000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=35000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=40000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=45000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
  
  p=50000
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  system.time(fitAdaLasso(d))
  system.time(fitPLSS(d))
  system.time(fitStablity(d))
}

rho = .5
Mtr = matrix(c(1,rho,rho,1),2,2)
m1 = rbind(Mtr,Mtr)

