#to check GIC for correlated variables
library(glmnet)
library(MASS)
library(lars)
library(hdi)
require(stabs)
library(c060)

#set.seed(10102016)
#q = .001
n=200
p=1000
err.sigma = 3
s = 20 #number of active features
true.beta = c(rep(3,s),rep(0,p-s))
S = c(1:s)

########generate orthonormal data
cov.mat = diag(1,p)
generate.data.orth = function(p, n, true.beta, err.sigma)
{
  x.train= mvrnorm(n=n, mu = rep(0,p), Sigma = cov.mat)
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train)
  return(result)
}


#with different number if observations
runAll = function(n)
{
  d = generate.data.orth(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  fitAdaLasso(d)
  fitPLSS(d)
  fitStablity(d)
}