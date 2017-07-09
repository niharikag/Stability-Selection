#Post Lasso Stability Selection (with Adaptive weights)
library(glmnet)
library(MASS)
library(lars)
library(hdi)
require(stabs)
library(c060)

set.seed(1)
n=200
p=1000
err.sigma = 1
s = 4 #number of active features
true.beta = c(rep(3,s),rep(0,p-s))
S = c(1:s)
cov.mat = diag(1,p)
cov.mat[1,5]=cov.mat[5,1]=.251
cov.mat[2,5]=cov.mat[5,2]=.251
cov.mat[3,5]=cov.mat[5,3]=.251
cov.mat[4,5]=cov.mat[5,4]=.251


########generate orthonormal data
generate.data = function(p, n, true.beta, err.sigma)
{
  x.train= mvrnorm(n=n, mu = rep(0,p), Sigma = cov.mat)
  #x.train= mvrnorm(n=300, mu = rep(0,p), Sigma = cov.mat, empirical = TRUE)
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train)
  return(result)
}

#with different number of observations
runAll = function(n)
{
  set.seed(0)
  d = generate.data(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  fitPLSS(d)
  fitStablity(d)
  
}



