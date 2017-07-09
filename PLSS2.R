#Post Lasso Stability Selection (with Adaptive weights)
library(glmnet)
library(MASS)
library(lars)
library(hdi)
require(stabs)
library(c060)

#set.seed(10102016)
n=200
p=1000
err.sigma = 1
s = 5 #number of active features
true.beta = c(rep(1,s),rep(0,p-s))
S = c(1:s)
cov.mat = diag(1,p)
cov.mat[1,6]=cov.mat[6,1]=.21
cov.mat[2,6]=cov.mat[6,2]=.21
cov.mat[3,6]=cov.mat[6,3]=.21
cov.mat[4,6]=cov.mat[6,4]=.21
cov.mat[5,6]=cov.mat[6,5]=.21

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
  d = generate.data(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  fitPLSS(d)
  fitStablity(d)
  
}



