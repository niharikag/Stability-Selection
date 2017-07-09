#to check GIC for correlated variables
library(glmnet)
library(MASS)
library(lars)
library(hdi)
require(stabs)

#set.seed(10102016)
#q = .001
n=400
p=1000
err.sigma = 3
s = 20 #number of active features
true.beta = c(rep(3,s),rep(0,p-s))
S = c(1:s)

########generate orthonormal data
cov.mat.cor = diag(1,p)
rho = 0.5
for(i in 1:p)
  for(j in 1:p)
  {
    pw = abs(i-j)
    cov.mat = (rho)^pw
    
  }
generate.data.cor = function(p, n, true.beta, err.sigma)
{
  x.train= mvrnorm(n=n, mu = rep(0,p), Sigma = cov.mat.cor)
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train)
  return(result)
}


#with different number if observations
runAll = function(n)
{
  d = generate.data.cor(p, n=n, true.beta, err.sigma)
  fitLasso(d)
  fitAdaLasso(d)
  fitPLSS(d)
  fitStablity(d)
}