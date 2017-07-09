#to check GIC for correlated variables
library(glmnet)
library(MASS)
library(lars)
library(hdi)
require(stabs)

set.seed(10102016)
q = .001
n=50
p=100
err.sigma = 1
s = 5 #number of active features
true.beta = c(rep(2,s),rep(0,p-s))
n.valid = 1000
S0 = c(1:s)
########generate orthonormal data
generate.data = function(p, n, true.beta, err.sigma)
{
  x.train=matrix(0,n,p)
  x.valid =matrix(0,n.valid,p)
  
  #independent features
  for(i in 1:p)
  {
    x.train[,i] = rnorm(n)
    x.valid[,i] = rnorm(n.valid)
  }
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  y.valid = drop (x.valid %*% true.beta) + err.sigma*rnorm (n.valid)
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

d = generate.data(p, n, true.beta, err.sigma)

fitLasso = function(d)
{
    fit.lasso = cv.glmnet(x, y , alpha = 1, intercept = FALSE, standardize = FALSE)
    coef.lasso = coef(fit.lasso)
    beta.lasso = coef.lasso[-1]
    which(beta.lasso!=0)
}

fitPLSS = function(d)
{
  x = scale(d$x.train)
  y = scale(d$y.train, scale = FALSE)
  
  
  lambda.min = 1.5*sqrt(log(p)/n)
  
  res.lasso = glmnet(x, y, alpha = 1, intercept = FALSE, standardize = FALSE)
  
  lasso.coef =  coef(res.lasso, s=lambda.min, mode="lambda")
  lasso.beta = lasso.coef[-1]
  s.lasso = which(lasso.beta!=0)
  s.lasso
  x.red = x[,s.lasso]
  dim(x.red)
  
  fit.stab <- stabsel(x.red, y, fitfun = glmnet.lasso, cutoff = 0.6, PFER = 4)
  fit.stab$selected
  plot(fit.stab, main = "Lasso")  
}

fitStablity = function()
{
  x = scale(d$x.train)
  y = scale(d$y.train, scale = FALSE)
  
  fit.stab <- stabsel(x, y, fitfun = glmnet.lasso, cutoff = 0.6, PFER = 4)
  print(fit.stab$selected)
}

{
fit.stab <- hdi(x.red, y, method = "stability", threshold = 0.5, EV=1)
fit.stab$selected
fit.stab$freq

A = diag(1,4)
rho = -0.4
B = rep(rho,4)
C = cbind(A,B)
B = c(rep(rho,4),1)
C = rbind(C,B)
C
eigen(C)
}

