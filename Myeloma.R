#PLSS for Real Data set- tian
library(glmnet)
library(MASS)
library(hdi)
require(stabs)
library(c060)
library(datamicroarray)

data.tian = data('tian', package = 'datamicroarray')
dim(tian$x)
#x = scale(tian$x) #4088 gene predictor
#y = scale(tian$y, scale = FALSE)
x = tian$x #12625 gene predictor
y = tian$y

#Set the active predictors as follows:
# 20 predictors most correlated with the response
n = dim(x)[1]
p = 1000





tianLasso = funtion()
{
res = cv.glmnet(x, y, alpha = 1, intercept = FALSE, standardize = FALSE)
fit = coef(res, res$lambda.min, mode="lambda")
beta.lasso = fit[-1]
print("LASSO")
s.lasso = which(beta.lasso!=0)
TP = length(intersect(s.lasso, S.tian))
FP = length(s.lasso) - TP
print(TP)
print(FP)
}

#fit Ada Lasso, weights are computed using Ada Lasso, 
tianAdaLasso = function()
{
  fit.lasso = cv.glmnet(x, y , alpha = 1, intercept = FALSE, standardize = FALSE)
  coef.lasso = coef(fit.lasso)
  beta.lasso = coef.lasso[-1]
  print("tian LASSO")
  s.lasso = which(beta.lasso!=0)
  TP = length(intersect(s.lasso, S.tian))
  FP = length(s.lasso) - TP
  print(TP)
  print(FP)
  
  print(s.lasso)
  x.red = x[,s.lasso]
  dim(x.red)
  
  weights = abs(beta.lasso[s.lasso])
  x.red = t(t(x.red)*weights)
  
  fit.adalasso = cv.glmnet(x.red, y , alpha = 1, intercept = FALSE, standardize = FALSE)
  coef.adalasso = coef(fit.adalasso)
  beta.adalasso = coef.adalasso[-1]
  print("tian ADA LASSO")
  s.adalasso = s.lasso[which(beta.adalasso!=0)]
  TP = length(intersect(s.adalasso, S.tian))
  FP = length(s.adalasso) - TP
  print(TP)
  print(FP)
}

tianPLSS = function(d)
{
  lambda.min = sqrt(log(p)/n)
  
  res.lasso = glmnet(x, y, alpha = 1, intercept = FALSE, standardize = FALSE)
  
  lasso.coef =  coef(res.lasso, s=lambda.min, mode="lambda")
  lasso.beta = lasso.coef[-1]
  s.lasso = which(lasso.beta!=0)
  print(s.lasso)
  x.red = x[,s.lasso]
  dim(x.red)
  
  weights = abs(lasso.beta[s.lasso])
  x.red = t(t(x.red)*weights)
  
  #fit.stab = hdi(x.red, y, method="stability", B=100, EV=5)
  #fit.stab$freq
  
  spath = stabpath( y, x.red, size=0.5, weakness=1)
  st = stabsel(spath)
  print(st$stable) #number of selected variables
  
  fit.stab = plot(spath, error=0.5, type="pcer", pi_thr=0.9, xvar=c("lambda"),
                  col.all="black", col.sel="red")
  
  print("tian PLSS")  
  s.plss = s.lasso[fit.stab$stable]
  TP = length(intersect(s.plss, S.tian))
  FP = length(s.plss) - TP
  print(TP)
  print(FP)
  
}

#fit stability selection using Lasso
tianStablity = function()
{
  spath = stabpath( y, x, size=0.5, weakness=1)
  st = stabsel(spath)
  print(st$stable) #number of selected variables
  fit.stab = plot(spath, error=0.5, type="pcer", pi_thr=0.9, xvar=c("lambda"),
                  col.all="black", col.sel="red")
  
  s.stab = fit.stab$stable
  TP = length(intersect(s.stab, S.tian))
  FP = length(s.stab) - TP
  print(TP)
  print(FP)
  
}
sig = var(tian$x);
top1000 = function()
{
  dig = diag(sig)
  sdiag = sort(dig, decreasing = TRUE)
  top100 = sdiag[1000]
  pickCol = which(dig >= top100)
  x.red = as.matrix(tian$x[,pickCol])
  return(x.red)
}

p = 1000
err.sigma = 3
S.tian =sample(200,5)
S = S.tian

gen.data.tian = function()
{
  
  x.train= top1000()
  x.train = scale(x.train)
  true.beta = rep(0, p)
  true.beta[S.tian] = 3
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train)
  return(result)
}

runAll = function(n)
{
  #S.tian = c(1, 3, 15, 76, 82) #sample(150,5)
  S.tian = sample(500,10)
  s=10
  S = S.tian
  d = gen.data.tian()
  fitLasso(d)
  fitAdaLasso(d)
  fitPLSS(d)
  fitStablity(d)
}

