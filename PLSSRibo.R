#PLSS for Real Data set- Riboflavin
library(glmnet)
library(MASS)
library(lars)
library(hdi)
require(stabs)
library(c060)

rib = data(riboflavin)
dim(riboflavin)
#x = scale(riboflavin$x) #4088 gene predictor
#y = scale(riboflavin$y, scale = FALSE)
x = riboflavin$x #4088 gene predictor
y = riboflavin$y

#Set the active predictors as follows:
# 20 predictors most correlated with the response
#n = dim(x)[1]
#p = dim(x)[2]





riboLasso = funtion()
{
  res = cv.glmnet(x, y, alpha = 1, intercept = FALSE, standardize = FALSE)
  fit = coef(res, res$lambda.min, mode="lambda")
  beta.lasso = fit[-1]
  print("LASSO")
  s.lasso = which(beta.lasso!=0)
  TP = length(intersect(s.lasso, S.ribo))
  FP = length(s.lasso) - TP
  print(TP)
  print(FP)
}

#fit Ada Lasso, weights are computed using Ada Lasso, 
riboAdaLasso = function()
{
  fit.lasso = cv.glmnet(x, y , alpha = 1, intercept = FALSE, standardize = FALSE)
  coef.lasso = coef(fit.lasso)
  beta.lasso = coef.lasso[-1]
  print("Ribo LASSO")
  s.lasso = which(beta.lasso!=0)
  TP = length(intersect(s.lasso, S.ribo))
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
  print("Ribo ADA LASSO")
  s.adalasso = s.lasso[which(beta.adalasso!=0)]
  TP = length(intersect(s.adalasso, S.ribo))
  FP = length(s.adalasso) - TP
  print(TP)
  print(FP)
}

riboPLSS = function(d)
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
  
  print("Ribo PLSS")  
  s.plss = s.lasso[fit.stab$stable]
  TP = length(intersect(s.plss, S.ribo))
  FP = length(s.plss) - TP
  print(TP)
  print(FP)
  
}

#fit stability selection using Lasso
riboStablity = function()
{
  spath = stabpath( y, x, size=0.5, weakness=1)
  st = stabsel(spath)
  print(st$stable) #number of selected variables
  fit.stab = plot(spath, error=0.5, type="pcer", pi_thr=0.9, xvar=c("lambda"),
                  col.all="black", col.sel="red")
  
  s.stab = fit.stab$stable
  TP = length(intersect(s.stab, S.ribo))
  FP = length(s.stab) - TP
  print(TP)
  print(FP)
  
}

top1000 = function()
{
  sig = var(x);
  dig = diag(sig)
  sdiag = order(dig, decreasing = TRUE)
  top100 = sdiag[1000]
  x.100 = as.matrix.data.frame(x[,top100])
  dim(x.100)
}

n = 71
p = 500
err.sigma = 1
S.ribo =sample(200,5)
S = S.ribo

gen.data.ribo = function()
{
  corVec = cor(y,x)
  corVec = as.vector(corVec)
  st = sort(abs(corVec), decreasing = TRUE)
  thr = st[p]
  temp = which(abs(corVec)>=thr)
  #S.ribo =sample(500,10)
  
  x.train= x[,temp]
  true.beta = rep(0, p)
  true.beta[S.ribo] = 3
  
  y.train = drop (x.train %*% true.beta) + err.sigma*rnorm (n)
  
  result = list(x.train = x.train,  y.train = y.train)
  return(result)
}

runAll = function(n)
{
  #S.ribo = c(1, 3, 15, 76, 82) #sample(150,5)
  S.ribo = sample(20,10)
  S = S.ribo
  d = gen.data.ribo()
  fitLasso(d)
  fitAdaLasso(d)
  fitPLSS(d)
  fitStablity(d)
}

