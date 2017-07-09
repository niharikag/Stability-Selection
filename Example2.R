##EXample 2A from the paper, Correlated variables in regression:
# clustering and sparse estimation

library(glmnet)
library(c060)
library(gglasso)
library(ClustOfVar)
library(SGL)
library(parcor)
library(psych)
library(hdi)


set.seed (2221)
n=100
p=100
err.sigma = 15
s0 = seq(1:20) #true active set
#seq1 = seq(2.1,4,.1)
seq1 = seq(.1,2,.1)
lambda = seq(.1,.8,.02)
k=10

true.beta <- c(rep(1,20),rep(0,p-20))
true.beta[1:20] = sample(seq1,20, replace=FALSE)    
s0Comp = which(true.beta==0) #noise feature index

varMat = matrix(1,10,10)
for(i in 1:10)
  for(j in 1:10)
  {
    if(i!=j)
      varMat[i,j] = .9
  }
#create block diagonal matrix
blockMat = superMatrix(rep(list(varMat),10))

########generate data
generate.data2 = function(p,n,true.beta,err.sigma)
{
  x.train = mvrnorm(n,rep(0,p),blockMat)
  y.train = drop (x.train %*% true.beta) + rnorm (n, 0, err.sigma)
  
  x.valid = mvrnorm(n,rep(0,p), blockMat)
  y.valid = drop (x.valid %*% true.beta) + rnorm (n, 0, err.sigma)
  
  
  result = list(x.train = x.train,  y.train = y.train, x.valid = x.valid, y.valid=y.valid)
  return(result)
}

#the following function performs cluster representative Lasso
clustGroupLasso = function(d)
{
  x.train = d$x.train
  x.valid = d$x.valid
  y.train = d$y.train
  y.valid = d$y.valid
  
  p = ncol(x.train)
  
  colNames =  paste("X", 1:p, sep = "")
  
  tr = hclustvar(x.train)
  plot(tr)
  
  rect.hclust(tr, k=k, border="red")
  part <- cutreevar(tr,k)
  clust = part$cluster
  
  #fit group lasso penalized least squares
  fit1 <- gglasso(x.train,y.train,group=clust,loss="ls",lambda=lambda )
  
  y.estimate = predict.gglasso(fit1, x.valid,lambda = lambda)
  
  gridLen = length(lambda)
  
  #compute residuals
  y.residuals = y.valid - y.estimate
  MSE = rep(0,gridLen)
  for(i in 1:gridLen){
    MSE[i] = mean(y.residuals[,i]^2)
  }
  
  #minimum Prediction Error
  indx = which.min(MSE)
  predErr = MSE[indx]
  fit2 = coef(fit1,fit1$lambda[indx]) 
  
  result = list(predErr=predErr, fit=fit2)  
  return(result)
}

#the following function performs cluster representative Lasso
clustRepLasso = function(d, k=10)
{
  x.train = d$x.train
  x.valid = d$x.valid
  y.train = d$y.train
  y.valid = d$y.valid
  
  tr = hclustvar(x.train)
  plot(tr)
  
  rect.hclust(tr, k=k, border="red")
  part <- cutreevar(tr,k)
  clust = part$cluster
  #get cluster representatives
  x.clust.train =  matrix(0,n,k)
  x.clust.valid =  matrix(0,n,k)
  
  for(i in 1:n)
    for(j in 1:k){
      clus.indx = which(clust==j)
      x.clust.train[i,j] =  mean(x.train[i,clus.indx])
      x.clust.valid[i,j] =  mean(x.valid[i,clus.indx])
    }
  
  d1 = list(x.train = x.clust.train, y.train = y.train, 
            x.valid = x.clust.valid, y.valid = y.valid)
  
  res = performLasso(d1)
  
  fit = res$fit[-1]
  selClust = which(fit!=0) #selected clusters
  selIndx = which(clust %in% selClust)
  print(selIndx)
  #result = list(predErr=predErr, fit=fit2)  
  return(selIndx)
}

#the following function performs cluster representative Lasso
StabSelCRL = function(d, k=10)
{
  x.train = d$x.train
  x.valid = d$x.valid
  y.train = d$y.train
  y.valid = d$y.valid
  
  tr = hclustvar(x.train)
  plot(tr)
  
  rect.hclust(tr, k=k, border="red")
  part <- cutreevar(tr,k)
  clust = part$cluster
  #get cluster representatives
  x.clust.train =  matrix(0,n,k)
  x.clust.valid =  matrix(0,n,k)
  
  for(i in 1:n)
    for(j in 1:k){
      clus.indx = which(clust==j)
      x.clust.train[i,j] =  mean(x.train[i,clus.indx])
      x.clust.valid[i,j] =  mean(x.valid[i,clus.indx])
    }
  
  #res = performStabSelection(x.clust.train, y.train)
  fit.stab <- hdi(x.clust.train, y.train, method = "stability", EV = 2)
  fit.stab
  #print(res$stable)
  #result = list(predErr=predErr, fit=fit2)  
  return(fit.stab)  
}

######perform Lasso
performLasso = function(d){
  x.train = d$x.train
  x.valid = d$x.valid
  y.train = d$y.train
  y.valid = d$y.valid
  
  #fit Lasso
  #fit1 = cv.glmnet(x.train, y.train)
  fit1 = glmnet(x.train, y.train,lambda = lambda, alpha = 1)
  y.estimate = predict.glmnet(fit1, x.valid,lambda=lambda)
  
  gridLen = length(lambda)
  
  #compute residuals
  y.residuals = y.valid - y.estimate
  MSE = rep(0,gridLen)
  for(i in 1:gridLen){
    MSE[i] = mean(y.residuals[,i]^2)
  }
  
  #mean Prediction Error
  indx = which.min(MSE)
  predErr = MSE[indx]
  fit = coef(fit1,fit1$lambda[indx]) 
  result = list(predErr=predErr, fit=fit)
  return(result)
}

######perform Lasso
performStabSelection = function(x, y){
  #res <- stabpath(y.train,x.clust.train)
  res <- stabpath(y,x, weakness=1,mc.cores=2)
  st = stabsel(res)
  print(st$stable) #number of selected variables
  plot(res)
  
  #predErr = MSE[indx]
  #fit = coef(fit1,fit1$lambda[indx]) 
  #result = list(predErr=predErr, fit=fit)
  return(st)
}

compareMethods121 = function(p,n,true.beta,err.sigma)
{
  Lasso = NULL
  crl  = NULL 
  cgl  = NULL
  
  runLasso = TRUE
  runCGL = TRUE
  runCRL = TRUE
  lambda = seq(.1,.8,.01)
  
  for(k in 1:50)
  {
    d = generate.data2(p,n,true.beta, err.sigma)
    
    if(runLasso)
      lassoResult = performLasso(d)
    
    if(runCRL)
      crlResult = clustRepLasso(d)
    
    if(runCGL)
      cglResult = clustGroupLasso(d)
    
    
    if(is.null(Lasso))
    {
      if(runLasso)
        Lasso = lassoResult
      if(runCRL)
        crl  = crlResult 
      if(runCGL)
        cgl  = cglResult
      true.beta1 =true.beta
      true.beta2 = true.beta
      true.beta3 = true.beta
    }
    else 
    {
      if(runLasso)
        if(lassoResult$predErr < Lasso$predErr ){
          Lasso = lassoResult
          true.beta1 = true.beta
        }
      
      if(runCRL)
        if(crlResult$predErr < crl$predErr ){
          crl = crlResult
          true.beta2 = true.beta
        }
      if(runCGL)
        if(cglResult$predErr < cgl$predErr ){
          cgl  = cglResult
          true.beta3 = true.beta
        }
    }
    print(k)
  }
  
  print("Lasso")
  #Find TP and FP
  if(runLasso){
    cf = Lasso$fit[-1]
    #std1 = computeStd(Lasso$predErr,true.beta1,cf)
    set.predict = which(cf!=0) #predicted set
    print(set.predict)
  }
  
  print("CRL")
  if(runCRL){
    #cf = crl$fit[-1]
    #std2 = computeStd(crl$predErr,true.beta2, cf)
    #set.predict = which(cf!=0) #predicted set
    print(crl$fit)
  }
  
  print("CGL") 
  if(runCGL){
    cf = cgl$fit[-1]
    #std3 = computeStd(cgl$predErr,true.beta3, cf)
    set.predict = which(cf!=0) #predicted set
    print(set.predict)
  }
  
  #run stability selection and stability selction using CRL
  d = generate.data2(p, n, err.sigma, err.sigma)
  
  print("Stability Selection")
  fit.stab <- hdi(d$x.train, d$y.train, method = "stability", EV = 2)
  print(fit.stab$selected)
  
  fit.stab = StabSelCRL(d)
  print("Stability Selection using CRL") 
  print(fit.stab$selected)
}
d = generate.data2(p, n, true.beta, err.sigma)
res = StabSelCRL(d)
res$selected
res = performStabSelection(d$x.train,d$y.train)
res = clustRepLasso(d)
res = performLasso(d)

fit = res$fit[-1]
which(fit!=0)
