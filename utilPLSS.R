#PLSS utility file

#fit Lasso, using 10-fold cross varidation
fitLasso = function(d)
{
  x = scale(d$x.train)
  y = scale(d$y.train, scale = FALSE)
  
  fit.lasso = cv.glmnet(x, y , alpha = 1, intercept = FALSE, standardize = FALSE)
  coef.lasso = coef(fit.lasso)
  beta.lasso = coef.lasso[-1]
  print("LASSO")
  s.lasso = which(beta.lasso!=0)
  print(s.lasso)
  TP = length(intersect(s.lasso,S))
  FP = (length(s.lasso) - TP)/
  print(TP/s)
  print(FP/length(s.lasso))
}

#fit Ada Lasso, weights are computed using Ada Lasso, 
fitAdaLasso = function(d)
{
  x = scale(d$x.train)
  y = scale(d$y.train, scale = FALSE)
  
  fit.lasso = cv.glmnet(x, y , alpha = 1, intercept = FALSE, standardize = FALSE)
  coef.lasso = coef(fit.lasso)
  beta.lasso = coef.lasso[-1]
  print("LASSO")
  s.lasso = which(beta.lasso!=0)
  TP = length(intersect(s.lasso,S))
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
  print("ADA LASSO")
  s.adalasso = s.lasso[which(beta.adalasso!=0)]
  TP = length(intersect(s.adalasso,S))
  FP = length(s.adalasso) - TP
  print(TP/s)
  print(FP/length(s.adalasso))
}

#fit PLSS, Lasso with min Beta, then Stability with weighted Lasso
fitPLSS = function(d)
{
  x = scale(d$x.train)
  y = scale(d$y.train, scale = FALSE)
  
  #lambda.min = err.sigma*sqrt(log(p)/n)
  lambda.min = err.sigma*sqrt(log(p)/n)
  
  res.lasso = glmnet(x, y, alpha = 1, intercept = FALSE, standardize = FALSE)
  
  lasso.coef =  coef(res.lasso, s=lambda.min, mode="lambda")
  lasso.beta = lasso.coef[-1]
  s.lasso = which(lasso.beta!=0)
  print(s.lasso)
  x.red = x[,s.lasso]
  dim(x.red)
  
  weights = round(abs(lasso.beta[s.lasso]),2)
  x.red = t(t(x.red)*weights)
  
  fit.stab = hdi(x.red, y, method="stability", B=100, EV=2)
  fit.stab$freq
  
  spath = stabpath( y, x.red, size=0.5, weakness=1)
  st = stabsel(spath)
  print(st$stable) #number of selected variables
  
  fit.stab = plot(spath, error=0.4, type="pcer", pi_thr=0.9, xvar=c("lambda"),
       col.all="black", col.sel="red")#, penpath=FALSE)
  
  print("PLSS")  
  s.plss = s.lasso[fit.stab$stable]
  TP = length(intersect(s.plss,S))
  FP = length(s.plss) - TP
  print(TP/s)
  print(FP/length(s.plss))
  
}

#fit stability selection using Lasso
fitStablity = function(d)
{
  x = scale(d$x.train)
  y = scale(d$y.train, scale = FALSE)
  
  #fit.stab = hdi(x, y, method="stability", B=100, EV=15)
  #fit.stab$freq
  
  spath = stabpath( y, x, size=0.5, weakness=1)
  st = stabsel(spath)
  print(st$stable) #number of selected variables
  fit.stab = plot(spath, error=0.001, type="pcer", pi_thr=0.9, xvar=c("lambda"),
                  col.all="black", col.sel="red") #, penpath=FALSE)
  
  s.stab = fit.stab$stable
  TP = length(intersect(s.stab,S))
  FP = length(s.stab) - TP
  print(TP/s)
  print(FP/length(s.stab) )
  
}

#The same as PLSS, just trying something different
iterate = function(iter = 10){
  for(i in 1:iter)
  {
    d = generate.data(p, n=n, true.beta, err.sigma)
    x = scale(d$x.train)
    y = scale(d$y.train, scale = FALSE)
    
    
    #lambda.min = err.sigma*sqrt(log(p)/n)
    lambda.min = err.sigma*sqrt(log(p)/n)
    
    res.lasso = glmnet(x, y, alpha = 1, intercept = FALSE, standardize = FALSE)
    
    lasso.coef =  coef(res.lasso, s=lambda.min, mode="lambda")
    lasso.beta = lasso.coef[-1]
    s.lasso = which(lasso.beta!=0)
    print(s.lasso)
    x.red = x[,s.lasso]
    dim(x.red)
    
    weights = abs(lasso.beta[s.lasso])
    x.red = t(t(x.red)*weights)
    
    
    spath = stabpath( y, x.red, size=0.5, weakness=1, mc.cores=1, family="gaussian")
    st = stabsel(spath)
    print(st$stable) #number of selected variables
    #plot(spath)
    plot(spath, error=0.5, type="pcer", pi_thr=0.9, xvar=c("lambda"),
         col.all="black", col.sel="red", penpath=FALSE)
    
  }
}

#some un-used code, may use later
xyz = function(){
  d = generate.data(p, n=n, true.beta, err.sigma)
  x = scale(d$x.train)
  y = scale(d$y.train, scale = FALSE)
  
  spath = stabpath( y, x, weakness=.5)
  st = stabsel(spath)
  print(st$stable) #number of selected variables
  #plot(spath)
  plot(spath, error=0.05, type="pfer", pi_thr=0.9, xvar=c("lambda"),
       col.all="black", col.sel="red")
  
    
    spath = stabpath( y, x.red, size=0.5, weakness=1)
    st = stabsel(spath)
    print(st$stable) #number of selected variables
    #plot(spath)
    plot(spath, error=0.5, type="pcer", pi_thr=0.6, xvar=c("lambda"),
         col.all="black", col.sel="red",penpath=FALSE)
    
    
    fit.stab = hdi(x.red, y, method="stability", B=100, EV=10)
    #fit.stab = hdi(x.red, y)
    fit.stab$freq
    
    
    fit.stab <- stabsel(x.red, y, fitfun = glmnet.lasso, cutoff = 0.7, PFER = 1)
    print("PLSS")
    print(fit.stab$selected)
}