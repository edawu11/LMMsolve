#####Method of Moments####
#####input####
##y:n*1 vector
##X:n*p design matrix
##Z:n*c covariate matrix
##ifeff: if use the unbiased estimator of the trace of K^2
##B:if ifeff is TURE, then set B (e.g. B=100)
#####output####
##estimator: a list:sb2,se2,omega

LMM_MoM = function(y,X,Z,ifintercept=F,ifK2=T,B=100,seed=2023){
  ini_time = Sys.time()
  X = as.matrix(X)
  Z = as.matrix(Z)
  
  set.seed(seed)
  p = ncol(X)
  n = length(y)
  
  
  if(ifintercept==F){
    q = ncol(Z)
  }
  else{
    Z = cbind(1,Z)
    q = ncol(Z)
  }
  y = scale(y,center = T,scale = F) # y is centered, however, no difference?
  X = scale(X,center = T,scale = F) # each column of X is centered
  Xsd = sqrt(colMeans(X^2))
  X = t(t(X)/Xsd)/sqrt(p)
  
  K = X%*%t(X)
  V = diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)
  
  VK = V%*%K
  Vy = V%*%y
  trVK = sum(diag(VK))
  
  if(ifK2==T){
    trVK2 = sum(VK^2)
  }
  else{
    Zd = matrix(rnorm(n*B),n,B)
    trVK2 = sum((VK%*%Zd)^2)/B
  }
  
  # solve the normal equations
  VyK = t(Vy)%*%K%*%Vy
  Vyy = t(Vy)%*%y
  coeffmat = matrix(c(trVK2,trVK,trVK,n-q),2,2)
  # print(coeffmat)
  sigma2vec = solve(coeffmat,c(VyK,Vyy))
  sigma2_beta = sigma2vec[1]
  sigma2_e = sigma2vec[2]
  
  # solve w, use weighted least square
  U = sigma2_beta*K + sigma2_e*diag(n)
  U_inv = solve(U)
  omega = solve(t(Z)%*%U_inv%*%Z)%*%t(Z)%*%U_inv%*%y
  
  result = list()
  result[["sb2"]] = sigma2_beta
  result[["se2"]] = sigma2_e
  result[["omega"]] = omega
  result[["usetime"]] = Sys.time()-ini_time
  return(result)
}

