#####Method of Moments####
#####input####
##y:n*1 vector
##X:n*p design matrix
##Z:n*c covariate matrix
##ifintercept: if add 1 vector to Z
##ifK2: if use the unbiased estimator of the trace of K^2
##B:if ifK2 is FALSE, then set B (e.g. B=100)
##seed
#####output####
##estimator: a list:sb2,se2,omega,cov_sig2,h2,cov_h2

LMM_MoM = function(y,X,Z,ifintercept=F,ifK2=T,B=100,seed=2023){
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
  A = matrix(c(trVK2,trVK,trVK,n-q),2,2)
  # print(A)
  rc = c(VyK,Vyy)
  A_inv = solve(A)
  sigma2 = A_inv%*%rc
  
  sb2 = sigma2[1]
  se2 = sigma2[2]
  
  # solve w, use weighted least square
  U = sb2*K + se2*diag(n)
  U_inv = solve(U)
  omega = solve(t(Z)%*%U_inv%*%Z)%*%t(Z)%*%U_inv%*%y
  
  # compute the variance matrix of sigma2 using sandwich estimator
  
  cov_B = matrix(0,2,2)
  VSig = sb2*VK + se2*V
  VKVSig = VK%*%VSig
  
  cov_B[1,1] = 2*sum(VKVSig^2)
  cov_B[2,2] = 2*sum(VSig^2)
  cov_B[1,2] = 2*sum(diag(VKVSig%*%VSig))
  cov_B[2,1] = cov_B[1,2]
  
  cov_sig2 = A_inv%*%cov_B%*%A_inv
  
  # compute the variance matrix of h2 using delta method
  
  sigma2_total = sum(sigma2)
  h2 = sb2/sigma2_total
  
  dh2 = c(se2/(sigma2_total^2),-sb2/(sigma2_total^2))
  
  cov_h2 = as.numeric(t(dh2)%*%cov_sig2%*%dh2) 
  
  result = list(sb2 =sb2,
                se2 = se2,
                omega = omega,
                cov_sig2 = cov_sig2,
                h2 = h2,
                cov_h2 = cov_h2)
  return(result)
}

