#####MM Algorithm####
#####input####
##y:n*1 vector
##X:n*p design matrix
##Z:n*c covariate matrix
##ifeff: if use the unbiased estimator of the trace of K^2
##B:if ifeff is TURE, then set B (e.g. B=100)
#####output####
##estimator: a list:sb2,se2,omega

LMM_MM = function(y,X,Z,ifintercept=F,ifK2=T,B=100,seed=2023){
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
  
  # initialize 
  sigma2_beta = 0.2
  sigma2_e = 1
  allL = c()
  Sig = sigma2_beta*K+sigma2_e*diag(n)
  invSig = solve(Sig)
  
  for(t in 1:100){
    
    ZinvSig = t(Z) %*% invSig
    omega = solve(ZinvSig%*%Z)%*%ZinvSig%*%y
    res = y - Z%*%omega
    
    sigma2_beta = sigma2_beta*sqrt(t(res)%*%invSig%*%K%*%invSig%*%res/sum(invSig*K))
    
    sigma2_e = sigma2_e*sqrt(t(res)%*%invSig^2%*%res/sum(invSig*K))
    
    sigma2_beta = as.numeric(sigma2_beta)
    sigma2_e = as.numeric(sigma2_e)
    
    Sig = sigma2_beta*K+sigma2_e*diag(n)
    invSig = solve(Sig)
    ZinvSig = t(Z) %*% invSig
    L_sig = -0.5*log(det(Sig)) - 0.5 *t(res) %*%invSig %*%res
    allL = c(allL,L_sig)
  }
  
  
  
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
  
  sigma2_beta = sigma2[1]
  sigma2_e = sigma2[2]
  
  # solve w, use weighted least square
  U = sigma2_beta*K + sigma2_e*diag(n)
  U_inv = solve(U)
  omega = solve(t(Z)%*%U_inv%*%Z)%*%t(Z)%*%U_inv%*%y
  
  # compute the variance matrix of sigma2 using sandwich estimator
  
  cov_B = matrix(0,2,2)
  VSig = sigma2_beta*VK + sigma2_e*V
  VKVSig = VK%*%VSig
  
  cov_B[1,1] = 2*sum(VKVSig^2)
  cov_B[2,2] = 2*sum(VSig^2)
  cov_B[1,2] = 2*sum(diag(VKVSig%*%VSig))
  cov_B[2,1] = cov_B[1,2]
  
  cov_sigma2 = A_inv%*%cov_B%*%A_inv
  
  # compute the variance matrix of h2 using delta method
  
  sigma2_total = sum(sigma2)
  h2 = sigma2_beta/sigma2_total
  
  dh2 = c(sigma2_e/(sigma2_total^2),-sigma2_beta/(sigma2_total^2))
  
  cov_h2 = t(dh2)%*%cov_sigma2%*%dh2
  
  result = list()
  result[["sb2"]] = sigma2_beta
  result[["se2"]] = sigma2_e
  result[["omega"]] = omega
  result[["cov_sigma2"]] = cov_sigma2
  result[["h2"]] = h2
  result[["cov_h2"]] = cov_h2
  return(result)
}
