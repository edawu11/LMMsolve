#####MM Algorithm####
#####input####
##y:n*1 vector
##X:n*p design matrix
##Z:n*c covariate matrix
##ifeff: if use the unbiased estimator of the trace of K^2
##B:if ifeff is TURE, then set B (e.g. B=100)
#####output####
##estimator: a list:sb2,se2,omega

LMM_MM = function(y,X,Z,maxiter,tol=1e-6,sb2,se2,ifintercept=F,ifed=F,seed=2023){
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
  
  allL = c()
  allsb2 = sb2
  allse2 = se2
  
  if(ifed==T){
    # eigenK = eigen(K)
    # D = diag(eigenK$values)
    # U = eigenK$vectors
    # 
    # UZ = t(U)%*%Z
    # Uy = t(U)%*%y
    # d = solve(D*sb2+se2*diag(n))
    # for(t in 1:maxiter){
    # 
    #   UZd = t(UZ)%*%d
    #   omega = solve(UZd%*%UZ)%*%UZd%*%Uy
    # 
    #   res = Uy-UZ%*%omega
    # 
    #   sb2 = sb2*sqrt(t(res)%*%d%*%D%*%d%*%res/sum(d*D))
    # 
    #   se2 = se2*sqrt(t(res)%*%d%*%d%*%res/sum(d))
    # 
    #   sb2 = as.numeric(sb2)
    #   se2 = as.numeric(se2)
    # 
    #   d = solve(D*sb2+se2*diag(n))
    #   L_sig = 0.5*log(det(d)) - 0.5 *t(res) %*%d%*%res
    #   allL = c(allL,as.numeric(L_sig))
    # 
    #   if(t>1){
    #     if(abs(allL[t] - allL[t-1]) < tol){
    #       print(paste("The iters is:",t,sep=""))
    #       break
    #     }
    #   }
    # }
    eigenK = eigen(K)
    d = eigenK$values #vector
    U = eigenK$vectors #matrix
    
    UtZ = t(U)%*%Z
    Uty = t(U)%*%y
    invd = 1/(d*sb2+se2)
    
    for(t in 1:maxiter){
      
      tmp = t(UtZ*invd)
      # omega = solve(tmp%*%UZ) %*% tmp%*%Uty
      omega = solve(tmp%*%UtZ,tmp%*%Uty) 
      res = Uty-UtZ%*%omega
      
      sb2 = sb2*sqrt(sum(res^2*invd^2*d)/sum(invd*d))
      allsb2 = c(allsb2,sb2)
      
      se2 = se2*sqrt(sum(res^2*invd^2/sum(invd)))
      allse2 = c(allse2,se2)
      
      invd = 1/(d*sb2+se2)
      L_sig = 0.5*sum(log(invd)) - 0.5 *sum(res^2*invd)- n/2*log(2*pi)
      
      allL = c(allL,L_sig)
      if(t>1){
        if(abs(allL[t] - allL[t-1]) < tol){
          break
        }
      }
    }
    
    invSig = U%*%(1/(d*sb2+se2)*t(U)) 
  }
  else{
    Sig = sb2*K+se2*diag(n)
    invSig = solve(Sig)
    
    for(t in 1:maxiter){
      
      ZinvSig = t(Z) %*% invSig
      omega = solve(ZinvSig%*%Z)%*%ZinvSig%*%y
      res = y - Z%*%omega
      
      sb2 = sb2*sqrt(t(res)%*%invSig%*%K%*%invSig%*%res/sum(invSig*K))
      
      se2 = se2*sqrt(t(res)%*%invSig%*%invSig%*%res/sum(diag(invSig)))
      
      sb2 = as.numeric(sb2)
      se2 = as.numeric(se2)
      
      Sig = sb2*K+se2*diag(n)
      invSig = solve(Sig)
      ZinvSig = t(Z) %*% invSig
      L_sig = 0.5*log(det(Sig)) - 0.5 *t(res) %*%invSig %*%res- n/2*log(2*pi)
      
      allL = c(allL,L_sig)
      
      if(t>1){
        if((allL[t] - allL[t-1]) < tol){
          break
        }
      }
    }
  }
  
  #the covariance matrix
  
  invSigK = invSig%*%K
  FIM = matrix(0,2,2)
  FIM[1,1] = sum(invSigK^2)/2
  FIM[2,2] = sum(invSig^2)/2
  FIM[1,2] = sum(invSig*invSigK)/2
  FIM[2,1] = FIM[1,2]
  cov_sigma2 = solve(FIM)
  
  #delta method
  sigma2_total = sum(sb2+se2)
  h2 = sb2/sigma2_total
  
  dh2 = c(se2/(sigma2_total^2),-sb2/(sigma2_total^2))
  
  cov_h2 = as.numeric(t(dh2)%*%cov_sigma2%*%dh2) 
  
  return(list(sb2=sb2,se2=se2,allL=allL,allsb2=allsb2,allse2 = allse2,
              cov_sigma2 = cov_sigma2, cov_h2 = cov_h2))
}

