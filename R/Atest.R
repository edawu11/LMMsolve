#####Asymptoic test####
#####input####
##thetahat: estimator
##S: the covariance matrix of the estimator
#####output####
##p_value

Atest = function(thetahat,theta0,S){
  k = length(thetahat)
  S_inv = solve(S)
  statistic = t(thetahat-theta0) %*% S_inv %*% (thetahat-theta0)
  p_value = pchisq(q=as.numeric(statistic),df =k,lower.tail = FALSE)
  return(p_value)
}
