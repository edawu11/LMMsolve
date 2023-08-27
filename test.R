set.seed(100)
n <- 1000
p <- 50
q = 30
sb2 <- 0.1
se2 <- 1
X <- matrix(rnorm(n*p),n,p)
X <- scale(X)/sqrt(p)
Z = matrix(rnorm(n*q),n,q)
w <- c(rnorm(p,0,sqrt(sb2)))
y0 <- X%*%w
alpha = rnorm(q)
y <- Z%*%alpha + y0 + sqrt(se2)*rnorm(n)

# ini_time = Sys.time()
result1 = linRegMoM(X,y,Z,approx_se=F,rv_approx=F,n_rv=20,seed=2023)
# Sys.time() - ini_time
# c(result1$sb2,result1$se2)
result1$covSig
result2 = LMM_MoM(y,X,Z,ifintercept=T,ifK2=T,B=100,seed=2023)
result2$sb2
result2$se2

#######
eigenK <- eigen(K)
eVal <- eigenK$values
eVec <- eigenK$vectors
mat = t(eVec)%*%diag(n)%*%eVec

eigenI = eigen(diag(3))
eigenI$values

fit_MM=linRegMM(X=X,y=y,Z=Z,tol=1e-6,sb2=0.2,se2=1,maxIter = 500,verbose=F)
fit_MM$sb2
fit_MM$se2
fit_MM$allsb2
myMM = LMM_MM(y,X,Z,maxiter=500,tol=1e-6,sigma2_beta=0.2,sigma2_e=1,ifintercept=T,ifed=T,seed=100)
myMM$allsb2

a = sum(fit_MM$K)
aa = sum(myMM$K)

