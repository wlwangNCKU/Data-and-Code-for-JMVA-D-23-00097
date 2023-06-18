library(mvtnorm)
library(tmvtnorm)
library(cubature)
library(invgamma)
library(Bessel)

source(paste(WD.PATH, 'function/mni.r', sep=''))
source(paste(WD.PATH, 'function/tmni.r', sep=''))
source(paste(WD.PATH, 'function/TMNImoment.r', sep=''))
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))

# True parameters
mu=c(1,2,3)
Sigma=matrix(0.5*c(9,4,2,4,2,1,2,1,2), 3,3)
nu = 4.0001       # MVT
nu1 = 2           # MSL
alpha1 = 2; beta1 = 1 # MVG
alpha2 = 1; beta2 = 2 # MDE
vech.S = vech.posi(3)

n = c(10, 50, 100, seq(0, 10000, 200)[-1], seq(10000, 25000, 500)[-1], 
    seq(25000, 40000, 1000)[-1], seq(40000, 50000, 2500)[-1], 
    seq(50000, 100000, 5000)[-1])
Rep = length(n)

set.seed(12)
### MVT ###
SM.TMVT = matrix(0, nrow=Rep, ncol=9)
MY1 = TMNI.moment(mu, Sigma, nu, distr='MVT', a.low=a.low, a.upp=a.upp)
par.TMVT = c(MY1$EY, MY1$CovY[vech.S])
for(r in 1: Rep){
  cat('Rep = ', r, 'TMVT with sample size n = ', n[r], '\n')
  y1 = rtmni(n=n[r], mu, Sigma, nu, distr='MVT', lower=a.low, upper=a.upp)
  SM.TMVT[r, ] = c(colMeans(y1), cov(y1)[vech.S])
  write(c(SM.TMVT[r, ]), paste(WD.PATH1, 'mTMVT.txt',sep=""), ncol=9, append=T)
}

### MSL ###
SM.TMSL = matrix(0, nrow=Rep, ncol=9)
MY2 = TMNI.moment(mu, Sigma, nu=nu1, distr='MSL', a.low=a.low, a.upp=a.upp)
par.TMSL = c(MY2$EY, MY2$CovY[vech.S])
for(r in 1: Rep){
  cat('Rep = ', r, 'TMSL with sample size n = ', n[r], '\n')
  y2 = rtmni(n=n[r], mu, Sigma, nu=nu1, distr='MSL', lower=a.low, upper=a.upp)
  SM.TMSL[r, ] = c(colMeans(y2), cov(y2)[vech.S])
  write(c(SM.TMSL[r, ]), paste(WD.PATH1, 'mTMSL.txt',sep=""), ncol=9, append=T)
}

### MCN ###
SM.TMCN = matrix(0, nrow=Rep, ncol=9)
MY3 = TMNI.moment(mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a.low, a.upp=a.upp)
par.TMCN = c(MY3$EY, MY3$CovY[vech.S])
for(r in 1: Rep){
  cat('Rep = ', r, 'TMCN with sample size n = ', n[r], '\n')
  y3 = rtmni(n=n[r], mu, Sigma, nu=0.25, rho=0.2, distr='MCN', lower=a.low, upper=a.upp)
  SM.TMCN[r, ] = c(colMeans(y3), cov(y3)[vech.S])
  write(c(SM.TMCN[r, ]), paste(WD.PATH1, 'mTMCN.txt',sep=""), ncol=9, append=T)
}

### MVG ### 
SM.TMVG = matrix(0, nrow=Rep, ncol=9)
MY4 = TMNI.moment(mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG', a.low=a.low, a.upp=a.upp)
par.TMVG = c(MY4$EY, MY4$CovY[vech.S])
for(r in 1: Rep){
  cat('Rep = ', r, 'TMVG with sample size n = ', n[r], '\n')
  y4 = rtmni(n=n[r], mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG', lower=a.low, upper=a.upp)
  SM.TMVG[r, ] = c(colMeans(y4), cov(y4)[vech.S])
  write(c(SM.TMVG[r, ]), paste(WD.PATH1, 'mTMVG.txt',sep=""), ncol=9, append=T)
}

### MDE (Laplace) 
SM.TMDE = matrix(0, nrow=Rep, ncol=9)
MY5 = TMNI.moment(mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG', a.low=a.low, a.upp=a.upp)
par.TMDE = c(MY5$EY, MY5$CovY[vech.S])
for(r in 1: Rep){
  cat('Rep = ', r, 'TMDE with sample size n = ', n[r], '\n')
  y5 = rtmni(n=n[r], mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG', lower=a.low, upper=a.upp)
  SM.TMDE[r, ] = c(colMeans(y5), cov(y5)[vech.S])
  write(c(SM.TMDE[r, ]), paste(WD.PATH1, 'mTMDE.txt',sep=""), ncol=9, append=T)
}
