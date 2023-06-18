library(mvtnorm)
library(tmvtnorm)
library(cubature)
library(invgamma)
library(Bessel)

source(paste(WD.PATH, 'function/mni.r', sep=''))
source(paste(WD.PATH, 'function/tmni.r', sep=''))
source(paste(WD.PATH, 'function/TMNImoment.r', sep=''))
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))

rname = c('TMVT', 'TMSL', 'TMCN', 'TMVG', 'TMDE')
cname = c('mu1', 'mu2', 'mu3', 'sig11', 'sig21', 'sig22', 'sig31', 'sig32', 'sig33')

# True parameters
mu=c(1,2,3)
Sigma=matrix(0.5*c(9,4,2,4,2,1,2,1,2), 3,3)
nu = 4.0001       # MVT
nu1 = 2           # MSL
alpha1 = 2; beta1 = 1 # MVG
alpha2 = 1; beta2 = 2 # MDE
vech.S = vech.posi(3)

# doubly
ubd = 1
a.low=rep(-ubd, 3); a.upp=rep(ubd, 3)
### MVT ###
MY1 = TMNI.moment(mu, Sigma, nu, distr='MVT', a.low=a.low, a.upp=a.upp)
par.TMVT3 = c(MY1$EY, MY1$CovY[vech.S])
### MSL ###
MY2 = TMNI.moment(mu, Sigma, nu=nu1, distr='MSL', a.low=a.low, a.upp=a.upp)
par.TMSL3 = c(MY2$EY, MY2$CovY[vech.S])
### MCN ###
MY3 = TMNI.moment(mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a.low, a.upp=a.upp)
par.TMCN3 = c(MY3$EY, MY3$CovY[vech.S])
### MVG ### 
MY4 = TMNI.moment(mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG', a.low=a.low, a.upp=a.upp)
par.TMVG3 = c(MY4$EY, MY4$CovY[vech.S])
### MDE (Laplace) 
MY5 = TMNI.moment(mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG', a.low=a.low, a.upp=a.upp)
par.TMDE3 = c(MY5$EY, MY5$CovY[vech.S])

Table2 = rbind(par.TMVT3, par.TMSL3, par.TMCN3, par.TMVG3, par.TMDE3)
rownames(Table2) = rname
colnames(Table2) = cname
print(round(Table2, 4))
write.csv(Table2, paste(WD.PATH,'results/Table2.csv',sep=""), row.names = TRUE)

# right
ubd = 1
a.low=rep(-Inf, 3); a.upp=rep(ubd, 3)
### MVT ###
MY1 = TMNI.moment(mu, Sigma, nu, distr='MVT', a.low=a.low, a.upp=a.upp)
par.TMVT1 = c(MY1$EY, MY1$CovY[vech.S])
### MSL ###
MY2 = TMNI.moment(mu, Sigma, nu=nu1, distr='MSL', a.low=a.low, a.upp=a.upp)
par.TMSL1 = c(MY2$EY, MY2$CovY[vech.S])
### MCN ###
MY3 = TMNI.moment(mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a.low, a.upp=a.upp)
par.TMCN1 = c(MY3$EY, MY3$CovY[vech.S])
### MVG ### 
MY4 = TMNI.moment(mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG', a.low=a.low, a.upp=a.upp)
par.TMVG1 = c(MY4$EY, MY4$CovY[vech.S])
### MDE (Laplace) 
MY5 = TMNI.moment(mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG', a.low=a.low, a.upp=a.upp)
par.TMDE1 = c(MY5$EY, MY5$CovY[vech.S])

TableS1 = rbind(par.TMVT1, par.TMSL1, par.TMCN1, par.TMVG1, par.TMDE1)

# left
ubd = 1
a.low=rep(-ubd, 3); a.upp=rep(Inf, 3)                         # left
### MVT ###
MY1 = TMNI.moment(mu, Sigma, nu, distr='MVT', a.low=a.low, a.upp=a.upp)
par.TMVT2 = c(MY1$EY, MY1$CovY[vech.S])
### MSL ###
MY2 = TMNI.moment(mu, Sigma, nu=nu1, distr='MSL', a.low=a.low, a.upp=a.upp)
par.TMSL2 = c(MY2$EY, MY2$CovY[vech.S])
### MCN ###
MY3 = TMNI.moment(mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a.low, a.upp=a.upp)
par.TMCN2 = c(MY3$EY, MY3$CovY[vech.S])
### MVG ### 
MY4 = TMNI.moment(mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG', a.low=a.low, a.upp=a.upp)
par.TMVG2 = c(MY4$EY, MY4$CovY[vech.S])
### MDE (Laplace) 
MY5 = TMNI.moment(mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG', a.low=a.low, a.upp=a.upp)
par.TMDE2 = c(MY5$EY, MY5$CovY[vech.S])

TableS1 = rbind(TableS1, par.TMVT2, par.TMSL2, par.TMCN2, par.TMVG2, par.TMDE2)

rownames(TableS1) = rep(rname, 2)
colnames(TableS1) = cname
print(round(TableS1, 4))
write.csv(TableS1, paste(WD.PATH,'results/TableS1.csv',sep=""), row.names = TRUE)

save.image(paste(WD.PATH, 'data/mMNI.RData', sep=''))
