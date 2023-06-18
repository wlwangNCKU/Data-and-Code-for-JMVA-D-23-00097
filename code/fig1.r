library(mvtnorm)
library(tmvtnorm)
library(cubature)
library(Bessel)

source(paste(WD.PATH, 'function/mni.r', sep=''))
source(paste(WD.PATH, 'function/tmni.r', sep=''))
source(paste(WD.PATH, 'function/TMNImoment.r', sep=''))

mu = 0
Sigma = 1
ubd = 1
a1.low=-Inf; a1.upp=Inf
a2.low=-Inf; a2.upp=ubd
a3.low=-ubd; a3.upp=Inf
a4.low=-ubd; a4.upp=ubd

x = seq(-10, 10, length=1000)
m = length(x)

yMVN1 = yMVN2 = yMVN3 = yMVN4 =
yMVT1 = yMVT2 = yMVT3 = yMVT4 = numeric(m)
for(i in 1: m){
# MVN
  yMVN1[i] = dtmvnorm(x[i], mean=mu, sigma=Sigma)
  yMVN2[i] = dtmvnorm(x[i], mean=mu, sigma=Sigma, lower=a2.low, upper=a2.upp)
  yMVN3[i] = dtmvnorm(x[i], mean=mu, sigma=Sigma, lower=a3.low, upper=a3.upp)
  yMVN4[i] = dtmvnorm(x[i], mean=mu, sigma=Sigma, lower=a4.low, upper=a4.upp)
# MVT
  yMVT1[i] = dmni(x[i], mu, Sigma, nu=4, distr='MVT')
  yMVT2[i] = dtmni(x[i], mu, Sigma, nu=4, distr='MVT', a.low=a2.low, a.upp=a2.upp)
  yMVT3[i] = dtmni(x[i], mu, Sigma, nu=4, distr='MVT', a.low=a3.low, a.upp=a3.upp)
  yMVT4[i] = dtmni(x[i], mu, Sigma, nu=4, distr='MVT', a.low=a4.low, a.upp=a4.upp)
}

# MSL
yMSL1 = yMSL2 = yMSL3 = yMSL4 = numeric(m)
for(i in 1: m){
  yMSL1[i] = dmni(x[i], mu, Sigma, nu=2, distr='MSL')
  yMSL2[i] = dtmni(x[i], mu, Sigma, nu=2, distr='MSL', a.low=a2.low, a.upp=a2.upp)
  yMSL3[i] = dtmni(x[i], mu, Sigma, nu=2, distr='MSL', a.low=a3.low, a.upp=a3.upp)
  yMSL4[i] = dtmni(x[i], mu, Sigma, nu=2, distr='MSL', a.low=a4.low, a.upp=a4.upp)
}

# MCN
yMCN1 = yMCN2 = yMCN3 = yMCN4 = numeric(m)
for(i in 1: m){
  yMCN1[i] = dmni(x[i], mu, Sigma, nu=0.25, rho=0.2, distr='MCN')
  yMCN2[i] = dtmni(x[i], mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a2.low, a.upp=a2.upp)
  yMCN3[i] = dtmni(x[i], mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a3.low, a.upp=a3.upp)
  yMCN4[i] = dtmni(x[i], mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a4.low, a.upp=a4.upp)
}

# MVG
yMVG1 = yMVG2 = yMVG3 = yMVG4 = numeric(m)
for(i in 1: m){
  yMVG1[i] = dmni(x[i], mu, Sigma, alpha=2, beta=1, distr='MVG')
  yMVG2[i] = dtmni(x[i], mu, Sigma, alpha=2, beta=1, distr='MVG', a.low=a2.low, a.upp=a2.upp)
  yMVG3[i] = dtmni(x[i], mu, Sigma, alpha=2, beta=1, distr='MVG', a.low=a3.low, a.upp=a3.upp)
  yMVG4[i] = dtmni(x[i], mu, Sigma, alpha=2, beta=1, distr='MVG', a.low=a4.low, a.upp=a4.upp)
}

# MDE
yMDE1 = yMDE2 = yMDE3 = yMDE4 = numeric(m)
for(i in 1: m){
  yMDE1[i] = dmni(x[i], mu, Sigma, alpha=1, beta=2, distr='MVG')
  yMDE2[i] = dtmni(x[i], mu, Sigma, alpha=1, beta=2, distr='MVG', a.low=a2.low, a.upp=a2.upp)
  yMDE3[i] = dtmni(x[i], mu, Sigma, alpha=1, beta=2, distr='MVG', a.low=a3.low, a.upp=a3.upp)
  yMDE4[i] = dtmni(x[i], mu, Sigma, alpha=1, beta=2, distr='MVG', a.low=a4.low, a.upp=a4.upp)
}

a1 = max(c(yMVN1, yMVT1, yMSL1, yMCN1, yMVG1, yMDE1))
a2 = max(c(yMVN2, yMVT2, yMSL2, yMCN2, yMVG2, yMDE2))
a3 = max(c(yMVN3, yMVT3, yMSL3, yMCN3, yMVG3, yMDE3))
a4 = max(c(yMVN4, yMVT4, yMSL4, yMCN4, yMVG4, yMDE4))

#win.graph(width=25, height=20)
postscript(paste(WD.PATH, 'results/fig1.eps', sep=''), width=25, height=20)
par(mfrow=c(2,2), mar=c(2,4,2.5,0.5))
plot(x, yMVN1, type='l', lty=1, col=1, xlim=c(-6.5, 6.5), ylim=c(0,a1), main='(a) No Truncation', font.main=2, cex.main=1.5, las=1, ylab='Density')
lines(x, yMVT1, lty=2, col=2, lwd=1.5)
lines(x, yMSL1, lty=3, col=3, lwd=3)
lines(x, yMCN1, lty=4, col=4, lwd=1.5)
lines(x, yMVG1, lty=5, col=5)
lines(x, yMDE1, lty=6, col=6)
legend("topleft", c('Normal', "Student's t",'Slash', 'Contaminated Normal', 'Variance-gamma', 'Double Exponential'), lty=1:6, col=1:6, bty='n', lwd=c(1, 1.2, 2, 1.2, 1, 1), cex=1)

plot(x, yMVN2, type='l', lty=1, col=1, xlim=c(-7, 1.2), ylim=c(0,a2), main='(b) Right Truncation', font.main=2, cex.main=1.5, las=1, ylab='Density')
lines(x, yMVT2, lty=2, col=2, lwd=1.5)
lines(x, yMSL2, lty=3, col=3, lwd=3)
lines(x, yMCN2, lty=4, col=4, lwd=1.5)
lines(x, yMVG2, lty=5, col=5)
lines(x, yMDE2, lty=6, col=6)
legend("topleft", c('Normal', "Student's t",'Slash', 'Contaminated Normal', 'Variance-gamma', 'Double Exponential'), lty=1:6, col=1:6, bty='n', lwd=c(1, 1.2, 2, 1.2, 1, 1), cex=1)

plot(x, yMVN3, type='l', lty=1, col=1, xlim=c(-1.2, 7), ylim=c(0,a3), main='(c) Left Truncation', font.main=2, cex.main=1.5, las=1, ylab='Density')
lines(x, yMVT3, lty=2, col=2, lwd=1.5)
lines(x, yMSL3, lty=3, col=3, lwd=3)
lines(x, yMCN3, lty=4, col=4, lwd=1.5)
lines(x, yMVG3, lty=5, col=5)
lines(x, yMDE3, lty=6, col=6)
legend("topright", c('Normal', "Student's t",'Slash', 'Contaminated Normal', 'Variance-gamma', 'Double Exponential'), lty=1:6, col=1:6, bty='n', lwd=c(1, 1.2, 2, 1.2, 1, 1), cex=1)

plot(x, yMVN4, type='l', lty=1, col=1, xlim=c(-1.2, 1.2), ylim=c(0,a4), main='(d) Double Truncation', font.main=2, cex.main=1.5, las=1, ylab='Density')
lines(x, yMVT4, lty=2, col=2, lwd=1.5)
lines(x, yMSL4, lty=3, col=3, lwd=3)
lines(x, yMCN4, lty=4, col=4, lwd=1.5)
lines(x, yMVG4, lty=5, col=5)
lines(x, yMDE4, lty=6, col=6)
legend("topright", c('Normal', "Student's t",'Slash', 'Contaminated Normal', 'Variance-gamma', 'Double Exponential'), lty=1:6, col=1:6, bty='n', lwd=c(1, 1.2, 2, 1.2, 1, 1), cex=1)
dev.off()
