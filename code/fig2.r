library(mvtnorm)
library(tmvtnorm)
library(cubature)
library(Bessel)
library(plot3D)

source(paste(WD.PATH, 'function/mni.r', sep=''))
source(paste(WD.PATH, 'function/tmni.r', sep=''))
source(paste(WD.PATH, 'function/TMNImoment.r', sep=''))

# True parameters
mu = c(0,0)
Sigma = 0.5*matrix(c(2,1,1,4), 2,2)
nu = 4            # MVT
nu1 = 2           # MSL
alpha1 = 2; beta1 = 1 # MVG
alpha2 = 1; beta2 = 2 # MDE

# Truncations
ubd = 1
a1.low=rep(-Inf, 2); a1.upp=rep(Inf, 2)
a2.low=c(-Inf,-Inf); a2.upp=c(ubd,ubd)
a3.low=c(-ubd,-ubd); a3.upp=c(Inf,Inf)
a4.low=c(-ubd,-ubd); a4.upp=c(ubd,ubd)

x1 = x2 = seq(-10, 10, length=200)
m = length(x1)
yMVN1 = yMVN2 = yMVN3 = yMVN4 = 
yMVT1 = yMVT2 = yMVT3 = yMVT4 = matrix(0, m, m)
for(i in 1: m)for(j in 1: m){
  x12 = c(x1[i],x2[j])
# MVN
  yMVN1[i, j] = dtmvnorm(x12, mean=mu, sigma=Sigma)
  yMVN2[i, j] = dtmvnorm(x12, mean=mu, sigma=Sigma, lower=a2.low, upper=a2.upp)
  yMVN3[i, j] = dtmvnorm(x12, mean=mu, sigma=Sigma, lower=a3.low, upper=a3.upp)
  yMVN4[i, j] = dtmvnorm(x12, mean=mu, sigma=Sigma, lower=a4.low, upper=a4.upp)
# MVT
  yMVT1[i, j] = dmni(x12, mu, Sigma, nu, distr='MVT')
  yMVT2[i, j] = dtmni(x12, mu, Sigma, nu, distr='MVT', a.low=a2.low, a.upp=a2.upp)
  yMVT3[i, j] = dtmni(x12, mu, Sigma, nu, distr='MVT', a.low=a3.low, a.upp=a3.upp)
  yMVT4[i, j] = dtmni(x12, mu, Sigma, nu, distr='MVT', a.low=a4.low, a.upp=a4.upp)
}

# MSL
yMSL1 = yMSL2 = yMSL3 = yMSL4 = matrix(0, m, m)
for(i in 1: m)for(j in 1: m){
  x12 = c(x1[i],x2[j])
  yMSL1[i, j] = dmni(x12, mu, Sigma, nu1, distr='MSL')
  yMSL2[i, j] = dtmni(x12, mu, Sigma, nu1, distr='MSL', a.low=a2.low, a.upp=a2.upp)
  yMSL3[i, j] = dtmni(x12, mu, Sigma, nu1, distr='MSL', a.low=a3.low, a.upp=a3.upp)
  yMSL4[i, j] = dtmni(x12, mu, Sigma, nu1, distr='MSL', a.low=a4.low, a.upp=a4.upp)
}

# MCN
yMCN1 = yMCN2 = yMCN3 = yMCN4 = matrix(0, m, m)
for(i in 1: m)for(j in 1: m){
  x12 = c(x1[i],x2[j])
  yMCN1[i, j] = dmni(x12, mu, Sigma, nu=0.25, rho=0.2, distr='MCN')
  yMCN2[i, j] = dtmni(x12, mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a2.low, a.upp=a2.upp)
  yMCN3[i, j] = dtmni(x12, mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a3.low, a.upp=a3.upp)
  yMCN4[i, j] = dtmni(x12, mu, Sigma, nu=0.25, rho=0.2, distr='MCN', a.low=a4.low, a.upp=a4.upp)
}

# MVG
yMVG1 = yMVG2 = yMVG3 = yMVG4 = matrix(0, m, m)
for(i in 1: m)for(j in 1: m){
  x12 = c(x1[i],x2[j])
  yMVG1[i, j] = dmni(x12, mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG')
  yMVG2[i, j] = dtmni(x12, mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG', a.low=a2.low, a.upp=a2.upp)
  yMVG3[i, j] = dtmni(x12, mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG', a.low=a3.low, a.upp=a3.upp)
  yMVG4[i, j] = dtmni(x12, mu, Sigma, alpha=alpha1, beta=beta1, distr='MVG', a.low=a4.low, a.upp=a4.upp)
}

# MDE
yMDE1 = yMDE2 = yMDE3 = yMDE4 = matrix(0, m, m)
for(i in 1: m)for(j in 1: m){
  x12 = c(x1[i],x2[j])
  yMDE1[i, j] = dmni(x12, mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG')
  yMDE2[i, j] = dtmni(x12, mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG', a.low=a2.low, a.upp=a2.upp)
  yMDE3[i, j] = dtmni(x12, mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG', a.low=a3.low, a.upp=a3.upp)
  yMDE4[i, j] = dtmni(x12, mu, Sigma, alpha=alpha2, beta=beta2, distr='MVG', a.low=a4.low, a.upp=a4.upp)
}
save.image(paste(WD.PATH, 'data/dMNI.RData', sep=''))

load(paste(WD.PATH, 'data/dMNI.RData', sep=''))
#win.graph(width=20, height=50)
postscript(paste(WD.PATH, 'results/fig2.eps', sep=''), width=7, height=60)
a <- seq(0, 0.01, length = 5)
b1 <- seq(0.01, max(yMDE1), length = 30)
b2 <- seq(0.01, max(yMDE2), length = 30)
b3 <- seq(0.01, max(yMDE3), length = 30)
b4 <- seq(0.01, max(yMDE4), length = 30)

layout(matrix(c(0,1:34), 7, 5, byrow=TRUE), widths = c(1.5, rep(5,4)), heights = c(2, rep(5,6)))
par(mar=c(2,2,0,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(paste(A[0],'=(',-infinity,',',infinity,')')), 1, line=0, cex=1, font=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(paste(A[1],'=(',-infinity,',',1,')')), 1, line=0, cex=1, font=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(paste(A[2],'=(',-1,',',infinity,')')), 1, line=0, cex=1, font=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(paste(A[3],'=(',-1,',',1,')')), 1, line=0, cex=1, font=2)

# MVN
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MVN', 2, line=-1, cex=1.2, font=2)
contour2D(yMVN1, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b1), xlim = c(-5, 5), ylim = c(-5, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVN2, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b2), xlim = c(-5, 1.2), ylim = c(-5, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVN3, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b3), xlim = c(-1.2, 5), ylim = c(-1.2, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVN4, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b4), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.8, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)

# MVT
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MVT', 2, line=-1, cex=1.2, font=2)
contour2D(yMVT1, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b1), xlim = c(-5, 5), ylim = c(-5, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVT2, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b2), xlim = c(-5, 1.2), ylim = c(-5, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVT3, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b3), xlim = c(-1.2, 5), ylim = c(-1.2, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVT4, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b4), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.8, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)

# MSL
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MSL', 2, line=-1, cex=1.2, font=2)
contour2D(yMSL1, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b1), xlim = c(-5, 5), ylim = c(-5, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMSL2, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b2), xlim = c(-5, 1.2), ylim = c(-5, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMSL3, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b3), xlim = c(-1.2, 5), ylim = c(-1.2, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMSL4, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b4), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.8, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)

# MCN
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MCN', 2, line=-1, cex=1.2, font=2)
contour2D(yMCN1, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b1), xlim = c(-5, 5), ylim = c(-5, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMCN2, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b2), xlim = c(-5, 1.2), ylim = c(-5, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMCN3, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b3), xlim = c(-1.2, 5), ylim = c(-1.2, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMCN4, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b4), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.8, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)

# MVG
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MVG', 2, line=-1, cex=1.2, font=2)
contour2D(yMVG1, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b1), xlim = c(-5, 5), ylim = c(-5, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVG2, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b2), xlim = c(-5, 1.2), ylim = c(-5, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVG3, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b3), xlim = c(-1.2, 5), ylim = c(-1.2, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMVG4, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b4), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.8, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)

# MDE
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MDE', 2, line=-1, cex=1.2, font=2)
contour2D(yMDE1, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b1), xlim = c(-5, 5), ylim = c(-5, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMDE2, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b2), xlim = c(-5, 1.2), ylim = c(-5, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMDE3, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b3), xlim = c(-1.2, 5), ylim = c(-1.2, 5), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.9, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
contour2D(yMDE4, x1, x2, lty=1, drawlabels=FALSE , levels=c(a,b4), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), colkey = FALSE, las=1, 
           xlab = "", ylab = "", cex.axis = 0.8, tcl = -0.5, mgp = c(3, 0.8, 0), lwd = 0.8)
dev.off()
