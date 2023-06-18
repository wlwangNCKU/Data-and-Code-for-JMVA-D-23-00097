# generate random sample from truncated multivariate normal independent distribution
rtmni = function(n, mu, Sigma=diag(length(mu)), nu=NULL, rho=NULL, alpha=NULL, beta=NULL, distr=c('MVN', 'MVT', 'MSL', 'MCN', 'MVG'), lower=rep(-Inf,length(mu)), upper=rep(Inf, length(mu)))
{
  distr = distr[1]
  p = length(mu)
  if(p == 1){
    Sd = sqrt(Sigma)
  } else{
    Lam = diag(sqrt(diag(Sigma)))
    Lam.inv = diag(1/sqrt(diag(Sigma)))
    R = Lam.inv %*%Sigma%*% Lam.inv
  }
  #  cat('Generate random samples of size', n, 'from the', p,'-variate', distr, 'distribution', '\n')
  Y = matrix(NA, nrow=n, ncol=p)
  Tau = numeric(n)
  iter = 1
  while(iter <= n){
    if(distr=='MVN') tau = 1
    if(distr=='MVT') tau = rgamma(1, shape=nu/2, rate=nu/2)
    if(distr=='MSL') tau = rbeta(1, shape1=nu, shape2=1)
    if(distr=='MCN'){
      a = runif(1)
      tau = ifelse(a<=nu, rho, 1)
    }
    if(distr=='MVG') tau = rinvgamma(1, shape=alpha, scale=beta)
    tau.sq = sqrt(tau)
    if(p==1){
      yi = mu + 1/tau.sq * Sd * rnorm(1, mean=0, sd=1)
    } else  yi = mu + 1/tau.sq * Lam %*% c(rmvnorm(1, mean=rep(0,p), sigma=R))
    if((sum(yi>=lower)+sum(yi<=upper))==(2*p)){
      Y[iter, ] = yi
      Tau[iter] = tau
      iter = iter + 1
    } else next
  }
  return(Y)
}

# PDF of truncated multivariate normal-independent distribution
dtmni = function(x, mu, Sigma=diag(length(mu)), nu=NULL, rho=NULL, alpha=NULL, beta=NULL, distr=c('MVN', 'MVT', 'MSL', 'MCN', 'MVG'), a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{
  distr = distr[1]
  p = length(mu)
  if((sum(x>=a.low)+sum(x<=a.upp))!=(2*p)){
    Tden=0
  } else {
    if(distr=='MVN'){
      cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, distr='MVN')
      Tden = cdf.inv * dmni(x, mu, Sigma, nu, distr='MVN')
    }
    if(distr=='MVT'){
      cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MVT')
      Tden = cdf.inv * dmni(x, mu, Sigma, nu, distr='MVT')
    }
    if(distr=='MSL'){
      cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MSL')
      Tden = cdf.inv * dmni(x, mu, Sigma, nu, distr='MSL')
    }
    if(distr=='MCN'){
      cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu=nu, rho=rho, distr='MCN')
      Tden = cdf.inv * dmni(x, mu, Sigma, nu=nu, rho=rho, distr='MCN')
    }
    if(distr=='MVG'){
      cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, alpha=alpha, beta=beta, distr='MVG')
      Tden = cdf.inv * dmni(x, mu, Sigma, alpha=alpha, beta=beta, distr='MVG')
    }}
  return(Tden)
}

ptmni = function(lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)), mu, Sigma=diag(length(mu)), nu=NULL, rho=NULL, alpha=NULL, beta=NULL, distr=c('MVN', 'MVT', 'MSL', 'MCN', 'MVG'), a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  distr = distr[1]
  p = length(mu)
  if(p == 1) Sd = sqrt(Sigma)
  for(i in 1: p){
    if(lower[i]<a.low[i]) lower[i]=a.low[i]
    if(upper[i]>a.upp[i]) upper[i]=a.upp[i]
  }
  if(distr=='MVN'){
    cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, distr='MVN')
    if(p==1) Tcdf = cdf.inv * (pnorm(upper, mean=mu, sd=Sd) - pnorm(lower, mean=mu, sd=Sd))
    else Tcdf = cdf.inv * pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma, 4), algorithm = GB)[1]
  }
  if(distr=='MVT'){
    cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MVT')
    if(p==1) Tcdf = cdf.inv * (pt((upper-mu)/sqrt(Sigma), df=nu) - pt((lower-mu)/sqrt(Sigma), df=nu))
    else Tcdf = cdf.inv * pmvt(lower=lower, upper=upper, delta=mu, df=nu, sigma=Sigma)[1]
  }
  if(distr=='MSL'){
    cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MSL')
    Tcdf = cdf.inv * cubintegrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')$integral
  }
  if(distr=='MCN'){
    cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu=nu, rho=rho, distr='MCN')
    Tcdf = cdf.inv * cubintegrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, rho=rho, distr='MCN')$integral
  }
  if(distr=='MVG'){
    cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, alpha=alpha, beta=beta, distr='MVG')
    Tcdf = cdf.inv * cubintegrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, alpha=alpha, beta=beta, distr='MVG')$integral
  }
  return(Tcdf)
}
