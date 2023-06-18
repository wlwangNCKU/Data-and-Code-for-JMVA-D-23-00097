# PDF of multivariate normal-independent distribution
dmni = function(x, mu, Sigma=diag(length(mu)), nu=NULL, rho=NULL, alpha=NULL, beta=NULL, distr=c('MVN','MVT','MSL','MCN','MVG'))
{
  IGfn = function(a, b) pgamma(1, shape=a, rate=b) * gamma(a) / b^a
  distr = distr[1]
  p = length(mu)
  if(distr=='MVN'){
    if(p==1){
      den = dnorm(x, mu, sd=sqrt(Sigma))
    } else den = dmvnorm(x, mean=mu, sigma=round(Sigma,5), log=F)
  }
  if(distr=='MVT'){
    if(p==1){
      den = dt((x-mu)/sqrt(Sigma), df=nu)/sqrt(Sigma)
    } else den = dmvt(x, delta=mu, sigma=Sigma, df=nu, log=F)
  }
  if(distr=='MSL'){
    if(p==1){
      den = c(nu / ((2*pi)^(p/2)*sqrt(Sigma))) * IGfn(a=p/2+nu, b=c(((x-mu)^2/c(Sigma))/2))
      idx = which(x==mu)
      den[idx] = c(2*nu / ((2*pi)^(p/2)*(p+2*nu)*sqrt(Sigma)))
    } else{
      den = nu / ((2*pi)^(p/2)*sqrt(det(Sigma))) * IGfn(a=p/2+nu, b=c(t(x-mu)%*%solve(Sigma)%*%(x-mu))/2)
      if(sum(x==mu)==p) den = 2*nu / ((2*pi)^(p/2)*(p+2*nu)*sqrt(det(Sigma)))
    }}
  if(distr=='MCN'){
    if(p==1){
      Delta = c((x-mu)^2/c(Sigma))
      den = c(1/((2*pi)^(p/2)*sqrt(Sigma))) * (nu*rho^(p/2)*exp(-rho*Delta/2) + (1-nu)*exp(-Delta/2))
    } else{
      Delta = c(t(x-mu)%*%solve(Sigma)%*%(x-mu))
      den = 1/((2*pi)^(p/2)*sqrt(det(Sigma))) * (nu*rho^(p/2)*exp(-rho*Delta/2) + (1-nu)*exp(-Delta/2))
    }}
  if(distr=='MVG'){
    if(p==1){
      Delta = c((x-mu)^2/c(Sigma))
      den = c(2/((2*pi)^(p/2)*sqrt(Sigma)*gamma(alpha)*2^alpha)) * (sqrt(2/beta))^(p/2+alpha) * Delta^((alpha-p/2)/2) * BesselK(z=sqrt(2*Delta/beta), nu=(alpha-p/2))
      idx = which(x==mu)
      den[idx] = c(gamma(alpha-p/2)/((2*pi)^(p/2)*sqrt(Sigma)*gamma(alpha)*beta^(p/2)))
    } else{
      Delta = c(t(x-mu)%*%solve(Sigma)%*%(x-mu))
      den = 2/((2*pi)^(p/2)*sqrt(det(Sigma))*gamma(alpha)*2^alpha) * (sqrt(2/beta))^(p/2+alpha) * Delta^((alpha-p/2)/2) * BesselK(z=sqrt(2*Delta/beta), nu=(alpha-p/2))
      if(sum(x==mu)==p){
        if((alpha-p/2)==0){
          den = gamma(1e-5)/((2*pi)^(p/2)*sqrt(det(Sigma))*gamma(alpha)*beta^(p/2))
        } else den = gamma(alpha-p/2)/((2*pi)^(p/2)*sqrt(det(Sigma))*gamma(alpha)*beta^(p/2))
      }
    }}
  return(den)
}

# CDF of multivariate normal-independent distribution
pmni = function(lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)), mu, Sigma=diag(length(mu)), nu=NULL, rho=NULL, alpha=NULL, beta=NULL, distr=c('MVN','MVT','MSL','MCN','MVG'))
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  distr = distr[1]
  p = length(mu)
  if(distr=='MVN'){
    if(p==1){
      Sd = sqrt(Sigma)
      cdf = pnorm(upper, mean=mu, sd=Sd) - pnorm(lower, mean=mu, sd=Sd)
    } else cdf = pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma, 4), algorithm = GB)[1]
  }
  if(distr=='MVT'){
    if(p==1){
      Sd = sqrt(Sigma)
      cdf = pt((upper-mu)/Sd, df=nu) - pt((lower-mu)/Sd, df=nu)
    } else cdf = pmvt(lower=lower, upper=upper, delta=mu, df=nu, sigma=Sigma)[1]
  }
  if(distr=='MSL'){
    if(p==1) cdf = integrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')$value
    else cdf = cubintegrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, distr='MSL')$integral
  }
  if(distr=='MCN'){
    if(p==1) cdf = integrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, rho=rho, distr='MCN')$value
    cdf = cubintegrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, rho=rho, distr='MCN')$integral
  }
  if(distr=='MVG'){
    if(p==1){
      cdf = integrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, alpha=alpha, beta=beta, distr='MVG')$value
    } else{
      cdf = cubintegrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, alpha=alpha, beta=beta, distr='MVG')$integral
    }}
  return(cdf)
}
