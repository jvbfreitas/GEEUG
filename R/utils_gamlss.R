dUG<-function(y, mu=0.5, sigma=1, log = FALSE){
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n",""))
  a = mu^(1/sigma)
  d = a/(1-a)
  fy = (d^(sigma)/gamma(sigma))*y^(d-1)*log(1/y)^(sigma-1)
  if(log){return(log(fy))}else{return(fy)}
}
library(cubature)
pUG<- function(q, mu=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("y must be between 0 and 1", "\n",""))
  value = 0
  for(i in 1:length(q)){
    if(lower.tail){value[i] = integrate(function(x) dUG(x,mu[i],sigma[i]), 0, q[i])$value}else{value[i] = 1-integrate(function(x) dUG(x,mu[i],sigma[i]), 0, q[i])$value}
  }
  if(log.p){return(log(value))}else{return(value)}
}

qUG<-function(p, mu=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if(any(p <= 0)|any(p >= 1)) stop(paste("p must be between 0 and 1","\n",""))
  if(log.p){p = exp(p)}
  value = 0
  for(i in 1:length(p)){
    if(p[i]==1){
      value[i] = 0.99999
    }
    if(p[i]==0){
      value[i] = 0.00001
    }
    if(p[i]>0 & p[i]<1){
      if(lower.tail){value[i] = uniroot(function(x) pUG(x,mu[i],sigma[i])-p[i], lower = .Machine$double.eps,
                                        upper = 1-.Machine$double.eps)$root}else{
                                          value[i] = uniroot(function(x) (1-pUG(x,mu[i],sigma[i]))-p[i], lower = .Machine$double.eps,
                                                             upper = 1-.Machine$double.eps)$root
                                        }
    }
  }
  return(value)
}

rUG<-function(n,mu=0.5,sigma=1){
  if(any(sigma <= 0)){stop(paste("sigma must be positive","\n"," "))}
  if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if (any(n <= 0))
    stop(paste("n must be a positive integer", "\n",""))
  r<-runif(n)
  value = qUG(r,mu,sigma)
  return(value)
}

UG<-function (mu.link = "logit", sigma.link = "log")
{
  mstats <- checklink("mu.link", "unitg", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "unitg", substitute(sigma.link),
                      c("log", "sqrt", "1/mu^2", "identity"))
  structure(list(family = c("BET", "unitg"), parameters = list(mu = TRUE,sigma = TRUE),
                 nopar = 2,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 dldm = function(y,mu, sigma) {
                   a = mu^(1/sigma)
                   d = a/(1-a)
                   ys = log(y)
                   mus = -sigma/d
                   dldm <- d*(1+d)*(ys-mus)/(mu*sigma)
                   dldm
                 },
                 d2ldm2 = function(mu, sigma) {
                   a = mu^(1/sigma)
                   d = a/(1-a)
                   d2ldm2 <- -d*(1+d)/(mu^(2+1/sigma)*sigma)
                   d2ldm2
                 },
                 dldd = function(y, mu, sigma) {
                   a = mu^(1/sigma)
                   d = a/(1-a)
                   ys = log(y)
                   yss = log(-log(y))
                   mus = -sigma/d
                   muss = digamma(sigma)-log(d)
                   dldd <- (yss-muss)-d*(1+d)*log(mu)*(ys-mus)/(sigma^2)
                   dldd
                 },
                 d2ldd2 = function(mu, sigma) {
                   a = mu^(1/sigma)
                   d = a/(1-a)
                   d2ldd2 <- -(d/(sigma^2))*log(mu)*(1+d)*(1/d+log(mu)/sigma+log(mu)/(d*sigma))-trigamma(sigma)-log(mu)*(1+d)/(sigma^2)
                   d2ldd2
                 },
                 d2ldmdd = function(mu, sigma) {
                   a = mu^(1/sigma)
                   d = a/(1-a)
                   d2ldmdd <- (d*(1+d)/(mu*sigma))*(1/d+log(mu)/sigma+log(mu)/(sigma*d))
                   d2ldmdd
                 },
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * dUG(y,mu, sigma, log = TRUE),
                 rqres = expression(rqres(pfun = "pUG", type = "Continuous", y = y, mu = mu, sigma = sigma)),
                 mu.initial = expression({mu <- (y + mean(y))/2}),
                 sigma.initial = expression({sigma <- rep(0.5, length(y))}),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1),
                 mean = function(mu, sigma) mu),
            class = c("gamlss.family", "family"))
}

