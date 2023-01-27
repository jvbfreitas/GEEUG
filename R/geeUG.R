geeUG = function(formula.mu, formula.phi, data, id, tol = 0.001,
                 maxiter = 25, corstr = "independence", linkmu = "logit",
                 linkphi = "log", silence = FALSE){

  if(all(c("independence", "unstructured", "exchangeable", "AR-1") != corstr)){
    stop("the correlation structure is not defined")
  }
  if(all(c("logit", "identity","cloglog","probit") != linkmu)){
    stop("the link function is not defined")
  }
  data = data
  formula.mu = as.formula(formula.mu)
  nformula = all.vars(formula.mu)
  fnames = 0
  jaux = 1
  listaux = list(NULL)
  for(i in 1:length(nformula)){
    if(is.factor(data[,nformula[i]])){
      fnames[jaux] = nformula[i]
      listaux[[jaux]] = "contr.treatment"
      jaux = jaux+1
    }
  }
  call = match.call()
  if(jaux>1){
    names(listaux) = fnames
    X = as.matrix(model.matrix(formula.mu, data = data, contrasts = listaux)) # Matriz de especificação
  }else{
    X = as.matrix(model.matrix(formula.mu, data = data)) # Matriz de especificação
  }
  formula.phi = as.formula(formula.phi)
  nformula = all.vars(formula.phi)
  fnames = 0
  jaux = 1
  listaux = list(NULL)
  for(i in 1:length(nformula)){
    if(is.factor(data[,nformula[i]])){
      fnames[jaux] = nformula[i]
      listaux[[jaux]] = "contr.treatment"
      jaux = jaux+1
    }
  }
  call = match.call()
  if(jaux>1){
    names(listaux) = fnames
    Z = as.matrix(model.matrix(formula.phi, data = data, contrasts = listaux)) # Matriz de especificação
  }else{
    Z = as.matrix(model.matrix(formula.phi, data = data)) # Matriz de especificação
  }
  qq = ncol(Z)
  p = ncol(X) # Número de parâmetros
  y = model.frame(formula.mu, data = data)[,1] # Variável resposta
  t = as.vector(table(id)) # Número de repetições
  n = length(table(id)) # Número de unidades experimentais
  N = nrow(X)
  warn = getOption("warn")
  options(warn=-1)
  if(linkmu == "identity"){
    if(jaux>1){
      if(linkphi == "log"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "identity",
                                                                         sigma.link = "log"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "identity"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "identity",
                                                                         sigma.link = "identity"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "invsquare"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "identity",
                                                                         sigma.link = "1/mu^2"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "sqrt"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "identity",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
    }else{
      if(linkphi == "log"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "identity",
                                                                         sigma.link = "log"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "identity"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "identity",
                                                                         sigma.link = "identity"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "invsquare"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "identity",
                                                                         sigma.link = "1/mu^2"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "sqrt"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "identity",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data)
      }
    }
  }
  if(linkmu == "logit"){
    if(jaux>1){
      if(linkphi == "log"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "logit",
                                                                         sigma.link = "log"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "identity"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "logit",
                                                                         sigma.link = "identity"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "invsquare"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "logit",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "sqrt"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "logit",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
    }else{
      if(linkphi == "log"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "logit",
                                                                         sigma.link = "log"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "identity"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "logit",
                                                                         sigma.link = "identity"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "invsquare"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "logit",
                                                                         sigma.link = "1/mu^2"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "sqrt"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "logit",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data)
      }
    }
  }
  if(linkmu == "cloglog"){
    if(jaux>1){
      if(linkphi == "log"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "cloglog",
                                                                         sigma.link = "log"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "identity"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "cloglog",
                                                                         sigma.link = "identity"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "invsquare"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "cloglog",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "sqrt"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "cloglog",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
    }else{
      if(linkphi == "log"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "cloglog",
                                                                         sigma.link = "log"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "identity"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "cloglog",
                                                                         sigma.link = "identity"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "invsquare"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "cloglog",
                                                                         sigma.link = "1/mu^2"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "sqrt"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "cloglog",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data)
      }
    }
  }
  if(linkmu == "probit"){
    if(jaux>1){
      if(linkphi == "log"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "probit",
                                                                         sigma.link = "log"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "identity"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "probit",
                                                                         sigma.link = "identity"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "invsquare"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "probit",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
      if(linkphi == "sqrt"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "probit",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data,contrasts = listaux)
      }
    }else{
      if(linkphi == "log"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "probit",
                                                                         sigma.link = "log"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "identity"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "probit",
                                                                         sigma.link = "identity"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "invsquare"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "probit",
                                                                         sigma.link = "1/mu^2"),
                      trace = FALSE, data = data)
      }
      if(linkphi == "sqrt"){
        mod0 = gamlss(formula.mu,sigma.formula = formula.phi,family = UG(mu.link = "probit",
                                                                         sigma.link = "sqrt"),
                      trace = FALSE, data = data)
      }
    }
  }
  options(warn=warn)
  beta = mod0$mu.coefficients # Chute inicial para beta
  eta = X%*%beta
  if(linkmu == "identity"){
    G = diag(1,N,N) # G para a ligação logarítimica
  }
  if(linkmu == "logit"){
    mu = as.vector(exp(eta)/(1+exp(eta)))
    G = diag(as.vector(exp(eta)/((1+exp(eta))^2)),N,N) # G para a ligação logarítimica
  }
  if(linkmu == "cloglog"){
    mu = as.vector(1-exp(-exp(eta)))
    G = diag(as.vector(exp(eta-exp(eta))),N,N) # G para a ligação logarítimica
  }
  if(linkmu == "probit"){
    mu = as.vector(pnorm(eta))
    G = diag(as.vector(dnorm(eta)),N,N) # G para a ligação logarítimica
  }
  ni = mod0$sigma.coefficients
  delta = Z%*%ni
  cont = 1
  contmax = 0
  Q = as.matrix(bdiag(X,Z))
  repeat{
    if(silence == FALSE){
      print(cont)
      print(beta)
      print(ni)
    }
    eta = X%*%beta
    if(linkmu == "identity"){
      mu = as.vector(eta) # mi para a ligação logarítmica
    }
    if(linkmu == "logit"){
      mu = as.vector(exp(eta)/(1+exp(eta)))
    }
    if(linkmu == "cloglog"){
      mu = as.vector(1-exp(-exp(eta)))
    }
    if(linkmu == "probit"){
      mu = as.vector(pnorm(eta))
    }
    delta = Z%*%ni
    if (linkphi=="log") 		phi = as.vector(exp(delta))
    if (linkphi=="identity") 	phi = as.vector(delta)
    if (linkphi=="invsquare") 	phi = as.vector(1/delta^2)
    if (linkphi=="sqrt") 	phi = as.vector(delta^2)
    #Calculo do d
    a = mu^(1/phi)
    d = a/(1-a)

    # y estrela
    ys = log(y)
    yss = log(-log(y))

    # mi estrela
    mus = -phi/d
    muss = digamma(phi)-log(d)

    # Variância de b_ij
    vmus = d*(1+d)/(mu^(2+1/phi)*phi)
    vphis = (d/(phi^2))*log(mu)*(1+d)*(1/d+log(mu)/phi+log(mu)/(d*phi))+trigamma(phi)+log(mu)*(1+d)/(phi^2)

    # Vetor u_i
    u = d*(1+d)*(ys-mus)/(mu*phi)

    # Vetor v_i
    v = (yss-muss)-d*(1+d)*log(mu)*(ys-mus)/(phi^2)

    #Matrizes utilizadas para o cálculo da equação de estimação
    if(linkmu == "identity"){
      G = diag(1,N,N) # G para a ligação identidade
    }
    if(linkmu == "logit"){
      G = diag(as.vector(exp(eta)/((1+exp(eta))^2)),N,N) # G para a ligação logit
    }
    if(linkmu == "cloglog"){
      G = diag(as.vector(exp(eta-exp(eta))),N,N) # G para a ligação cloglog
    }
    if(linkmu == "probit"){
      G = diag(as.vector(dnorm(eta)),N,N) # G para a ligação probit
    }
    if (linkphi=="log") 		Fi = diag(phi)
    if (linkphi=="identity") 	Fi = diag(N)
    if (linkphi=="invsquare") 	Fi = diag(as.vector(-2/delta^3),N)
    if (linkphi=="sqrt") 	Fi = diag(as.vector(2*delta),N)
    A = diag(as.vector(vmus))
    M = diag(as.vector(vphis))
    c = -d*(1+d)*(1/d+log(mu)/phi+log(mu)/(phi*d))/(mu*phi)
    C = diag(c)
    Lambda1 = G%*%A
    Lambda2 = Fi%*%M
    Lambda12 = Fi%*%C
    Lambda21 = G%*%C
    eu2 = split(d*(1+d)/(mu^(2+1/phi)*phi),id)

    uc = split(u,id)
    if(corstr == "unstructured"){
      Rg = matrix(0,max(t),max(t))
      cnum = den1 = den2 = 0
      for(j in 1:(max(t))){
        for(k in j:(max(t))){
          for(i in 1:n){
            if(is.na(uc[[i]][j])||is.na(uc[[i]][k])){
              cnum = cnum
            }
            else{
              cnum = cnum + (uc[[i]][j])*(uc[[i]][k])
              den1 = den1 + (uc[[i]][j])^2
              den2 = den2 + (uc[[i]][k])^2
            }
          }
          Rg[j,k] = cnum/(sqrt(den1)*sqrt(den2))
          Rg[k,j] = Rg[j,k]
        }
      }
      diag(Rg) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = Rg[1:t[i],1:t[i]]
      }
      Rm = bdiag(R)
    }
    if(corstr == "AR-1"){
      cnum = den1 = den2 = 0
      for(i in 1:n){
        if(t[i]==1){
          cnum = cnum
          den1 = den1
          den2 = den2
        }else{
          for(j in 1:(max(t)-1)){
            cnum = cnum + uc[[i]][j]*uc[[i]][j+1]
            den1 = den1 + uc[[i]][j]^2
            den2 = den2 + uc[[i]][j+1]^2
          }
        }
      }
      alpha = cnum/sqrt(den1*den2)
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(0,t[i],t[i])
        for(j in 1:t[i]){
          for(l in 1:t[i]){
            R[[i]][j,l] = alpha^(abs(j-l))
          }
        }
      }
      # Matriz de correlação AR-1
      Rm = as.matrix(bdiag(R))
      R=R[[1]]
    }
    if(corstr == "exchangeable"){
      cnum = den =  0
      for(i in 1:n){
        for(j in 1:(t[i])){
          for(k in 1:t[i]){
            if(is.na(uc[[i]][j])||is.na(uc[[i]][k])){
              cnum = cnum
              den = den
            }else{
              if(j>k){
                cnum = cnum + (2/(t[i]*(t[i]-1)))*uc[[i]][j]*uc[[i]][k]
              }else{
                cnum = cnum
              }
            }
          }
          den = den + (1/t[i])*uc[[i]][j]^2
        }
      }
      alpha = cnum/den
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(0,t[i],t[i])
        for(j in 1:t[i]){
          for(l in 1:t[i]){
            R[[i]][j,l] = alpha
          }
        }
        diag(R[[i]]) = 1
      }
      # Matriz de correlação AR-1
      Rm = as.matrix(bdiag(R))
      R=R[[1]]
      diag(Rm) = 1
    }


    vc = split(v,id)
    Rg = matrix(0,max(t),max(t))
    cnum = den1 = den2 = 0
    for(j in 1:(max(t))){
      for(k in j:(max(t))){
        for(i in 1:n){
          if(is.na(vc[[i]][j])||is.na(vc[[i]][k])){
            cnum = cnum
          }
          else{
            cnum = cnum + (vc[[i]][j])*(vc[[i]][k])
            den1 = den1 + (vc[[i]][j])^2
            den2 = den2 + (vc[[i]][k])^2
          }
        }
        Rg[j,k] = cnum/(sqrt(den1)*sqrt(den2))
        Rg[k,j] = Rg[j,k]
      }
    }
    diag(Rg) = 1
    R2 = list(NULL)
    for(i in 1:n){
      R2[[i]] = Rg[1:t[i],1:t[i]]
    }
    Rm2 = bdiag(R2)

    sa = diag(sqrt(vmus),N)
    sm = diag(sqrt(vphis),N)
    Omega = sa%*%Rm%*%sa
    Omega2 = sm%*%Rm2%*%sm
    Lambda = rbind(cbind(Lambda1,Lambda12),cbind(Lambda21,Lambda2))
    indaux = cumsum(c(0,t))
    sN = bdiag(solve(Omega), solve(Omega2))
    W = Lambda%*%sN%*%t(Lambda)
    z = c(eta,delta) + solve(Lambda)%*%c(u,v)

    thetan = solve(as.matrix(t(Q)%*%W%*%Q))%*%(t(Q)%*%W%*%z)
    beta1 = thetan[1:p]
    ni1 = thetan[(p+1):(p+qq)]
    dif = sum(abs(beta1-beta)) + sum(abs(ni1-ni))
    if(dif<=(2*tol)){
      beta = beta1
      ni = ni1
      cat("The algorithm converged")
      converg = 1
      break
    }

    if(cont == maxiter){
      cat("Maximum number of iterations reached")
      beta = beta1
      ni = ni1
      contmax = 1
      converg = 0
      break
    }
    beta = beta1
    ni = ni1
    cont = cont + 1
  }
  S = -t(Q)%*%W%*%Q
  invOmega = solve(Omega)
  VarBeta = solve(S)%*%t(Q)%*%Lambda%*%sN%*%c(u,v)%*%t(c(u,v))%*%sN%*%t(Lambda)%*%Q%*%solve(S)
  SEbeta = sqrt(VarBeta[col(VarBeta)==row(VarBeta)])

  fit = list()
  fit$call = call
  class(fit) = "geeUG"
  fit$title = "geeUG:  UNIT GAMMA GENERALIZED ESTIMATING EQUATIONS"
  fit$model = list()
  fit$model$link.mu = linkmu
  fit$model$link.phi = linkphi
  fit$model$corstr = corstr
  fit$call = call
  fit$formula.mu = formula.mu
  fit$formula.phi = formula.phi
  fit$nclusters = n
  fit$clusters = t
  fit$nobs = N
  fit$contmax = contmax
  fit$iterations = cont
  fit$coefficients = c(beta,ni)
  fit$mu.coefficients = beta
  fit$phi.coefficients = ni
  eta = as.vector(X %*% fit$mu.coefficients)
  delta = as.vector(Z %*% fit$phi.coefficients)
  fit$mu.linear.predictors = eta
  fit$phi.linear.predictors = delta
  mu = as.vector(mu)
  phi = as.vector(phi)
  fit$mu.fitted.values = mu
  fit$phi.fitted.values = phi
  fit$family = "Unit Gamma"
  fit$y = as.vector(y)
  fit$id = as.vector(id)
  fit$max.id = max(t)
  fit$working.correlation = R
  fit$robust.variance = VarBeta
  fit$robust.se = SEbeta
  if(corstr == "unstructured"){
    fit$alpha = fit$working.correlation[upper.tri(fit$working.correlation)]
  }
  if(corstr == "AR-1"||corstr == "exchangeable"){
    fit$alpha = alpha
  }
  fit$comp$X = X
  fit$comp$Z = Z
  fit$comp$Q = Q
  fit$comp$W = W
  fit$comp$u = u
  fit$comp$v = v
  fit$comp$Lambda = Lambda
  fit$comp$G = G
  fit$comp$Rm = Rm
  fit$comp$Rm2 = Rm2
  fit$comp$sa = sa
  fit$comp$sm = sm
  fit$comp$Omega = Omega
  fit$comp$Omega2 = Omega2
  fit$comp$M = M
  fit$comp$tol = tol
  fit$comp$maxiter = maxiter
  resq = 0
  for(i in 1:N){
    resq[i] = qnorm(pUG(y[i],mu[i],phi[i]))
  }
  fit$residuals = resq
  auxsumm = summary(fit)
  print.summary.geeUG(auxsumm)
  return(fit)
}

print.geeUG = function(x, digits = NULL, quote = FALSE, prefix = "", ...){
  if(is.null(digits)) digits = options()$digits else options(digits =
                                                               digits)
  cat("\n", x$title)
  cat("n"," Model:\n")
  cat(" Link (mean):                     ", x$model$link.mu, "\n")
  cat(" Link (precision):                     ", x$model$link.phi, "\n")
  cat(" Correlation Structure:    ", x$model$corstr, "\n")
  cat("\nCall:\n")
  dput(x$call)
  ys = matrix(rep(matrix(x$id, ncol = 1), 5), ncol = 5)
  ys[, 2] = x$y
  ys[, 3] = x$mu.linear.predictors
  ys[, 4] = x$mu.fitted.values
  ys[, 5] = x$residuals
  dimnames(ys) = list(1:length(x$y), c("ID", "Y", "LP", "fitted",
                                       "Residual")) #       cat("\nFitted Values:\n")
  cat("\nNumber of observations : ", x$nobs, "\n")
  cat("\nMaximum cluster size   : ", x$max.id, "\n")
  cat("\n\nCoefficients (mean):\n")
  print(t(x$mu.coefficients), digits = digits)
  cat("\n\nCoefficients (precision):\n")
  print(t(x$phi.coefficients), digits = digits)
  cat("\nNumber of Iterations: ", x$iterations)
  invisible(x)
}

summary.geeUG = function(object,...){
  value = list()
  class(value) = "summary.geeUG"
  value$call = object$call
  value$corstr = object$model$corstr
  value$link.mu = object$model$link.mu
  value$link.phi = object$model$link.phi
  value$alpha = object$alpha
  value$nclusters = object$nclusters
  value$max.id = object$max.id
  valorpmu = 0
  waldmu = 0
  for(i in 1:length(object$mu.coefficients)){
    waldmu[i] = wald.test(object$robust.variance, object$coefficients, Terms = i)$result$chi2[1]
    valorpmu[i] = wald.test(object$robust.variance, object$coefficients, Terms = i)$result$chi2[3]
  }
  valorpphi = 0
  waldphi = 0
  for(i in (length(object$mu.coefficients)+1):length(object$coefficients)){
    waldphi[i-length(object$mu.coefficients)] = wald.test(object$robust.variance, object$coefficients, Terms = i)$result$chi2[1]
    valorpphi[i-length(object$mu.coefficients)] = wald.test(object$robust.variance, object$coefficients, Terms = i)$result$chi2[3]
  }
  value$coefficients.mu = data.frame(object$mu.coefficients,
                                     object$robust.se[1:length(object$mu.coefficients)],
                                     waldmu, valorpmu)
  value$coefficients.phi = data.frame(object$phi.coefficients,
                                      object$robust.se[(length(object$mu.coefficients) + 1):length(object$coefficients)],
                                      waldphi, valorpphi)
  colnames(value$coefficients.mu) = c("Estimate","Robust.Std.err", "Wald", "Pr(>|W|)")
  colnames(value$coefficients.phi) = c("Estimate","Robust.Std.err", "Wald", "Pr(>|W|)")
  return(value)
}

print.summary.geeUG = function(x, digits = NULL, quote = FALSE, prefix = "", ...){
  if(is.null(digits)) digits = options()$digits else options(digits = digits)
  cat("\nCall:\n")
  dput(x$call)
  cat("\nCoefficients (mean):\n")
  printCoefmat(as.matrix(x$coefficients.mu), digits = digits)
  cat("\nCoefficients (precision):\n")
  printCoefmat(as.matrix(x$coefficients.phi), digits = digits)
  cat("\nNumber of clusters : ", x$nclusters, "\n")
  cat("\nMaximum cluster size   : ", x$max.id, "\n")
  invisible(x)
}

plot.geeUG = function(x,...){
  rp = x$residuals
  qplot(1:length(rp),rp) + ylab("Quantile residuals") + xlab("Index")
}
