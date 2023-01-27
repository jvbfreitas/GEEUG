LI.geeUG = function(model, pert = c("case-weight", "response")){
  p = length(model$mu.coefficients)
  qq = length(model$phi.coefficients)
  S <- -t(model$comp$Q)%*%model$comp$W%*%model$comp$Q
  Sbb <- Sgg <- matrix(0,(p+qq),(p+qq))
  Sgg[(p+1):(p+qq),(p+1):(p+qq)] <- as.matrix(solve(S[(p+1):(p+qq),(p+1):(p+qq)]))
  Sbb[1:p,1:p] <- as.matrix(solve(S[1:p,1:p]))
  N = model$nobs
  if(pert == "case-weight"){
    #Local Influence - Case-Weights - Using Sensibility Matrix
    dv = c(model$comp$u,model$comp$v)
    LI <- -diag(dv)%*%solve(model$comp$Lambda)%*%model$comp$W%*%model$comp$Q%*%solve(S)%*%t(model$comp$Q)%*%model$comp$W%*%solve(model$comp$Lambda)%*%diag(dv)
    LIb <- -diag(dv)%*%solve(model$comp$Lambda)%*%model$comp$W%*%model$comp$Q%*%(solve(S)-Sgg)%*%t(model$comp$Q)%*%model$comp$W%*%solve(model$comp$Lambda)%*%diag(dv)
    LIg <- -diag(dv)%*%solve(model$comp$Lambda)%*%model$comp$W%*%model$comp$Q%*%(solve(S)-Sbb)%*%t(model$comp$Q)%*%model$comp$W%*%solve(model$comp$Lambda)%*%diag(dv)
    elib = eigen(LIb, symmetric=T)
    elig = eigen(LIg, symmetric=T)
    strb = sqrt(sum(elib$values^2))
    Di = 0
    for(i in 1:N){
      eij = rep(0,2*N)
      eij[i] = 1
      Di[i] = t(eij)%*%LIb%*%eij/strb
    }
    strg = sqrt(sum(elig$values^2))
    Di2 = 0
    for(i in 1:N){
      eij = rep(0,2*N)
      eij[i+N] = 1
      Di2[i] = t(eij)%*%LIg%*%eij/strg
    }
  }

  if(pert == "response"){
    am = model$mu.fitted.values^(1/model$phi.fitted.values)
    dm = am/(1-am)
    mi = model$mu.fitted.values
    phi = model$phi.fitted.values
    sy = sqrt(mi*((1/((2-mi^(1/phi))^phi))-mi))
    duw = dm*(1+dm)*sy/(model$mu.fitted.values*model$phi.fitted.values*(model$y))
    dvw = (sy/model$y)*(1/log(model$y)-dm*(1+dm)*log(model$y)/(model$phi.fitted.values^2))
    B = c(duw,dvw)
    LI = -diag(B)%*%solve(model$comp$Lambda)%*%model$comp$W%*%model$comp$Q%*%solve(S)%*%t(model$comp$Q)%*%model$comp$W%*%solve(model$comp$Lambda)%*%diag(B)
    LIb = -diag(B)%*%solve(model$comp$Lambda)%*%model$comp$W%*%model$comp$Q%*%(solve(S)-Sgg)%*%t(model$comp$Q)%*%model$comp$W%*%solve(model$comp$Lambda)%*%diag(B)
    LIg = -diag(B)%*%solve(model$comp$Lambda)%*%model$comp$W%*%model$comp$Q%*%(solve(S)-Sbb)%*%t(model$comp$Q)%*%model$comp$W%*%solve(model$comp$Lambda)%*%diag(B)
    elib = eigen(LIb, symmetric=T)
    elig = eigen(LIg, symmetric=T)
    strb = sqrt(sum(elib$values^2))
    Di = 0
    for(i in 1:N){
      eij = rep(0,2*N)
      eij[i] = 1
      Di[i] = t(eij)%*%LIb%*%eij/strb
    }
    strg = sqrt(sum(elig$values^2))
    Di2 = 0
    for(i in 1:N){
      eij = rep(0,2*N)
      eij[i+N] = 1
      Di2[i] = t(eij)%*%LIg%*%eij/strg
    }
  }
  dret = cbind(Di,Di2)
  colnames(dret) = c("mean", "precision")
  return(dret)
}

