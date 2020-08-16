###-------------------------------------------------------------###
###    Global-local shrinkage prior (GLSP) for count data       ###
###-------------------------------------------------------------###
library(GIGrvg)
library(MCMCpack)
library(nleqslv)



###  GLSP given adjustment term   ###
# Y: vector of response counts 
# Eta: vector of adjustment factors 
# prior: EH or IG or PG
# HP: hyperparameter of gamma prior

GLSP.count <- function(Y, Eta=NULL, prior="EH", mc=3000, burn=500, HP=c(1,1)){
  m <- length(Y)    # number of observations 
  MC <- mc + burn    # length of MCMC
  if( is.null(Eta) ){ Eta <- rep(1 ,m) }
  
  ## MCMC sample box
  Lam.pos <- matrix(NA, MC, m)
  U.pos <- matrix(NA, MC, m)
  Beta.pos <- c()
  Alpha.pos <- c()
  Gam.pos <- c()
  
  ## initial values
  Lam <- Nu <- Y + 0.1
  U <- rep(1, m)
  if( prior=="EH" ){  V <- W <- rep(1, m)  }
  Beta <- Alpha <- Gam <- 1
  
  ## MCMC iterations
  for(r in 1:MC){
    # Lambda
    Lam <- rgamma(m, Y+Alpha, Eta+Beta/U)
    Lam.pos[r,] <- Lam
    
    ## U and Gam (EH prior)
    if( prior=="EH" ){
      # U
      V <- rgamma(m, 1+Gam, 1+log(1+U))
      W <- rgamma(m, 1+V, 1+U)
      for(i in 1:m){
        U[i] <- rgig(1, lambda=1-Alpha, chi=2*Beta*Lam[i], psi=2*W[i])
      }
      U.pos[r,] <- U
      # Gam
      ss <- sum( log(1+log(1+U)) )
      Gam <- rgamma(1, HP[1]+m, HP[2]+ss)
      Gam.pos[r] <- Gam
    }
    
    ## U and Gam (IG prior)
    if( prior=="IG" ){
      # U
      U <- rinvgamma(m, Alpha+Gam, Beta*Lam+Gam)
      # Gam
      bb <- 0.1
      Gam.new <- min(max(Gam+bb*rnorm(1), 0.001), 150)
      L1 <- m*Gam*log(Gam) - m*log(gamma(Gam)) - Gam*sum(log(U)) - Gam*sum(1/U)
      L2 <- m*Gam.new*log(Gam.new) - m*log(gamma(Gam.new)) - Gam.new*sum(log(U)) - Gam.new*sum(1/U)
      pp <- min(exp(L2-L1), 1)
      Gam <- Gam + rbinom(1, 1, pp)*(Gam.new-Gam)
      Gam.pos[r] <- Gam
    }
    
    # Beta
    Beta <- rgamma(1, HP[1]+m*Alpha, HP[2]+sum(Lam/U))
    Beta.pos[r] <- Beta
    
    # alpha
    Alpha <- rgamma(1, HP[1]+sum(Nu), HP[2]+sum(log(1+Eta*U/Beta)))
    Alpha.pos[r] <- Alpha
    
    # Nu
    for(i in 1:m){
      if(Y[i]==0){ Nu[i]=0 }
      if(Y[i]>0){
        pp <- Alpha/((1:Y[i])-1+Alpha)
        Nu[i] <- sum( rbinom(length(pp), 1, pp) )
      }
    }
  }
  
  # Summary
  om <- 1:burn
  Lam.pos <- Lam.pos[-om,]
  U.pos <- U.pos[-om,]
  Beta.pos <- Beta.pos[-om]
  Alpha.pos <- Alpha.pos[-om]
  Gam.pos <- Gam.pos[-om]
  Res <- list(lam=Lam.pos, u=U.pos, beta=Beta.pos, alpha=Alpha.pos, gam=Gam.pos)
  if( prior=="PG" ){ Res <- list(lam=Lam.pos, beta=Beta.pos, alpha=Alpha.pos) }
  return(Res)
}






###  GLSP with regression   ###
# Y: vector of response counts 
# X: matrix of covaraites 
# offset: 
# prior: EH or IG or PG
# HP: hyperparameter of gamma prior

GLSP.count.reg <- function(Y, X, offset=NULL, prior="EH", mc=3000, burn=500, HP=c(1,1)){
  m <- length(Y)    # number of observations 
  p <- dim(X)[2]
  MC <- mc + burn    # length of MCMC
  if( is.null(offset) ){ offset <- rep(1, m) }
  Om <- (1/100)*diag(p)    # precision for priors in regression coeffieicnts
  
  ## MCMC sample box
  Lam.pos <- matrix(NA, MC, m)
  U.pos <- matrix(NA, MC, m)
  Beta.pos <- c()
  Alpha.pos <- c()
  Gam.pos <- c()
  Reg.pos <- matrix(NA, MC, p)
  
  
  ## initial values
  Lam <- Nu <- Y + 0.1
  U <- rep(1, m)
  if( prior=="EH" ){  V <- W <- rep(1, m)  }
  Beta <- Alpha <- Gam <- 1
  Reg <- coef( glm(Y~X, offset=offset, family="poisson") )[-1]
  
  ## MCMC iterations
  for(r in 1:MC){
    # Regression part
    Q <- function(b){  t(X)%*%(Y - Lam*exp(offset + X%*%b)) }
    hReg <- nleqslv(Reg, Q)$x   # mode
    hmu <- as.vector( Lam*exp(offset + X%*%hReg) )
    mS <- t(X*hmu)%*%X
    A1 <- solve( mS + Om )
    A2 <- as.vector( A1%*%( mS%*%hReg ) )
    Reg.prop <- mvrnorm(1, A2, A1)  # proporsal 
    
    T1 <- t(Y)%*%X%*%(Reg.prop-Reg) - sum( Lam*as.vector(exp(offset + X%*%Reg.prop) - exp(offset + X%*%Reg)) )
    T2 <- 0.5*( t(Reg.prop-hReg)%*%mS%*%(Reg.prop-hReg) - t(Reg-hReg)%*%mS%*%(Reg-hReg) )  
    log.ratio <- T1 + T2
    pp <- min(exp(log.ratio), 1)
    ch <- rbinom(1, 1, pp)
    Reg <- as.vector( Reg+ch*(Reg.prop-Reg) )
    Eta <- exp(offset + as.vector(X%*%Reg))    # adjustment term
    Reg.pos[r,] <- Reg
    
    # Lambda
    Lam <- rgamma(m, Y+Alpha, Eta+Beta/U)
    Lam.pos[r,] <- Lam
    
    ## U and Gam (EH prior)
    if( prior=="EH" ){
      # U
      V <- rgamma(m, 1+Gam, 1+log(1+U))
      W <- rgamma(m, 1+V, 1+U)
      for(i in 1:m){
        U[i] <- rgig(1, lambda=1-Alpha, chi=2*Beta*Lam[i], psi=2*W[i])
      }
      U.pos[r,] <- U
      # Gam
      ss <- sum( log(1+log(1+U)) )
      Gam <- rgamma(1, HP[1]+m, HP[2]+ss)
      Gam.pos[r] <- Gam
    }
    
    ## U and Gam (IG prior)
    if( prior=="IG" ){
      # U
      U <- rinvgamma(m, Alpha+Gam, Beta*Lam+Gam)
      # Gam
      bb <- 0.1
      Gam.new <- min(max(Gam+bb*rnorm(1), 0.001), 150)
      L1 <- m*Gam*log(Gam) - m*log(gamma(Gam)) - Gam*sum(log(U)) - Gam*sum(1/U)
      L2 <- m*Gam.new*log(Gam.new) - m*log(gamma(Gam.new)) - Gam.new*sum(log(U)) - Gam.new*sum(1/U)
      pp <- min(exp(L2-L1), 1)
      Gam <- Gam + rbinom(1, 1, pp)*(Gam.new-Gam)
      Gam.pos[r] <- Gam
    }
    
    # Beta
    Beta <- rgamma(1, HP[1]+m*Alpha, HP[2]+sum(Lam/U))
    Beta.pos[r] <- Beta
    
    # alpha
    Alpha <- rgamma(1, HP[1]+sum(Nu), HP[2]+sum(log(1+Eta*U/Beta)))
    Alpha.pos[r] <- Alpha
    
    # Nu
    for(i in 1:m){
      if(Y[i]==0){ Nu[i]=0 }
      if(Y[i]>0){
        pp <- Alpha/((1:Y[i])-1+Alpha)
        Nu[i] <- sum( rbinom(length(pp), 1, pp) )
      }
    }
  }
  
  # Summary
  om <- 1:burn
  Lam.pos <- Lam.pos[-om,]
  U.pos <- U.pos[-om,]
  Beta.pos <- Beta.pos[-om]
  Alpha.pos <- Alpha.pos[-om]
  Gam.pos <- Gam.pos[-om]
  Reg.pos <- Reg.pos[-om,]
  Res <- list(lam=Lam.pos, u=U.pos, beta=Beta.pos, alpha=Alpha.pos, gam=Gam.pos, reg=Reg.pos)
  if( prior=="PG" ){ Res <- list(lam=Lam.pos, beta=Beta.pos, alpha=Alpha.pos, reg=Reg.pos) }
  return(Res)
}


