library(GIGrvg)
library(MCMCpack)


###  Function implementing the proposed methods
GLSP.PG=function(Y,Eta,cc=1,mc=3000,burn=500,prior="EH"){
  m=length(Y)
  # matrix for posterior samples
  Lam.pos=matrix(NA,mc,m)
  Beta.pos=rep(NA,mc)
  Alpha.pos=rep(NA,mc)
  Nu.pos=matrix(NA,mc,m)
  Gam.pos=rep(NA,mc)
  U.pos=matrix(NA,mc,m)
  if(prior=="EH"){
    V.pos=matrix(NA,mc,m)
    W.pos=matrix(NA,mc,m)
    TT.pos=matrix(NA,mc,m)
  }

  # initial values
  Lam=Y
  U=V=W=TT=rep(1,m)
  Beta=Alpha=Gam=1
  Nu=Y
  
  # MCMC
  for(r in 1:mc){
    ### sampling common parameters ###
    # Lambda
    Lam=rgamma(m,Y+Alpha,Eta+Beta/U)
    Lam.pos[r,]=Lam
    # Beta
    Beta=rgamma(1,cc+m*Alpha,cc+sum(Lam/U))
    Beta.pos[r]=Beta
    # alpha
    Alpha=rgamma(1,cc+sum(Nu),cc+sum(log(1+Eta*U/Beta)))
    Alpha.pos[r]=Alpha
    # Nu
    for(i in 1:m){
      if(Y[i]==0){ Nu[i]=0 }
      if(Y[i]>0){
        pp=Alpha/((1:Y[i])-1+Alpha)
        d=rbinom(length(pp),1,pp)
        Nu[i]=sum(d)
      }
    }
    Nu.pos[r,]=Nu
    
    ###  EH prior  ###
    if(prior=="EH"){
      # U
      a1=1-Alpha 
      a2=2*Beta*Lam
      a3=2*(V+TT)    
      for(i in 1:m){ U[i]=rgig(1,a1,a2[i],a3[i]) }
      U.pos[r,]=U
      # V
      V=rgamma(m,1,1+U)
      V.pos[r,]=V
      # W
      W=rgamma(m,1+Gam,1+log(1+U))
      W.pos[r,]=W
      # TT
      TT=rgamma(m,W,1+U)
      TT.pos[r,]=TT
      # Gamma
      ss=sum(log(1+log(1+U)))
      Gam=rgamma(1,cc+m,cc+ss)+0.001
      Gam.pos[r]=Gam
    }
    
    ###  IG prior  ###
    if(prior=="IG"){
      # U
      U=rinvgamma(m,Alpha+Gam,Beta*Lam+Gam)
      for(i in 1:m){ U[i]=max(U[i],0.001) }
      U.pos[r,]=U
      # Gamma
      bb=0.05
      prop=max(Gam+bb*rnorm(1),0.001)
      L1=m*Gam*log(Gam)-m*log(gamma(Gam))-Gam*sum(log(U.pos[r,]))-Gam*sum(1/U.pos[r,])
      L2=m*prop*log(prop)-m*log(gamma(prop))-prop*sum(log(U.pos[r,]))-prop*sum(1/U.pos[r,])
      pp=min(exp(L2-L1),1)
      Gam=Gam+rbinom(1,1,pp)*(prop-Gam)
      Gam.pos[r]=Gam
    }
  }
  om=1:burn
  Res=list(Lam.pos[-om,],Alpha.pos[-om],Beta.pos[-om])
  names(Res)=c("lam","alpha","beta")
  return(Res)
}








###   Function implementing the conventional Poisson-gamma model
PG=function(Y,Eta,cc=1,mc=3000,burn=500){
  m=length(Y)
  # matrix for posterior samples
  Lam.pos=matrix(NA,mc,m)
  Beta.pos=rep(NA,mc)
  Alpha.pos=rep(NA,mc)
  Nu.pos=matrix(NA,mc,m)
  
  # Initial values
  Lam=Nu=Y
  Beta=Alpha=1
  
  # MCMC
  for(r in 1:mc){
    # Lambda
    Lam=rgamma(m,Y+Alpha,Eta+Beta)
    Lam.pos[r,]=Lam
    # Beta
    Beta=rgamma(1,cc+m*Alpha,cc+sum(Lam))
    Beta.pos[r]=Beta
    # alpha
    Alpha=rgamma(1,cc+sum(Nu),cc+sum(log(1+Eta/Beta)))
    Alpha.pos[r]=Alpha
    # Nu
    for(i in 1:m){
      if(Y[i]==0){ Nu[i]=0 }
      if(Y[i]>0){ 
        pp=Alpha/((1:Y[i])-1+Alpha)
        Nu[i]=sum(rbinom(length(pp),1,pp))
      }
    }
    Nu.pos[r,]=Nu
  }
  om=1:burn
  Res=list(Lam.pos[-om,],Alpha.pos[-om],Beta.pos[-om])
  names(Res)=c("lam","alpha","beta")
  return(Res)
}



