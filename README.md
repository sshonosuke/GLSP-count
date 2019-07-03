# GLSP-count: Global-local shrinkage priors for analyzing sequence of counts
This package implements a new shrinkage method based on global-local shrinkage priors for Poisson likelihood via Markov Chain Monte Carlo algorithm, as proposed by the following paper.

Hamura, Y., Irie, K. and Sugasawa, S. (2019). On Global-local Shrinkage Priors for Count Data. 
https://arxiv.org/abs/1907.01333

Functions are implemented in *GLSP-function.R*.

Install R functions and an example dataset.
```{r}
load("Example-Data.RData") 
source("GLSP-function.R")
```

Apply the proposed method.
- `Y`: vector of observed counts
- `Eta`: vector of structural part (adjustment factor)
- `cc`: common hyperparameter in inverse gamma priors for scale parameters (default: 1)
- `mc`: length of MCMC run (default: 3000)
- `burn`: number of burn-in period (default: 500)
- `prior`: prior for local parameters: EH (extremely heavy tailed) or IG (inverse gamma)
```{r}
set.seed(1)
fit1=GLSP.PG(Y,Eta,prior="EH")
fit2=GLSP.PG(Y,Eta,prior="IG")
```


Apply the conventional Poisson-gamma model.
```{r}
fit3=PG(Y,Eta)
```

Compute and compare posterior means.
```{r}
Est1=apply(fit1$lam,2,mean)
Est2=apply(fit2$lam,2,mean)
Est3=apply(fit3$lam,2,mean)
plot(Est1,Est3); abline(0,1)
```

Compute squared errors.
```{r}
Est=cbind(Est1,Est2,Est3,Y)
apply(((Est-la)^2)[out,],2,mean)
apply(((Est-la)^2)[-out,],2,mean)
```

Compute and show credible intervals.
```{r}
quant=function(x){ quantile(x, prob=c(0.025,0.975)) }
ID=c(1,2,3,11,12,13)  # selected areas
L=length(ID)

CI=array(NA,c(3,2,L))
CI[1,,]=apply(fit1$lam,2,quant)[,ID]
CI[2,,]=apply(fit2$lam,2,quant)[,ID]
CI[3,,]=apply(fit3$lam,2,quant)[,ID]
dimnames(CI)[[1]]=c("EH","IG","PG")
dimnames(CI)[[3]]=ID

Est=cbind(Est1,Est2,Est3)[ID,]

ran=range(CI)
color=c(1,2,4)
plot(NA,ylim=ran,xlim=c(1-0.5,L+0.5),ylab="lambda",xlab="area",xaxt="n",main="Credible interval")
axis(1,at=1:L,label=ID)
for(s in 1:L){
  ss=seq(s-0.2,s+0.2,length=3)
  for(k in 1:3){
    lines(c(ss[k],ss[k]),CI[k,,s],lty=1,col=color[k],lwd=1.2)
    lines(c(ss[k]-0.05,ss[k]+0.05),CI[k,c(1,1),s],lty=1,col=color[k])
    points(ss[k],Est[s,k],pch=18)
    lines(c(ss[k]-0.05,ss[k]+0.05),CI[k,c(2,2),s],lty=1,col=color[k])
  }
}
legend("topright",c("EH","IG","PG","Posterior mean"),pch=c(NA,NA,NA,18),lty=c(1,1,1,NA),col=c(1,2,4,1),ncol=2,lwd=1.2)
```

