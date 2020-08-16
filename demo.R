## load simulated data 
load("Example.RData")

## load functions 
source("GLSP-function.R")


## fitting models 
fit1 <- GLSP.count(Y, Eta, mc=3000, burn=500, prior="IG")
fit2 <- GLSP.count(Y, Eta, mc=3000, burn=500, prior="EH")
fit3 <- GLSP.count(Y, Eta, mc=3000, burn=500, prior="PG")


## posterior mean 
est1 <- apply(fit1$lam, 2, mean)
est2 <- apply(fit2$lam, 2, mean)
est3 <- apply(fit3$lam, 2, mean)

plot(est1, est3)
points(est2, est3, col=2)
abline(0,1)
legend("topleft", c("IG vs PG", "EH vs PG"), lty=1, col=1:2)


# squared error
Est <- cbind(est1, est2, est3, Y/Eta)
apply(((Est-la)^2)[1:out,], 2, mean)      # outlier
apply(((Est-la)^2)[-(1:out),], 2, mean)   # non-outlier


## credible interval
quant <- function(x){ quantile(x, prob=c(0.025,0.975)) }
ID <- 1:10  # selected areas
L <- length(ID)

CI <- array(NA,c(3,2,L))
CI[1,,] <- apply(fit1$lam, 2, quant)[,ID]
CI[2,,] <- apply(fit2$lam, 2, quant)[,ID]
CI[3,,] <- apply(fit3$lam, 2, quant)[,ID]
dimnames(CI)[[1]] <- c("EH","IG","PG")
dimnames(CI)[[3]] <- ID


# Plot
ran <- range(CI)
color <- c(1,2,4)
plot(NA, ylim=ran, xlim=c(1-0.5,L+0.5), ylab="lambda", xlab="area", xaxt="n", main="Credible interval")
axis(1, at=1:L, label=ID)
for(s in 1:L){
  ss <- seq(s-0.2, s+0.2, length=3)
  for(k in 1:3){
    lines(c(ss[k],ss[k]), CI[k,,s], lty=1, col=color[k], lwd=1.2)
    lines(c(ss[k]-0.05,ss[k]+0.05), CI[k,c(1,1),s], lty=1, col=color[k])
    points(ss[k], Est[ID[s], k], pch=18)
    lines(c(ss[k]-0.05,ss[k]+0.05), CI[k,c(2,2),s], lty=1, col=color[k])
  }
}
legend("topright",c("EH", "IG", "PG", "Posterior mean"), pch=c(NA,NA,NA,18), lty=c(1,1,1,NA), col=c(1,2,4,1), ncol=2, lwd=1.2)
