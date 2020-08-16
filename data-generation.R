### Code for generating "Example.RData"  ###

set.seed(123)

# setting
sc <- 1   # scenario number (1-4)
m <- 200  # number of observations
om <- 0.1  # outlier ratio 
out <- round(m*om)


## generation of true lambda 
if( sc==1 ){ 
  la <- rgamma(m, 2, 2)
  la[1:out] <- rgamma(out, 14, 2)
}
if( sc==2 ){
  ch <- rbinom(m, 1, 0.25)
  la <- ch*1 + (1-ch)*rgamma(m, 2, 2)
  la[1:out] <- rgamma(out, 14, 2)
}
if( sc==3 ){
  ch <- rbinom(m, 1, 0.5)
  la <- ch*1 + (1-ch)*rgamma(m, 2, 2)
  la[1:out] <- rgamma(out, 14, 2)
}
if( sc==4 ){
  la <- runif(m, 0, 2)
  la[1:out] <- 4 + abs( rt(out, 3) )
}

## generation of observations
Eta <- runif(m, 1, 5)      # adjustment factor (known) 
Y <- rpois(m, Eta*la)    # observations 


## save data
save(Y, Eta, la, out, file="Example.RData")