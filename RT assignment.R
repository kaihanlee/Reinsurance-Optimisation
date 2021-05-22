# Part 1
# declare variables
mu <- 1
lam <- 1/mu
theta <- 0.2
thetaR <- 0.3

alpha <- c(0.5, 0.8, 0.9)
U <- c(15, 20, 25)

# simulate Poisson claim rate
simulate.claim.until <- function(lambda, final.time){
  i<-1 # counter for number of observed points
  #Note we generate the time of the first point beyond the interval
  claim.times <- 0
  while(claim.times[i]<= final.time){
    inter<-rexp(1,lambda) #inter-claim time
    claim.times[i+1]<- claim.times[i]+inter
    i<-i+1
  }
  # Discard last claim since it will be beyond final.time
  claim.times<-claim.times[-length(claim.times)]
  claim.times
}
final.time <- 365   # duration

# simulate ruin processes
ruinsim <- function(alpha, U){
  
  N <- 10000  # number of simulations
  ruinsim <- c()  # declare empty list for number of ruins
  
  for(n in 1:N){
    Ut = c()
    claim.times <- simulate.claim.until(lam, final.time)  # simulate claim times
    claim_amount <- rexp(length(claim.times),1)   # simulate claim amounts
    sum_claims <- c(0, cumsum(claim_amount))   # sum up claims
    
    # simulate individual ruin processes
    for (i in 2:length(claim.times)){
      c <- (1+theta)*lam - (1+thetaR)*lam*(1-alpha)   # insurance premium
      Ut[1] = U + c*(claim.times[1])     # surplus at time 1
      Ut[i] = U + c*(claim.times[i]) - sum_claims[i]*alpha   # surplus at time i
    }
    
    ruinsim[n] <- sum(Ut < 0)   # check occurrence of U(t)<0
  }
  
  ruin <- sum(ruinsim>0)/N; ruin    # P(Ruin) = # of ruins / Total simulations
  return(ruin)
}

# create empty matrix
ruinmat <- matrix(rep(0,9),nrow=3,ncol=3)
rownames(ruinmat) <- c("a=0.5","a=0.8","a=0.9")
colnames(ruinmat) <- c("u=15","u=20","u=25")

# fill in matrix
for(i in 1:3){
  for(j in 1:3){
    ruinmat[i,j] <- ruinsim(alpha[i],U[j])
  }
}

ruinmat   # display probabilities of ruin (Different alpha and U)

alpha=1
# create empty matrix
ruinmat1 <- matrix(rep(0,6),nrow=2,ncol=3)
rownames(ruinmat1) <- c("approx","exact")
colnames(ruinmat1) <- c("u=15","u=20","u=25")

# fill in matrix
for(g in 1:3){
  ruinmat1[1,g] = ruinsim(1,U[g])
  ruinmat1[2,g] = (1 / (1 + theta)) * exp(-alpha * theta * U[g] / (1 + theta))
}

ruinmat1   # display probabilities of ruin (Approx. vs Exact)

# plot graph
plot(U, ruinmat[1,], main="Probability of Ruin, Psi(U) against Initial Surplus U", 
     xlab="Initial Surplus U", ylab="Prob. of Ruin, Psi(U)", type="l", lty=1, 
     col=1, ylim=c(0,0.07), xlim=c(14,26))
lines(U,ruinmat[2,], lty=2, col=2)
lines(U,ruinmat[3,], lty=3, col=3)
lines(U,ruinmat1[1,], lty=4, col=4)
legend(x=23, y=0.07, cex=0.7, c("a=0.5","a=0.8","a=0.9","a=1.0"), lty=c(1,2,3,4), col=c(1,2,3,4))

# Part 2

# Use Inverse CDF method to simulate a Pareto distribution
rpareto <- function(n, palpha, plambda){
  uu <- runif(n, min=0, max=1)
  palpha <- palpha
  plambda <- plambda
  
  rpa <- plambda*uu^(-1/palpha)-plambda
  return (rpa)
}

# simulate ruin processes
ruinsim2 <- function(a, b, U){
  
  N <- 10000     # number of simulations
  ruinsim <- c()
  m1 <- b/(a-1)     # 1st moment of Pareto distribution
  
  for(n in 1:N){
    Ut = c()  # declare empty list for number of ruins
    claim.times <- simulate.claim.until(lam, final.time)  # simulate claim times
    claim_amount <- rpareto(length(claim.times), a, b)   # simulate claim amounts
    sum_claims <- c(0, cumsum(claim_amount))   # sum up claims
    
    # simulate individual ruin processes
    for (i in 2:length(claim.times)){
      c <- (1+theta)*lam - (1+thetaR)*lam*(1-alpha)   # insurance premium
      Ut[1] = U + c*(claim.times[1])     # surplus at time 1
      Ut[i] = U + c*(claim.times[i]) - sum_claims[i]*alpha   # surplus at time i
    }
    
    ruinsim[n] <- sum(Ut < 0)   # check occurrence of U(t)<0
  }
  
  ruin <- sum(ruinsim>0)/N; ruin    # P(Ruin) = # of ruins / Total simulations
  return(ruin)
}

# create empty matrix
ruinmat2 <- matrix(rep(0,3),nrow=1,ncol=3)
colnames(ruinmat2) <- c("u=15","u=20","u=25")

# fill in matrix
for(i in 1:3){
  ruinmat2[1,i] <- ruinsim2(5,4,U[i])
}

ruinmat2   # display probabilities of ruin (Pareto iid claims)
