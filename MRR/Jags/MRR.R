
# FB 15/02/2019 -- my version with only one live stage but a possibility for carcasses to remain. 

rm(list=ls())
library("R2jags")

## function
source('simul.R')

# -------------------------------------------------
s <- 0.8 #survival probability
eta <-0.1 # probability of carcasse staying uneaten/undecomposed 
r <- 0.2 # dead recovery probability
p <- 0.1 # live resighting/recapture probability

n.occasions <- 50
n.states <- 3
n.obs <- 3
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(10, n.occasions)	# Releases in study area

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      s, 1-s, 0,
      0,   eta, 1-eta,
      0,    0,   1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      p, 0, 1-p,
      0, r, 1-r,
      0, 0, 1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

image(CH)

# Compute date of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed!
# 1 = alive and in study are, 2 = recovered dead, 3 = not seen or recovered
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 3

# Parameters:
# s: survival probability (no distinction between apparent and true to simplify)
# r: recovery probability
# eta <-0.1 # probability of carcasse staying uneaten/undecomposed 
# p: recapture/resighting probability
# -------------------------------------------------
# States (S):
# 1 alive i
# 2 recently dead (and possibly recovered)
# 3 dead (absorbing)
# Observations (O):
# 1 seen alive
# 2 recovered dead
# 3 neither seen nor recovered
# -------------------------------------------------

sink("MRR.jags")
cat("
model {

    # Priors and constraints
    for (t in 1:(n.occasions-1)){
    s[t] <- mean.s
    eta[t] <- mean.eta
    r[t] <- mean.r
    p[t] <- mean.p
    }
    mean.s ~ dunif(0, 1)     # Prior for mean survival
    mean.eta ~ dunif(0,0.5)     # Prior for mean carcass survival prob
    mean.r ~ dunif(0, 0.5)     # Prior for mean recovery
    mean.p ~ dunif(0, 0.5)     # Prior for mean recapture
    
    # Define state-transition and observation matrices 	
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- s[t]
    ps[1,i,t,2] <- (1-s[t])
    ps[1,i,t,3] <- 0
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- eta[t]
    ps[2,i,t,3] <- (1-eta[t])
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0
    ps[3,i,t,3] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- p[t]
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 1-p[t]
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- r[t]
    po[2,i,t,3] <- 1-r[t]
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- 1
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } #t
    } #i
    }
    ",fill = TRUE)
sink()


# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])

# Initial values

ld.init <- function(ch, f){
  
  for (i in 1:nrow(ch)){ # for all individuals
   x = ch[i,]
   F1 = f[i] #first one
   L1 = max(which(x==1))#last one
   F2 = min(which(x==2))#first two
   L2 = max(which(x==2))#last two
   
   if(!is.finite(F2)) #if there are only ones, not twos
   {
     ch[i,L1:F1]<-1
   }
   else #if there are both ones and twos
   {
     ch[i,F1:L1]<-1
     ch[i,L1:F2]<-1 #if there are zeroes between 1s and 2s
     ch[i,F2:L2]<-2
   }
   
   ch[i,1:F1]<-NA # the states up to the first capture are NAs
   
   LM = max(L1,L2)
   if(LM<ncol(ch)){
     ch[i,LM+1]<-2
     if (LM<(ncol(ch)-1)){ch[i,(LM+2):ncol(ch)]<-3 }
     #after last true state 2 (recently dead) the next real state is 3 (dead)
     #but after last true state 1, next true state is necessarily 2! Hence the 2 value
   } 
   
  } #end of loop on individuals
   return(ch)
}
  

## check initial states
z = ld.init(CH, f)
for (i in 1:10)
{
  print(z[i,])
  print(CH[i,])
}

inits <- function(){list(mean.s = runif(1, 0, 1), mean.eta = runif(1, 0, 0.5), 
                         mean.p = runif(1, 0, 0.5), mean.r = runif(1, 0, 0.5), z = ld.init(CH, f))}  
## essayons sans spécifier z -> node inconsistent with parents (as could be expected...)
#inits <- function(){list(mean.s = runif(1, 0, 1), mean.eta = runif(1, 0, 0.5), mean.p = runif(1, 0, 0.5), mean.r = runif(1, 0, 0.5))}  

# Parameters monitored
parameters <- c("mean.s", "mean.eta", "mean.r", "mean.p")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 80 min)
lifedead <- jags(jags.data, inits, parameters, "MRR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(lifedead, digit = 3)
