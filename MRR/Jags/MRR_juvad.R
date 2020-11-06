# -------------------- Mark-recovery-recapture models for gyrs -----------------------------------#
# FB 15/02/2019 -- my version with only one live stage but a possibility for carcasses to remain. 
# FB 21/02/2019 -- adding a juvenile and adult stage, or rather adult juvenile differences in the survival probs
# An individual only stays juvenile for two years and then moves probabilistically up to the adult class

rm(list=ls())
library("R2jags")

## function
source('simul.R')

# -------------------------------------------------
surv <- c(0.25,0.8) #survival probability for juveniles and adults
eta <-0.1 # probability of carcasse staying uneaten/undecomposed 
r <- 0.2 # dead recovery probability
p <- 0.1 # live resighting/recapture probability

n.occasions <- 50 #10 for speed
n.states <- 3
n.obs <- 3
### Matrix of marked individuals 
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(10, n.occasions)	# Releases in study area - 10 juveniles a year

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time

# 1. State process matrix -- modified
totrel <- sum(marked)*(n.occasions-1)
# that doesn't sound right
totind<-sum(marked)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totind, n.occasions-1))

# First year of capture / year of marking
# nb of individuals marked as young per year
fj <- as.numeric(gl(n.occasions-1,10))
fmarking<-fj

# Matrix giving the stage of all individuals that are marked
xj <- matrix(NA, ncol = n.occasions-1, nrow = length(fj))
for (i in 1:length(fj)){
  for (t in fj[i]:(n.occasions-1)){
    #if (fj[i]>(n.occasions-1)) next
    xj[i,t] <- 2
    maxyearjuv = min(n.occasions-1,fj[i]+1)
    xj[i,fj[i]:maxyearjuv] <- 1   
  } #t
} #i
xj[is.na(xj)]=2

for (i in 1:length(fj)){ 
  #before it was up to totrel, but totrel is sum(marked)*(n.occasions-1) which is more than the 
  #number of individuals
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      surv[xj[i,t]], 1-surv[xj[i,t]], 0,
      0,   eta, 1-eta,
      0,    0,   1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totind, n.occasions-1))
for (i in 1:totind){
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
    s[t,1] <- mean.s[1]
    s[t,2] <- mean.s[2]
    eta[t] <- mean.eta
    r[t] <- mean.r
    p[t] <- mean.p
    }
    mean.s[1] ~ dunif(0, 1)     # Prior for mean juvenile survival
    mean.s[2] ~ dunif(0, 1)     # Prior for mean adult survival

    mean.eta ~ dunif(0,0.5)     # Prior for mean carcass survival prob
    mean.r ~ dunif(0, 0.5)     # Prior for mean recovery
    mean.p ~ dunif(0, 0.5)     # Prior for mean recapture
    
    # Define state-transition and observation matrices 	
    for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
    ps[1,i,t,1] <- s[t,xj[i,t]]
    ps[1,i,t,2] <- (1-s[t,xj[i,t]])
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
jags.data <- list(y = rCH, f = f,  xj=xj, n.occasions = dim(rCH)[2],nind = dim(rCH)[1])

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

inits <- function(){list(mean.s = runif(2, 0, 1), mean.eta = runif(1, 0, 0.5), 
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
# 
# Inference for Bugs model at "MRR.jags", fit using jags,
# 3 chains, each with 4000 iterations (first 1000 discarded), n.thin = 3
# n.sims = 3000 iterations saved
# mu.vect sd.vect    2.5%     25%     50%     75%   97.5%  Rhat n.eff
# mean.eta    0.072   0.040   0.010   0.041   0.067   0.096   0.160 1.004   590
# mean.p      0.049   0.013   0.029   0.040   0.048   0.057   0.079 1.019   110
# mean.r      0.217   0.020   0.178   0.203   0.217   0.230   0.256 1.001  3000
# mean.s[1]   0.337   0.035   0.267   0.313   0.339   0.362   0.403 1.013   210
# mean.s[2]   0.817   0.039   0.741   0.790   0.817   0.845   0.891 1.286    11
# deviance  700.223  10.439 682.053 693.029 699.493 706.733 723.113 1.028    76
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 53.1 and DIC = 753.3
# DIC is an estimate of expected predictive error (lower deviance is better).
