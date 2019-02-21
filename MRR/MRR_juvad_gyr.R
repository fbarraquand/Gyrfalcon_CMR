# -------------------- Mark-recovery-recapture models for gyrs -----------------------------------#
# FB 15/02/2019 -- my version with only one live stage but a possibility for carcasses to remain. 
# FB 21/02/2019 -- adding a juvenile and adult stage, or rather adult juvenile differences in the survival probs
# An individual only stays juvenile for two years and then moves probabilistically up to the adult class
# Applied to the gyr data rather than simulated data

rm(list=ls())
library("R2jags")

## function
source('simul.R')

# -------------------------------------------------
# Loading data

CH=read.csv(file = "../data/CaptureHistories.csv")
CH[,1]<-NULL
CH=as.matrix(CH)
image(CH)

get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)
#read.csv(file = "data/fmarking.csv") would give only the dates

# Parameters used before for the simulation, for memory
# surv <- c(0.25,0.8) #survival probability for juveniles and adults
# eta <-0.1 # probability of carcasse staying uneaten/undecomposed 
# r <- 0.2 # dead recovery probability
# p <- 0.1 # live resighting/recapture probability

n.occasions <- ncol(CH)
n.states <- 3
n.obs <- 3

stagemarking = read.csv(file ="../data/stagemarking.csv")
stagemarking = stagemarking[,2]
table(stagemarking) #OK

# Matrix giving the stage of all individuals that are marked
xj <- matrix(NA, ncol = n.occasions-1, nrow = length(f))
for (i in 1:length(f)){
  for (t in f[i]:(n.occasions-1)){
    #if (fj[i]>(n.occasions-1)) next
    xj[i,t] <- 2
    if (stagemarking[i]==1){
    maxyearjuv = min(n.occasions-1,f[i]+1)
    xj[i,f[i]:maxyearjuv] <- 1   
    }
  } #t
} #i
image(xj)
nrow(xj[stagemarking==2,]) #to check the individuals tagged as adults are still there

xj[is.na(xj)]=2

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
rCH=CH
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
# Inference for Bugs model at "MRR.jags", fit using jags,
# 3 chains, each with 4000 iterations (first 1000 discarded), n.thin = 3
# n.sims = 3000 iterations saved
#             mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# mean.eta     0.011   0.011    0.000    0.003    0.007    0.015    0.040 1.007   310
# mean.p       0.039   0.006    0.028    0.034    0.038    0.042    0.052 1.013   180
# mean.r       0.138   0.009    0.122    0.132    0.138    0.144    0.156 1.001  3000
# mean.s[1]    0.291   0.017    0.257    0.280    0.292    0.303    0.326 1.138    20
# mean.s[2]    0.803   0.021    0.761    0.790    0.803    0.818    0.843 1.039    57
# deviance  1728.141   6.982 1716.264 1722.931 1727.822 1732.926 1742.715 1.171    17
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 21.4 and DIC = 1749.5
# DIC is an estimate of expected predictive error (lower deviance is better).

pdf(file="TraceDens_ni4000_juvad_gyr.pdf", width = 6, height = 8)
plot(as.mcmc(lifedead))
dev.off()

postsamples=cbind(lifedead$BUGSoutput$sims.list$mean.eta,
                  lifedead$BUGSoutput$sims.list$mean.p,
                  lifedead$BUGSoutput$sims.list$mean.r,
                  lifedead$BUGSoutput$sims.list$mean.s)
png(file="PairPosteriorPlot_juvad_gyr.png", width = 1200, height = 1200,res=300)
pairs(postsamples,c("eta","p","r","s_j","s_a"))
dev.off()

### Convergence does not seem exceptional with 4000, let's try more.. 


# MCMC settings
nc <- 3 #number of chains
nb <- 4000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-10000
nt <- 10 # “thinning”

# Call JAGS from R (BRT 80 min)
lifedead <- jags(jags.data, inits, parameters, "MRR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(lifedead, digit = 3)


