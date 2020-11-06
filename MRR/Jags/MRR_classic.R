### FB 21/02/2019 -- version written but not fitted by Kéry and Schaub
### see Kéry & Schaub BPA book 

rm(list=ls())
library("R2jags")

## function
source('simul.R')

# 9.5. Joint analysis of capture-recapture and mark-recovery data
# 9.5.1. Model description
# 9.5.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals 
s <- 0.8
F <- 0.6
r <- 0.1
p <- 0.5
n.occasions <- 50  
n.states <- 4
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
totrel <- sum(marked)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      s*F, s*(1-F), 1-s, 0,
      0,   s,       1-s, 0,
      0,   0,       0,   1,
      0,   0,       0,   1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      p, 0, 1-p,
      0, 0, 1,
      0, r, 1-r,
      0, 0, 1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute date of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed!
# 1 = alive and in study are, 2 = recovered dead, 3 = not seen or recovered
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 3


# 9.5.3. Analysis of the model
# Specify model in BUGS language
sink("lifedead.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# s: true survival probability
# F: fidelity probability
# r: recovery probability
# p: recapture/resighting probability
# -------------------------------------------------
# States (S):
# 1 alive in study area
# 2 alive outside study area
# 3 recently dead and recovered
# 4 recently dead, but not recovered, or dead (absorbing)
# Observations (O):
# 1 seen alive
# 2 recovered dead
# 3 neither seen nor recovered
# -------------------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   s[t] <- mean.s
   F[t] <- mean.f
   r[t] <- mean.r
   p[t] <- mean.p
   }
mean.s ~ dunif(0, 1)     # Prior for mean survival
mean.f ~ dunif(0, 1)     # Prior for mean fidelity
mean.r ~ dunif(0, 1)     # Prior for mean recovery
mean.p ~ dunif(0, 1)     # Prior for mean recapture

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- s[t]*F[t]
      ps[1,i,t,2] <- s[t]*(1-F[t])
      ps[1,i,t,3] <- (1-s[t])
      ps[1,i,t,4] <- 0 
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- s[t]
      ps[2,i,t,3] <- (1-s[t])
      ps[2,i,t,4] <- 0
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- 1
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- p[t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-p[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- 0
      po[2,i,t,3] <- 1
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- r[t]
      po[3,i,t,3] <- 1-r[t]
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 1
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
# FB: my version for initial values of states

ld.init <- function(ch, f){ 
  ## 21/02/2019: Modified because there are four true (latent)states now
  for (i in 1:nrow(ch)){ # for all individuals
   x = ch[i,]
   F1 = f[i] #first one
   L1 = max(which(x==1))#last one in the true capture history
   F2 = min(which(x==2))#first two in the true capture history
   L2 = max(which(x==2))#last two in the true capture history
   
   if(!is.finite(F2)) #if there are only ones, not twos
   {
     ch[i,L1:F1]<-1 #we assume the individual is in the population
   } else #if there are both ones and twos
   {
     ch[i,F1:L1]<-1
     ch[i,L1:F2]<-1 #if there are zeroes between 1s and 2s
     ch[i,F2:L2]<-3 #we set the true state to 3, recently dead
   }
   
   ch[i,1:F1]<-NA # the states up to the first capture are NAs
   
   LM = max(L1,L2)
   if(LM<ncol(ch)){
     #ch[i,LM+1]<-2 ## Problematic for this model with eta=0 [see other code]
     if (LM==L1){ch[i,LM+1]<-3}else{ch[i,LM+1]<-4} #only one true state 2 in the sequence here
     if (LM<(ncol(ch)-1)){ch[i,(LM+2):ncol(ch)]<-4 }
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


inits <- function(){list(mean.s = runif(1, 0, 1), mean.f = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.r = runif(1, 0, 1), z = ld.init(rCH, f))}  

# Parameters monitored
parameters <- c("mean.s", "mean.f", "mean.r", "mean.p")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R 
lifedead <- jags(jags.data, inits, parameters, "lifedead.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(lifedead, digit = 3)
# Note that convergence can sometimes be tricky
plot(as.mcmc(lifedead))

# Keep in mind
# s <- 0.8
# F <- 0.6
# r <- 0.1
# p <- 0.5

# With dunif(0,1) priors for all parameters
# Inference for Bugs model at "lifedead.jags", fit using jags,
# 3 chains, each with 4000 iterations (first 1000 discarded), n.thin = 3
# n.sims = 3000 iterations saved
#            mu.vect sd.vect    2.5%     25%     50%     75%   97.5%  Rhat n.eff
# mean.f     0.915   0.018   0.877   0.903   0.915   0.927   0.948 1.001  3000
# mean.p     0.703   0.031   0.639   0.683   0.704   0.725   0.760 1.002  1600
# mean.r     0.091   0.013   0.067   0.082   0.090   0.099   0.119 1.001  3000
# mean.s     0.475   0.016   0.444   0.464   0.475   0.486   0.508 1.001  2200
# deviance 631.615  14.079 607.015 621.593 630.870 640.479 661.706 1.005   690
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 98.9 and DIC = 730.5
# DIC is an estimate of expected predictive error (lower deviance is better).

#Plot correlations by pairs? Might be substantial correlations in posteriors there. 

### Plot pair posterior densities
postsamples=cbind(lifedead$BUGSoutput$sims.list$mean.f,
                  lifedead$BUGSoutput$sims.list$mean.p,
                  lifedead$BUGSoutput$sims.list$mean.r,
                  lifedead$BUGSoutput$sims.list$mean.s)
png(file="PairPosteriorPlot_MRRclassic.png", width = 1200, height = 1200,res=300)
pairs(postsamples,c("f","p","r","s"))
dev.off()
# May be that JAGS does not update correctly the latent states
# And that's why we have a very high site fidelity here
# There could identifiability problems, 
# partly solved by the reparameterization proposed by Kéry & Schaud

### Not clear I can do the stuff with known states suggested in the BPA book here

