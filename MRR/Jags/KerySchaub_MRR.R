### Kéry & Schaub BPA code for the MRR model (special version)

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
n.occasions <- 10  
n.states <- 4
n.obs <- 3
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)	# Releases in study area

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
      ps[1,i,t,3] <- (1-s[t])*r[t]
      ps[1,i,t,4] <- (1-s[t])*(1-r[t])
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- s[t]
      ps[2,i,t,3] <- (1-s[t])*r[t]
      ps[2,i,t,4] <- (1-s[t])*(1-r[t])
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
      po[3,i,t,2] <- 1
      po[3,i,t,3] <- 0
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
# In order to run the model, we must provide sound initial values for the true latent state z. The difficulty is that observed states do not alsways correspond to the true latent state. For example, in our observation the state 2 refers to an individuals whose ring has been reported, while the true state 2 refers to an individuals that is avlive, but outside the study area. Consequently, we cannot use the observed states as the initial values for the true state in JAGS (in BUGS this works well). The function known.ld provides initial values for the state z for our model. The important things are i) that the observations correspond to the true state (i.e. all "2" are converted into "3"), ii) that states after a "3" are all "4", iii) that all non-observations between "1" become "1", and iv) that all remaining originally "3" after the first observation become "1".

ld.init <- function(ch, f){
   ch[ch==3] <- NA
   v2 <- which(ch==2, arr.ind = T)
   ch[v2] <- 3
   for (i in 1:nrow(v2)){
      ifelse(v2[i,2]!=ncol(ch), ch[v2[i,1], (v2[i,2]+1):ncol(ch)] <- 4, next)}
   for (i in 1:nrow(ch)){
      m <- max(which(ch[i,]==1))
      ch[i,f[i]:m] <- 1
      }
   for (i in 1:nrow(v2)){
      u1 <- min(which(ch[v2[i,1],]==1))
      ch[v2[i],u1:(v2[i,2]-1)] <- 1
      }
   for (i in 1:nrow(ch)){
      for (j in f[i]:ncol(ch)){
         if(is.na(ch[i,j])==1) ch[i,j] <- 1
         }
      }
   return(ch)
   }

inits <- function(){list(mean.s = runif(1, 0, 1), mean.f = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.r = runif(1, 0, 1), z = ld.init(rCH, f))}  

# Parameters monitored
parameters <- c("mean.s", "mean.f", "mean.r", "mean.p")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 80 min)
lifedead <- jags(jags.data, inits, parameters, "lifedead.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(lifedead, digit = 3)

# Note that convergence is hard to get, much longer chains or more informative priors would be necessary to get convergence quicker.


# Add-in for JAGS
# Since we have created a matrix with initial values for the true state z, we can use part of this information as data (see chapter 7.3.1) which can help with convergence and computing time). 
# Here we give those initial values that are based on an actual observation. 
# Since the first observation is deterministic, it must be excluded. 
# The following code constructs the data matrix:

ch <- rCH
ch[ch==3] <- NA
z.known <- ld.init(rCH, f)
z.known[is.na(ch)] <- NA
for (i in 1:nrow(ch)){
   z.known[i,f[i]] <- NA
   }

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = z.known)

inits <- function(){list(mean.s = runif(1, 0, 1), mean.f = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.r = runif(1, 0, 1), z = ld.init(rCH, f))}  

# Parameters monitored
parameters <- c("mean.s", "mean.f", "mean.r", "mean.p")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 80 min)
lifedead <- jags(jags.data, inits, parameters, "lifedead.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(lifedead, digit = 3)
