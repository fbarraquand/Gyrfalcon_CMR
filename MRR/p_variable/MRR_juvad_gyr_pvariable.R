# -------------------- Mark-recovery-recapture models for gyrs -----------------------------------#
# FB 15/02/2019 -- my version with only one live stage but a possibility for carcasses to remain. 
# FB 21/02/2019 -- adding a juvenile and adult stage, or rather adult juvenile differences in the survival probs
# An individual only stays juvenile for two years and then moves probabilistically up to the adult class
# Applied to the gyr data rather than simulated data
# FB 04/05/2019 -- implements variable recapture probability, which is probably U-shaped given the number of live recaptures over time

rm(list=ls())
library("R2jags")

## function
source('../simul.R')

# -------------------------------------------------
# Loading data

CH=read.csv(file = "../../data/CaptureHistories.csv")
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

stagemarking = read.csv(file ="../../data/stagemarking.csv")
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
    reltime[t] <-t/(n.occasions-1)
    p[t] <- (1-gamma)*mean.p+gamma*4*(reltime[t]^2-reltime[t]+0.25)
    }
    mean.s[1] ~ dunif(0, 1)     # Prior for mean juvenile survival
    mean.s[2] ~ dunif(0, 1)     # Prior for mean adult survival

    mean.eta ~ dunif(0,0.5)    # Prior for mean carcass survival prob
    mean.r ~ dunif(0, 0.5)     # Prior for mean recovery
    mean.p ~ dunif(0, 0.5)     # Prior for mean recapture
    gamma ~ dunif(0,1)         # Prior for the weight given to the U-shaped probability
    
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

inits <- function(){list(mean.s = runif(2, 0, 1), mean.eta = runif(1, 0, 0.5), gamma = runif(1,0.2,0.8),
                         mean.p = runif(1, 0, 0.5), mean.r = runif(1, 0, 0.5), z = ld.init(CH, f))}  
## essayons sans spécifier z -> node inconsistent with parents (as could be expected...)
#inits <- function(){list(mean.s = runif(1, 0, 1), mean.eta = runif(1, 0, 0.5), mean.p = runif(1, 0, 0.5), mean.r = runif(1, 0, 0.5))}  

# Parameters monitored
parameters <- c("mean.s", "mean.eta", "mean.r", "mean.p","gamma")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 80 min)
lifedead <- jags(jags.data, inits, parameters, "MRR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(lifedead, digit = 3)

########## Old results without variable p ##########################################
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

#             mu.vect  2.5%     97.5%  Rhat n.eff
# mean.eta     0.011   0.000    0.040 1.007   310
# mean.p       0.039   0.028    0.052 1.013   180
# mean.r       0.138   0.122    0.156 1.001  3000
# mean.s[1]    0.291   0.257    0.326 1.138    20
# mean.s[2]    0.803   0.761    0.843 1.039    57

########## New results with variable p #################################################

# Inference for Bugs model at "MRR.jags", fit using jags,
# 3 chains, each with 4000 iterations (first 1000 discarded), n.thin = 3
# n.sims = 3000 iterations saved
# mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# gamma        0.186   0.033    0.129    0.163    0.183    0.206    0.263 1.022    99
# mean.eta     0.010   0.009    0.000    0.003    0.007    0.014    0.034 1.005   400
# mean.p       0.004   0.004    0.000    0.001    0.003    0.005    0.013 1.010   210
# mean.r       0.138   0.009    0.122    0.132    0.138    0.144    0.155 1.001  3000
# mean.s[1]    0.283   0.023    0.238    0.268    0.283    0.298    0.331 1.435     8
# mean.s[2]    0.791   0.018    0.754    0.779    0.791    0.803    0.825 1.150    18
# deviance  1670.165   8.782 1651.485 1664.639 1670.010 1676.287 1686.770 1.128    20
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 34.6 and DIC = 1704.8
# DIC is an estimate of expected predictive error (lower deviance is better).

########################################################################################

pdf(file="TraceDens_ni4000_juvad_gyr.pdf", width = 6, height = 8)
plot(as.mcmc(lifedead))
dev.off()

postsamples=cbind(lifedead$BUGSoutput$sims.list$mean.eta,
                  lifedead$BUGSoutput$sims.list$mean.p,
                  lifedead$BUGSoutput$sims.list$gamma,
                  lifedead$BUGSoutput$sims.list$mean.r,
                  lifedead$BUGSoutput$sims.list$mean.s)
png(file="PairPosteriorPlot_juvad_gyr.png", width = 1200, height = 1200,res=300)
pairs(postsamples,c("eta","p","gamma","r","s_j","s_a"))
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

############# Old results without p variable ##############################
# Inference for Bugs model at "MRR.jags", fit using jags,
# 3 chains, each with 10000 iterations (first 4000 discarded), n.thin = 10
# n.sims = 1800 iterations saved
# mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# mean.eta     0.010   0.010    0.000    0.003    0.007    0.014    0.036 1.032    70
# mean.p       0.033   0.005    0.024    0.029    0.033    0.037    0.045 1.149    18
# mean.r       0.139   0.009    0.123    0.133    0.138    0.145    0.156 1.005   460
# mean.s[1]    0.318   0.020    0.279    0.303    0.318    0.333    0.359 1.407     9
# mean.s[2]    0.813   0.018    0.776    0.801    0.814    0.825    0.847 1.014   150
# deviance  1740.633   7.481 1725.727 1735.359 1741.373 1745.606 1755.261 1.415     9
############################################################################

############ New results with p variable ###################################
# Inference for Bugs model at "MRR.jags", fit using jags,
# 3 chains, each with 10000 iterations (first 4000 discarded), n.thin = 10
# n.sims = 1800 iterations saved
# mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# gamma        0.146   0.028    0.095    0.127    0.145    0.164    0.208 1.029    96
# mean.eta     0.009   0.010    0.000    0.003    0.006    0.013    0.036 1.003  1100
# mean.p       0.003   0.003    0.000    0.001    0.002    0.004    0.011 1.004  1000
# mean.r       0.139   0.009    0.122    0.133    0.139    0.144    0.156 1.004   470
# mean.s[1]    0.314   0.020    0.276    0.300    0.313    0.327    0.353 1.074    52
# mean.s[2]    0.812   0.019    0.776    0.799    0.812    0.825    0.849 1.334    10
# deviance  1691.011   9.310 1671.344 1685.155 1691.051 1696.744 1710.081 1.100    31
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 40.5 and DIC = 1731.5
# DIC is an estimate of expected predictive error (lower deviance is better).

#############################################################################

pdf(file="TraceDens_ni10000_juvad_gyr.pdf", width = 6, height = 8)
plot(as.mcmc(lifedead))
dev.off()

# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
ni<-34000
nt <- 10 # “thinning”

# Call JAGS from R (BRT 80 min)
lifedead <- jags(jags.data, inits, parameters, "MRR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(lifedead, digit = 3)

pdf(file="TraceDens_ni20000_juvad_gyr.pdf", width = 6, height = 8)
plot(as.mcmc(lifedead))
dev.off()

#################### Old results p constant ########################################################
# Inference for Bugs model at "MRR.jags", fit using jags,
# 3 chains, each with 34000 iterations (first 24000 discarded), n.thin = 10
# n.sims = 3000 iterations saved
# mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# mean.eta     0.009   0.009    0.000    0.003    0.006    0.013    0.033 1.006  1300
# mean.p       0.026   0.004    0.019    0.023    0.026    0.029    0.036 1.100    25
# mean.r       0.140   0.009    0.123    0.134    0.140    0.146    0.158 1.003   790
# mean.s[1]    0.362   0.019    0.324    0.349    0.363    0.376    0.397 1.546     7
# mean.s[2]    0.825   0.015    0.796    0.815    0.826    0.835    0.853 1.052    50
# deviance  1758.349   7.440 1745.239 1752.969 1758.524 1763.088 1773.663 1.422     8

### I did not store the latent states... if I want to compare to true histories I need to do this
### I would need to run this even longer perhaps (need to recode this in STAN?)
######################################################################################################

################## New results p variable ############################################################
# Inference for Bugs model at "MRR.jags", fit using jags,
# 3 chains, each with 34000 iterations (first 24000 discarded), n.thin = 10
# n.sims = 3000 iterations saved
# mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# gamma        0.116   0.023    0.076    0.100    0.114    0.130    0.167 1.244    13
# mean.eta     0.009   0.009    0.000    0.003    0.006    0.013    0.031 1.002  2200
# mean.p       0.002   0.002    0.000    0.001    0.002    0.003    0.008 1.011   200
# mean.r       0.139   0.009    0.122    0.134    0.139    0.145    0.157 1.001  2700
# mean.s[1]    0.357   0.024    0.313    0.340    0.356    0.374    0.405 1.910     5
# mean.s[2]    0.825   0.014    0.796    0.816    0.825    0.835    0.853 1.007   830
# deviance  1711.315  10.295 1689.308 1704.927 1711.798 1718.271 1730.169 1.895     5
######################################################################################################

# MCMC settings
nc <- 3 #number of chains
nb <- 34000 # “burn in”
ni <- 54000
nt <- 10 # “thinning”

parameters <- c("mean.s", "mean.eta", "mean.r", "mean.p","z")
#https://stackoverflow.com/questions/16723036/strange-jags-parallel-error-avoiding-lazy-evaluation-in-function-call
#mrr <- do.call(jags.parallel, 
#           list(jags.data, inits, parameters, "MRR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd()))
mrr <- jags(jags.data, inits, parameters, "MRR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
head(print(mrr, digit = 3))

jags.sum<-mrr$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_MRR.txt")

### check out autojags if need be. 
# https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/autojags

### plot densities
library(mcmcplots)
denplot(mrr, c("mean.s", "mean.eta", "mean.r", "mean.p"))
### Trace plots
traplot(mrr, c("mean.s", "mean.eta", "mean.r", "mean.p"))

postsamples=cbind(mrr$BUGSoutput$sims.list$mean.eta,
                  mrr$BUGSoutput$sims.list$mean.p,
                  mrr$BUGSoutput$sims.list$mean.r,
                  mrr$BUGSoutput$sims.list$mean.s)

head(z)
timespan=1:ncol(z)
zmean=mrr$BUGSoutput$mean$z
head(zmean)
### Plotting the capture histories
library('RColorBrewer')
library('fields')
cols = brewer.pal(6,"RdBu")
rf <- colorRampPalette(cols)   # make colors
# voir aussi  http://www.sthda.com/french/wiki/couleurs-dans-r
image.plot(1:nrow(z),timespan,zmean,col=cols)

### Let's take a couple of individuals

CH[44,]
CH[100,]

z[44,]
z[100,]

zmean[44,]
zmean[100,]

## Perhaps have a look at the probability distributions even 
C1 = mrr$BUGSoutput$sims.array[,1,'mean.r']
C1 = mrr$BUGSoutput$sims.array[,1,'z'] # does not work
# and denplot(mrr, c("z")) does not work either

#First row of the array for the first chain
mrr$BUGSoutput$sims.array[1,1,1:10]
mrr$BUGSoutput$sims.array[1,1,100:110]
### They do not seem to be in the right order!!

### Recomputing properly the estimated latent state array
ptm <- proc.time()
zchain1 = array(NA,c(2000,nind,length(timespan)))
chain1 = mrr$BUGSoutput$sims.array[,1,]
mycolnames = colnames(chain1)
for (i in 1:nind){
  for (t in f[i]:length(timespan)){
    zname = paste("z[",i,",",t,"]",sep="")
    zindic<-grepl(zname,mycolnames,fixed=TRUE)
    #sum(zindic) # 1 if one is true, for checking
    zchain1[,i,t] = chain1[,zindic]
  }
} # long but at least we know it is correct...
proc.time()-ptm #

mstate1=apply(zchain1,2:3,function(x) mean(x,na.rm=TRUE))
image.plot(1:nrow(z),timespan,mstate1,col=cols)

#### Best to do that for the other chains too
#### And average over all

ptm <- proc.time()
zchain2 = array(NA,c(2000,nind,length(timespan)))
chain2 = mrr$BUGSoutput$sims.array[,2,]
zchain3 = array(NA,c(2000,nind,length(timespan)))
chain3 = mrr$BUGSoutput$sims.array[,3,]
mycolnames2 = colnames(chain2)
mycolnames3 = colnames(chain3)
for (i in 1:nind){
  for (t in f[i]:length(timespan)){
    zname = paste("z[",i,",",t,"]",sep="")
    zindic2<-grepl(zname,mycolnames2,fixed=TRUE)
    zindic3<-grepl(zname,mycolnames3,fixed=TRUE)
    #sum(zindic) # 1 if one is true, for checking
    zchain2[,i,t] = chain2[,zindic2]
    zchain3[,i,t] = chain3[,zindic3]
  }
} # long but at least we know it is correct...
proc.time()-ptm #

mstate2=apply(zchain2,2:3,function(x) mean(x,na.rm=TRUE))
image.plot(1:nrow(z),timespan,mstate2,col=cols)

mstate3=apply(zchain3,2:3,function(x) mean(x,na.rm=TRUE))
image.plot(1:nrow(z),timespan,mstate3,col=cols)

mstate=(mstate1+mstate2+mstate3)/3
image.plot(1:nrow(z),timespan,mstate,col=cols)
### I could also plot the states for which we know that Pr(1) = 0. 
### That's another way to represent the uncertainty

mstate_adults = mstate[stagemarking==2,]
image.plot(1:nrow(mstate_adults),timespan,mstate_adults,col=cols)

mstate_juvs = mstate[stagemarking==1,]
image.plot(1:nrow(mstate_juvs),timespan,mstate_juvs,col=cols)

### So that's quite nice, a number of juveniles keep on surviving in this model.

#################################################
### Let's do that for individual 44
i=44
for (t in f[i]:length(timespan)){
  zname = paste("z[",i,",",t,"]",sep="")
  zindic<-grepl(zname,colnames(chain1),fixed=TRUE)
  #sum(zindic) # 1 if one is true, for checking
  zchain1[,i,t] = chain1[,zindic]
}
# Examining the first chain for individual 44
table(zchain1[,44,]) ## this does not fully work
zchain1[1,44,23]
zchain1[1,44,28] # just to check

### Mean state recomputed
apply(zchain1[,44,],2,mean)
### Basic JAGS output for mean state
zmean[44,]

z[44,]
### Following is not working...
apply(zchain1[,44,],2, function(x) barplot(table(x)))
#########################################################


### stagemarking will give which individuals are taken as adults vs juveniles

write.csv(mstate,"meanLatentStates.csv")

cols = brewer.pal(15,"YlOrRd")
rf <- colorRampPalette(cols)   # make colors
# voir aussi  http://www.sthda.com/french/wiki/couleurs-dans-r

freq1_chain1=apply(zchain1,2:3,function(x) sum(x==1)/length(x))
image.plot(1:nrow(z),timespan,freq1_chain1,col=cols,xlab="Individuals",ylab="Time",main="Pr(alive)")

freq1_chain2=apply(zchain2,2:3,function(x) sum(x==1)/length(x))
image.plot(1:nrow(z),timespan,freq1_chain2,col=cols,xlab="Individuals",ylab="Time",main="Pr(alive)")

freq1_chain3=apply(zchain3,2:3,function(x) sum(x==1)/length(x))
image.plot(1:nrow(z),timespan,freq1_chain3,col=cols,xlab="Individuals",ylab="Time",main="Pr(alive)")


freq1 = (freq1_chain1+freq1_chain2+freq1_chain3)/3
image.plot(1:nrow(z),1973:2018,freq1,col=cols,xlab="Individuals",ylab="Time",main="Pr(alive)")

write.csv(freq1,"PrAlive.csv")

#######################################################
### Ideas for faster code for latent states computing
########################################################

### 1. recup ttes les valeurs znames pour tout i pour tout t
### 2. Sort sur i et t
### 3. re-order dans la matrice

strsplit(mycolnames,split=c("z[","]"),fixed=T)
strsplit(mycolnames,split=c("z["),fixed=T)

tmp1=strsplit(mycolnames,split="z[",fixed=T)
tmp2=lapply(tmp1,function(x) strsplit(x[2],",",fixed=T))

### Other ideas
### 1. Parse column names (vector form)
### 2. reorder
