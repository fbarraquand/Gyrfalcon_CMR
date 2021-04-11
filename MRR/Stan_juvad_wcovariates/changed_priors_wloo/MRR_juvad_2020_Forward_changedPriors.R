## Stan code for Juvenile-Adult MRR model designed for the Icelandic Gyrfalcon
## FB 25/03/2020 -- edited for covariates 20/04/2020
## with modification of priors 11/03/2021
## computes PSIS-LOO 13/03/2021
# .libPaths("/home/frederic/myRpackages/")

library(ggplot2)
library(xtable)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)
library(bridgesampling)
library(loo)

## function
source('../../simul.R')
source('../../figlabel.R')

# -------------------------------------------------
# Loading data

CH=read.csv(file = "../../../data/CaptureHistories1981til2020.csv")
CH[,1]<-NULL
CH=as.matrix(CH)
image(CH)
ncol(CH) #39

get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)
n.occasions <- ncol(CH)
n.states <- 3
n.obs <- 3

stagemarking = read.csv(file ="../../../data/stagemarking1981til2020.csv")
stagemarking = stagemarking[,2]
table(stagemarking) #OK

# Matrix giving the stage of all individuals that are marked
xj <- matrix(NA, ncol = n.occasions-1, nrow = length(f))
for (i in 1:length(f)){
  for (t in f[i]:(n.occasions-1)){
    if (f[i]>(n.occasions-1)) next
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

# Removing the elements that verify f[i]<n_occasions
indices=which(f<n.occasions)
f=f[indices]
CH = CH[indices,]
xj = xj[indices,]
n.occasions=length(f)

# Get ptarmigan data
ptarmigan = read.csv(file = "../../../data/Ptarmigan_Data.csv")
mean(ptarmigan$Mean.density) #for the priors
nrow(ptarmigan)

##########################################################
## Collating other covariates -- adding in weather. 
##########################################################

#standardize previous variable for an easier comparison of coefficients

prey_abund = (ptarmigan$Mean.density - mean(ptarmigan$Mean.density))/sd(ptarmigan$Mean.density)

#spatially averaged weather
av.weather = read.csv(file = "../../../data/weather/average_weatherNEIceland.csv")
#last two variables computed over winter of the current for month 1-2-3 and previous for month 10-11-12 years. 
temp=lograin=rep(0,39)
for (year in 1981:2019){
  #we now temporally average
  winter_temp = av.weather$temp[(av.weather$year==(year-1))&(av.weather$month %in%c(10,11,12))]
  winter_temp = c(winter_temp,av.weather$temp[(av.weather$year==(year))&(av.weather$month %in%c(1,2,3))])
  temp[year-1981+1] = mean(winter_temp)
  #for precipitation
  winter_rain = av.weather$logRainfall[(av.weather$year==(year-1))&(av.weather$month %in%c(10,11,12))]
  winter_rain = c(winter_rain,av.weather$logRainfall[(av.weather$year==(year))&(av.weather$month %in%c(1,2,3))])
  lograin[year-1981+1] = mean(winter_rain)
}

temp = (temp - mean(temp))/sd(temp)
lograin = (lograin - mean(lograin))/sd(lograin)

### Correlations? 
covar = data.frame(prey_abund,temp,lograin)
plot(covar)
par(pch=20)
tmax=nrow(covar)
plot(1:tmax,prey_abund,type="o",ylim=c(-3,3))
lines(1:tmax,temp,type="o",col="red")
lines(1:tmax,lograin,type="o",col="blue")
### No apparent correlation. 

  
# Bundle data
#stan.data.R <- list(y = CH, f = f,  xj=xj, n_occasions = dim(CH)[2],nind = dim(CH)[1], prey_abund = prey_abund, temp = temp,lograin=lograin)
stan.data.R <- list(y = CH, f = f,  xj=xj, n_occasions = dim(CH)[2],nind = dim(CH)[1], prey_abund = prey_abund) #first check

######################################################################################################
### First re-running the model with the prey-only covariate (Model B) with improved parameterization
######################################################################################################

## Parameters monitored
params <- c("mean_s2", "mean_eta", "mean_r", "mean_p","beta","mu_juvsurv","log_lik")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

inits <- function(){list(mean_s2 = runif(1, 0.2, 1), mean_eta = runif(1, 0.05, 0.4), 
                         mean_p = runif(1, 0.05, 0.4), mean_r = runif(1, 0.05, 0.4))}  #z = ld.init(CH, f)

mrr_juvad2 <- stan("mrr_w1covar.stan",
                  data = stan.data.R , init = inits, pars = params,
                  chains = nc, iter = ni, warmup = nb, thin = nt,
                  seed = 1)
print(mrr_juvad2, digits = 3,max=100) #similar to previous
# Inference for Stan model: mrr_w1covar.
# 3 chains, each with iter=2000; warmup=1000; thin=1; 
# post-warmup draws per chain=1000, total post-warmup draws=3000.
# 
# mean se_mean    sd      2.5%       25%       50%       75%     97.5% n_eff  Rhat
# mean_s2           0.828   0.000 0.020     0.787     0.815     0.829     0.842     0.865  3413 1.000
# mean_eta          0.009   0.000 0.009     0.000     0.003     0.006     0.012     0.031  3247 1.000
# mean_r            0.137   0.000 0.009     0.120     0.131     0.137     0.143     0.155  3381 0.999
# mean_p            0.021   0.000 0.005     0.013     0.017     0.020     0.023     0.030  2822 1.000
# beta              0.077   0.002 0.097    -0.111     0.009     0.077     0.144     0.272  3642 0.999
# mu_juvsurv       -0.476   0.002 0.129    -0.733    -0.560    -0.476    -0.388    -0.226  3076 0.999
# log_lik[1]       -7.609   0.004 0.216    -8.044    -7.749    -7.603    -7.463    -7.191  2780 1.000
# log_lik[2]       -0.173   0.000 0.011    -0.194    -0.180    -0.173    -0.165    -0.152  3636 0.999
# log_lik[3]       -0.173   0.000 0.011    -0.194    -0.180    -0.173    -0.165    -0.152  3636 0.999
# log_lik[4]       -0.173   0.000 0.011    -0.194    -0.180    -0.173    -0.165    -0.152  3636 0.999
# [ reached getOption("max.print") -- omitted 1727 rows ]

log_lik_2 <- extract_log_lik(mrr_juvad2, merge_chains = FALSE)
loo_juvad2<-loo(mrr_juvad2, save_psis = TRUE)
print(loo_juvad2)
########################## Old one ###################################
# Computed from 3000 by 1791 log-likelihood matrix
# Estimate    SE
# elpd_loo  -1335.0  69.9
# p_loo         6.0   0.5
# looic      2670.0 139.8
# ------
#   Monte Carlo SE of elpd_loo is NA.
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     1730  96.6%   2297      
# (0.5, 0.7]   (ok)          0   0.0%   <NA>      
#   (0.7, 1]   (bad)         0   0.0%   <NA>      
#   (1, Inf)   (very bad)   61   3.4%   1500      
# See help('pareto-k-diagnostic') for details.
### Can it be due to the fact that I had to se log-lik = 0 for some individuals? 
### Quite likely since there is cross-validation! 

######################## New one ###################################
### estimates removing the ones for which I had log-lik = 0
# Computed from 3000 by 1730 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -1335.0  69.7
# p_loo         6.0   0.5
# looic      2670.0 139.3
# ------
#   Monte Carlo SE of elpd_loo is 0.0.
# 
# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.

params2 <- extract(mrr_juvad2)
save(params2, file ="param_logisticv2_chains.RData")

## Real prey data ptarmigan$Mean.density
min(prey_abund) # -1.439433
max(prey_abund) #  2.580465
minp = min(prey_abund)-1
maxp = max(prey_abund)+1
prey_abund_vec=seq(from = minp, to = maxp, length.out = 200)

nsamples = length(params2$beta)
F=matrix(0,nrow  = nsamples,ncol=length(prey_abund_vec))
m=matrix(0,nrow = 1,ncol=length(prey_abund_vec))
q=matrix(0,nrow = 2,ncol=length(prey_abund_vec))

for (p in 1:length(prey_abund_vec)){
  F[,p] = 1/(1+exp(- (params2$mu_juvsurv + params2$beta*prey_abund_vec[p])))
  m[p] = mean(F[,p])
  q[,p] = quantile(F[,p],c(0.05,0.95))
}
pdf(file = "FigLogisticPtarmigan_v2.pdf",width = 5,height=8)
par(mfrow=c(2,1))
plot(prey_abund_vec,m,col="red",type="l",lwd=2,xlab="",ylab="Predicted juvenile survival",xlim=c(minp,maxp),ylim=c(0.15,0.6))
lines(prey_abund_vec,q[1,],col="black")
lines(prey_abund_vec,q[2,],col="black")
fig_label("A")
### This is a bit puzzling. On average all these values are above the estimate without the covariate. 
d <- density(prey_abund)
plot(d, xlab = "Ptarmigan abundance (stdized) ",ylab="Kernel density of abundance",main="",xlim=c(minp,maxp),ylim=c(0,0.75))
polygon(d, col="red", border="black") 
fig_label("B")
dev.off()

#################################################################################
### Model C with 3 covariates -- adding weather (temperature and precipitation)
#################################################################################

# Bundle data
stan.data.R <- list(y = CH, f = f,  xj=xj, n_occasions = dim(CH)[2],nind = dim(CH)[1], prey_abund = prey_abund, temp = temp,lograin=lograin)

## Parameters monitored
params <- c("mean_s2", "mean_eta", "mean_r", "mean_p","beta","mu_juvsurv","log_lik")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

inits <- function(){list(mean_s2 = runif(1, 0.2, 1), mean_eta = runif(1, 0.05, 0.4), 
                         mean_p = runif(1, 0.05, 0.4), mean_r = runif(1, 0.05, 0.4))} 

mrr_juvad3 <- stan("mrr_w3covar.stan",
                   data = stan.data.R , init = inits, pars = params,
                   chains = nc, iter = ni, warmup = nb, thin = nt,
                   seed = 1)
print(mrr_juvad3, digits = 3)

loo_juvad3<-loo(mrr_juvad3, save_psis = TRUE)
print(loo_juvad3)
# Computed from 3000 by 1730 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -1334.1  69.7
# p_loo         8.3   0.7
# looic      2668.3 139.4
# ------
#   Monte Carlo SE of elpd_loo is 0.1.
# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.

comp <- loo_compare(loo_juvad2, loo_juvad3)
print(comp, digits = 2)
#        elpd_diff se_diff
# model2  0.00      0.00  
# model1 -0.86      2.61  

# Following https://avehtari.github.io/modelselection/CV-FAQ.html
# It seems that the differences are small (and very variable)

params3 <- extract(mrr_juvad3)
save(params3, file ="param_logistic3covar_chains.RData")
# Inference for Stan model: mrr_w3covar.
# 3 chains, each with iter=2000; warmup=1000; thin=1; 
# post-warmup draws per chain=1000, total post-warmup draws=3000.
# 
# mean se_mean    sd      2.5%       25%       50%       75%     97.5% n_eff  Rhat
# mean_s2        0.833   0.000 0.018     0.797     0.822     0.834     0.846     0.869  3787 1.000
# mean_eta       0.009   0.000 0.009     0.000     0.002     0.006     0.012     0.033  4368 1.000
# mean_r         0.139   0.000 0.009     0.122     0.133     0.139     0.145     0.157  4255 1.000
# mean_p         0.029   0.000 0.005     0.020     0.026     0.029     0.032     0.040  2635 1.000
# beta[1]        0.276   0.003 0.134     0.014     0.189     0.275     0.359     0.548  1765 1.000
# beta[2]        0.231   0.003 0.120    -0.005     0.152     0.228     0.311     0.470  2038 1.000
# beta[3]        0.270   0.003 0.126     0.020     0.188     0.270     0.351     0.520  1962 1.001
# mu_juvsurv    -0.504   0.002 0.113    -0.730    -0.578    -0.503    -0.429    -0.278  2557 1.000
# lp__       -1364.566   0.059 2.106 -1369.553 -1365.811 -1364.224 -1362.973 -1361.521  1296 1.002
# 
# Samples were drawn using NUTS(diag_e) at Thu Mar 11 13:07:30 2021.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at 
#                                                                   convergence, Rhat=1).
### For extraction of the results
xtable(summary(mrr_juvad3)$summary,digits=3)

######################################################################
### Violin plot of the coefficients for model C (if we need them)
###################################################################### 
## extracting data 
# beta=data.frame(params3$beta)
# colnames(beta)=c("beta1","beta2","beta3")
# library("tidyr")
# beta %>% pivot_longer(names_to = "coefficient", values_to = "value") ## putting into the right format for ggplot2
# Violin plot basics
# http://www.sthda.com/french/wiki/ggplot2-violin-plot-guide-de-demarrage-rapide-logiciel-r-et-visualisation-de-donnees
# https://www.benjaminackerman.com/post/2019-03-08-equation_labels/
# https://www.r-graph-gallery.com/95-violin-plot-with-ggplot2.html
### (All this data-wrangling for so little is getting ridiculous)
# Red for prey, orange for temp, blue for lograin
# Use for legend c(expression(beta[1]),expression(beta[2]),expression(beta[3]))
# p <- ggplot(beta, aes(y=value,x=coefficient)) +geom_violin()

#############################################################################
### Using bridgesampling to compare models through Bayes factors
#############################################################################
### See e.g. https://www.r-bloggers.com/2019/05/bayesian-modeling-using-stan-a-case-study/
### Also https://www.jstatsoft.org/article/view/v092i10 for the theory and implementation

# Bayes factor for model B over model C
bridge_juvad2 = bridge_sampler(mrr_juvad2)
bridge_juvad3 = bridge_sampler(mrr_juvad3)
BF = bf(bridge_juvad2,bridge_juvad3) # or bayes_factor()? 
BF # Estimated Bayes factor in favor of bridge_juvad2 over bridge_juvad3: 3.55905 ### WTF that's different from before, why? 
# Check errors
print(error_measures(bridge_juvad2)$percentage) #1%
print(error_measures(bridge_juvad3)$percentage) #1%
# compute posterior model probabilities (assuming equal prior model probabilities)
post1 <- post_prob(bridge_juvad2,bridge_juvad3)
print(post1)

### Let's do this with different options to ensure robustness
## Bayes factor for model B over model C
bridge_juvadB = bridge_sampler(mrr_juvad2,method = "warp3")
bridge_juvadC = bridge_sampler(mrr_juvad3,method = "warp3")
BFbis = bf(bridge_juvadB,bridge_juvadC) # 
BFbis # Estimated Bayes factor in favor of bridge_juvadB over bridge_juvadC: 3.55798

#### ---------------------------------------------------------------------------#####
#### Comparison to "null" model - Model A with only juvenile vs adult survival 
#### ---------------------------------------------------------------------------#####

# Bundle data (needed again? probably not) 
stan.data.R <- list(y = CH, f = f,  xj=xj, n_occasions = dim(CH)[2],nind = dim(CH)[1])
inits <- function(){list(mean_s = runif(2, 0.2, 1), mean_eta = runif(1, 0.05, 0.4), 
                         mean_p = runif(1, 0.05, 0.4), mean_r = runif(1, 0.05, 0.4))}  
## Parameters monitored
params <- c("mean_s", "mean_eta", "mean_r", "mean_p","log_lik")
## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3
## Call Stan from R
mrr_juvad0 <- stan("mrr0.stan",
                  data = stan.data.R , init = inits, pars = params,
                  chains = nc, iter = ni, warmup = nb, thin = nt,
                  seed = 1)
print(mrr_juvad0, digits = 3)
# Inference for Stan model: mrr0.
# 3 chains, each with iter=2000; warmup=1000; thin=1; 
# post-warmup draws per chain=1000, total post-warmup draws=3000.
# 
# mean se_mean    sd      2.5%       25%       50%       75%     97.5% n_eff  Rhat
# mean_s[1]     0.384   0.001 0.029     0.330     0.364     0.384     0.404     0.442  2858 1.000
# mean_s[2]     0.829   0.000 0.020     0.788     0.815     0.829     0.843     0.866  3422 1.000
# mean_eta      0.008   0.000 0.008     0.000     0.003     0.006     0.012     0.030  3587 1.001
# mean_r        0.137   0.000 0.009     0.120     0.131     0.137     0.143     0.154  3566 1.000
# mean_p        0.021   0.000 0.004     0.013     0.018     0.021     0.024     0.031  2812 0.999
# lp__      -1351.983   0.041 1.584 -1355.784 -1352.829 -1351.653 -1350.805 -1349.870  1499 1.000
# 
# Samples were drawn using NUTS(diag_e) at Thu Mar 11 18:37:24 2021.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at 
#                                                                   convergence, Rhat=1).

# Bayes factor of model B over model A
bridge_juvad2 = bridge_sampler(mrr_juvad2) ## (we redo-it just to check)
bridge_juvad0 = bridge_sampler(mrr_juvad0)
BF20 = bf(bridge_juvad2,bridge_juvad0) 
BF20 ### Estimated Bayes factor in favor of bridge_juvad2 over bridge_juvad0: 0.11289 instead of 0.00911 = 10^{-2}

# Bayes factor of model A over model B
BF02 = bf(bridge_juvad0,bridge_juvad2) # just to check computation
BF02 ## 8.85789 close to 10 still decisive -- combination to prior odds ? 

# Bayes factor of model A over model C
BF03 = bf(bridge_juvad0,bridge_juvad3) # 31.10027 (checking this was consistent)
BF03

loo_juvad0<-loo(mrr_juvad0, save_psis = TRUE)
# print(loo_juvad0)
# Computed from 3000 by 1730 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -1334.0  69.6
# p_loo         4.6   0.4
# looic      2668.1 139.2
# ------
#   Monte Carlo SE of elpd_loo is 0.0.
# 
# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.
comp <- loo_compare(loo_juvad0, loo_juvad2,loo_juvad3)
print(comp, digits = 2)

################### Some thoughts on possible prior and posterior odds ###############################
## Let say Pr(A)/Pr(B) = 0.1/0.9 = 0.11 a priori so that posterior odds = 0.11 x 9.02 approx 1

######### A reality check -- Model S
# We consider a model where adult and juvenile survival is equal -- obviously, a stupid model

## Parameters monitored
params <- c("mean_s", "mean_eta", "mean_r", "mean_p","log_lik")
inits <- function(){list(mean_s = runif(1, 0.2, 1), mean_eta = runif(1, 0.05, 0.4), 
                         mean_p = runif(1, 0.05, 0.4), mean_r = runif(1, 0.05, 0.4))}  
## Call Stan from R
mrr_0 <- stan("mrr00.stan",   ## it has only been changed so that mean_s is now a scalar not a vector. 
                   data = stan.data.R , init = inits, pars = params,
                   chains = nc, iter = ni, warmup = nb, thin = nt,
                   seed = 1)
print(mrr_0, digits = 3)

bridge_juvadS = bridge_sampler(mrr_0)  # S for "stupid" model since we know that adult and juvenile survival probs are different
BF0S = bf(bridge_juvad0,bridge_juvadS) # comparison to the null juvenile adult model
BF0S
### Estimated Bayes factor in favor of bridge_juvad0 over bridge_juvadS: 3853537887457871237128126464.00000 = 10^27 // same conclusion. 

loo_00<-loo(mrr_0, save_psis = TRUE)
# print(loo_juvad0)
# Computed from 3000 by 1730 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -1334.0  69.6
# p_loo         4.6   0.4
# looic      2668.1 139.2
# ------
#   Monte Carlo SE of elpd_loo is 0.0.
# 
# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.
comp <- loo_compare(loo_00,loo_juvad0, loo_juvad2,loo_juvad3) ## rankings not really different. 
print(comp, digits = 2)

