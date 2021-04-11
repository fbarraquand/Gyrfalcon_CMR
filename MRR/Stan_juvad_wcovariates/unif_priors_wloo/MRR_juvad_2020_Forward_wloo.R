## Stan code for Juvenile-Adult MRR model designed for the Icelandic Gyrfalcon
## FB 25/03/2020 -- edited for covariates 20/04/2020
## with modification of priors 11/03/2021
## computes PSIS-LOO 13/03/2021
## back to old priors 29/03/2021
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
# mean_s2           0.828   0.000 0.021     0.784     0.815     0.829     0.843     0.866  3832 1.000
# mean_eta          0.009   0.000 0.009     0.000     0.003     0.006     0.012     0.031  3268 1.000
# mean_r            0.137   0.000 0.009     0.120     0.131     0.137     0.143     0.155  3796 0.999
# mean_p            0.021   0.000 0.004     0.013     0.018     0.020     0.024     0.031  3240 1.001
# beta              0.079   0.002 0.094    -0.110     0.016     0.080     0.142     0.268  3661 1.000
# mu_juvsurv       -0.491   0.002 0.126    -0.739    -0.575    -0.490    -0.406    -0.240  3397 1.000

#log_lik_2 <- extract_log_lik(mrr_juvad2, merge_chains = FALSE)
loo_juvad2<-loo(mrr_juvad2, save_psis = TRUE)
print(loo_juvad2)
# Computed from 3000 by 1730 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -1334.9  69.7
# p_loo         5.8   0.5
# looic      2669.8 139.4
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
print(mrr_juvad3, digits = 3,max=100)
# Inference for Stan model: mrr_w3covar.
# 3 chains, each with iter=2000; warmup=1000; thin=1; 
# post-warmup draws per chain=1000, total post-warmup draws=3000.
# 
# mean se_mean    sd      2.5%       25%       50%       75%     97.5% n_eff  Rhat
# mean_s2           0.831   0.000 0.020     0.791     0.818     0.832     0.845     0.869  3309 1.000
# mean_eta          0.009   0.000 0.008     0.000     0.002     0.006     0.012     0.032  3948 0.999
# mean_r            0.138   0.000 0.009     0.122     0.133     0.138     0.144     0.156  3740 1.000
# mean_p            0.021   0.000 0.005     0.013     0.017     0.020     0.023     0.030  2664 1.000
# beta[1]           0.339   0.003 0.146     0.057     0.236     0.337     0.440     0.628  1873 1.000
# beta[2]           0.253   0.003 0.131    -0.004     0.161     0.251     0.344     0.511  2460 1.000
# beta[3]           0.286   0.003 0.135     0.033     0.191     0.285     0.377     0.560  2124 1.000
# mu_juvsurv       -0.484   0.002 0.127    -0.735    -0.571    -0.483    -0.398    -0.236  2876 1.000

loo_juvad3<-loo(mrr_juvad3, save_psis = TRUE)
print(loo_juvad3)
# Computed from 3000 by 1730 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -1334.2  69.7
# p_loo         8.3   0.7
# looic      2668.5 139.3
# ------
#   Monte Carlo SE of elpd_loo is 0.1.
# 
# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.

comp <- loo_compare(loo_juvad2, loo_juvad3)
print(comp, digits = 2)
# elpd_diff se_diff
# model2  0.00      0.00  
# model1 -0.65      2.62  

# Following https://avehtari.github.io/modelselection/CV-FAQ.html
# It seems that the differences are small (and very variable)

params3 <- extract(mrr_juvad3)
save(params3, file ="param_logistic3covar_chains.RData")
xtable(summary(mrr_juvad3)$summary,digits=3)

######################################################################
### Violin plot of the coefficients for model C (if we need them)
###################################################################### 
## extracting data 
beta=data.frame(params3$beta)
colnames(beta)=c("beta1","beta2","beta3")
## we do by hand what pivot_longer() does
# library("tidyr")
# beta %>% pivot_longer(names_to = "coefficient", values_to = "value") ## putting into the right format for ggplot2
kd = nrow(beta)
beta_long = c(beta$beta1,beta$beta2,beta$beta3)
coeff = c(rep("beta1",kd),rep("beta2",kd),rep("beta3",kd))
beta_long = data.frame(coeff,beta_long)
names(beta_long) = c("coefficient","value")
# Violin plot basics
# http://www.sthda.com/french/wiki/ggplot2-violin-plot-guide-de-demarrage-rapide-logiciel-r-et-visualisation-de-donnees
# https://www.benjaminackerman.com/post/2019-03-08-equation_labels/
# https://www.r-graph-gallery.com/95-violin-plot-with-ggplot2.html
### (All this data-wrangling for so little is getting ridiculous)
# Red for prey, orange for temp, blue for lograin
# Use for legend
pdf(file = "ViolinPlotCoeffs_modelC.pdf",width = 8,height=6)
xlabels = c(expression(beta[1]),expression(beta[2]),expression(beta[3]))
p <- ggplot(beta_long, aes(y=value,x=coefficient, fill=coefficient)) + geom_violin(trim=FALSE) + theme_bw() + 
  stat_summary(
  mapping = aes(x = coefficient, y = value),
  fun.min = function(z) { quantile(z,0.025) },
  fun.max = function(z) { quantile(z,0.975) },
  fun = mean) + 
scale_x_discrete(labels = xlabels) + scale_fill_manual(values=c("#FF0000", "#FFA500", "#00BFFF"))
plot(p)
dev.off()
#############################################################################
### Using bridgesampling to compare models through Bayes factors
#############################################################################
### See e.g. https://www.r-bloggers.com/2019/05/bayesian-modeling-using-stan-a-case-study/
### Also https://www.jstatsoft.org/article/view/v092i10 for the theory and implementation

# Bayes factor for model B over model C
bridge_juvad2 = bridge_sampler(mrr_juvad2)
bridge_juvad3 = bridge_sampler(mrr_juvad3)
BF = bf(bridge_juvad2,bridge_juvad3) # or bayes_factor()? 
BF # Estimated Bayes factor in favor of bridge_juvad2 over bridge_juvad3: 3.53911 ### Lower than without removing the last year
# Check errors
print(error_measures(bridge_juvad2)$percentage) #1%
print(error_measures(bridge_juvad3)$percentage) #1%
# compute posterior model probabilities (assuming equal prior model probabilities)
post1 <- post_prob(bridge_juvad2,bridge_juvad3)
print(post1) #0.7796923     0.2203077 

### Let's do this with different options to ensure robustness
## Bayes factor for model B over model C
bridge_juvadB = bridge_sampler(mrr_juvad2,method = "warp3")
bridge_juvadC = bridge_sampler(mrr_juvad3,method = "warp3")
BFbis = bf(bridge_juvadB,bridge_juvadC) # 
BFbis # Estimated Bayes factor in favor of bridge_juvadB over bridge_juvadC: 3.54027

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
print(mrr_juvad0, digits = 3,max=100)
# Inference for Stan model: mrr0.
# 3 chains, each with iter=2000; warmup=1000; thin=1; 
# post-warmup draws per chain=1000, total post-warmup draws=3000.
# 
# mean se_mean    sd      2.5%       25%       50%       75%     97.5% n_eff  Rhat
# mean_s[1]         0.386   0.001 0.028     0.333     0.366     0.384     0.404     0.441  2902 1.000
# mean_s[2]         0.830   0.000 0.020     0.790     0.816     0.831     0.844     0.868  3599 1.000
# mean_eta          0.009   0.000 0.009     0.000     0.003     0.006     0.013     0.033  2935 1.000
# mean_r            0.137   0.000 0.008     0.121     0.131     0.137     0.143     0.154  2865 1.000
# mean_p            0.021   0.000 0.005     0.013     0.018     0.021     0.024     0.031  2594 1.000

# Bayes factor of model B over model A
bridge_juvad2 = bridge_sampler(mrr_juvad2) ## (we redo-it just to check)
bridge_juvad0 = bridge_sampler(mrr_juvad0)
BF20 = bf(bridge_juvad2,bridge_juvad0) 
BF20 ### Estimated Bayes factor in favor of bridge_juvad2 over bridge_juvad0: 0.04492 (before, 0.11289, and 0.00911 = 10^{-2})

# Bayes factor of model A over model B
BF02 = bf(bridge_juvad0,bridge_juvad2) # just to check computation
BF02 ## 22.26222 close to 10 still decisive -- combination to prior odds ? 

# Bayes factor of model A over model C
BF03 = bf(bridge_juvad0,bridge_juvad3) # 78.24 (checking this was consistent)
BF03

loo_juvad0<-loo(mrr_juvad0, save_psis = TRUE)
print(loo_juvad0)
# Computed from 3000 by 1730 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -1334.1  69.5
# p_loo         4.6   0.4
# looic      2668.1 139.1
# ------
#   Monte Carlo SE of elpd_loo is 0.0.
# 
# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.
comp <- loo_compare(loo_juvad0, loo_juvad2,loo_juvad3)
print(comp, digits = 2)
# elpd_diff se_diff
# model1  0.00      0.00  
# model3 -0.16      2.78  
# model2 -0.81      0.92  

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
### Estimated Bayes factor in favor of bridge_juvad0 over bridge_juvadS: 1696739028080432700425502720.00000 = 10^27 // same conclusion. 

loo_00<-loo(mrr_0, save_psis = TRUE)
print(loo_00)
# Computed from 3000 by 1730 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -1400.6  73.8
# p_loo         7.4   0.8
# looic      2801.1 147.6
# ------
#   Monte Carlo SE of elpd_loo is 0.1.
# 
# All Pareto k estimates are good (k < 0.5).
# See help('pareto-k-diagnostic') for details.
comp <- loo_compare(loo_00,loo_juvad0, loo_juvad2,loo_juvad3) ## rankings not really different. 
# print(comp, digits = 2)
# elpd_diff se_diff
# model2   0.00      0.00 
# model4  -0.16      2.78 
# model3  -0.81      0.92 
# model1 -66.48     10.59 

