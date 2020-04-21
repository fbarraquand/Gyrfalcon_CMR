## Stan code for Juvenile-Adult MRR model designed for the Icelandic Gyrfalcon
## FB 25/03/2020 -- edited for covariates 20/04/2020
# .libPaths("/home/frederic/myRpackages/")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## function
source('../simul.R')
source('../figlabel.R')

# -------------------------------------------------
# Loading data

CH=read.csv(file = "../../data/CaptureHistories1981til2020.csv")
CH[,1]<-NULL
CH=as.matrix(CH)
image(CH)
ncol(CH) #39

get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)
n.occasions <- ncol(CH)
n.states <- 3
n.obs <- 3

stagemarking = read.csv(file ="../../data/stagemarking1981til2020.csv")
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

# Get ptarmigan data
ptarmigan = read.csv(file = "../../data/Ptarmigan_Data.csv")
mean(ptarmigan$Mean.density) #for the priors
nrow(ptarmigan)

# Bundle data
stan.data.R <- list(y = CH, f = f,  xj=xj, n_occasions = dim(CH)[2],nind = dim(CH)[1], prey_abund = ptarmigan$Mean.density)


inits <- function(){list(mean_s2 = runif(1, 0.2, 1), mean_eta = runif(1, 0.05, 0.4), 
                         mean_p = runif(1, 0.05, 0.4), mean_r = runif(1, 0.05, 0.4))}  #z = ld.init(CH, f)
## Looks like I don't need to initialize the states...


## Parameters monitored
params <- c("mean_s2", "mean_eta", "mean_r", "mean_p","gamma_prey","mu_prey")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 2

## Call Stan from R
#mrr_juvad <- stan("mrr.stan",
#                 data = stan.data.R , init = inits, pars = params,
#                 chains = nc, iter = ni, warmup = nb, thin = nt,
#                 seed = 1,
#                 open_progress = TRUE) #does not work on my desktop PC for some reason. 
#print(mrr_juvad, digits = 3)

## Call Stan from R -- that finally works. 
mrr_juvad <- stan("mrr.stan",
                  data = stan.data.R , init = inits, pars = params,
                  chains = nc, iter = ni, warmup = nb, thin = nt,
                  seed = 1)
print(mrr_juvad, digits = 3)
#saveRDS(mrr_juvad, "fit1.rds") #save(mrr_juvad, file ="fit1.RData")



## Other way of calling it
#mrr_juvad <-stan_model("mrr.stan",verbose=TRUE)
#fit.mrr_juvad <- sampling(mrr_juvad, data = stan.data.R, iter = 1000, init = inits,chains = 2, cores = 2)
#fit.mrr_juvad <- sampling(mrr_juvad, data = stan.data.R, iter = 1000, init = inits,chains = 2)

### plot Juvenile survival = f(ptarmigan density?)? 
params <- extract(mrr_juvad)
save(params, file ="param_chains.RData")

## Real prey data ptarmigan$Mean.density
minp = 0 #min(ptarmigan$Mean.density)
maxp = 15 #max(ptarmigan$Mean.density)
prey_abund=seq(from = minp, to = maxp, length.out = 200)

nsamples = length(params$gamma_prey)
F=matrix(0,nrow  = nsamples,ncol=length(prey_abund))
m=matrix(0,nrow = 1,ncol=length(prey_abund))
q=matrix(0,nrow = 2,ncol=length(prey_abund))

for (p in 1:length(prey_abund)){
  F[,p] = 1/(1+exp(-params$gamma_prey*(prey_abund[p] - params$mu_prey) ))
  m[p] = mean(F[,p])
  q[,p] = quantile(F[,p],c(0.05,0.95))
}
pdf(file = "FigLogisticPtarmigan.pdf",width = 5,height=8)
par(mfrow=c(2,1))
plot(prey_abund,m,col="red",type="l",lwd=2,xlab="",ylab="Predator juvenile survival")
lines(prey_abund,q[1,],col="black")
lines(prey_abund,q[2,],col="black")
fig_label("A")
### This is a bit puzzling. On average all these values are above the estimate without the covariate. 
d <- density(ptarmigan$Mean.density)
plot(d, xlab = "Ptarmigan abundance",ylab="Kernel density of abundance",main="")
polygon(d, col="red", border="black") 
fig_label("B")
dev.off()

# or dlogis(ptarmigan$Mean.density,mu_prey,1/gamma_prey)); 

## Doing it with three covariates -- adding in weather. 

#standardize previous variables

prey_abund = (ptarmigan$Mean.density - mean(ptarmigan$Mean.density))/sd(ptarmigan$Mean.density)

#spatially averaged weather
av.weather = read.csv(file = "../../data/weather/average_weatherNEIceland.csv")
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
  
# Bundle data
#stan.data.R <- list(y = CH, f = f,  xj=xj, n_occasions = dim(CH)[2],nind = dim(CH)[1], prey_abund = prey_abund, temp = temp,lograin=lograin)
stan.data.R <- list(y = CH, f = f,  xj=xj, n_occasions = dim(CH)[2],nind = dim(CH)[1], prey_abund = prey_abund) #first check


## Parameters monitored
params <- c("mean_s2", "mean_eta", "mean_r", "mean_p","beta","mu_juvsurv")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 2


inits <- function(){list(mean_s2 = runif(1, 0.2, 1), mean_eta = runif(1, 0.05, 0.4), 
                         mean_p = runif(1, 0.05, 0.4), mean_r = runif(1, 0.05, 0.4))}  #z = ld.init(CH, f)
## Looks like I don't need to initialize the states...

mrr_juvad2 <- stan("mrr_w3covar.stan",
                  data = stan.data.R , init = inits, pars = params,
                  chains = nc, iter = ni, warmup = nb, thin = nt,
                  seed = 1)
print(mrr_juvad2, digits = 3)

params2 <- extract(mrr_juvad2)
save(params2, file ="param_logisticv2_chains.RData")

## Real prey data ptarmigan$Mean.density
minp = min(prey_abund)
maxp = max(prey_abund)
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
plot(prey_abund_vec,m,col="red",type="l",lwd=2,xlab="",ylab="Predator juvenile survival",xlim=c(minp,maxp),ylim=c(0,1))
lines(prey_abund_vec,q[1,],col="black")
lines(prey_abund_vec,q[2,],col="black")
fig_label("A")
### This is a bit puzzling. On average all these values are above the estimate without the covariate. 
d <- density(prey_abund)
plot(d, xlab = "Ptarmigan abundance (stdized) ",ylab="Kernel density of abundance",main="",xlim=c(minp,maxp),ylim=c(0,1))
polygon(d, col="red", border="black") 
fig_label("B")
dev.off()


