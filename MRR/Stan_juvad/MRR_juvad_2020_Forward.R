## Stan code for Juvenile-Adult MRR model designed for the Icelandic Gyrfalcon
## FB 25/03/2020 
.libPaths("/home/frederic/myRpackages/")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## function
source('../simul.R')

# -------------------------------------------------
# Loading data

CH=read.csv(file = "../../data/CaptureHistories2020.csv")
CH[,1]<-NULL
CH=as.matrix(CH)
image(CH)

get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)
n.occasions <- ncol(CH)
n.states <- 3
n.obs <- 3

stagemarking = read.csv(file ="../../data/stagemarking2020.csv")
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



# Bundle data
stan.data.R <- list(y = CH, f = f,  xj=xj, n_occasions = dim(CH)[2],nind = dim(CH)[1])


inits <- function(){list(mean_s = runif(2, 0.2, 1), mean_eta = runif(1, 0.05, 0.4), 
                         mean_p = runif(1, 0.05, 0.4), mean_r = runif(1, 0.05, 0.4))}  #z = ld.init(CH, f)
## Looks like I don't need to initialize the states...


## Parameters monitored
params <- c("mean_s", "mean_eta", "mean_r", "mean_p")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 2

## Call Stan from R
mrr_juvad <- stan("mrr.stan",
                 data = stan.data.R , init = inits, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt,
                 seed = 1,
                 open_progress = TRUE)
print(mrr_juvad, digits = 3)



## Other way of calling it
mrr_juvad <-stan_model("mrr.stan",verbose=TRUE)
fit.mrr_juvad <- sampling(mrr_juvad, data = stan.data.R, iter = 1000, init = inits,chains = 2, cores = 2)



